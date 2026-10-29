#!/usr/bin/env python
"""Run a list of jobs with bounded concurrency, per-job timeouts, RSS watchdog and per-job logs (DESIGN v3 §4, §12).

Jobs file: JSON list (or JSON-lines) of {"name": str, "cmd": [argv...] | "shell string", "n": int (largest budget
the job evaluates; optional), "timeout_s": int (optional), "env": {..} (optional)}.

Concurrency policy (config/protocol.yaml `compute`):
  * small jobs (n <= small_job_max_n): up to small_job_concurrency at once, small_job_threads threads each, while
    the summed RSS of running jobs stays <= rss_cap_gb_total;
  * other jobs: at most max_concurrent_default running jobs in total, threads_default threads;
  * n >= alone_min_n: runs alone (RSS cap rss_cap_gb_alone).
Optional MEMORY GATE (opt-in, `--min-free-gb H`; D10, 2026-10-29): replaces the class limits above. A job starts only
if fewer than --gate-max-jobs (default 6) run, at most --gate-max-big (default 2) jobs with n >= alone_min_n run, and
the machine's unused memory (macOS vm_stat free + speculative pages; Linux MemAvailable) minus the memory running jobs
have not yet claimed (their `rss_gb` hint minus current RSS) is >= the job's `rss_gb` hint + H. With nothing running a
job always starts (progress). Jobs may carry "rss_gb" (expected peak RSS; default rss_cap_gb_per_job) and
"rss_cap_gb" (kill cap overriding the class cap). Without the flag, behaviour is unchanged.
Threads are passed as OMP/OPENBLAS/MKL/VECLIB/NUMEXPR_NUM_THREADS and HSL_THREADS. A job exceeding its timeout or
its RSS cap is killed (process group). stdout+stderr -> logs/<name>.log; one JSON line per job ->
logs/launch_<timestamp>.jsonl. Exit code = number of failed jobs (capped at 255).

    .venv/bin/python scripts/launch.py jobs.json [--max-concurrent 2] [--dry-run]
        [--min-free-gb 2 --gate-max-jobs 6 --gate-max-big 2]
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import argparse
import json
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

from hsl import config, paths

LOGS = Path(os.environ.get("HSL_LOGS_DIR") or paths.PROJECT_ROOT / "logs")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS", "HSL_THREADS")


def load_jobs(path: Path) -> list[dict]:
    text = Path(path).read_text().strip()
    jobs = json.loads(text) if text.startswith("[") else [json.loads(ln) for ln in text.splitlines() if ln.strip()]
    names = [j["name"] for j in jobs]
    if len(set(names)) != len(names):
        raise ValueError("job names must be unique")
    return jobs


def job_class(job: dict, c: dict) -> str:
    n = job.get("n")
    if n is not None and int(n) >= int(c["alone_min_n"]):
        return "alone"
    if n is not None and int(n) <= int(c["small_job_max_n"]):
        return "small"
    return "large"


def threads_for(cls: str, c: dict) -> int:
    return int(c["small_job_threads"] if cls == "small" else c["threads_default"])


def rss_cap(cls: str, c: dict) -> float:
    return float(c["rss_cap_gb_alone"] if cls == "alone" else c["rss_cap_gb_per_job"])


def tree_rss_gb(pid: int) -> float:
    """RSS of a process and its descendants (ps)."""
    try:
        out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True).stdout
    except OSError:
        return 0.0
    kids, rss = {}, {}
    for ln in out.splitlines():
        p = ln.split()
        if len(p) == 3:
            kids.setdefault(int(p[1]), []).append(int(p[0]))
            rss[int(p[0])] = int(p[2])
    tot, stack = 0, [pid]
    while stack:
        q = stack.pop()
        tot += rss.get(q, 0)
        stack.extend(kids.get(q, []))
    return tot / 2**20


def can_start(cls: str, running: list[dict], c: dict, max_conc: int | None) -> bool:
    if not running:
        return True
    if cls == "alone" or any(r["cls"] == "alone" for r in running):
        return False
    total_rss = sum(r.get("rss", 0.0) for r in running)
    if cls == "small" and all(r["cls"] == "small" for r in running):
        limit = int(c["small_job_concurrency"])
        if total_rss + rss_cap("small", c) / 2 > float(c["rss_cap_gb_total"]):
            return False
    else:
        limit = int(c["max_concurrent_default"])
    if max_conc is not None:
        limit = min(limit, int(max_conc))
    return len(running) < limit


def unused_gb() -> float | None:
    """Unused physical memory (GB): macOS vm_stat (free + speculative pages); Linux /proc/meminfo MemAvailable."""
    try:
        out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=10).stdout
        ps = int(re.search(r"page size of (\d+) bytes", out).group(1))
        pages = {}
        for ln in out.splitlines()[1:]:
            k, _, v = ln.partition(":")
            v = v.strip().rstrip(".")
            if v.isdigit():
                pages[k.strip()] = int(v)
        return (pages.get("Pages free", 0) + pages.get("Pages speculative", 0)) * ps / 2**30
    except (OSError, AttributeError, ValueError, subprocess.SubprocessError):
        pass
    try:
        for ln in Path("/proc/meminfo").read_text().splitlines():
            if ln.startswith("MemAvailable:"):
                return int(ln.split()[1]) / 2**20
    except OSError:
        pass
    return None


def job_hint(job: dict, c: dict) -> float:
    return float(job.get("rss_gb", c["rss_cap_gb_per_job"]))


def gate_ok(job: dict, running: list[dict], c: dict, headroom: float, max_jobs: int, max_big: int,
            free_fn=unused_gb) -> bool:
    """Memory gate (see module docstring)."""
    if not running:
        return True
    if len(running) >= int(max_jobs):
        return False
    big = lambda j: j.get("n") is not None and int(j["n"]) >= int(c["alone_min_n"])  # noqa: E731
    if big(job) and sum(big(r) for r in running) >= int(max_big):
        return False
    free = free_fn()
    if free is None:
        return False
    pending = sum(max(0.0, float(r.get("rss_hint", 0.0)) - float(r.get("rss", 0.0))) for r in running)
    return free - pending >= job_hint(job, c) + float(headroom)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("jobs", type=Path)
    ap.add_argument("--max-concurrent", type=int, default=None, help="further cap on concurrency")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--min-free-gb", type=float, default=None,
                    help="enable the memory gate with this headroom (GB) above each job's rss_gb hint")
    ap.add_argument("--gate-max-jobs", type=int, default=6)
    ap.add_argument("--gate-max-big", type=int, default=2)
    a = ap.parse_args(argv)
    c = config.protocol()["compute"]
    jobs = load_jobs(a.jobs)
    LOGS.mkdir(exist_ok=True)
    summary = LOGS / f"launch_{time.strftime('%Y%m%d_%H%M%S')}.jsonl"
    queue = list(jobs)
    running: list[dict] = []
    failed = 0
    poll = float(c["poll_s"])
    while queue or running:
        def startable(job):
            if a.min_free_gb is None:
                return can_start(job_class(job, c), running, c, a.max_concurrent)
            lim = a.gate_max_jobs if a.max_concurrent is None else min(a.gate_max_jobs, a.max_concurrent)
            return gate_ok(job, running, c, a.min_free_gb, lim, a.gate_max_big)

        while queue and startable(queue[0]):
            job = queue.pop(0)
            cls = job_class(job, c)
            th = threads_for(cls, c)
            env = {**os.environ, **{k: str(th) for k in THREAD_VARS},
                   **{k: str(v) for k, v in job.get("env", {}).items()}}  # job env overrides thread caps
            if not job.get("env", {}).get("HSL_FINAL"):
                env.pop("HSL_FINAL", None)              # final mode only when a job asks for it explicitly
            cmd = job["cmd"] if isinstance(job["cmd"], list) else shlex.split(job["cmd"])
            log = LOGS / f"{job['name']}.log"
            rec = {"name": job["name"], "cmd": cmd, "cls": cls, "threads": th, "log": os.path.relpath(
                log, paths.PROJECT_ROOT), "timeout_s": int(job.get("timeout_s", c["job_timeout_s"])),
                   "rss_cap_gb": float(job.get("rss_cap_gb", rss_cap(cls, c))), "n": job.get("n"),
                   "rss_hint": job_hint(job, c)}
            print(f"[launch] start {job['name']} ({cls}, {th} threads)", flush=True)
            if a.dry_run:
                rec.update(status="dry-run")
                with open(summary, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                continue
            fh = open(log, "w")
            proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=paths.PROJECT_ROOT,
                                    start_new_session=True)
            running.append({**rec, "proc": proc, "fh": fh, "t0": time.time(), "rss": 0.0, "peak_rss_gb": 0.0})
        if not running:
            continue
        time.sleep(poll)
        for r in list(running):
            p = r["proc"]
            r["rss"] = tree_rss_gb(p.pid)
            r["peak_rss_gb"] = max(r["peak_rss_gb"], r["rss"])
            reason = None
            if p.poll() is None:
                if time.time() - r["t0"] > r["timeout_s"]:
                    reason = "timeout"
                elif r["rss"] > r["rss_cap_gb"]:
                    reason = "rss_cap"
                if reason:
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
                else:
                    continue
            r["fh"].close()
            status = reason or ("ok" if p.returncode == 0 else f"exit_{p.returncode}")
            failed += status != "ok"
            out = {k: v for k, v in r.items() if k not in ("proc", "fh", "rss", "t0")}
            out.update(status=status, wall_s=round(time.time() - r["t0"], 1), returncode=p.returncode)
            with open(summary, "a") as f:
                f.write(json.dumps(out) + "\n")
            print(f"[launch] {r['name']}: {status} ({out['wall_s']} s, peak {out['peak_rss_gb']:.2f} GB)", flush=True)
            running.remove(r)
    return min(failed, 255)


if __name__ == "__main__":
    sys.exit(main())
