#!/usr/bin/env python
"""Runtime / peak-memory benchmark of upstream methods vs number of screened queries. NO ACCURACY: scores are
discarded; no label of any target (test query) is read or computed.

Screened sets: nested random draws from the DEV POOL only (config/split/dev_pool.txt): one permutation of the sorted
dev pool with default_rng(runtime_benchmark_draw_seed) (config/seeds.yaml); the first n ORFs are screened.
Targets: the test queries (config/split/test_queries.txt) — only their NAMES enter the view (ScreenView never holds
target rows).
Each (method, n) runs in a fresh subprocess (clean peak RSS = the child's ru_maxrss), 3 BLAS/OMP threads; a
watchdog kills a job whose RSS exceeds --mem-cap-gb (default 6) or that exceeds --timeout-s.

Outputs: results/runtime_benchmark/runtime.json (+ runtime.md table). Usage:
    .venv/bin/python scripts/02_runtime_benchmark.py                 # full grid
    .venv/bin/python scripts/02_runtime_benchmark.py --worker METHOD N OUT.json
"""
from __future__ import annotations

import hsl  # noqa: F401  (YEAST_SLs write guard + paths; must be first)

import argparse
import json
import os
import platform
import resource
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import yaml

from hsl import paths

OUT = paths.output_root() / "results" / "runtime_benchmark"
SIZES = (10, 30, 100, 300, 650, 1300, 2600)
METHODS = ("stack_lr", "stack_lr_wr2", "m1_hgb")
THREADS = 3


def _seed() -> int:
    return int(yaml.safe_load((paths.PROJECT_ROOT / "config" / "seeds.yaml").read_text())["runtime_benchmark_draw_seed"])


def _list(name: str) -> list[str]:
    return (paths.PROJECT_ROOT / "config" / "split" / name).read_text().split()


def screened_set(n: int) -> np.ndarray:
    dev = np.sort(np.asarray(_list("dev_pool.txt"), dtype=str))
    perm = np.random.default_rng(_seed()).permutation(len(dev))
    return np.sort(dev[perm[:n]])


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


def worker(method: str, n: int, out: Path) -> None:
    rec = {"method": method, "n_screened": n, "status": "error"}
    t0 = time.perf_counter()
    try:
        from hsl import upstream as U
        b = U.load_s1()
        rec["t_load_bundle_s"] = time.perf_counter() - t0
        scr = screened_set(n)
        test = _list("test_queries.txt")
        assert not set(scr) & set(test) and set(scr) <= set(_list("dev_pool.txt"))
        t1 = time.perf_counter()
        view = U.make_view(b, scr, test)
        rec["t_view_s"] = time.perf_counter() - t1
        rec["n_targets"], rec["n_arrays"] = int(view.n_targets), int(view.n_arrays)
        # training-side facts only (screened rows): needed to understand small-n failures
        rec["n_screened_on_array"] = int((view.train_col >= 0).sum())
        rec["n_SL_in_screened_rows"] = int(np.nansum(view.Y))
        S, info = U.run_method(view, method, n_threads=THREADS)
        rec["t_fit_predict_s"] = info["fit_predict_s"]
        rec["scores_finite"] = info["finite"]
        rec["score_shape"] = list(S.shape)
        rec["n_training_cells"] = info["method_info"].get("n_training_cells")
        del S
        rec["status"] = "ok"
    except Exception as e:  # noqa: BLE001  (record failures: whether a method runs at small n is a result)
        rec["error"] = f"{type(e).__name__}: {e}"
        rec["traceback"] = traceback.format_exc()[-3000:]
    rec["t_total_s"] = time.perf_counter() - t0
    rec["peak_rss_gb"] = peak_rss_gb()
    out.write_text(json.dumps(rec, indent=1))


def _rss_gb(pid: int) -> float:
    r = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    return int(r) / 2**20 if r else 0.0


def run_job(method: str, n: int, mem_cap: float, timeout: float) -> dict:
    tmp = OUT / "jobs" / f"{method}_n{n}.json"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, **{k: str(THREADS) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                                          "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")})
    t0 = time.time()
    p = subprocess.Popen([sys.executable, __file__, "--worker", method, str(n), str(tmp)], env=env,
                         cwd=paths.PROJECT_ROOT)
    max_seen, killed = 0.0, None
    while p.poll() is None:
        time.sleep(1.0)
        rss = _rss_gb(p.pid)
        max_seen = max(max_seen, rss)
        if rss > mem_cap:
            killed = f"killed: RSS {rss:.2f} GB > cap {mem_cap} GB"
        elif time.time() - t0 > timeout:
            killed = f"killed: wall time > {timeout:.0f} s"
        if killed:
            p.kill()
            p.wait()
            break
    if killed or not tmp.exists():
        rec = {"method": method, "n_screened": n, "status": "killed" if killed else "crashed",
               "error": killed or f"exit code {p.returncode}", "t_total_s": time.time() - t0}
    else:
        rec = json.loads(tmp.read_text())
    rec["watchdog_max_rss_gb"] = max_seen
    rec["wall_s_incl_startup"] = time.time() - t0
    return rec


def render_md(res: dict) -> str:
    rows = res["runs"]
    lines = ["# Runtime benchmark (no accuracy)", "",
             f"Generated by `scripts/02_runtime_benchmark.py` on {res['hardware']['machine']} "
             f"({res['hardware']['cpu_count']} CPUs), {THREADS} BLAS/OMP threads, one fresh process per cell. "
             "Screened sets = nested random draws from the dev pool; targets = the 656 test queries "
             "(names only; no target label is read). fit+predict = the method's `fit_predict` only (bundle load "
             "~1–2 s and view construction excluded); peak RSS = the child process's maximum resident set.", "",
             "| method | n screened | status | fit+predict (s) | peak RSS (GB) | SLs in screened rows |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        fp = f"{r['t_fit_predict_s']:.1f}" if "t_fit_predict_s" in r else "—"
        rss = r.get("peak_rss_gb", r.get("watchdog_max_rss_gb"))
        st = r["status"] if r["status"] == "ok" else f"**{r['status']}**: {r.get('error', '')[:120]}"
        lines.append(f"| {r['method']} | {r['n_screened']} | {st} | {fp} | {rss:.2f} | {r.get('n_SL_in_screened_rows', '—')} |")
    if res.get("skipped"):
        lines += ["", "Skipped: " + "; ".join(res["skipped"])]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", nargs=3, metavar=("METHOD", "N", "OUT"))
    ap.add_argument("--methods", nargs="*", default=list(METHODS))
    ap.add_argument("--sizes", nargs="*", type=int, default=list(SIZES))
    ap.add_argument("--mem-cap-gb", type=float, default=6.0)
    ap.add_argument("--timeout-s", type=float, default=3600)
    ap.add_argument("--slow-skip-s", type=float, default=1500,
                    help="skip larger sizes of a method once one size took longer than this (fit+predict)")
    a = ap.parse_args(argv)
    if a.worker:
        worker(a.worker[0], int(a.worker[1]), Path(a.worker[2]))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"description": "runtime and peak memory only; no accuracy computed",
           "draw_seed": _seed(), "threads": THREADS, "mem_cap_gb": a.mem_cap_gb,
           "targets": "config/split/test_queries.txt (656 test queries; names only)",
           "screened": "first n of default_rng(draw_seed).permutation(sorted dev pool) (nested)",
           "hardware": {"machine": platform.machine(), "system": f"{platform.system()} {platform.release()}",
                        "cpu_count": os.cpu_count(), "python": platform.python_version()},
           "started": time.strftime("%Y-%m-%dT%H:%M:%S"), "runs": [], "skipped": []}
    for m in a.methods:
        stop = None
        for n in a.sizes:
            if stop:
                res["skipped"].append(f"{m} n={n} ({stop})")
                continue
            rec = run_job(m, n, a.mem_cap_gb, a.timeout_s)
            print(json.dumps({k: rec.get(k) for k in ("method", "n_screened", "status", "t_fit_predict_s",
                                                      "peak_rss_gb", "error")}), flush=True)
            res["runs"].append({k: v for k, v in rec.items() if k != "traceback"} |
                               ({"traceback_tail": rec["traceback"][-800:]} if "traceback" in rec else {}))
            if rec.get("t_fit_predict_s", 0) > a.slow_skip_s:
                stop = f"previous size took {rec['t_fit_predict_s']:.0f} s > {a.slow_skip_s:.0f} s"
            elif rec["status"] == "killed":
                stop = f"previous size {rec['error']}"
            res["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
            (OUT / "runtime.json").write_text(json.dumps(res, indent=1))
            (OUT / "runtime.md").write_text(render_md(res))


if __name__ == "__main__":
    main()
