#!/usr/bin/env python
"""ONE-SHOT final evaluation on the held-out test set T (DESIGN v3.1 §8, §10) — and its dry run on a dev split.

Plan, arms, R, seeds and grids: config/finalists.yaml (hsl.final). This script is the ONLY place that loads the
unblinded bundle (views.final_bundle, final mode, inside launch.py jobs carrying HSL_FINAL=1).

    check   --mode final|dryrun           evaluate every refusal; no side effects (exit 1 if a final refusal fails)
    start   --mode M [--launch] [--max-concurrent k]
            REFUSES unless config/finalists.yaml is committed, the git tree is clean (ignored outputs allowed),
            config/frozen.yaml hashes match, and <out>/.lock is absent (dry run: only the lock is enforced, the rest
            is reported). Then writes the lock (git hash, plan hash), run/plan.json, run/jobs.json for launch.py and
            a ledger line; --launch execs launch.py (one PID).
    work    --mode M --job NAME            one job (called by launch.py): selection saved first (or loaded, never
            recomputed), then each not-yet-done unit: predictions saved, then one row. `_labels` writes the E_fixed
            labels + label-only manifest for the score step.
    resume  --mode M [--launch]           crash re-run: lock must exist, HEAD == locked commit, clean tree, frozen
            hashes still match; a ledger line; run/jobs_resume_<t>.json with the incomplete jobs only.
    status  --mode M                       units done / total by arm role
    score   --mode M [--partial]           -> scripts/score_final.py (separate step; no bundle, no HSL_FINAL)
    cost                                   projected job-hours / wall-clock of the real final from measured runtimes
Working files: <out>/run/ (gitignored: selections, parts, predictions, labels, boot cache, logs, ledger).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path

from hsl import final as F, views

ROOT = F.ROOT


def get_bundle(P: dict):
    """Final: the unblinded bundle (requires HSL_FINAL=1 + clean tree, enforced by views); dry run: dev bundle."""
    if P["mode"] == "final":
        F.guard_final(P)
        return views.final_bundle()
    return views.dev_bundle()


def _print_checks(checks) -> None:
    for c in checks:
        print(f"  [{'ok' if c.ok else 'REFUSE'}] {c.name}: {c.detail}")


def cmd_check(a) -> int:
    cfg, P = F.mode_params(a.mode)
    checks = F.preflight(ROOT, P["lock"])
    print(f"preflight ({a.mode}; refusals {'enforced' if P['enforce'] else 'reported only, lock enforced'}):")
    _print_checks(checks)
    try:
        F.check_global_knobs(cfg)
        print("  [ok] global selector knobs match finalists.yaml")
    except F.Refusal as e:
        print(f"  [REFUSE] {e}")
        return 1
    bad = [c for c in checks if not c.ok and (P["enforce"] or c.name == "lock_absent")]
    return 1 if bad else 0


def start(mode: str, root: Path = ROOT) -> dict:
    """Refusals -> lock -> plan -> job file. Returns the lock record (raises F.Refusal)."""
    cfg, P = F.mode_params(mode, root)
    checks = F.preflight(root, P["lock"])
    F.check_global_knobs(cfg)
    hard = [c for c in checks if not c.ok and (P["enforce"] or c.name == "lock_absent")]
    if hard:
        raise F.Refusal("; ".join(f"{c.name}: {c.detail}" for c in hard))
    plan = F.build_plan(cfg, P)
    gs = views.git_state()
    lock = {"git_head": gs["head"], "git_dirty_at_start": gs["dirty"], "mode": mode, "split": P["split"],
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "plan_sha256": plan["sha256"],
            "n_units": plan["n_units"], "n_jobs": len(plan["jobs"]) + 1,
            "finalists_sha256": F.sha256_file(root / F.FINALISTS),
            "frozen_sha256": F.sha256_file(root / F.FROZEN) if (root / F.FROZEN).exists() else None,
            "preflight": [c.__dict__ for c in checks], "host": platform.node(), "python": sys.version.split()[0]}
    F.write_lock(P["lock"], lock)                       # O_EXCL: never overwritten
    P["run"].mkdir(parents=True, exist_ok=True)
    (P["run"] / "plan.json").write_text(json.dumps(plan, indent=0) + "\n")
    jobs = F.job_entries(plan, P)
    (P["run"] / "jobs.json").write_text(json.dumps(jobs, indent=1) + "\n")
    F.ledger_event(P, "start", git_head=gs["head"], plan_sha256=plan["sha256"], n_units=plan["n_units"],
                   n_jobs=len(jobs))
    return lock


def _launch_cmd(P: dict, jobs_file: Path, max_conc) -> tuple[list, dict]:
    cmd = [F.PY, "scripts/launch.py", str(jobs_file.relative_to(ROOT))]
    if max_conc:
        cmd += ["--max-concurrent", str(max_conc)]
    env = {**os.environ, "HSL_LOGS_DIR": str(P["run"] / "logs")}   # logs inside the ignored run dir (clean tree)
    env.pop("HSL_FINAL", None)                                       # only final jobs carry it (job env)
    return cmd, env


def _exec_or_print(P, jobs_file, a) -> int:
    cmd, env = _launch_cmd(P, jobs_file, a.max_concurrent)
    print("launch: HSL_LOGS_DIR=" + env["HSL_LOGS_DIR"] + " " + " ".join(cmd))
    if a.launch:
        os.chdir(ROOT)
        os.execvpe(cmd[0], cmd, env)                  # one PID (LESSONS 8)
    return 0


def cmd_start(a) -> int:
    cfg, P = F.mode_params(a.mode)
    try:
        lock = start(a.mode)
    except F.Refusal as e:
        print(f"REFUSED: {e}")
        return 2
    print(f"lock written: {P['lock']} (git {lock['git_head'][:10]}, {lock['n_units']} units, {lock['n_jobs']} jobs)")
    return _exec_or_print(P, P["run"] / "jobs.json", a)


def cmd_work(a) -> int:
    cfg, P = F.mode_params(a.mode)
    if F.read_lock(P["lock"]) is None:
        print(f"REFUSED: no lock at {P['lock']}")
        return 2
    plan = F.load_plan(P)
    try:
        F.guard_final(P, plan)
        F.check_global_knobs(cfg)
        b = get_bundle(P)
        if a.job == F.LABELS_JOB:
            man = F.write_labels(P, b)
            F.ledger_event(P, "labels", eset=man["eset"], eset_digest=man["eset_digest"])
            print(f"[final] labels written ({man['counts']['n_cells']} cells)")   # counts only
            return 0
        F.work(P, plan, a.job, b)
    except F.Refusal as e:
        print(f"REFUSED: {e}")
        return 2
    except Exception as e:                                # noqa: BLE001  (recorded, then re-raised)
        F.ledger_event(P, "job_fail", job=a.job, error=f"{type(e).__name__}: {e}")
        raise
    return 0


def resume(mode: str, root: Path = ROOT) -> tuple[dict, Path, list]:
    cfg, P = F.mode_params(mode, root)
    lk = F.read_lock(P["lock"])
    if lk is None:
        raise F.Refusal(f"no lock at {P['lock']}: nothing to resume (use start)")
    if P["enforce"]:
        bad = [c for c in (F.check_finalists_committed(root), F.check_tree_clean(root), F.check_frozen(root))
               if not c.ok]
        if bad:
            raise F.Refusal("; ".join(f"{c.name}: {c.detail}" for c in bad))
        if views.git_state()["head"] != lk["git_head"]:
            raise F.Refusal("HEAD differs from the locked commit")
    plan = F.load_plan(P)
    if plan["sha256"] != lk["plan_sha256"]:
        raise F.Refusal("run/plan.json differs from the locked plan")
    prog = F.progress(P, plan)
    todo = [j for j, (d, t) in prog.items() if d < t]
    labels_missing = not (P["run"] / "labels.npz").exists()
    names = set(todo) | ({F.LABELS_JOB} if labels_missing else set())
    jobs = F.job_entries(plan, P, names)
    jf = P["run"] / f"jobs_resume_{time.strftime('%Y%m%d_%H%M%S')}.json"
    jf.write_text(json.dumps(jobs, indent=1) + "\n")
    F.ledger_event(P, "resume", n_jobs_incomplete=len(todo), labels_missing=labels_missing,
                   units_done=sum(d for d, _ in prog.values()), units_total=sum(t for _, t in prog.values()),
                   jobs_file=str(jf.relative_to(root)))
    return P, jf, jobs


def cmd_resume(a) -> int:
    try:
        P, jf, jobs = resume(a.mode)
    except F.Refusal as e:
        print(f"REFUSED: {e}")
        return 2
    print(f"resume: {len(jobs)} incomplete jobs -> {jf}")
    return _exec_or_print(P, jf, a)


def cmd_status(a) -> int:
    cfg, P = F.mode_params(a.mode)
    plan = F.load_plan(P)
    prog = F.progress(P, plan)
    by = {}
    roles = {cid: d["role"] for cid, d in plan["arms"].items()}
    for j in plan["jobs"]:
        d, t = prog[j["name"]]
        k = (roles[j["cid"]], j["group"])
        x = by.setdefault(k, [0, 0, 0])
        x[0] += d
        x[1] += t
        x[2] += d == t
    for (role, g), (d, t, jd) in sorted(by.items()):
        print(f"{role:10s} {g:5s} units {d:5d}/{t:5d}")
    print(f"total units {sum(d for d, _ in prog.values())}/{sum(t for _, t in prog.values())}; labels "
          f"{'yes' if (P['run'] / 'labels.npz').exists() else 'no'}")
    return 0


def cmd_score(a) -> int:
    import runpy
    sys.argv = ["score_final.py", "--mode", a.mode] + (["--partial"] if a.partial else [])
    runpy.run_path(str(ROOT / "scripts" / "score_final.py"), run_name="__main__")
    return 0


def cmd_cost(a) -> int:
    import runpy
    sys.argv = ["score_final.py", "--mode", "final", "--cost"]
    runpy.run_path(str(ROOT / "scripts" / "score_final.py"), run_name="__main__")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("check", "start", "work", "resume", "status", "score"):
        s = sub.add_parser(nm)
        s.add_argument("--mode", required=True, choices=list(F.MODES))
        if nm in ("start", "resume"):
            s.add_argument("--launch", action="store_true")
            s.add_argument("--max-concurrent", type=int, default=None)
        if nm == "work":
            s.add_argument("--job", required=True)
        if nm == "score":
            s.add_argument("--partial", action="store_true")
    sub.add_parser("cost")
    a = ap.parse_args(argv)
    return {"check": cmd_check, "start": cmd_start, "work": cmd_work, "resume": cmd_resume, "status": cmd_status,
            "score": cmd_score, "cost": cmd_cost}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
