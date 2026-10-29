#!/usr/bin/env python
"""ONE-SHOT final evaluation on the held-out test set T (DESIGN v3.1 §8, §10) — and its dry run on a dev split.

Plan, arms, R, seeds and grids: config/finalists.yaml (hsl.final). This script is the ONLY place that loads the
unblinded bundle (final mode, inside launch.py jobs carrying HSL_FINAL=1).

    check   --mode final|dryrun           evaluate every refusal; no side effects (exit 1 if a final refusal fails)
    start   --mode M [--launch] [--max-concurrent k]
            REFUSES unless config/finalists.yaml is committed, the git tree is clean (ignored outputs allowed),
            config/frozen.yaml hashes match, <out>/.lock is absent, the YEAST_SLs read-only check passes and the
            upstream data checksums match (dry run: only the lock is enforced, the rest is reported). Then writes
            run/plan.json.tmp, the lock (O_EXCL: git hash, plan hash, frozen.yaml hash), run/plan.json,
            run/jobs.json and a ledger line; --launch execs launch.py (one PID, PID lock run/launcher.pid, output
            -> run/launch_<t>.out).
    work    --mode M --job NAME            one job (called by launch.py): selection saved first (or loaded, never
            recomputed), then each not-yet-done unit: predictions saved, then one row. `_labels` writes the E_fixed
            labels + label-only manifest for the score step.
    resume  --mode M [--launch]           crash re-run: lock must exist, plan == locked plan, frozen hashes still
            match, no live launcher holds run/launcher.pid; a ledger line; run/jobs_resume_<t>.json with the
            incomplete jobs only.
    status  --mode M                       units done / total by arm role
    score   --mode M [--partial]           -> scripts/score_final.py (separate step; no bundle, no HSL_FINAL)
    calibrate --pred P --n N               RSS/time calibration on the dev re-split (~1,050 targets; DEV ONLY)
    cost                                   projected job-hours / wall-clock of the real final from measured runtimes
Working files: <out>/run/ (gitignored: selections, parts, predictions, labels, boot cache, logs, ledger).

Per-job integrity (DESIGN deviation (iv), review F1): every final job checks HSL_FINAL=1, the lock, the plan hash and
every frozen-file hash (frozen.yaml: DESIGN §8 list + extended group + the clusters cache), not a whole-tree-clean
HEAD, so non-frozen notes may be committed during the run. `start` keeps the full §10 check.
RECOVERY
  * a job refused "frozen files changed": restore them (`git checkout <locked sha> -- <files>`), then `resume`.
  * jobs crashed / machine restarted: `resume --mode final --launch` (only incomplete jobs; selections are loaded,
    never recomputed; a truncated last row is dropped and its unit re-run).
  * `start` died after the lock but before run/plan.json: `resume` promotes run/plan.json.tmp if its hash equals the
    locked plan hash; otherwise stop and ask the PI (never delete the lock).
  * resume refuses "live launcher": a launcher is still running; wait for it or stop it first.
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
    """Final: the unblinded bundle after the per-job integrity check (F.guard_final). On a clean tree this is the
    frozen loader itself; on a tree with non-frozen edits (deviation (iv)) the same loader and assertions minus
    views.final_bundle's whole-tree-clean check. Dry run: dev bundle."""
    if P["mode"] != "final":
        return views.dev_bundle()
    F.guard_final(P)
    if not views.git_state()["dirty"]:
        return views.final_bundle()
    return _final_bundle_frozen_integrity(P)


def _final_bundle_frozen_integrity(P: dict):
    """= views.final_bundle() minus `git_state()["dirty"]` (DESIGN deviation (iv)); guard_final has verified the
    lock, plan and frozen hashes (views.py, guard.py and the ysl loader included). HSL_FINAL=1 opens the label gate."""
    if os.environ.get("HSL_FINAL") != "1":
        raise PermissionError("unblinded load requires HSL_FINAL=1")
    from hsl import paths
    from ysl.data.dense import load_bundle
    b = load_bundle(paths.costanzo_dataset(), unblind=True, cache=False)
    views._assert_full_names(b)
    b.meta = dict(b.meta, hsl_view="final_bundle", loaded_under="deviation_iv_frozen_integrity")
    F.ledger_event(P, "bundle_loaded_dirty_tree", git=views.git_state())
    return b.freeze()


def _print_checks(checks) -> None:
    for c in checks:
        print(f"  [{'ok' if c.ok else 'REFUSE'}] {c.name}: {c.detail}")


def cmd_check(a) -> int:
    cfg, P = F.mode_params(a.mode)
    checks = F.preflight(ROOT, P["lock"], deep=a.deep)
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


def start(mode: str, root: Path = ROOT, deep: bool = True) -> dict:
    """Refusals -> lock -> plan -> job file. Returns the lock record (raises F.Refusal)."""
    cfg, P = F.mode_params(mode, root)
    checks = F.preflight(root, P["lock"], deep=deep)
    F.check_global_knobs(cfg)
    hard = [c for c in checks if not c.ok and (P["enforce"] or c.name == "lock_absent")]
    if hard:
        raise F.Refusal("; ".join(f"{c.name}: {c.detail}" for c in hard))
    plan = F.build_plan(cfg, P)
    P["run"].mkdir(parents=True, exist_ok=True)
    (P["run"] / "plan.json.tmp").write_text(json.dumps(plan, indent=0) + "\n")   # before the lock (recovery)
    gs = views.git_state()
    lock = {"git_head": gs["head"], "git_dirty_at_start": gs["dirty"], "mode": mode, "split": P["split"],
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "plan_sha256": plan["sha256"],
            "n_units": plan["n_units"], "n_jobs": len(plan["jobs"]) + 1,
            "finalists_sha256": F.sha256_file(root / F.FINALISTS),
            "frozen_sha256": F.sha256_file(root / F.FROZEN) if (root / F.FROZEN).exists() else None,
            "preflight": [c.__dict__ for c in checks], "host": platform.node(), "python": sys.version.split()[0]}
    F.write_lock(P["lock"], lock)                       # O_EXCL: never overwritten
    os.replace(P["run"] / "plan.json.tmp", P["run"] / "plan.json")
    jobs = F.job_entries(plan, P)
    (P["run"] / "jobs.json").write_text(json.dumps(jobs, indent=1) + "\n")
    F.ledger_event(P, "start", git_head=gs["head"], plan_sha256=plan["sha256"], n_units=plan["n_units"],
                   n_jobs=len(jobs))
    return lock


def _launch_cmd(P: dict, jobs_file: Path, max_conc) -> tuple[list, dict]:
    cmd = [F.PY, "scripts/launch.py", str(jobs_file.relative_to(ROOT))]
    g = (P.get("compute") or {}).get("gate")
    if g:                                                            # D10 memory gate (opt-in in launch.py)
        cmd += ["--min-free-gb", str(g["headroom_gb"]), "--gate-max-jobs", str(g["max_jobs"]),
                "--gate-max-big", str(g["max_big"])]
    if max_conc:
        cmd += ["--max-concurrent", str(max_conc)]
    cmd += ["--pid-lock", str((P["run"] / "launcher.pid").relative_to(ROOT))]
    env = {**os.environ, "HSL_LOGS_DIR": str(P["run"] / "logs")}   # logs inside the ignored run dir (clean tree)
    env.pop("HSL_FINAL", None)                                       # only final jobs carry it (job env)
    return cmd, env


def _exec_or_print(P, jobs_file, a) -> int:
    cmd, env = _launch_cmd(P, jobs_file, a.max_concurrent)
    print("launch: HSL_LOGS_DIR=" + env["HSL_LOGS_DIR"] + " " + " ".join(cmd))
    if a.launch:
        os.chdir(ROOT)
        out = P["run"] / f"launch_{time.strftime('%Y%m%d_%H%M%S')}.out"   # launcher output inside the run dir
        print(f"launcher output -> {out}", flush=True)
        fd = os.open(out, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
        os.dup2(fd, 1)
        os.dup2(fd, 2)
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
    pidf = P["run"] / "launcher.pid"
    if pidf.exists():
        try:
            pid = int(pidf.read_text().strip() or 0)
        except ValueError:
            pid = 0
        if pid and F._pid_alive(pid):
            raise F.Refusal(f"a live launcher (pid {pid}) holds {pidf}: wait for it or stop it first")
    if P["enforce"]:
        c = F.check_frozen(root)                       # deviation (iv): frozen hashes, not a clean HEAD
        if not c.ok:
            raise F.Refusal(f"frozen files changed: {c.detail} (restore them, then resume)")
    if not (P["run"] / "plan.json").exists() and (P["run"] / "plan.json.tmp").exists():
        tmp = json.loads((P["run"] / "plan.json.tmp").read_text())
        if tmp["sha256"] != lk["plan_sha256"]:
            raise F.Refusal("run/plan.json missing and plan.json.tmp does not match the lock: ask the PI")
        os.replace(P["run"] / "plan.json.tmp", P["run"] / "plan.json")
        (P["run"] / "jobs.json").write_text(json.dumps(F.job_entries(tmp, P), indent=1) + "\n")
        F.ledger_event(P, "plan_promoted_from_tmp")
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
                   jobs_file=os.path.relpath(jf, root))
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


def cmd_calibrate(a) -> int:
    """RSS / time calibration for the larger target grid of T (PI 2026-10-29, DEV ONLY): one unit per process on the
    re-split check split (targets = inner folds 4+5, ~1,050 queries; candidates = folds 1-3), screened set = first n
    genes of a seeded permutation of its pool (dry-run trajectory seed 0). Nothing saved but one JSON row:
    peak RSS (ps-polled), fit+predict seconds -> results/final_dryrun/rss_calibration.jsonl."""
    import numpy as np
    from hsl import harness as H, p1
    cfg, Pd = F.mode_params("dryrun")
    sp = views.get_split(a.split)
    if sp.kind != "dev" or a.split == "final":
        raise SystemExit("calibration runs on a dev split only")
    n = min(int(a.n), len(sp.pool))
    rs = H.RandomSelector()
    sel_ss, ev_ss = H.stream_seeds(F.traj_seed(Pd, 0))
    rs.reset(sp.pool, np.random.default_rng(sel_ss))
    method, params = Pd["preds"][a.pred]
    with p1.PeakRSS() as pk:
        res = H.run_predictor(H.RunSpec(a.split, tuple(rs._perm[:n]), ("E_fixed",), ()), method, params,
                              {"model_seed": H.model_seed_of(ev_ss), "traj_index": 0}, save=False, ledger=False,
                              tag="rss_calibration")
    row = {"pred": a.pred, "method": method, "n": n, "n_requested": int(a.n), "split": a.split,
           "n_targets": int(len(sp.targets)), "peak_rss_gb": pk.peak, "fit_predict_s": res.manifest["runtime_s"]["fit_predict_s"],
           "total_s": res.manifest["runtime_s"]["total_s"], "n_threads": res.manifest["n_threads"],
           "t": time.strftime("%Y-%m-%dT%H:%M:%S")}
    p1.append_row(ROOT / "results" / "final_dryrun" / "rss_calibration.jsonl", row)   # no accuracy recorded
    print(json.dumps(row))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("check", "start", "work", "resume", "status", "score"):
        s = sub.add_parser(nm)
        s.add_argument("--mode", required=True, choices=list(F.MODES))
        if nm == "check":
            s.add_argument("--deep", action="store_true", help="also YEAST_SLs read-only + upstream checksums")
        if nm in ("start", "resume"):
            s.add_argument("--launch", action="store_true")
            s.add_argument("--max-concurrent", type=int, default=None)
        if nm == "work":
            s.add_argument("--job", required=True)
        if nm == "score":
            s.add_argument("--partial", action="store_true")
    sub.add_parser("cost")
    cb = sub.add_parser("calibrate")
    cb.add_argument("--pred", required=True)
    cb.add_argument("--n", type=int, required=True)
    cb.add_argument("--split", default="resplit", choices=["resplit", "dev1", "dev2", "dev3"])
    a = ap.parse_args(argv)
    return {"check": cmd_check, "start": cmd_start, "work": cmd_work, "resume": cmd_resume, "status": cmd_status,
            "score": cmd_score, "cost": cmd_cost, "calibrate": cmd_calibrate}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
