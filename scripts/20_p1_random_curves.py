#!/usr/bin/env python
"""P1 dev random curves (DESIGN v3.1 §5, §7.1, §8; DEV DATA ONLY).

For dev split k ∈ {1,2,3}: R random permutations of the candidate pool (R0 trajectories, seeds
SeedSequence(dev_traj_base_seed).spawn(R)[r]; identical to harness.Trajectory + RandomSelector), budgets = budgets.yaml
random_grid truncated at the pool size plus the pool maximum; predictors stack_lr_wr2 (P*), stack_lr, B1_array_rate,
B2_product; m1_hgb (genome dropped) on split 1 only, 5 permutations, budgets m1_hgb_q1_grid (2627 -> pool max).
Every run: harness.run_predictor (ledger line; float32 E_fixed_dev cell predictions in results/predictions/p1/).
Rows -> results/p1/parts/<job>.jsonl; `merge` -> results/p1/curves.csv.

    run    --split 1 --perms 0-9 --budgets 10,30 --predictors stack_lr_wr2,stack_lr [--upstream-at 660] --job NAME
    jobs   --out jobs.json [--splits 1,2,3] [--perms 0-9] [--segments small,large,alone] [--m1] [--pilot]
    merge  (all parts -> results/p1/curves.csv; upstream-protocol rows -> results/p1/upstream_protocol.csv)
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from hsl import config, p1, views

PARTS = p1.P1_DIR / "parts"
BASE_PREDICTORS = ("stack_lr_wr2", "stack_lr", "B1_array_rate", "B2_product")
PY = ".venv/bin/python"


def _ints(s: str) -> list[int]:
    out = []
    for tok in str(s).split(","):
        if "-" in tok:
            a, b = tok.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif tok.strip():
            out.append(int(tok))
    return out


def _done(part: Path) -> set:
    return {(r["split"], r["perm"], r["n"], r["predictor"], r.get("kind", "curve")) for r in p1.read_rows([part])}


def cmd_run(a) -> int:
    p1.use_p1_prediction_dir()
    split = f"dev{a.split}"
    part = PARTS / f"{a.job}.jsonl"
    done = _done(part)
    preds = [m for m in a.predictors.split(",") if m]
    up_methods = [m for m in a.upstream_methods.split(",") if m]
    pool_n = len(views.get_split(split).pool)
    budgets = [pool_n if tok == "max" else int(tok) for tok in a.budgets.split(",") if tok.strip()] \
        if a.budgets else p1.dev_grid(pool_n)
    if max(budgets) > pool_n:
        raise ValueError(f"budget above the pool size {pool_n}")
    fails = 0
    for r in _ints(a.perms):
        t = p1.random_traj(split, r)
        p1.write_traj_log(t)
        for n in sorted(budgets):
            for m in preds:
                if (split, r, n, m, "curve") in done:
                    continue
                t0 = time.time()
                try:
                    row = p1.run_cell(t, n, m)
                    row["kind"] = "curve"
                    p1.append_row(part, row)
                    print(f"[p1] {split} r{r} n{n} {m}: AUPRC {row['auprc']:.4f} fit {row['fit_predict_s']:.1f}s "
                          f"peak {row['peak_rss_fit_gb']:.2f} GB", flush=True)
                except Exception as e:  # noqa: BLE001  (ledger line already written by the harness)
                    fails += 1
                    print(f"[p1] FAIL {split} r{r} n{n} {m}: {type(e).__name__}: {e} ({time.time() - t0:.0f}s)",
                          flush=True)
            for um in (up_methods if a.upstream_at and n == a.upstream_at else []):
                if (split, r, n, um, "upstream") in done:
                    continue
                row = p1.run_upstream_protocol(t, n, um)
                row["kind"] = "upstream"
                p1.append_row(part, row)
                print(f"[p1] upstream-protocol {um} {split} r{r} n{n}: E_upstream_dev AUPRC "
                      f"{row['E_upstream_dev.auprc']:.4f} vs E_fixed_dev {row['E_fixed_dev.auprc']:.4f}", flush=True)
    return min(fails, 255)


def _job(name, split, perms, budgets, preds, n_max, extra=()) -> dict:
    cmd = [PY, "scripts/20_p1_random_curves.py", "run", "--split", str(split), "--perms", perms,
           "--budgets", ",".join(str(b) for b in budgets), "--predictors", ",".join(preds),
           "--job", name, *extra]
    return {"name": name, "cmd": cmd, "n": int(n_max)}


def segments(pool_n: int, grid) -> dict:
    c = config.protocol()["compute"]
    small = [n for n in grid if n <= int(c["small_job_max_n"])]
    alone = [n for n in grid if n >= int(c["alone_min_n"])]
    large = [n for n in grid if n not in small and n not in alone]
    return {"small": small, "large": large, "alone": alone}


def cmd_jobs(a) -> int:
    jobs = {"small": [], "large": [], "alone": []}
    up = int(config.budgets()["q2_checkpoints"][-1])           # 660 ≈ 20% of 3,283
    for k in _ints(a.splits):
        split = f"dev{k}"
        pool_n = len(views.get_split(split).pool)
        seg = segments(pool_n, p1.dev_grid(pool_n))
        for s in a.segments.split(","):
            if not seg[s]:
                continue
            for r in _ints(a.perms):
                extra = ["--upstream-at", str(up)] if (k == 1 and s == "small" and up in seg[s]) else []
                jobs[s].append(_job(f"p1_{split}_r{r}_{s}", k, str(r), seg[s], BASE_PREDICTORS, max(seg[s]), extra))
        if a.m1 and k == 1:
            seg = segments(pool_n, p1.m1_grid(pool_n))
            R = int(config.budgets()["m1_hgb_q1_trajectories"])
            for s in ("small", "large", "alone"):
                for r in range(R):
                    if seg[s]:
                        jobs[s].append(_job(f"p1_{split}_r{r}_m1_{s}", k, str(r), seg[s], ("m1_hgb",), max(seg[s])))
    out = jobs["small"] + jobs["large"] + jobs["alone"]
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"{len(out)} jobs -> {a.out} ({', '.join(f'{k}: {len(v)}' for k, v in jobs.items())})")
    return 0


def cmd_merge(a) -> int:
    rows = p1.read_rows(sorted(PARTS.glob("*.jsonl")))
    cur = pd.DataFrame([r for r in rows if r.get("kind", "curve") == "curve"])
    ups = pd.DataFrame([r for r in rows if r.get("kind") == "upstream"])
    if len(cur):
        cur = cur.drop_duplicates(["split", "perm", "n", "predictor"], keep="last")
        cur = cur.sort_values(["split", "predictor", "perm", "n"])
        cur.to_csv(p1.P1_DIR / "curves.csv", index=False)
    if len(ups):
        ups.drop_duplicates(["split", "perm", "n", "predictor"], keep="last").sort_values(
            ["split", "predictor", "perm", "n"]).to_csv(p1.P1_DIR / "upstream_protocol.csv", index=False)
    print(f"curves: {len(cur)} rows; upstream-protocol: {len(ups)} rows")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--split", type=int, required=True)
    r.add_argument("--perms", required=True)
    r.add_argument("--budgets", default="")
    r.add_argument("--predictors", default=",".join(BASE_PREDICTORS))
    r.add_argument("--upstream-at", type=int, default=0)
    r.add_argument("--upstream-methods", default=p1.PSTAR)
    r.add_argument("--job", required=True)
    j = sub.add_parser("jobs")
    j.add_argument("--out", required=True)
    j.add_argument("--splits", default="1,2,3")
    j.add_argument("--perms", default="0-9")
    j.add_argument("--segments", default="small,large,alone")
    j.add_argument("--m1", action="store_true")
    sub.add_parser("merge")
    a = ap.parse_args(argv)
    return {"run": cmd_run, "jobs": cmd_jobs, "merge": cmd_merge}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
