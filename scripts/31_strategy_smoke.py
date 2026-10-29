#!/usr/bin/env python
"""Strategy SMOKE trajectories (pipeline proof only — NOT a comparison, NOT for decisions): one trajectory per
strategy arm (R1, H, H-static, F, F-static, C, M) on dev split 1 to n = 100 with predictor stack_lr evaluated at
n = 30, 60, 100 on E_fixed_dev. Trajectory seed = SeedSequence(strategy_smoke_seed).spawn(1)[0] (disjoint from the
dev/final trajectory bases). ≤ 2 threads. Output: results/strategies/smoke_trajectories.{json,md}; trajectory logs and
LEDGER lines (arm = smoke_<strategy>) via the harness; per-round selector diagnostics next to each trajectory log.

    .venv/bin/python scripts/31_strategy_smoke.py
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "2")
os.environ.setdefault("HSL_THREADS", "2")

import hsl  # noqa: E402,F401  (guard + label gate)

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

from hsl import config, harness as H, paths, strategies  # noqa: E402

CHECKPOINTS = (30, 60, 100)
PREDICTORS = (("stack_lr", {}),)


def main():
    seed = np.random.SeedSequence(config.seed("strategy_smoke_seed")).spawn(1)[0]
    out = {"label": "SMOKE — pipeline proof only, not for decisions", "split": "dev1", "predictor": "stack_lr",
           "checkpoints": list(CHECKPOINTS), "arms": {}}
    on_array = set(_on_array())
    for name in strategies.ARMS:
        t0 = time.perf_counter()
        sel = strategies.make(name)
        log = strategies.run_trajectory(H.Trajectory(sel, "dev1", seed, 0, arm=f"smoke_{name}",
                                                     checkpoints=CHECKPOINTS, max_n=max(CHECKPOINTS),
                                                     predictors=PREDICTORS))
        sel_s = [r["select_s"] for r in log["rounds"]]
        out["arms"][name] = {
            "traj_id": log["traj_id"],
            "auprc": {n: d["stack_lr"]["metrics"]["E_fixed_dev"]["auprc"] for n, d in log["checkpoints"].items()},
            "select_s_mean": float(np.mean(sel_s)), "select_s_max": float(np.max(sel_s)),
            "phases": [d.get("phase") for d in sel.diagnostics],
            "selector_diagnostics_file": log.get("selector_diagnostics_file"),
            "n_offarray_picked": int(sum(1 for g in log["selected"] if g not in on_array)),
            "wall_s": time.perf_counter() - t0}
        print(name, json.dumps(out["arms"][name]["auprc"]), f"sel {np.mean(sel_s):.2f}s/round", flush=True)
    d = paths.output_root() / "results" / "strategies"
    d.mkdir(parents=True, exist_ok=True)
    (d / "smoke_trajectories.json").write_text(json.dumps(out, indent=1) + "\n")
    lines = ["# Strategy smoke trajectories (dev split 1) — SMOKE, NOT FOR DECISIONS", "",
             "One trajectory per arm, stack_lr, E_fixed_dev pooled AUPRC. Proves the pipeline end to end "
             "(SelectorView -> selector -> harness -> predictor -> metrics); single trajectories carry no comparative "
             "information (no random arm here by design; comparisons are P4, after the P* freeze).", "",
             "| arm | " + " | ".join(f"AUPRC@{n}" for n in CHECKPOINTS) + " | select s/round (mean / max) | off-array picked |",
             "|---|" + "---|" * (len(CHECKPOINTS) + 2)]
    for name, a in out["arms"].items():
        lines.append(f"| {name} | " + " | ".join(f"{a['auprc'][str(n)]:.4f}" for n in CHECKPOINTS)
                     + f" | {a['select_s_mean']:.2f} / {a['select_s_max']:.2f} | {a['n_offarray_picked']} |")
    (d / "smoke_trajectories.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _on_array():
    from hsl import views
    return views.dev_bundle().array_orfs


if __name__ == "__main__":
    main()
