#!/usr/bin/env python
"""P0 dev smoke run (DEV DATA ONLY): dev split 1, R0 random selector, 1 trajectory (dev trajectory seed r = 0),
predictors stack_lr and stack_lr_wr2 evaluated at n = 30 and 100 on E_fixed_dev. Exercises the whole harness
(Trajectory -> run_predictor -> predictions saved -> metrics -> LEDGER). Output: results/smoke/dev_smoke.json.

    .venv/bin/python scripts/10_dev_smoke.py
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import json
import time

from hsl import config, harness as H, paths

CHECKPOINTS = (30, 100)
PREDICTORS = (("stack_lr", {}), ("stack_lr_wr2", {}))


def main() -> dict:
    t0 = time.perf_counter()
    seed = config.trajectory_seeds("dev", 1)[0]
    log = H.Trajectory(H.RandomSelector(), "dev1", seed, 0, arm="R0_random_smoke", checkpoints=CHECKPOINTS,
                       max_n=max(CHECKPOINTS), predictors=PREDICTORS).run()
    rows = []
    for n, d in log["checkpoints"].items():
        for method, r in d.items():
            m = r["metrics"]["E_fixed_dev"]
            rows.append({"n": int(n), "method": method, "auprc": m["auprc"], "lift": m["lift"],
                         "macro_auprc": m["macro_auprc"], "auroc": m["auroc"], "auprc_AO": m["auprc_AO"],
                         "auprc_VV": m["auprc_VV"], "fit_predict_s": r["runtime_s"].get("fit_predict_s"),
                         "total_s": r["runtime_s"]["total_s"], "run_id": r["run_id"]})
    any_m = next(iter(next(iter(log["checkpoints"].values())).values()))["metrics"]["E_fixed_dev"]
    out = {"split": "dev1", "traj_id": log["traj_id"], "E_fixed_dev": {k: any_m[k] for k in
                                                                       ("n_cells", "n_pos", "prevalence",
                                                                        "prev_AO", "prev_VV")},
           "results": rows, "wall_s": time.perf_counter() - t0, "peak_rss_gb": log["peak_rss_gb"],
           "selection_s_total": sum(r["select_s"] for r in log["rounds"])}
    d = paths.output_root() / "results" / "smoke"
    d.mkdir(parents=True, exist_ok=True)
    (d / "dev_smoke.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
