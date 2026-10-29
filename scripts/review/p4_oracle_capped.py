#!/usr/bin/env python
"""P4 review (independent, DEV ONLY, evaluation-only diagnostic): H's own selection machinery (randomised top-k greedy
+ per-batch cluster cap, same seeds) driven by the TRUE query SL degree instead of q̂ (hsl.oracles.true_query_degree,
dev labels; never a strategy). Separates 'knowing the degree' from 'how the hubs are picked' (P3's oracle_hub takes
the strict global top-n with no cap and no randomisation). P* at n = 160; nothing saved, no LEDGER line.
Output: results/review_p4/oracle_capped.jsonl.
"""
from __future__ import annotations

import hsl  # noqa: F401

import json
import sys
from pathlib import Path

import numpy as np

from hsl import harness as H, oracles, q2dev, views
from hsl.strategies.hub import Hub

OUT = Path(__file__).resolve().parents[2] / "results" / "review_p4"


class CappedOracleHub(Hub):
    name = "oracle_hub_capped"

    def __init__(self, split, **kw):
        super().__init__(static=False, **kw)
        self._split = split

    def scores(self, view):
        deg = np.asarray(oracles.true_query_degree(self._split, np.asarray(view.candidates)), dtype=np.float64)
        return np.log1p(deg), "oracle_true_degree"


def main(argv=None) -> int:
    splits = (argv or sys.argv[1:] or ["dev1,dev2,dev3"])[0].split(",")
    b = views.dev_bundle()
    on = set(b.orf_to_col)
    OUT.mkdir(parents=True, exist_ok=True)
    for split in splits:
        sp = q2dev.dev_split(split)
        qd = dict(zip(sp.pool, views.query_sl_counts(b, sp.pool)))
        for r in range(3):
            sel = CappedOracleHub(split)
            t = H.Trajectory(sel, split, q2dev.traj_seed(r), r, arm="oracle_hub_capped", checkpoints=(160,),
                             max_n=160, predictors=q2dev.PREDICTORS, save=False, ledger=False, write_log=False)
            log = t.run()
            S = log["selected"][:160]
            m = log["checkpoints"]["160"][q2dev.PSTAR]["metrics"]["E_fixed_dev"]
            rec = {"split": split, "traj": r, "n": 160, "auprc": m["auprc"], "auprc_AO": m.get("auprc_AO"),
                   "auprc_VV": m.get("auprc_VV"), "mean_qdeg": float(np.mean([qd[g] for g in S])),
                   "offarray": float(np.mean([g not in on for g in S])), "model_seed": log["model_seed"]}
            with open(OUT / "oracle_capped.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
