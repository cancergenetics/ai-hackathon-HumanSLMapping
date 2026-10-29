#!/usr/bin/env python
"""P4 review (independent, DEV ONLY): is the adaptive-selection gain at n = 160 specific to P*?

Refits a non-P* predictor on the STORED selected sets (first n genes) of finalist trajectories, with the trajectory's
own model seed (CRN with the random arm), and compares with the P1 random rows of the same predictor (same seeds,
same E_fixed_dev cells): Welch per split with the shared 'q2_dev' query weights (hsl.eval.welch).
Predictors: stack_lr_wr2 (pre-P2 P*, DESIGN §4 'unmodified predictor'), B1_array_rate (array-side SL rate only: a
pure array-propensity channel), m1_hgb (genome dropped; dev1 only, where P1 has 5 random rows).
Nothing is saved to the prediction store and no LEDGER line is written. Output: results/review_p4/altpred.jsonl.

    python scripts/review/p4_altpred.py --predictor stack_lr_wr2 --arms H,C,F_r16 --splits dev1,dev2,dev3 --trajs 0-2
"""
from __future__ import annotations

import hsl  # noqa: F401

import argparse
import json
from pathlib import Path

import numpy as np

import hsl.predictors  # noqa: F401
from hsl import eval as E, harness as H, p1, paths, q2dev

OUT = Path(__file__).resolve().parents[2] / "results" / "review_p4"


def rand_rows(pred, split, n):
    out = {}
    for x in p1.read_rows(sorted((p1.P1_DIR / "parts").glob("*.jsonl"))):
        if x.get("predictor") == pred and x["split"] == split and x["n"] == n and x.get("kind", "curve") == "curve" \
                and x.get("pred_file") and (paths.PROJECT_ROOT / x["pred_file"]).exists():
            out[int(x["perm"])] = x
    return [out[k] for k in sorted(out)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictor", required=True)
    ap.add_argument("--arms", default="H,C,F_r16")
    ap.add_argument("--splits", default="dev1,dev2,dev3")
    ap.add_argument("--trajs", default="0-2")
    ap.add_argument("--n", type=int, default=160)
    a = ap.parse_args(argv)
    lo, hi = (int(t) for t in a.trajs.split("-"))
    prm = p1.PREDICTORS[a.predictor]
    OUT.mkdir(parents=True, exist_ok=True)
    for split in a.splits.split(","):
        es = q2dev._eset(split)
        W = q2dev.weights(split)
        rr = rand_rows(a.predictor, split, a.n)
        Ar = np.column_stack([q2dev.boot_vector(x["pred_file"], split, OUT / "boot")[1:] for x in rr])
        er = np.array([q2dev.boot_vector(x["pred_file"], split, OUT / "boot")[0] for x in rr])
        for arm in a.arms.split(","):
            As, ess = [], []
            for r in range(lo, hi + 1):
                row = next(x for x in p1.read_rows([q2dev.RESULTS / "p4" / "parts" / f"{arm}_{split}_r{r}.jsonl"]))
                sel = json.loads((H.TRAJ_DIR / split / f"{row['traj_id']}.json").read_text())["selected"][: a.n]
                res = H.run_predictor(H.RunSpec(split, tuple(sel), ("E_fixed",), ()), a.predictor, prm,
                                      {"model_seed": int(row["model_seed"]), "traj_index": r}, save=False,
                                      ledger=False)
                s = res.scores["E_fixed"]
                e = E.WeightedAP(es.y, s, es.rows, W.shape[1])
                As.append(np.concatenate([e.ap(W[i:i + 250]) for i in range(0, len(W), 250)]))
                ess.append(float(e.ap(np.ones(W.shape[1]))[0]))
                print(split, arm, r, a.predictor, round(ess[-1], 5), round(res.manifest["runtime_s"]["fit_predict_s"],
                                                                          1), flush=True)
            w = E.welch(np.column_stack(As), np.array(ess), Ar, er)
            rec = {"predictor": a.predictor, "split": split, "arm": arm, "n": a.n, "R_s": len(ess), "R_r": len(er),
                   "s_mean": float(np.mean(ess)), "r_mean": float(er.mean()), "delta": w["delta"],
                   "rel": w["delta"] / float(er.mean()), "lo": w["ci95"][0], "hi": w["ci95"][1],
                   "p": w["p_one_sided"], "s_each": ess}
            with open(OUT / "altpred.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
