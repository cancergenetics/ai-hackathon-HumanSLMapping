#!/usr/bin/env python
"""P4 review (independent, DEV ONLY): recompute H vs random Δ(160) per split (Welch, DESIGN §6.4) and pooled gap
closure G for H from the stored prediction vectors, using only hsl.eval.WeightedAP / bootstrap_weights and hand-written
Welch / joint-bootstrap code (not q2dev.welch_at / contrast_suite / gap_closure_joint). Unit-weight AP is checked
against sklearn. Two weight sets: the shared 'q2_dev' matrix (as P4) and a fresh independent stratified matrix.
"""
from __future__ import annotations

import hsl  # noqa: F401

import json
from pathlib import Path

import numpy as np
from scipy import stats
from sklearn.metrics import average_precision_score

from hsl import eval as E, p1, q2dev, views

ARM = "H"
SPLITS = ("dev1", "dev2", "dev3")
B = 2000


def stored(arm, split, n):
    if arm == "R0":
        return sorted([x for x in q2dev.random_rows() if x["split"] == split and x["n"] == n], key=lambda x: x["traj"])
    out = []
    for f in sorted((q2dev.RESULTS / "p4" / "parts").glob(f"{arm}_{split}_r*.jsonl")):
        out += [x for x in p1.read_rows([f]) if x["n"] == n and x["config"] == arm]
    return sorted(out, key=lambda x: x["traj"])


def reps(rows, es, W, check_sk=False):
    A, est = [], []
    for x in rows:
        s = p1.load_scores(x["pred_file"])
        e = E.WeightedAP(es.y, s, es.rows, W.shape[1])
        u = float(e.ap(np.ones(W.shape[1]))[0])
        if check_sk:
            sk = average_precision_score(es.y.astype(int), s.astype(np.float64))
            assert abs(sk - u) < 1e-6, (sk, u)
        assert abs(u - float(x["auprc"])) < 1e-9, (u, x["auprc"])
        est.append(u)
        A.append(np.concatenate([e.ap(W[i:i + 250]) for i in range(0, len(W), 250)]))
    return np.column_stack(A), np.array(est)


def welch(As, es_, Ar, er):
    d = es_.mean() - er.mean()
    VQ = np.var(As.mean(1) - Ar.mean(1), ddof=1)
    a, b = np.var(es_, ddof=1) / len(es_), np.var(er, ddof=1) / len(er)
    se = np.sqrt(VQ + a + b)
    nu = se ** 4 / (a ** 2 / (len(es_) - 1) + b ** 2 / (len(er) - 1))
    c = stats.t.ppf(0.975, nu)
    return {"delta": d, "lo": d - c * se, "hi": d + c * se, "p": float(stats.t.sf(d / se, nu)), "V_Q": VQ,
            "se": se, "nu": nu, "s_sd": np.sqrt(np.var(es_, ddof=1)), "r_sd": np.sqrt(np.var(er, ddof=1))}


def main() -> int:
    out = {}
    pairs = {"q2_dev": [], "fresh": []}
    for split in SPLITS:
        es = q2dev._eset(split)
        sp = views.get_split(split)
        strata = views.query_strata(sp.targets, views.dev_bundle())
        Ws = {"q2_dev": p1.query_weights(split, "q2_dev"),
              "fresh": E.bootstrap_weights(strata, B, np.random.default_rng(424242 + int(split[3:])))}
        for wn, W in Ws.items():
            As, ess = reps(stored(ARM, split, 160), es, W, check_sk=(wn == "q2_dev"))
            Ar, er = reps(stored("R0", split, 160), es, W)
            Ar6, er6 = reps(stored("R0", split, 660), es, W)
            w = welch(As, ess, Ar, er)
            w["rel"] = w["delta"] / er.mean()
            w["G_split"] = (ess.mean() - er.mean()) / (er6.mean() - er.mean())
            out[f"{split}/{wn}"] = {k: float(v) for k, v in w.items()}
            pairs[wn].append((As, ess, Ar, er, Ar6, er6))
            print(split, wn, {k: round(float(v), 6) for k, v in w.items()}, flush=True)
    for wn, pr in pairs.items():
        num = sum(p[1].mean() - p[3].mean() for p in pr)
        den = sum(p[5].mean() - p[3].mean() for p in pr)
        rng = np.random.default_rng(777)
        nb = np.zeros(B)
        db = np.zeros(B)
        for As, ess, Ar, er, Ar6, er6 in pr:
            i_s = rng.integers(0, As.shape[1], size=(B, As.shape[1]))
            i_r = rng.integers(0, Ar.shape[1], size=(B, Ar.shape[1]))     # whole random trajectories, same draw
            ms = np.take_along_axis(As, i_s, 1).mean(1)                     # for 160 and 660
            mr = np.take_along_axis(Ar, i_r, 1).mean(1)
            mr6 = np.take_along_axis(Ar6, i_r, 1).mean(1)
            nb += ms - mr
            db += mr6 - mr
        G = nb / db
        out[f"G_pooled/{wn}"] = {"G": float(num / den), "lo": float(np.quantile(G, 0.025)),
                                 "hi": float(np.quantile(G, 0.975))}
        print("G pooled", wn, out[f"G_pooled/{wn}"], flush=True)
    d = Path(__file__).resolve().parents[2] / "results" / "review_p4"
    d.mkdir(parents=True, exist_ok=True)
    (d / "recompute_stats.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
