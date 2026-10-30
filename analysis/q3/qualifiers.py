"""'Where the gain is' qualifiers (final-results review): pooled AUPRC gain vs per-query macro-AUPRC gain, AUROC
difference, and (final) the gain by SL-count tertile and on E_TT, for the finalists vs random at n = 160 / 660.

    .venv/bin/python -m analysis.q3.qualifiers --source dev|final   -> <out>/gain_qualifiers.csv
dev: macro AUPRC / AUROC stored with the P1/P4 rows. final: recomputed from the saved predictions with the frozen
hsl.eval.summarise on the scorer's T labels; tertile / TT relative Δ read from results/final/contrasts.csv.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from analysis.q3 import common as C

BUDGETS = (160, 660)


def secondary(source: str) -> pd.DataFrame:
    rows = C.rows(source)
    rows = rows[rows.arm.isin(("R0",) + C.FINALISTS) & rows.n.isin(BUDGETS)]
    if source == "dev":
        p4 = pd.concat([C._p1_random_rows(), C._p4_rows()], ignore_index=True)
        p4 = p4[p4.arm.isin(("R0",) + C.FINALISTS) & p4.n.isin(BUDGETS)]
        return p4[["arm", "split", "traj", "n", "auprc", "macro_auprc", "auroc"]]
    from hsl import eval as E
    L = C.eval_labels(source, "final" if source == "final" else "dev1")
    out = []
    for _, x in rows.iterrows():
        m = E.summarise(L["y"], C.load_pred(x.pred_file), L["rows"])
        out.append({"arm": x.arm, "split": x.split, "traj": x.traj, "n": x.n, "auprc": x.auprc,
                    "macro_auprc": m["macro_auprc"], "auroc": m["auroc"]})
    return pd.DataFrame(out)


def analyse(source: str, out=None) -> pd.DataFrame:
    O = C.out_dir(source, out)
    S = secondary(source)
    M = S.groupby(["arm", "split", "n"]).mean(numeric_only=True).reset_index()
    rec = []
    for a in C.FINALISTS:
        for n in BUDGETS:
            per = []
            for sp in sorted(M.split.unique()):
                s = M[(M.arm == a) & (M.split == sp) & (M.n == n)]
                r = M[(M.arm == "R0") & (M.split == sp) & (M.n == n)]
                if len(s) and len(r):
                    s, r = s.iloc[0], r.iloc[0]
                    per.append({"pooled": s.auprc / r.auprc - 1, "macro": s.macro_auprc / r.macro_auprc - 1,
                                "auroc_diff": s.auroc - r.auroc})
            if not per:
                continue
            P = pd.DataFrame(per)
            rec.append({"arm": a, "n": n, "pooled_auprc_rel": P.pooled.mean(), "macro_auprc_rel": P.macro.mean(),
                        "auroc_diff": P.auroc_diff.mean(), "n_splits": len(P)})
    Q = pd.DataFrame(rec)
    if source != "dev":
        c = pd.read_csv(C.final_dir(source) / "contrasts.csv")
        c = c[(c.b == "R0") & (c.pred_a == "pstar") & (c.n_a == c.n_b) & c.a.isin(C.FINALISTS)]
        for st in ("TT", "onrows", "tert1", "tert2", "tert3"):
            x = c[c.set == st][["a", "n_a", "rel", "lo", "hi", "s_b"]]
            x = x.assign(**{f"{st}_rel": x.rel, f"{st}_ci": [f"[{lo:+.4f}, {hi:+.4f}]" for lo, hi in zip(x.lo, x.hi)]})
            Q = Q.merge(x[["a", "n_a", f"{st}_rel", f"{st}_ci"]].rename(columns={"a": "arm", "n_a": "n"}),
                        on=["arm", "n"], how="left")
    Q.to_csv(O / "gain_qualifiers.csv", index=False)
    return Q


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    print(analyse(a.source, a.out).round(3).to_string())


if __name__ == "__main__":
    main()
