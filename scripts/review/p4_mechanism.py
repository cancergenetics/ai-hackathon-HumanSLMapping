#!/usr/bin/env python
"""P4 review (independent, DEV ONLY, evaluation-only diagnostic): what do the screened sets of adaptive arms reveal
in the predictor-visible but selector-masked columns, compared with the hub oracle and random?

For every stored trajectory (P4 parts, P3 oracle_hub parts, P1 random permutations 0-9) at n in NS:
  sl_rows        strong SLs in the screened rows (all columns)
  sl_efix_cols   ... restricted to E_fixed_dev columns (array-only + on-array D_val columns): array-side propensity
                 information for the evaluation columns
  sl_vv_cols     ... restricted to on-array D_val columns (target genes' own array side = query-propensity proxy)
  sl_ao_cols     ... restricted to array-only columns
  vv_targets_hit on-array D_val targets with >= 1 SL revealed in their column
  rho_vv         Spearman over on-array D_val targets: revealed SLs in the target's column vs the target's true
                 E_fixed_dev row SL count (how informative the revealed target-column counts are)
  mean_qdeg      mean true query SL count of the screened genes; offarray share
plus the stored AUPRC (all / AO / VV / on- and off-array target rows). Output: results/review_p4/mechanism.csv.
"""
from __future__ import annotations

import hsl  # noqa: F401

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from hsl import harness as H, p1, q2dev, views

OUT = Path(__file__).resolve().parents[2] / "results" / "review_p4"
NS = (30, 100, 160, 330, 660)
ARMS = ("H", "C", "F_r16", "M_lam0.03", "H-static", "R1", "F-static_r8")


def main() -> int:
    b = views.dev_bundle()
    Y = np.asarray(b.Y)
    recs = []
    for split in ("dev1", "dev2", "dev3"):
        sp = q2dev.dev_split(split)
        es = q2dev._eset(split)
        efix_cols = np.unique(es.cols)
        vv_cols = np.unique(es.cols[np.asarray(es.strata["VV"], bool)])
        ao_cols = np.unique(es.cols[np.asarray(es.strata["AO"], bool)])
        on_targets = [t for t in sp.targets if t in b.orf_to_col]
        tcol = np.array([b.orf_to_col[t] for t in on_targets])
        # true E_fixed_dev row SL count of each on-array target
        ti = {t: i for i, t in enumerate(es.targets)}
        row_sl = np.bincount(es.rows, weights=es.y.astype(float), minlength=len(es.targets))
        t_rowsl = np.array([row_sl[ti[t]] for t in on_targets])
        qdeg = dict(zip(sp.pool, views.query_sl_counts(b, sp.pool)))
        on = set(b.orf_to_col)
        trajs = []                                       # (arm, traj, selected list, rows by n)
        for arm in ARMS:
            for f in sorted((q2dev.RESULTS / "p4" / "parts").glob(f"{arm}_{split}_r*.jsonl")):
                rows = p1.read_rows([f])
                sel = json.loads((H.TRAJ_DIR / split / f"{rows[0]['traj_id']}.json").read_text())["selected"]
                trajs.append((arm, rows[0]["traj"], sel, {x["n"]: x for x in rows}))
        for f in sorted((q2dev.RESULTS / "p3" / "parts").glob(f"oracle_hub_{split}_r*.jsonl")):
            rows = p1.read_rows([f])
            sel = json.loads((H.TRAJ_DIR / split / f"{rows[0]['traj_id']}.json").read_text())["selected"]
            trajs.append(("oracle_hub", rows[0]["traj"], sel, {x["n"]: x for x in rows}))
        rnd = [x for x in q2dev.random_rows() if x["split"] == split]
        for r in range(10):
            t = p1.random_traj(split, r)
            trajs.append(("R0", r, list(t.perm), {x["n"]: x for x in rnd if x["traj"] == r}))
        for arm, r, sel, byn in trajs:
            for n in NS:
                if n not in byn or len(sel) < n:
                    continue
                S = sel[:n]
                YS = Y[b.gene_rows(S)]
                sl = YS == 1
                rec = {"split": split, "arm": arm, "traj": int(r), "n": n,
                       "sl_rows": int(sl.sum()), "sl_efix_cols": int(sl[:, efix_cols].sum()),
                       "sl_vv_cols": int(sl[:, vv_cols].sum()), "sl_ao_cols": int(sl[:, ao_cols].sum()),
                       "vv_targets_hit": int((sl[:, tcol].sum(0) > 0).sum()),
                       "rho_vv": float(spearmanr(sl[:, tcol].sum(0), t_rowsl)[0]),
                       "mean_qdeg": float(np.mean([qdeg[g] for g in S])),
                       "offarray": float(np.mean([g not in on for g in S]))}
                x = byn[n]
                for k in ("auprc", "auprc_AO", "auprc_VV", "auprc_on_array_rows", "auprc_off_array_rows"):
                    rec[k] = x.get(k)
                recs.append(rec)
        print(split, len(recs), flush=True)
    df = pd.DataFrame(recs)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "mechanism.csv", index=False)
    g = df.groupby(["n", "split", "arm"]).mean(numeric_only=True).drop(columns="traj")
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_rows", 400):
        print(g.loc[160].round(4))
        print(g.loc[30].round(4))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
