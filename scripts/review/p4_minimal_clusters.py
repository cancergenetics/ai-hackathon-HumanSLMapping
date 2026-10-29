#!/usr/bin/env python
"""P4 review (independent, DEV ONLY): the HAP1-minimal scenario (<X>__min) still uses AP-MS through the batch-cap
clusters (hsl.features.clusters: k-means on coexpression + ESM + AP-MS kernel). Re-run <X>__min with clusters fitted
WITHOUT the AP-MS kernel (same k, embedding dim, seed, universe) and compare selections and P* AUPRC at n = 160 with
the stored <X>__min trajectories (same seeds). Nothing saved to the prediction store, no LEDGER line.
Output: results/review_p4/minimal_clusters.jsonl.

    python scripts/review/p4_minimal_clusters.py H dev1,dev2,dev3
"""
from __future__ import annotations

import hsl  # noqa: F401

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from hsl import config, features, harness as H, p1, q2dev, strategies

OUT = Path(__file__).resolve().parents[2] / "results" / "review_p4"


def clusters_without_apms() -> pd.Series:
    from sklearn.cluster import KMeans
    G = features._G()
    pr = config.protocol()["selector"]
    k, dim, seed = int(pr["kmeans_k"]), int(pr["cluster_embedding_dim"]), config.seed("cluster_seed")
    genes = features.cluster_universe()
    cur = G.current_names(genes)
    Kc = np.clip(np.nan_to_num(G.similarity("coexp", cur, cur, current=True)), -1, 1)
    e = G._esm_current()
    X = e.reindex(cur).to_numpy(dtype=np.float64)
    mu, sd = np.nanmean(e.to_numpy(np.float64), 0), np.nanstd(e.to_numpy(np.float64), 0)
    X = (X - mu) / np.where(sd > 0, sd, 1)
    X = X / np.linalg.norm(np.nan_to_num(X), axis=1, keepdims=True).clip(1e-12)
    K = (Kc + np.nan_to_num(X @ X.T)) / 2.0
    K = (K + K.T) / 2
    np.fill_diagonal(K, 1.0)
    n = len(K)
    J = np.eye(n) - 1.0 / n
    w, V = np.linalg.eigh(J @ K @ J)
    top = np.argsort(w)[::-1][:dim]
    Z = V[:, top] * np.sqrt(np.clip(w[top], 0, None))
    lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(Z).astype(np.int32)
    return pd.Series(lab, index=pd.Index(genes, name="orf"), name="cluster")


def main(argv=None) -> int:
    argv = argv or sys.argv[1:]
    arm, splits = argv[0], argv[1].split(",")
    old = features._cluster_table()
    new = clusters_without_apms()
    from sklearn.metrics import adjusted_rand_score
    ari = float(adjusted_rand_score(old.reindex(new.index).to_numpy(), new.to_numpy()))
    print("ARI(old, no-AP-MS clusters) =", round(ari, 3), flush=True)
    features._cluster_table = lambda: new                      # every view.clusters() call now uses the new table
    cid = f"{arm}__min"
    for split in splits:
        for r in range(3):
            stored_row = next(x for x in p1.read_rows([q2dev.RESULTS / "p4" / "parts" / f"{cid}_{split}_r{r}.jsonl"])
                              if x["n"] == 160)
            stored = json.loads((H.TRAJ_DIR / split / f"{stored_row['traj_id']}.json").read_text())["selected"][:160]
            spec = q2dev.parse(cid)
            sel = q2dev.make_selector(spec, split)
            try:
                t = H.Trajectory(sel, split, q2dev.traj_seed(r), r, arm=cid + "_noapmscl", checkpoints=(160,),
                                 max_n=160, predictors=q2dev.PREDICTORS, save=False, ledger=False, write_log=False)
                log = strategies.run_trajectory(t)
            finally:
                q2dev._unpatch_minimal()
            S = log["selected"][:160]
            m = log["checkpoints"]["160"][q2dev.PSTAR]["metrics"]["E_fixed_dev"]
            rec = {"arm": cid, "split": split, "traj": r, "ari": ari, "auprc_noapms_clusters": m["auprc"],
                   "auprc_stored_min": stored_row["auprc"], "overlap160": len(set(S) & set(stored)),
                   "first_diff": next((i for i, (x, y) in enumerate(zip(S, stored)) if x != y), None)}
            with open(OUT / "minimal_clusters.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
