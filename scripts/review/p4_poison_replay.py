#!/usr/bin/env python
"""P4 review (independent): replay the SELECTION part of stored dev trajectories under bundle poisonings.

DEV ONLY. No predictor is fitted. One poisoning per process (module-level caches can never carry a clean result into a
poisoned run). Output: one JSON line per run appended to results/review_p4/poison_replay.jsonl.

    python scripts/review/p4_poison_replay.py --config H --split dev1 --traj 0 --poison dval_rows --n 160

poisons:
  none            baseline replay (must reproduce the stored trajectory and its traj_id)
  dval_rows       (a) every strain row of D_val genes garbled (E, SP, Y incl. missingness)
  dval_cols       (b) columns of the split's excluded genes (D_val, and T columns present in the dev bundle) garbled
                      in EVERY row
  unscreened_rows (c) every row of pool genes NOT in the stored selected[:n] garbled (never screened in the replay
                      if selection is unchanged)
  dval_smf        (d) array SMF of D_val (and T) columns garbled
  pc_one_col      positive control: one visible column (the candidate column with most SLs among the stored first 30
                      screened rows) garbled in every row
  pc_vis_cols     positive control: every visible on-array pool-gene column garbled in every row
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard; must be first)

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np

from hsl import harness as H, q2dev, views

OUT = Path(__file__).resolve().parents[2] / "results" / "review_p4"


def garble(b, rows=None, cols=None, seed=20261029):
    """Copy of bundle b with E/SP/Y garbled on rows x cols (None = all). Y random {0,1} (missingness removed),
    E ~ N(-0.15, 0.3), SP ~ U(0, 6)."""
    rng = np.random.default_rng(seed)
    E, SP, Y = (np.array(getattr(b, k), copy=True) for k in ("E", "SP", "Y"))
    r = np.arange(E.shape[0]) if rows is None else np.flatnonzero(rows)
    c = np.arange(E.shape[1]) if cols is None else np.flatnonzero(cols)
    ix = np.ix_(r, c)
    shp = (len(r), len(c))
    Y[ix] = rng.choice(np.array([0, 1], dtype=Y.dtype), size=shp, p=[0.7, 0.3])
    E[ix] = rng.normal(-0.15, 0.3, size=shp).astype(E.dtype)
    SP[ix] = rng.uniform(0, 6, size=shp).astype(SP.dtype)
    return b.with_matrices(E=E, SP=SP, Y=Y)


def poisoned_bundle(kind: str, split, stored: list[str], n: int):
    b = views.dev_bundle()
    if kind == "none":
        return b, {}
    if kind == "dval_rows":
        rows = np.isin(b.strain_orf, split.targets)
        return garble(b, rows=rows), {"n_rows": int(rows.sum())}
    if kind == "dval_cols":
        cols = views.cols_matching(b, split.excluded)
        return garble(b, cols=cols), {"n_cols": int(cols.sum())}
    if kind == "unscreened_rows":
        unscreened = np.setdiff1d(split.pool, stored[:n])
        rows = np.isin(b.strain_orf, unscreened)
        return garble(b, rows=rows), {"n_rows": int(rows.sum())}
    if kind == "dval_smf":
        cols = views.cols_matching(b, split.excluded)
        smf = np.array(b.array_smf, copy=True)
        smf[cols] = np.random.default_rng(7).uniform(0.2, 1.3, size=int(cols.sum())).astype(smf.dtype)
        return dataclasses.replace(b, array_smf=smf), {"n_cols": int(cols.sum())}
    if kind in ("pc_one_col", "pc_vis_cols"):
        vis_pool = np.array([b.orf_to_col[g] for g in split.pool if g in b.orf_to_col])
        if kind == "pc_vis_cols":
            cols = np.zeros(b.n_arrays, dtype=bool)
            cols[vis_pool] = True
            return garble(b, cols=cols), {"n_cols": int(cols.sum())}
        Y30 = np.asarray(b.Y[b.gene_rows(stored[:30])])
        cand30 = vis_pool[~np.isin(np.asarray(b.array_orfs)[vis_pool], stored[:30])]
        j = int(cand30[np.argmax((Y30[:, cand30] == 1).sum(0))])
        cols = np.zeros(b.n_arrays, dtype=bool)
        cols[j] = True
        return garble(b, cols=cols), {"col": str(b.array_orfs[j]), "col_sl_first30": int((Y30[:, j] == 1).sum())}
    raise KeyError(kind)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--traj", type=int, default=0)
    ap.add_argument("--poison", required=True)
    ap.add_argument("--n", type=int, default=160)
    a = ap.parse_args(argv)
    sp = q2dev.dev_split(a.split)
    part = Path(q2dev.RESULTS / "p4" / "parts" / f"{a.config}_{a.split}_r{a.traj}.jsonl")
    tid = json.loads(part.read_text().splitlines()[0])["traj_id"]
    stored = json.loads((H.TRAJ_DIR / a.split / f"{tid}.json").read_text())["selected"]
    t0 = time.perf_counter()
    b, info = poisoned_bundle(a.poison, sp, stored, a.n)
    spec = q2dev.parse(a.config)
    sel = q2dev.make_selector(spec, a.split)
    traj = H.Trajectory(sel, a.split, q2dev.traj_seed(a.traj), a.traj, arm=a.config, batch_size=spec.batch_size,
                        checkpoints=(), max_n=a.n, predictors=(), save=False, ledger=False, write_log=False, bundle=b)
    log = traj.run()
    got = log["selected"][: a.n]
    ref = stored[: a.n]
    first_diff = next((i for i, (x, y) in enumerate(zip(got, ref)) if x != y), None)
    rec = {"config": a.config, "split": a.split, "traj": a.traj, "poison": a.poison, "n": a.n, **info,
           "identical": got == ref, "first_diff": first_diff, "set_overlap": len(set(got) & set(ref)),
           "traj_id_replay": log["traj_id"], "traj_id_stored": tid, "wall_s": round(time.perf_counter() - t0, 1),
           "peak_rss_gb": round(H.peak_rss_gb(), 2)}
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "poison_replay.jsonl", "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
