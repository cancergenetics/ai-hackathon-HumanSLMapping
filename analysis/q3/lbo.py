"""Leave-batch-out (LBO) cross-validation of P* on SCREENED rows only (DESIGN §7.4: online stopping estimator and
score thresholds). Labels used: rows of genes in the screened prefix S_n only; never targets, never T rows.

For a trajectory (ordered selection) and checkpoint n:
  S_n = first n selected genes; batches = consecutive blocks of `batch` genes in selection order (batch b -> fold
  b mod K, K = min(k_folds, #batches)). For each fold f: P* is refit on S_n minus fold f (same model seed as the
  trajectory) and predicts the fold's rows. CV cells = fold rows x columns that are (i) not screened in S_n (gene
  identity by current SGD name, so the reciprocal feature is never available, as in E_fixed), (ii) not a test-gene
  column (T is never used), (iii) not the row gene itself, (iv) measured. CV AUPRC = pooled tie-grouped AP over all
  folds' cells (frozen hsl.eval).

Jobs: `python -m analysis.q3.lbo run --source dev --arm R0 --split dev1 --traj 0 --checkpoints 30,60,...`
One fit at a time, BLAS threads 2. Rows appended to <out>/lbo/lbo_rows.jsonl; cell scores (float32) and labels
(int8) per (trajectory, n) saved to <out>/predictions/lbo/ (gitignored) for the threshold analysis.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import resource
import sys
import time
from pathlib import Path

from analysis.q3 import common as C

C.set_threads(2)

import numpy as np  # noqa: E402


def fold_of_batches(n: int, batch: int, k_folds: int) -> np.ndarray:
    """Fold id of each of the first n genes: batch index (i // batch) mod K, K = min(k_folds, #batches)."""
    nb = int(np.ceil(n / batch))
    K = max(2, min(int(k_folds), nb))
    return (np.arange(n) // batch) % K


def cv_cells(b, fold_genes, screened, excluded_cols):
    """(row index into sorted fold_genes, column index, label) of CV cells; see module docstring."""
    from hsl import views
    tg = np.sort(np.asarray(fold_genes, dtype=str))
    col_ok = ~views.cols_matching(b, screened) & ~views.cols_matching(b, excluded_cols)
    Y = b.Y[b.gene_rows(tg)]
    row_cur = views.current_of(b, tg)
    not_self = row_cur[:, None] != b.array_orfs_current[None, :]
    M = (Y >= 0) & col_ok[None, :] & not_self
    r, c = np.nonzero(M)
    return tg, r.astype(np.int32), c.astype(np.int32), Y[r, c].astype(np.int8)


def lbo_at(b, selected, n: int, model_seed: int, batch: int = 10, k_folds: int = 5, n_threads: int = 2,
           method: str | None = None, params: dict | None = None, excluded_cols=None, fit_fn=None) -> dict:
    """LBO-CV of the first n genes of `selected`. fit_fn(train, targets, n_folds) -> score matrix replaces P*
    (toy tests)."""
    from hsl import eval as E, upstream as U, views
    S = np.asarray(selected[:n], dtype=str)
    if len(set(S)) != len(S):
        raise ValueError("duplicate genes in selection")
    excluded_cols = views.test_queries() if excluded_cols is None else np.asarray(excluded_cols, dtype=str)
    if np.isin(S, excluded_cols).any():
        raise PermissionError("a screened gene is an excluded (test) gene")
    folds = fold_of_batches(n, batch, k_folds)
    ys, ss, rows_g, cols, fold_id, fits = [], [], [], [], [], []
    pr = C.protocol()["predictor"]
    for f in np.unique(folds):
        test = S[folds == f]
        train = S[folds != f]
        nf = int(min(int(pr["n_folds"]), max(len(train), 2)))
        t0 = time.perf_counter()
        if fit_fn is None:
            view = U.construct_view(b, train, test)
            Smat, _ = U.run_method(view, method, model_seed=int(model_seed), split_seed=0, n_threads=n_threads,
                                   n_folds=nf, chunk_targets=int(pr["chunk_targets"]), params=params)
        else:
            Smat = np.asarray(fit_fn(train, test, nf))
        dt = time.perf_counter() - t0
        if fit_fn is None:
            del view
        tg, r, c, y = cv_cells(b, test, S, excluded_cols)
        ys.append(y)
        ss.append(np.asarray(Smat, dtype=np.float32)[r, c])
        rows_g.append(tg[r])
        cols.append(c)
        fold_id.append(np.full(len(y), f, dtype=np.int16))
        del Smat
        gc.collect()
        fits.append({"fold": int(f), "n_train": int(len(train)), "n_test": int(len(test)), "fit_s": dt,
                     "n_cells": int(len(y)), "n_pos": int(y.sum())})
    y = np.concatenate(ys)
    s = np.concatenate(ss)
    out = {"n": int(n), "k": int(len(np.unique(folds))), "n_cells": int(len(y)), "n_pos": int(y.sum()),
           "prevalence": float(y.mean()) if len(y) else float("nan"),
           "cv_auprc": E.average_precision(y, s) if y.any() else float("nan"),
           "cv_auprc_fold_mean": float(np.nanmean([E.average_precision(a, z) if a.any() else np.nan
                                                  for a, z in zip(ys, ss)])),
           "fits": fits, "fit_s_total": float(sum(x["fit_s"] for x in fits))}
    out["cv_lift"] = out["cv_auprc"] / out["prevalence"] if out["prevalence"] else float("nan")
    out["_cells"] = {"y": y, "s": s, "row_gene": np.concatenate(rows_g), "col": np.concatenate(cols),
                     "fold": np.concatenate(fold_id)}
    return out


def _maxrss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


def rss_watchdog(cap_gb: float, poll_s: float = 0.5) -> None:
    """Kill this process (exit 3) if its current RSS exceeds cap_gb (PI limit: peak RSS <= 3 GB)."""
    import subprocess
    import threading

    def loop():
        pid = str(os.getpid())
        while True:
            try:
                kb = int(subprocess.run(["ps", "-o", "rss=", "-p", pid], capture_output=True, text=True).stdout)
                if kb / 2**20 > cap_gb:
                    sys.stderr.write(f"RSS watchdog: {kb / 2**20:.2f} GB > {cap_gb} GB, exiting\n")
                    os._exit(3)
            except ValueError:
                pass
            time.sleep(poll_s)

    threading.Thread(target=loop, daemon=True).start()


def run(args) -> None:
    rss_watchdog(float(args.rss_cap_gb) - 0.05)
    import yaml
    sys.path[:0] = [str(C.ROOT / "src"), str(C.ROOT / "third_party")]
    import hsl  # noqa: F401
    import hsl.predictors  # noqa: F401  (registers P*)
    pz = yaml.safe_load((C.ROOT / "config" / "predictor.yaml").read_text())["primary"]
    method, params = pz["method"], dict(pz["params"] or {})
    out = C.out_dir(args.source, args.out)
    pred_dir = out / "predictions" / "lbo"
    pred_dir.mkdir(parents=True, exist_ok=True)
    part = out / "lbo" / "lbo_rows.jsonl"
    part.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if part.exists():
        for ln in part.read_text().splitlines():
            if ln.strip():
                x = json.loads(ln)
                done.add((x["arm"], x["split"], x["traj"], x["n"]))
    sels = C.selections(args.source, arms=[args.arm])
    key = (args.arm, args.split, int(args.traj))
    if key not in sels:
        raise KeyError(f"no stored selection for {key}")
    sel = sels[key]
    b = C.dev_bundle()
    rss_cap = float(args.rss_cap_gb)
    for n in [int(x) for x in args.checkpoints.split(",")]:
        if (args.arm, args.split, int(args.traj), n) in done:
            continue
        if n > len(sel["selected"]):
            break
        t0 = time.perf_counter()
        r = lbo_at(b, sel["selected"], n, sel["model_seed"], batch=int(args.batch), k_folds=int(args.k),
                   n_threads=2, method=method, params=params)
        cells = r.pop("_cells")
        tag = f"{args.arm}_{args.split}_r{args.traj}_n{n}"
        np.savez_compressed(pred_dir / f"{tag}.npz", **cells)
        row = {"arm": args.arm, "split": args.split, "traj": int(args.traj), "n": n, "model_seed": sel["model_seed"],
               "method": method, "batch": int(args.batch), "k_folds": int(args.k), **r,
               "wall_s": time.perf_counter() - t0, "maxrss_gb": _maxrss_gb(), "cells_file":
               os.path.relpath(pred_dir / f"{tag}.npz", C.ROOT), "t_done": time.strftime("%Y-%m-%dT%H:%M:%S")}
        with open(part, "a") as fh:
            fh.write(json.dumps(row, default=float) + "\n")
        print(f"{tag}: cv_auprc={r['cv_auprc']:.4f} prev={r['prevalence']:.4f} fit_s={r['fit_s_total']:.0f} "
              f"rss={row['maxrss_gb']:.2f} GB", flush=True)
        if row["maxrss_gb"] > rss_cap:
            raise MemoryError(f"peak RSS {row['maxrss_gb']:.2f} GB > cap {rss_cap} GB; stopping")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--source", choices=C.SOURCES, default="dev")
    r.add_argument("--out", default=None)
    r.add_argument("--arm", required=True)
    r.add_argument("--split", required=True)
    r.add_argument("--traj", type=int, required=True)
    r.add_argument("--checkpoints", default="30,60,100,160,220,330,490,660")
    r.add_argument("--batch", type=int, default=10)
    r.add_argument("--k", type=int, default=5)
    r.add_argument("--rss-cap-gb", default=3.0)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        run(a)


if __name__ == "__main__":
    main()
