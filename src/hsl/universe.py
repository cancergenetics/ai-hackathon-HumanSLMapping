"""Gene sub-universes for the scaling-law experiment (DESIGN v3.1 §7.3; DEV DATA ONLY).

Universe units = CURRENT SGD names of the non-test genes ((S1 query ∪ array ORFs) \\ T, minus any name whose current
name equals a test gene's current name), so genes that are one current gene enter or leave together (identity rules
i–iii stay empty by construction). A draw keeps round(level · #units) units (SeedSequence(universe_subsample_seed,
spawn_key=(level_pct, draw))). The same genes leave the query side and the array side.

sub_bundle(U)    dev bundle restricted to rows D ∩ U and columns A ∩ U (T columns dropped: T is outside every
                 sub-universe); all_query_orfs = Q ∩ U, so array-only / on-array strata are those of the sub-world.
sub_split(...)   inner 5 folds REBUILT inside U (zero-SL stratum + pandas.qcut(nonzero SL counts over U's columns, 10,
                 duplicates='drop'); per stratum ORFs sorted, shuffled with one Generator consumed stratum by stratum,
                 round-robin over folds — the hsl.views.inner_folds rule), val = fold `val_fold`, candidates = the rest.
                 E_fixed_dev(U) = views.e_fixed(sub_bundle, sub_split) (frozen constructor, unchanged rules).
run_sub(...)     upstream method on the sub-universe view (ScreenView over the sub-bundle: every GI-derived statistic —
                 array rates, LOO, query propensity, kNN, tied low rank, 2-hop, kernel normalisation over training
                 genes, PSD projection, LR standardisation — is computed from the view and so recomputed inside U),
                 scored on E_fixed_dev(U); predictions saved (float32) before scoring; one ledger line.

What is NOT restricted (stated in the results): the upstream pair-feature store (coexpression Fisher z, AP-MS
neighbourhood cosine / edge score, ESM-2 standardised cosine, paralog identity, InterPro Jaccard) is built by the
upstream builders over each SOURCE's own gene set (all genes of each expression study; all AP-MS network genes,
"so degrees are not truncated"; all ESM-embedded proteins), independent of the S1 universe — rebuilding it on U with
the upstream code gives identical values for pairs inside U. These source-level normalisations cannot be restricted
without re-deriving the sources, so the sub-universes keep full-genome external context (results labelled a bound).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import dataclasses
import os
import time

import numpy as np
import pandas as pd

from hsl import config, eval as E, harness as H, views

LEVELS = (50, 25)
N_DRAWS = 3
N_FOLDS = 5


def units() -> pd.DataFrame:
    """One row per original ORF of the non-test universe: orf, current, is_query, is_array, in_dev_pool."""
    b = views.dev_bundle()
    T = views.test_queries()
    t_cur = set(views.current_of(b, T))
    orfs = np.setdiff1d(np.union1d(b.all_query_orfs, b.array_orfs), T)
    cur = views.current_of(b, orfs)
    keep = ~np.isin(cur, list(t_cur))
    orfs, cur = orfs[keep], cur[keep]
    return pd.DataFrame({"orf": orfs, "current": cur, "is_query": np.isin(orfs, b.all_query_orfs),
                         "is_array": np.isin(orfs, b.array_orfs), "in_dev_pool": np.isin(orfs, views.dev_pool())})


def draw_universe(level_pct: int, draw: int) -> np.ndarray:
    u = units()
    names = np.unique(u.current.to_numpy())
    ss = np.random.SeedSequence(config.seed("universe_subsample_seed"), spawn_key=(int(level_pct), int(draw)))
    rng = np.random.default_rng(ss)
    k = int(round(level_pct / 100 * len(names)))
    keep = np.sort(names[rng.choice(len(names), size=k, replace=False)])
    return np.sort(u.orf[u.current.isin(keep)].to_numpy())


def sub_bundle(U):
    b = views.dev_bundle()
    U = np.asarray(U, dtype=str)
    rows = np.isin(b.strain_orf, U)
    cols = np.isin(b.array_orfs, U)
    mats = {k: np.ascontiguousarray(getattr(b, k)[rows][:, cols]) for k in ("E", "SP", "Y", "DMF", "DMF_sd")}
    strain = {k: getattr(b, k)[rows] for k in ("strain_ids", "strain_orf", "strain_canonical", "strain_is_supp",
                                                "strain_temperature")}
    arr = {k: getattr(b, k)[cols] for k in ("array_strains", "array_orfs", "array_orfs_current", "array_smf")}
    qk = np.isin(b.all_query_orfs, U)
    sb = dataclasses.replace(b, **mats, **strain, **arr, all_query_orfs=b.all_query_orfs[qk],
                             all_query_orfs_current=b.all_query_orfs_current[qk],
                             meta=dict(b.meta, hsl_view="dev_sub_universe", n_universe=int(len(U))))
    for k in ("observed", "query_orfs", "canon_row", "orf_to_col", "current_name"):
        sb.__dict__.pop(k, None)                        # cached_property values of the parent must not leak
    return sb.freeze()


def sub_split(sb, level_pct: int, draw: int, val_fold: int = 1) -> views.Split:
    D = np.intersect1d(views.dev_pool(), sb.query_orfs)
    cnt = pd.Series((sb.Y[sb.gene_rows(D)] == 1).sum(1), index=D)
    st = pd.Series(0, index=D)
    nz = cnt[cnt > 0]
    st.loc[nz.index] = pd.qcut(nz, q=10, labels=False, duplicates="drop").astype(int) + 1
    ss = np.random.SeedSequence(config.seed("universe_subsample_seed"), spawn_key=(int(level_pct), int(draw), 1))
    rng = np.random.default_rng(ss)
    fold, off = {}, 0
    for h in sorted(st.unique()):
        orfs = np.sort(st.index[st == h].to_numpy())
        for i, o in enumerate(orfs[rng.permutation(len(orfs))]):
            fold[o] = (off + i) % N_FOLDS + 1
        off += len(orfs)
    cur = views.current_of(sb, D)
    for c in np.unique(cur):
        m = np.sort(D[cur == c])
        for o in m[1:]:
            fold[o] = fold[m[0]]
    val = np.sort(np.array([o for o in D if fold[o] == val_fold]))
    cand = np.setdiff1d(D, val)
    excl = np.union1d(val, views.test_queries())
    return views.Split(f"u{level_pct}d{draw}", "dev", val, cand, excl, cand)


def scaled_grid(pool_n: int, full_pool: int = 2101) -> list[int]:
    """Budgets at the same pool fractions as the full-universe subset grid (30 … 1640 of ~2,101) plus the sub-pool
    maximum; budgets < 10 dropped."""
    base = [30, 60, 100, 160, 220, 330, 490, 660, 980, 1310, 1640]
    g = sorted({int(round(n * pool_n / full_pool)) for n in base} | {int(pool_n)})
    return [n for n in g if 10 <= n <= pool_n]


def run_sub(sb, split: views.Split, es, S, method: str, model_seed: int, traj_index: int, parent: str,
            tag: str = "p1_universe", n_threads: int | None = None, save: bool = True) -> dict:
    from hsl import upstream as U
    pr = config.protocol()
    n_threads = int(n_threads or os.environ.get("HSL_THREADS") or pr["compute"]["threads_default"])
    S = np.sort(np.asarray(S, dtype=str))
    n_folds = int(min(int(pr["predictor"]["n_folds"]), max(len(S), 2)))
    prm = H._method_params(method, {}, len(S))
    seeds = {"model_seed": int(model_seed), "traj_index": int(traj_index), "split_seed": 0}
    cfg = {"split": split.name, "split_digest": split.digest(), "method": method, "params": prm, "seeds": seeds,
           "screened_sha": H._hash(S.tolist()), "n_folds": n_folds, "tag": tag,
           "universe_n": sb.meta.get("n_universe")}
    run_id = f"{split.name}_{method}_n{len(S)}_{H._hash(cfg)}"
    rec = {"id": run_id, "parent": parent, "config_hash": H._hash(cfg), "split": split.name, "method": method,
           "n_screened": int(len(S)), "outcome": "error", "metrics": {}}
    man = {"run_id": run_id, "tag": tag, "params": prm, "seeds": seeds, "n_folds_passed": n_folds,
           "n_threads": n_threads, "git": H._git(), "split_digest": split.digest(),
           "universe_n": sb.meta.get("n_universe"), "n_arrays": int(sb.n_arrays)}
    stages, t0 = {}, time.perf_counter()
    try:
        view = U.construct_view(sb, S, split.targets)
        Smat, _ = U.run_method(view, method, model_seed=int(model_seed), split_seed=0, n_threads=n_threads,
                               n_folds=n_folds, chunk_targets=int(pr["predictor"]["chunk_targets"]), params=prm)
        stages["fit_predict_s"] = time.perf_counter() - t0
        if not np.isfinite(Smat).all():
            raise FloatingPointError("non-finite scores")
        s = es.scores(Smat).astype(np.float32)
        if save:
            man["predictions"] = H.save_predictions(run_id, split, {"E_fixed": s}, {"E_fixed": es})
        m = E.summarise(es.y, s, es.rows)
        for st, mask in es.strata.items():
            m[f"auprc_{st}"] = E.average_precision(es.y[mask], s[mask]) if es.y[mask].any() else None
        rec.update(outcome="ok", metrics={es.name: m})
    finally:
        stages["total_s"] = time.perf_counter() - t0
        man["runtime_s"] = stages
        rec["manifest"] = man
        H.append_ledger(rec)
    return {"run_id": run_id, "metrics": m, "fit_predict_s": stages["fit_predict_s"],
            "pred_file": man.get("predictions", {}).get("E_fixed", {}).get("file")}
