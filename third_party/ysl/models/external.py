"""B5: external-only method (no GI information about the held-out query x).

score(x, a) = HistGradientBoosting log-odds on
  * all external pair features f(x, a) (src/ysl/features/pairs.py; NaN handled natively by HGB), and
  * array-side information about a only: log shrunk array positive rate over the screened rows (exact
    leave-one-query-out for the training rows, DESIGN §3 self-exclusion) and the array gene's SMF.
Nothing about x comes from the NxN matrix: no x array-side column / propensity, no reciprocal, no x SMF, no
on-array flag. Training cells = every observed cell of the view's screened rows; the pair features are external
and identical at train and test time, so no cross-fitting is needed beyond the LOO array rate.
Hyper-parameters are fixed a priori (not tuned on dev). `use_array_rate=False` gives the pure pair-feature variant.
"""
from __future__ import annotations

import numpy as np

from ysl.features import pairs as PF
from ysl.models import baselines as B
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

B5_DEFAULTS = {"rate_prior": B.DEFAULT_RATE_PRIOR, "use_array_rate": True, "max_iter": 200, "learning_rate": 0.05,
               "max_leaf_nodes": 31, "min_samples_leaf": 500, "l2_regularization": 1.0}


def _pair_stack(rows: np.ndarray, cols: np.ndarray, names) -> list[np.ndarray]:
    return [PF.block(n, rows, cols) for n in names]


def _array_cols(v: ScreenView, ctx: MethodContext, log_rate: np.ndarray) -> list[np.ndarray]:
    """Array-side columns [n_rows x n_a] (or broadcastable)."""
    out = [np.broadcast_to(v.array_smf[None, :], log_rate.shape).astype(np.float32)]
    if ctx.get("use_array_rate", True):
        out.append(log_rate.astype(np.float32))
    return out


def _flatten(cols: list[np.ndarray], keep=None) -> np.ndarray:
    X = np.stack([np.ascontiguousarray(c, dtype=np.float32).ravel() for c in cols], axis=1)
    return X if keep is None else X[keep]


@register("B5_external", defaults=B5_DEFAULTS, logit_scale=True,
          description="external-only HGB: pair features + array-side rate (LOO) and array SMF; nothing about x from GI")
def b5_external(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    from sklearn.ensemble import HistGradientBoostingClassifier
    names = PF.feature_names()
    a0 = ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)
    cidx = PF.index_of(v.array_orfs_current)
    # ---- training: screened rows, LOO array rate
    tr_idx = PF.index_of(v.train_orfs_current)
    y = v.Y.ravel()
    keep = ~np.isnan(y)
    loo_rate = np.log(B.loo_additive(v, a0)["array_rate"])
    X = _flatten(_pair_stack(tr_idx, cidx, names) + _array_cols(v, ctx, loo_rate), keep)
    yk = y[keep].astype(np.int8)
    clf = HistGradientBoostingClassifier(
        max_iter=ctx.get("max_iter"), learning_rate=ctx.get("learning_rate"), max_leaf_nodes=ctx.get("max_leaf_nodes"),
        min_samples_leaf=ctx.get("min_samples_leaf"), l2_regularization=ctx.get("l2_regularization"),
        early_stopping=False, random_state=int(ctx.model_seed))
    clf.fit(X, yk)
    ctx.info["n_training_cells"] = int(len(yk))
    ctx.info["features"] = list(names) + ["array_smf"] + (["log_array_rate"] if ctx.get("use_array_rate", True) else [])
    del X
    # ---- targets: full grid, chunked
    rate = np.log(B.array_rate(v, a0))[None, :]
    tg_idx = PF.index_of(v.target_orfs_current)
    out = np.empty((v.n_targets, v.n_arrays), dtype=np.float32)
    for i in range(0, v.n_targets, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, v.n_targets))
        n = sl.stop - sl.start
        Xc = _flatten(_pair_stack(tg_idx[sl], cidx, names) + _array_cols(v, ctx, np.broadcast_to(rate, (n, v.n_arrays))))
        out[sl] = clf.decision_function(Xc).reshape(n, v.n_arrays)
    return out


# ===================================================================================================== B5 variants
# Added after the first B5_external dev run (0.0410 < B1_array_rate 0.0430): the flexible HGB can recover
# query-level gene properties from pair features (e.g. abund_sum with abund_absdiff) and fit noisy per-query rates of
# the 657 screened queries. Two constrained forks (B5_external itself is left unchanged):
#   B5_external_es : same HGB, early stopping on a held-out 10% of TRAINING QUERIES (query-grouped X_val, DESIGN §6)
#   B5_external_lr : logistic regression; each pair feature (monotone transform) NaN -> 0 next to one availability
#                    flag per mask group; + log array rate (LOO in training) + array SMF (+ flag)
_LOG1P = {"apms_shared", "apms_edge", "dom_pfam_shared", "dom_ipr_shared", "gen_dist_kb"}


@register("B5_external_es", defaults=B5_DEFAULTS | {"max_iter": 500, "val_query_frac": 0.1, "n_iter_no_change": 10},
          logit_scale=True, description="B5_external + early stopping on held-out training queries (grouped)")
def b5_external_es(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    from sklearn.ensemble import HistGradientBoostingClassifier
    names = PF.feature_names()
    a0 = ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)
    cidx = PF.index_of(v.array_orfs_current)
    tr_idx = PF.index_of(v.train_orfs_current)
    loo_rate = np.log(B.loo_additive(v, a0)["array_rate"])
    X = _flatten(_pair_stack(tr_idx, cidx, names) + _array_cols(v, ctx, loo_rate))
    y = v.Y.ravel()
    val_q = ctx.rng("B5_es_val").random(v.n) < ctx.get("val_query_frac")
    is_val = np.repeat(val_q, v.n_arrays)
    obs = ~np.isnan(y)
    tr, va = obs & ~is_val, obs & is_val
    clf = HistGradientBoostingClassifier(
        max_iter=ctx.get("max_iter"), learning_rate=ctx.get("learning_rate"), max_leaf_nodes=ctx.get("max_leaf_nodes"),
        min_samples_leaf=ctx.get("min_samples_leaf"), l2_regularization=ctx.get("l2_regularization"),
        early_stopping=True, scoring="loss", n_iter_no_change=ctx.get("n_iter_no_change"),
        random_state=int(ctx.model_seed))
    clf.fit(X[tr], y[tr].astype(np.int8), X_val=X[va], y_val=y[va].astype(np.int8))
    del X
    ctx.info["n_iter"] = int(clf.n_iter_)
    ctx.info["n_training_cells"] = int(tr.sum())
    rate = np.log(B.array_rate(v, a0))[None, :]
    tg_idx = PF.index_of(v.target_orfs_current)
    out = np.empty((v.n_targets, v.n_arrays), dtype=np.float32)
    for i in range(0, v.n_targets, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, v.n_targets))
        n = sl.stop - sl.start
        Xc = _flatten(_pair_stack(tg_idx[sl], cidx, names) + _array_cols(v, ctx, np.broadcast_to(rate, (n, v.n_arrays))))
        out[sl] = clf.decision_function(Xc).reshape(n, v.n_arrays)
    return out


def _lr_design(rows: np.ndarray, cidx: np.ndarray, names, v: ScreenView, log_rate: np.ndarray) -> np.ndarray:
    cols, flags = [], {}
    for n in names:
        M = PF.block(n, rows, cidx)
        if n == "coex_mr":
            M = np.log(M)
        elif n in _LOG1P:
            M = np.log1p(np.maximum(M, 0))
        avail = ~np.isnan(M)
        flags.setdefault(PF.FEATURES[n].mask, avail)
        cols.append(np.where(avail, M, 0.0))
    cols += [f.astype(np.float32) for f in flags.values()]
    smf = np.broadcast_to(v.array_smf[None, :], log_rate.shape)
    cols += [np.where(np.isnan(smf), 0.0, smf), (~np.isnan(smf)).astype(np.float32), log_rate]
    return _flatten(cols)


@register("B5_external_lr", defaults={"rate_prior": B.DEFAULT_RATE_PRIOR, "C": 1.0, "max_iter": 300},
          logit_scale=True, description="external-only LR: pair features (+availability flags), log array rate, SMF")
def b5_external_lr(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    names = PF.feature_names()
    a0 = ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)
    cidx = PF.index_of(v.array_orfs_current)
    y = v.Y.ravel()
    keep = ~np.isnan(y)
    X = _lr_design(PF.index_of(v.train_orfs_current), cidx, names, v,
                   np.log(B.loo_additive(v, a0)["array_rate"]))[keep]
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(C=ctx.get("C"), max_iter=ctx.get("max_iter")).fit(sc.transform(X), y[keep].astype(np.int8))
    del X
    ctx.info["n_training_cells"] = int(keep.sum())
    rate = np.log(B.array_rate(v, a0))[None, :]
    tg_idx = PF.index_of(v.target_orfs_current)
    out = np.empty((v.n_targets, v.n_arrays), dtype=np.float32)
    for i in range(0, v.n_targets, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, v.n_targets))
        n = sl.stop - sl.start
        Xc = _lr_design(tg_idx[sl], cidx, names, v, np.broadcast_to(rate, (n, v.n_arrays)))
        out[sl] = clf.decision_function(sc.transform(Xc)).reshape(n, v.n_arrays)
    return out


# B5_external_cf: B5_external and B5_external_es both scored BELOW B1 (0.041 vs 0.043) although early stopping did not
# bind. Cause (diagnosed after those runs): with an exact leave-one-query-out array rate, positives of an array have a
# rate exactly 1/(n-1)-shifted below its negatives; a tree that can also identify the array (array_smf, exact rate
# level) decodes the training label from that shift — target-encoding leakage that cannot transfer to test rows.
# Fork: same HGB with the array rate CROSS-FITTED over K = ctx.n_folds query folds (rows of fold k get the rate over
# the other folds' rows: constant within an array column of a fold, so no per-row label shift).
def crossfit_array_rate(v: ScreenView, ctx: MethodContext, a0: float) -> np.ndarray:
    folds = ctx.folds(v)
    obs = v.obs.astype(np.float64)
    pos = np.nan_to_num(v.Y).astype(np.float64)
    k_all, n_all = pos.sum(0), obs.sum(0)
    out = np.empty((v.n, v.n_arrays), dtype=np.float64)
    for f in range(ctx.n_folds):
        r = folds == f
        k, n = k_all - pos[r].sum(0), n_all - obs[r].sum(0)
        out[r] = B.shrunk_rate(k, n, k.sum() / n.sum(), a0)[None, :]
    return out


@register("B5_external_cf", defaults=B5_DEFAULTS, logit_scale=True,
          description="B5_external with a K-fold cross-fitted (not LOO) array rate in training")
def b5_external_cf(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    from sklearn.ensemble import HistGradientBoostingClassifier
    names = PF.feature_names()
    a0 = ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)
    cidx = PF.index_of(v.array_orfs_current)
    y = v.Y.ravel()
    keep = ~np.isnan(y)
    X = _flatten(_pair_stack(PF.index_of(v.train_orfs_current), cidx, names)
                 + _array_cols(v, ctx, np.log(crossfit_array_rate(v, ctx, a0))), keep)
    clf = HistGradientBoostingClassifier(
        max_iter=ctx.get("max_iter"), learning_rate=ctx.get("learning_rate"), max_leaf_nodes=ctx.get("max_leaf_nodes"),
        min_samples_leaf=ctx.get("min_samples_leaf"), l2_regularization=ctx.get("l2_regularization"),
        early_stopping=False, random_state=int(ctx.model_seed))
    clf.fit(X, y[keep].astype(np.int8))
    del X
    ctx.info["n_training_cells"] = int(keep.sum())
    rate = np.log(B.array_rate(v, a0))[None, :]
    tg_idx = PF.index_of(v.target_orfs_current)
    out = np.empty((v.n_targets, v.n_arrays), dtype=np.float32)
    for i in range(0, v.n_targets, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, v.n_targets))
        n = sl.stop - sl.start
        Xc = _flatten(_pair_stack(tg_idx[sl], cidx, names) + _array_cols(v, ctx, np.broadcast_to(rate, (n, v.n_arrays))))
        out[sl] = clf.decision_function(Xc).reshape(n, v.n_arrays)
    return out
