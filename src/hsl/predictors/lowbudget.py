"""Low-budget variants of P* = stack_lr_wr2 (DESIGN v3.1 §4, §8 predictor-change gate, §9 P2).

`stack_lr_wr2x` is upstream `foldin2.stack_wr2_predict` (same folds, seeds, LOO additive statistics, within-row
features, standardisation, LR) with optional, individually switchable changes. With every switch off it reproduces
stack_lr_wr2 bit for bit (tests/test_p2_predictors.py). Switches (ctx.params):

  gi_n0  > 0     n-dependent GI-kernel weight in the multi-kernel fold-in `mk`: w_gi(view) = n_view/(n_view + gi_n0)
                 (mk__w_gi per view; with upstream's lam_mode 'gi' this is equivalent to up-weighting the external
                 kernels and the ridge by (n + n0)/n relative to the GI kernel).
  extra          list of extra LR column groups:
     "prop"      upstream query-propensity sub-model (propensity.propensity_predict on every fold view and the full
                 view; nested cross-fitted), columns prop_logit, prop_logit x off-array (as stack_lr_prop / m1_hgb).
     "dense"     dense degree statistics from the screened rows' epsilon (not only the strong-SL label), which are
                 much less noisy than label rates when few queries are screened:
                   arr_mild   log shrunk rate of mild negative epsilon (eps < -0.08) in array a's column
                   arr_negeps shrunk mean of -clip(eps) in array a's column
                   q_mild, q_negeps  the same statistics in the target's own array-side column (on-array targets;
                              0 for off-array targets, next to upstream's on_array flag)
                 Training rows use exact leave-one-query-out values (as upstream's loo_additive for array rate and
                 query propensity); test rows use the full view. Only screened rows' epsilon enters.
     "dense_arr" the array-side dense columns only (arr_mild, arr_negeps). Query-level columns are learned from
                 only n distinct training queries, which overfits at n = 30 (P2 interim: "dense" -3% at n = 30).
     "aprior"    empirical-Bayes array-degree prior: per view, a Poisson GLM of each array's strong-SL count over
                 the view's screened rows (exposure = observed cells) on a cubic spline of array SMF + upstream's
                 external gene scalars of the array gene (ysl.features.genes.gene_matrix; label-free; median-imputed,
                 standardised, with upstream's has_* availability flags). One column: arr_prior = log predicted rate
                 of array a (array side only, for the same reason as dense_arr).
                 Fold views fit their own GLM (their rows exclude the fold), the full view one for the test targets.
  wr_min_n > 0   below this number of screened queries the within-row features are dropped (stack_lr design).

Data access: only the ScreenView (screened rows + target names/array columns) and upstream's label-free external
stores — no genome coordinates, no other GI source, no query SMF, no off-array indicator beyond upstream's flag.

`feature_cache()` (context manager) memoises the expensive per-view features (kNN, within-row scores, propensity)
by view content, so several configurations evaluated on the same screened set share them; results are identical
with or without the cache.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard; must be first)

import contextlib
import hashlib
import json
import warnings

import numpy as np

from ysl.models import baselines as B
from ysl.models import completion as CM
from ysl.models import foldin2 as FD
from ysl.models import propensity as PR
from ysl.models.registry import METHODS, register
from ysl.models.view import MethodContext, ScreenView

STACK_BASE = CM.STACK_BASE
WR2_NAMES = tuple(FD.WR2_STACK_NAMES)
MILD_EPS = -0.08                 # "mild negative" epsilon threshold (Costanzo's intermediate cut-off)
NEGEPS_PRIOR_M = 5.0             # pseudo-observations shrinking the mean -eps towards the global mean
DENSE_COLS = ("arr_mild", "arr_negeps", "q_mild", "q_negeps")
APRIOR_COLS = ("arr_prior",)
APRIOR_ALPHA = 1e-3              # L2 penalty of the Poisson GLM (fixed a priori, not tuned)
APRIOR_KNOTS = 6

LOWB_DEFAULTS = B._LR_DEFAULTS | B.KNN_DEFAULTS | PR.PROP_DEFAULTS | {
    "wr": list(WR2_NAMES), "prop_fset": "all", "gi_n0": 0.0, "extra": [], "wr_min_n": 0}

# ============================================================================ feature cache
_CACHE: dict | None = None


@contextlib.contextmanager
def feature_cache():
    global _CACHE
    old, _CACHE = _CACHE, {}
    try:
        yield _CACHE
    finally:
        _CACHE = old


def _vkey(v: ScreenView) -> str:
    h = hashlib.sha1("|".join(v.train_orfs.tolist()).encode())
    h.update(b"#")
    h.update("|".join(v.target_orfs.tolist()).encode())
    return h.hexdigest()


def _memo(kind: str, v: ScreenView, spec, fn):
    if _CACHE is None:
        return fn()
    key = (kind, _vkey(v), json.dumps(spec, sort_keys=True, default=str))
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


# ============================================================================ feature functions (one per view)
def _knn_params(params: dict) -> dict:
    keys = {**B.KNN_DEFAULTS, "rate_prior": B.DEFAULT_RATE_PRIOR}
    return {k: params.get(k, d) for k, d in keys.items()}


def pair_features(v: ScreenView, params: dict, knn) -> dict:
    return _memo("pair", v, [_knn_params(params), list(knn)], lambda: B.pair_features(v, params, knn))


def _wr_view_params(v: ScreenView, params: dict) -> dict:
    p = dict(params)
    n0 = float(p.get("gi_n0") or 0.0)
    if n0 > 0:
        p["mk__w_gi"] = float(v.n) / (float(v.n) + n0)
    return p


def wr_features(v: ScreenView, params: dict, names, seed: int) -> dict:
    """foldin2.wr2_features computed feature by feature (each is an independent function of the view and seed)."""
    p = _wr_view_params(v, params)
    out = {}
    for nm in names:
        own = {k: val for k, val in p.items() if k.startswith(nm + "__")}
        out[nm] = _memo("wr", v, [nm, own, int(seed)], lambda nm=nm: FD.wr2_features(v, p, (nm,), seed)[nm])
    return out


def prop_feature(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    spec = {k: ctx.params.get(k) for k in ("prop_model", "prop_k", "prop_folds", "prop_hgb_iter", "prop_test_view",
                                           "prop_fset")} | {"seed": [ctx.model_seed, ctx.split_seed]}

    def fn():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return np.asarray(PR.propensity_predict(v, ctx), dtype=np.float64)
    return _memo("prop", v, spec, fn)


def _dense_parts(v: ScreenView):
    """Per-row contributions over v's training rows: observed, mild-negative indicator, -clip(eps)."""
    obs = v.obs
    E = np.where(obs, v.E, np.nan).astype(np.float64)
    mild = np.where(obs, (E < MILD_EPS).astype(np.float64), 0.0)
    ne = np.where(obs, -np.clip(np.nan_to_num(E), *B.EPS_CLIP), 0.0)
    return obs.astype(np.float64), mild, ne


def _dense_from_sums(N, M, S, Ntot, Mtot, Stot, a0: float):
    """Shrunk mild-rate (Beta prior with a0 pseudo-positives, mean = global rate) and shrunk mean -eps."""
    g_m = np.maximum(Mtot / np.maximum(Ntot, 1.0), 1e-6)
    g_e = Stot / np.maximum(Ntot, 1.0)
    mild = (M + a0) / (N + a0 / g_m)
    ne = (S + NEGEPS_PRIOR_M * g_e) / (N + NEGEPS_PRIOR_M)
    return np.log(mild), ne


def dense_test(v: ScreenView, a0: float) -> dict:
    """Array-side statistics [n_a] over v's training rows and the targets' query-side values [n_targets]."""
    O, Mi, Ne = _dense_parts(v)
    N, M, S = O.sum(0), Mi.sum(0), Ne.sum(0)
    am, ae = _dense_from_sums(N, M, S, N.sum(), M.sum(), S.sum(), a0)
    c = v.target_col
    qm = np.where(c >= 0, am[np.maximum(c, 0)], 0.0)
    qe = np.where(c >= 0, ae[np.maximum(c, 0)], 0.0)
    return {"arr_mild": am, "arr_negeps": ae, "q_mild": qm, "q_negeps": qe}


def dense_loo(v: ScreenView, a0: float) -> dict:
    """Exact leave-one-query-out versions for every TRAINING row t of v: array side [n x n_a] from the rows other
    than t; query side [n] = t's array-side column over the rows other than t (cf. baselines.loo_additive)."""
    O, Mi, Ne = _dense_parts(v)
    N, M, S = O.sum(0)[None, :] - O, Mi.sum(0)[None, :] - Mi, Ne.sum(0)[None, :] - Ne
    Nt, Mt, St = (O.sum() - O.sum(1))[:, None], (Mi.sum() - Mi.sum(1))[:, None], (Ne.sum() - Ne.sum(1))[:, None]
    am, ae = _dense_from_sums(N, M, S, Nt, Mt, St, a0)
    c = v.train_col
    r = np.arange(v.n)
    qm = np.where(c >= 0, am[r, np.maximum(c, 0)], 0.0)
    qe = np.where(c >= 0, ae[r, np.maximum(c, 0)], 0.0)
    return {"arr_mild": am, "arr_negeps": ae, "q_mild": qm, "q_negeps": qe}


def dense_design(d: dict, n_a: int, rows=slice(None), arr_rows=None) -> np.ndarray:
    """[n_rows*n_a x 4]; `arr_rows` selects rows of 2-D (LOO) array statistics (None: 1-D, broadcast)."""
    qm, qe = np.asarray(d["q_mild"])[rows], np.asarray(d["q_negeps"])[rows]
    n_r = len(qm)
    if arr_rows is None:
        am = np.broadcast_to(d["arr_mild"][None, :], (n_r, n_a))
        ae = np.broadcast_to(d["arr_negeps"][None, :], (n_r, n_a))
    else:
        am, ae = d["arr_mild"][arr_rows], d["arr_negeps"][arr_rows]
    cols = [am, ae, np.broadcast_to(qm[:, None], (n_r, n_a)), np.broadcast_to(qe[:, None], (n_r, n_a))]
    return np.stack([np.asarray(c, dtype=np.float32).ravel() for c in cols], axis=1)


_ADESIGN: dict = {}


def _array_design(v: ScreenView) -> np.ndarray:
    """Label-free design over the view's arrays: SMF spline + standardised external gene scalars (cached)."""
    key = hashlib.sha1("|".join(v.array_orfs_current.tolist()).encode()).hexdigest()
    if key not in _ADESIGN:
        from sklearn.preprocessing import SplineTransformer
        from ysl.features import genes as G
        smf = np.asarray(v.array_smf, dtype=np.float64)
        smf = np.where(np.isnan(smf), np.nanmedian(smf), smf)
        Sp = SplineTransformer(n_knots=APRIOR_KNOTS, degree=3).fit_transform(smf[:, None])
        X = G.gene_matrix(v.array_orfs_current, current=True).to_numpy(dtype=np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(X, 0)
        X = np.where(np.isnan(X), np.nan_to_num(med)[None, :], X)
        sd = X.std(0)
        X = (X - X.mean(0)) / np.where(sd > 0, sd, 1.0)
        _ADESIGN[key] = np.hstack([Sp, X])
    return _ADESIGN[key]


def array_prior(v: ScreenView) -> np.ndarray:
    """log predicted strong-SL rate per array from the view's screened rows (Poisson GLM, exposure = observed)."""
    def fn():
        from sklearn.linear_model import PoissonRegressor
        X = _array_design(v)
        k, n = v.array_counts()
        ok = n > 0
        m = PoissonRegressor(alpha=APRIOR_ALPHA, max_iter=1000)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(X[ok], k[ok] / n[ok], sample_weight=n[ok])
        return np.log(np.maximum(m.predict(X), 1e-6))
    return _memo("aprior", v, [APRIOR_ALPHA, APRIOR_KNOTS], fn)


def aprior_design(ap: np.ndarray, cols: np.ndarray, n_a: int, rows=slice(None)) -> np.ndarray:
    n_r = len(np.asarray(cols)[rows])
    return np.broadcast_to(ap[None, :], (n_r, n_a)).reshape(-1, 1).astype(np.float32)


# ============================================================================ the stack
def stack_predict(v: ScreenView, ctx: MethodContext, inject=None) -> np.ndarray:
    """stack_lr_wr2 with the switches in the module docstring. `inject(view, parent_rows) -> [n_rows*n_a x d]` adds
    columns (evaluation-only diagnostics in scripts/30_p2_*.py; never used by a registered method)."""
    p = ctx.params
    a0 = ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)
    names = tuple(ctx.get("wr", WR2_NAMES))
    if int(p.get("wr_min_n") or 0) > 0 and v.n < int(p["wr_min_n"]):
        names = ()
    extra = tuple(p.get("extra") or ())
    bad = set(extra) - {"prop", "dense", "dense_arr", "aprior"}
    if {"dense", "dense_arr"} <= set(extra):
        raise ValueError("'dense' already contains 'dense_arr'")
    dcols = DENSE_COLS if "dense" in extra else DENSE_COLS[:2] if "dense_arr" in extra else ()
    if bad:
        raise ValueError(f"unknown extra groups {sorted(bad)}")
    base = STACK_BASE
    knn = B._knn_needed(base)
    folds = ctx.folds(v)
    loo = B.loo_additive(v, a0)
    seed = int(ctx.rng("tied").integers(2**31))
    dl = dense_loo(v, a0) if dcols else None
    n_a = v.n_arrays
    Xs, ys = [], []
    for k in range(ctx.n_folds):
        rows = np.flatnonzero(folds == k)
        fv = v.without_queries(v.train_orfs[rows])
        f = B._loo_override(pair_features(fv, p, knn), loo, rows)
        parts = [B.design_matrix(f, base)]
        if names:
            fw, _ = CM._standardise_wr(wr_features(fv, p, names, seed), None)
            parts.append(CM._wr_design(fw, names))
        if "prop" in extra:
            parts.append(PR._prop_design(prop_feature(fv, ctx), f["cols"], n_a))
        if dl is not None:
            parts.append(dense_design({kk: (vv if vv.ndim == 2 else vv[rows]) for kk, vv in dl.items()}, n_a,
                                      arr_rows=rows)[:, :len(dcols)])
        if "aprior" in extra:
            parts.append(aprior_design(array_prior(fv), f["cols"], n_a))
        if inject is not None:
            parts.append(inject(fv, rows))
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        Xs.append(np.hstack(parts)[keep])
        ys.append(y[keep])
    X, y = np.concatenate(Xs), np.concatenate(ys).astype(np.int8)
    del Xs
    sc, clf = B.fit_lr(X, y, ctx)
    del X
    f = pair_features(v, p, knn)
    fw = CM._standardise_wr(wr_features(v, p, names, seed), None)[0] if names else None
    prop = prop_feature(v, ctx) if "prop" in extra else None
    dt = dense_test(v, a0) if dl is not None else None
    ap = array_prior(v) if "aprior" in extra else None
    n_t = v.n_targets
    S = np.empty((n_t, n_a), dtype=np.float32)
    for i in range(0, n_t, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, n_t))
        parts = [B.design_matrix(f, base, sl)]
        if names:
            parts.append(CM._wr_design(fw, names, sl))
        if prop is not None:
            parts.append(PR._prop_design(prop, f["cols"], n_a, sl))
        if dt is not None:
            parts.append(dense_design(dt, n_a, sl)[:, :len(dcols)])
        if ap is not None:
            parts.append(aprior_design(ap, f["cols"], n_a, sl))
        if inject is not None:
            parts.append(inject(v, None)[sl.start * n_a: sl.stop * n_a])
        S[sl] = clf.decision_function(sc.transform(np.hstack(parts))).reshape(-1, n_a)
    cols = list(base) + list(names) + (list(PR.PROP_COLS) if "prop" in extra else []) + \
        list(dcols) + (list(APRIOR_COLS) if "aprior" in extra else [])
    ctx.info["lr_coef"] = dict(zip(cols, np.round(clf.coef_[0][:len(cols)], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    ctx.info["wr_names"] = list(names)
    return S


def lowb_method(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return stack_predict(v, ctx)


if "stack_lr_wr2x" not in METHODS:
    register("stack_lr_wr2x", defaults=dict(LOWB_DEFAULTS), logit_scale=True,
             description="stack_lr_wr2 + optional low-budget switches (n-dependent GI weight gi_n0, extra groups "
                         "'prop'/'dense', wr_min_n); all off = stack_lr_wr2 bit for bit")(lowb_method)


# ============================================================================ frozen P* (P2 freeze, config/predictor.yaml)
# Passed the DESIGN §8 predictor-change gate on dev (results/p2/summary.md); chosen as the passing configuration with
# the largest mean relative gain. config/predictor.yaml records the same recipe (tests assert they agree).
PSTAR_NAME = "stack_lr_wr2_lb"
PSTAR_PARAMS = {"extra": ["dense_arr", "aprior"], "gi_n0": 300.0}
if PSTAR_NAME not in METHODS:
    register(PSTAR_NAME, defaults=LOWB_DEFAULTS | PSTAR_PARAMS, logit_scale=True,
             description="frozen P* (P2): stack_lr_wr2 + array-side dense degree statistics + empirical-Bayes array "
                         "prior (Poisson GLM on SMF spline + external gene scalars) + n-dependent GI weight "
                         "n/(n+300); = stack_lr_wr2x with PSTAR_PARAMS")(lowb_method)
