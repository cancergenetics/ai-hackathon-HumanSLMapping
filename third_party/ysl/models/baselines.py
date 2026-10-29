"""Baselines B0-B4 and the cross-fitted stacked LR (DESIGN.md §6, §9), ported from the pilot to the ScreenView API.

Every method is `fit_predict(view, ctx) -> float32 [n_targets x n_arrays]` (full grid, finite) and reads only the view.
Notation: x = held-out query ORF (target), a = array strain, t = training query ORF. Rates are shrunk towards the
global training rate g with `rate_prior` pseudo-positives: rate = (k + a0) / (n + a0 / g)   (Beta prior, mean g).

  array_rate(a)      : shrunk positive rate of array a over training rows                                      [B1]
  query_prop(x)      : shrunk positive rate of x AS AN ARRAY over training rows; NaN if x has no array strain
  reciprocal(x, a)   : eps / signed -log10 P / label of (a's ORF as training query, x's array strain); NaN if a's
                       ORF is not a training query, x has no array strain or the cell is unobserved             [B3]
  kNN, query side    : score(x,a) = (sum_t w(x,t) Y(t,a) + m rate_a) / (sum_t w(x,t) obs(t,a) + m), top-k t   [B4]
       'col'   : sim = Pearson(x's array column, t's array column) over training rows (t on the array)
       'mixed' : sim = Pearson(E(q, x), E(t, q)) over training queries q on the array (x's column vs t's row).
                 NOTE (DESIGN §8): gains on the reciprocal stratum come from its q = a term (an implicit reciprocal x
                 degree feature), not from profile similarity.
  kNN, array side    : score(x,a) = weighted avg over a's neighbours t (Pearson of array columns) of Y(t, x)
  Weights: w = max(sim, 0)^p on each row's top-k, normalised to sum 1 (pilot defaults k = 50, p = 1, m = 0.05).

Cross-fitted LRs (B2_lr, B3_lr, stack_lr): training features for fold k come from view.without_queries(fold k)
(K = 10), with exact leave-one-query-out values for the additive statistics (array rate, query propensity) when
params["additive"] == "loo" (DESIGN §6); test features from the full view — one feature function for both.
NaN features become 0 only inside the LR design matrix, always next to an availability flag.
"""
from __future__ import annotations

import numpy as np

from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

EPS_CLIP = (-1.0, 0.5)
SP_CLIP = 12.0
DEFAULT_RATE_PRIOR = 1.0
KNN_DEFAULTS = {"k": 50, "p": 1.0, "m": 0.05}


# ----------------------------------------------------------------------------- similarity
def masked_pearson(A: np.ndarray, B: np.ndarray, min_n: int = 10) -> np.ndarray:
    """Pairwise-complete Pearson between columns of A [n x p] and B [n x q] (NaN = missing) -> [p x q].
    Pairs with < min_n shared entries or zero variance get 0."""
    MA, MB = ~np.isnan(A), ~np.isnan(B)
    A0 = np.where(MA, A, 0.0).astype(np.float64)
    B0 = np.where(MB, B, 0.0).astype(np.float64)
    MA, MB = MA.astype(np.float64), MB.astype(np.float64)
    n = MA.T @ MB
    sa, sb = A0.T @ MB, MA.T @ B0
    with np.errstate(divide="ignore", invalid="ignore"):
        cov = A0.T @ B0 - sa * sb / n
        va = (A0 ** 2).T @ MB - sa ** 2 / n
        vb = MA.T @ (B0 ** 2) - sb ** 2 / n
        r = cov / np.sqrt(va * vb)
    bad = (n < min_n) | ~(va > 1e-12) | ~(vb > 1e-12) | ~np.isfinite(r)
    r[bad] = 0.0
    return r.astype(np.float32)


def topk_weights(sim: np.ndarray, k: int, p: float) -> np.ndarray:
    """Dense weight matrix: max(sim,0)^p on each row's top-k, zero elsewhere, rows normalised to sum 1.
    Ties at the k-th similarity are broken by candidate index (stable), independent of target order."""
    W = np.zeros_like(sim, dtype=np.float64)
    k = min(k, sim.shape[1])
    idx = np.argsort(-sim, axis=1, kind="stable")[:, :k]
    rows = np.arange(sim.shape[0])[:, None]
    W[rows, idx] = np.maximum(sim[rows, idx], 0.0) ** p
    s = W.sum(1, keepdims=True)
    np.divide(W, s, out=W, where=s > 0)
    return W


# ----------------------------------------------------------------------------- rates
def shrunk_rate(k, n, g: float, a0: float = DEFAULT_RATE_PRIOR):
    return (np.asarray(k, dtype=np.float64) + a0) / (np.asarray(n, dtype=np.float64) + a0 / g)


def global_rate(v: ScreenView) -> float:
    K, N = v.global_counts()
    return K / N


def array_rate(v: ScreenView, a0: float = DEFAULT_RATE_PRIOR) -> np.ndarray:
    k, n = v.array_counts()
    return shrunk_rate(k, n, global_rate(v), a0).astype(np.float32)


def query_prop(v: ScreenView, cols: np.ndarray, a0: float = DEFAULT_RATE_PRIOR) -> np.ndarray:
    """x's array-side rate over the view's training rows (NaN if x has no array strain)."""
    out = np.full(len(cols), np.nan, dtype=np.float32)
    on = cols >= 0
    k, n = v.array_counts()
    out[on] = shrunk_rate(k[cols[on]], n[cols[on]], global_rate(v), a0)
    return out


def loo_additive(v: ScreenView, a0: float = DEFAULT_RATE_PRIOR) -> dict:
    """Exact leave-one-query-out additive features for every TRAINING row t of v (self-exclusion, DESIGN §3):
    array_rate [n_train x n_arrays] and query_prop [n_train] computed from the training rows other than t."""
    kl, nl = v.loo_array_counts()
    Kl, Nl = v.loo_global_counts()
    g = (Kl / Nl)[:, None]
    ar = ((kl + a0) / (nl + a0 / g)).astype(np.float32)
    qp = np.full(v.n, np.nan, dtype=np.float32)
    on = np.flatnonzero(v.train_col >= 0)
    c = v.train_col[on]
    qp[on] = (kl[on, c] + a0) / (nl[on, c] + a0 / g[on, 0])
    return {"array_rate": ar, "query_prop": qp, "global_rate": (Kl / Nl).astype(np.float64)}


# ----------------------------------------------------------------------------- reciprocal
def reciprocal(v: ScreenView, cols: np.ndarray) -> dict:
    """eps, SP, label of (a's ORF as training query, x as array) -> [n_t x n_a], NaN where unavailable."""
    n_t, n_a = len(cols), v.n_arrays
    out = {k: np.full((n_t, n_a), np.nan, dtype=np.float32) for k in ("eps", "sp", "label")}
    on = np.flatnonzero(cols >= 0)
    rows = v.arr_tr
    ra = np.flatnonzero(rows >= 0)
    for k, M in (("eps", v.E), ("sp", v.SP), ("label", v.Y)):
        out[k][np.ix_(on, ra)] = M[rows[ra]][:, cols[on]].T
    return out


# ----------------------------------------------------------------------------- kNN
def _training_on_array(v: ScreenView) -> np.ndarray:
    return np.flatnonzero(v.train_col >= 0)


def query_side_similarity(v: ScreenView, cols_on: np.ndarray, kind: str) -> tuple[np.ndarray, np.ndarray]:
    q_on = _training_on_array(v)
    if kind == "col":
        return masked_pearson(v.E[:, cols_on], v.E[:, v.train_col[q_on]]), q_on
    if kind == "mixed":
        A = v.E[q_on][:, cols_on]                  # E(q, x) over training queries q on the array
        B = v.E[:, v.train_col[q_on]].T            # E(t, q) for every training t
        return masked_pearson(A, B), np.arange(v.n)
    raise ValueError(kind)


def array_side_similarity(v: ScreenView, kind: str) -> tuple[np.ndarray, np.ndarray]:
    q_on = _training_on_array(v)
    if kind == "col":
        return masked_pearson(v.E, v.E[:, v.train_col[q_on]]), q_on
    raise ValueError(kind)


def _weighted_avg(W: np.ndarray, Yc: np.ndarray, prior: np.ndarray, m: float) -> np.ndarray:
    O = (~np.isnan(Yc)).astype(np.float64)
    Y0 = np.nan_to_num(Yc).astype(np.float64)
    return ((W @ Y0 + m * prior) / (W @ O + m)).astype(np.float32)


def knn_query_side(v: ScreenView, cols: np.ndarray, kind: str, k: int, p: float, m: float, a0: float) -> np.ndarray:
    """[n_t x n_a]; NaN rows for targets without an array strain."""
    S = np.full((len(cols), v.n_arrays), np.nan, dtype=np.float32)
    on = np.flatnonzero(cols >= 0)
    if len(on):
        sim, cand = query_side_similarity(v, cols[on], kind)
        S[on] = _weighted_avg(topk_weights(sim, k, p), v.Y[cand], array_rate(v, a0)[None, :], m)
    return S


def knn_array_side(v: ScreenView, cols: np.ndarray, kind: str, k: int, p: float, m: float, a0: float) -> np.ndarray:
    S = np.full((len(cols), v.n_arrays), np.nan, dtype=np.float32)
    on = np.flatnonzero(cols >= 0)
    if len(on):
        sim, cand = array_side_similarity(v, kind)                  # [n_a x n_cand]
        Vx = v.Y[cand][:, cols[on]]                                  # Y(t, x)
        prior = query_prop(v, cols[on], a0)[None, :]
        S[on] = _weighted_avg(topk_weights(sim, k, p), Vx, prior, m).T
    return S


KNN_FEATURES = {  # fixed a priori (pilot defaults; re-selection on S1 dev is a later PI step)
    "knn_qs_col": ("query", "col"),
    "knn_qs_mixed": ("query", "mixed"),
    "knn_as_col": ("array", "col"),
}


# ----------------------------------------------------------------------------- one feature function
def pair_features(v: ScreenView, params: dict, knn: tuple[str, ...] = ()) -> dict:
    """Per-pair features for v's targets ([n_targets x n_arrays] float32, NaN = unavailable)."""
    a0 = params.get("rate_prior", DEFAULT_RATE_PRIOR)
    cols = v.target_col
    n_t, n_a = len(cols), v.n_arrays
    f = {"cols": cols, "global_rate": global_rate(v)}
    f["array_rate"] = np.broadcast_to(array_rate(v, a0)[None, :], (n_t, n_a))
    f["array_smf"] = np.broadcast_to(v.array_smf[None, :], (n_t, n_a))
    f["query_prop"] = np.broadcast_to(query_prop(v, cols, a0)[:, None], (n_t, n_a))
    r = reciprocal(v, cols)
    f["recip_eps"], f["recip_sp"], f["recip_label"] = r["eps"], r["sp"], r["label"]
    kp = {kk: params.get(kk, dv) for kk, dv in KNN_DEFAULTS.items()}
    for name in knn:
        side, kind = KNN_FEATURES[name]
        fn = knn_query_side if side == "query" else knn_array_side
        f[name] = fn(v, cols, kind, kp["k"], kp["p"], kp["m"], a0)
    return f


def _loo_override(f: dict, loo: dict, rows: np.ndarray) -> dict:
    """Replace additive features of fold targets (= training rows `rows` of the parent view) by exact LOO values."""
    f = dict(f)
    n_t, n_a = f["recip_eps"].shape
    f["array_rate"] = loo["array_rate"][rows]
    f["query_prop"] = np.broadcast_to(loo["query_prop"][rows][:, None], (n_t, n_a))
    return f


def design_matrix(f: dict, names: tuple[str, ...], rows=slice(None)) -> np.ndarray:
    """LR inputs [n_rows*n_a x d] for target rows `rows`; unavailable values -> 0 next to an availability flag."""
    cols_ = []
    sub = {k: (np.asarray(v)[rows] if isinstance(v, np.ndarray) and v.ndim >= 1 else v) for k, v in f.items()}
    n_t, n_a = sub["recip_eps"].shape
    on_arr = np.broadcast_to((sub["cols"] >= 0)[:, None], (n_t, n_a)).astype(np.float32)
    qp = np.where(np.isnan(sub["query_prop"]), 1.0, sub["query_prop"])
    rav = ~np.isnan(sub["recip_eps"])
    smf = sub["array_smf"]
    for nm in names:
        if nm == "log_array_rate":
            cols_.append(np.log(sub["array_rate"]))
        elif nm == "array_smf":
            cols_.append(np.where(np.isnan(smf), 0.0, smf))
        elif nm == "array_smf_avail":
            cols_.append((~np.isnan(smf)).astype(np.float32))
        elif nm == "log_query_prop":
            cols_.append(np.log(qp) * on_arr)
        elif nm == "on_array":
            cols_.append(on_arr)
        elif nm == "recip_avail":
            cols_.append(rav.astype(np.float32))
        elif nm == "recip_negeps":
            cols_.append(np.where(rav, -np.clip(sub["recip_eps"], *EPS_CLIP), 0.0))
        elif nm == "recip_label":
            cols_.append(np.where(rav, sub["recip_label"], 0.0))
        elif nm == "recip_negsp":
            cols_.append(np.where(rav, -np.clip(sub["recip_sp"], -SP_CLIP, SP_CLIP), 0.0))
        elif nm.startswith("log_knn_"):
            k = sub[nm[4:]]
            cols_.append(np.where(np.isnan(k), 0.0, np.log(np.where(np.isnan(k), 1.0, k))))
        else:
            raise ValueError(nm)
    return np.stack([np.asarray(c, dtype=np.float32).ravel() for c in cols_], axis=1)


LR_FEATURES = {
    "B2_lr": ("log_array_rate", "log_query_prop", "on_array"),
    "B3_lr": ("log_array_rate", "log_query_prop", "on_array", "recip_avail", "recip_negeps", "recip_label",
              "recip_negsp"),
    "stack_lr": ("log_array_rate", "log_query_prop", "on_array", "recip_avail", "recip_negeps", "recip_label",
                 "recip_negsp", "log_knn_qs_col", "log_knn_qs_mixed", "log_knn_as_col", "array_smf",
                 "array_smf_avail"),
}


def _knn_needed(names) -> tuple[str, ...]:
    return tuple(n[4:] for n in names if n.startswith("log_knn_"))


def crossfit_training_set(v: ScreenView, ctx: MethodContext, names: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Cross-fitted design matrix + labels over all observed training cells (K = ctx.n_folds query folds)."""
    folds = ctx.folds(v)
    knn = _knn_needed(names)
    loo = loo_additive(v, ctx.get("rate_prior", DEFAULT_RATE_PRIOR)) if ctx.get("additive", "loo") == "loo" else None
    X, ys = [], []
    for k in range(ctx.n_folds):
        rows = np.flatnonzero(folds == k)
        fv = v.without_queries(v.train_orfs[rows])
        f = fold_features(fv, ctx.params, knn)
        if loo is not None:
            f = _loo_override(f, loo, rows)
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        X.append(design_matrix(f, names)[keep])
        ys.append(y[keep])
    return np.concatenate(X), np.concatenate(ys).astype(np.int8)


def fold_features(fv: ScreenView, params: dict, knn: tuple[str, ...]) -> dict:
    """Features of a fold view's targets — the SAME function as test time (pair_features)."""
    return pair_features(fv, params, knn)


def fit_lr(X: np.ndarray, y: np.ndarray, ctx: MethodContext):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=ctx.get("max_iter", 500), C=ctx.get("C", 1.0)).fit(sc.transform(X), y)
    return sc, clf


def predict_lr(fitted, f: dict, names, chunk: int) -> np.ndarray:
    sc, clf = fitted
    n_t, n_a = f["recip_eps"].shape
    out = np.empty((n_t, n_a), dtype=np.float32)
    for i in range(0, n_t, chunk):
        sl = slice(i, min(i + chunk, n_t))
        out[sl] = clf.decision_function(sc.transform(design_matrix(f, names, sl))).reshape(-1, n_a)
    return out


def crossfit_lr_predict(v: ScreenView, ctx: MethodContext, names) -> np.ndarray:
    X, y = crossfit_training_set(v, ctx, names)
    fitted = fit_lr(X, y, ctx)
    del X
    f = pair_features(v, ctx.params, _knn_needed(names))
    S = predict_lr(fitted, f, names, ctx.chunk_targets)
    ctx.info["lr_coef"] = dict(zip(names, np.round(fitted[1].coef_[0], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    return S


# ----------------------------------------------------------------------------- registered methods
def _fill_off_array(S: np.ndarray, v: ScreenView, a0: float) -> np.ndarray:
    """Targets without an array strain get the array rate (pilot fallback)."""
    rate = array_rate(v, a0)
    return np.where(np.isnan(S), rate[None, :], S).astype(np.float32)


@register("B0_random", description="uniform random scores (seeded)")
def b0_random(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return ctx.rng("B0_random").random((v.n_targets, v.n_arrays), dtype=np.float32)


@register("B1_array_rate", defaults={"rate_prior": DEFAULT_RATE_PRIOR}, description="shrunk array degree rate")
def b1_array_rate(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    r = array_rate(v, ctx.get("rate_prior", DEFAULT_RATE_PRIOR))
    return np.ascontiguousarray(np.broadcast_to(r[None, :], (v.n_targets, v.n_arrays)), dtype=np.float32)


@register("B1p_array_smf", description="-SMF of the array gene (missing SMF -> median)")
def b1p_array_smf(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    smf = np.where(np.isnan(v.array_smf), np.nanmedian(v.array_smf), v.array_smf)
    return np.ascontiguousarray(np.broadcast_to(-smf[None, :], (v.n_targets, v.n_arrays)), dtype=np.float32)


@register("B2_product", defaults={"rate_prior": DEFAULT_RATE_PRIOR},
          description="array rate x query propensity (global rate for off-array x)")
def b2_product(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    a0 = ctx.get("rate_prior", DEFAULT_RATE_PRIOR)
    qp = query_prop(v, v.target_col, a0)
    qp = np.where(np.isnan(qp), global_rate(v), qp)
    return (qp[:, None] * array_rate(v, a0)[None, :]).astype(np.float32)


_LR_DEFAULTS = {"rate_prior": DEFAULT_RATE_PRIOR, "additive": "loo", "C": 1.0, "max_iter": 500}


@register("B2_lr", defaults=_LR_DEFAULTS, logit_scale=True,
          description="cross-fitted LR: log array rate, log query propensity, on-array flag")
def b2_lr(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return crossfit_lr_predict(v, ctx, LR_FEATURES["B2_lr"])


@register("B3_lr", defaults=_LR_DEFAULTS, logit_scale=True,
          description="B2_lr + reciprocal (availability, -eps, label, -signed log10 P)")
def b3_lr(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return crossfit_lr_predict(v, ctx, LR_FEATURES["B3_lr"])


@register("stack_lr", defaults=_LR_DEFAULTS | KNN_DEFAULTS, logit_scale=True,
          description="permanent comparator: B3_lr + 3 kNN scores (log) + array SMF; K=10 cross-fitted")
def stack_lr(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return crossfit_lr_predict(v, ctx, LR_FEATURES["stack_lr"])


def _knn_method(name: str):
    side, kind = KNN_FEATURES[name]

    def fn(v: ScreenView, ctx: MethodContext) -> np.ndarray:
        a0 = ctx.get("rate_prior", DEFAULT_RATE_PRIOR)
        kp = {kk: ctx.get(kk, dv) for kk, dv in KNN_DEFAULTS.items()}
        f = knn_query_side if side == "query" else knn_array_side
        return _fill_off_array(f(v, v.target_col, kind, kp["k"], kp["p"], kp["m"], a0), v, a0)
    return fn


for _n, _label in (("knn_qs_col", "B4_knn_qs_col"), ("knn_qs_mixed", "B4_knn_qs_mixed"),
                   ("knn_as_col", "B4_knn_as_col")):
    register(_label, defaults={"rate_prior": DEFAULT_RATE_PRIOR} | KNN_DEFAULTS,
             description=f"kNN {KNN_FEATURES[_n][0]}-side '{KNN_FEATURES[_n][1]}' (off-array x -> array rate)")(
        _knn_method(_n))
