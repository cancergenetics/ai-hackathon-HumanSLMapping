"""Query-propensity sub-model (DESIGN.md §6 B2 sub-model, §8 roadmap rule) and the stacked LR that uses it.

Target: for each screened (training) query t, R_t = logit((k_t + 0.5) / (n_t + 1)) — its smoothed query-row positive
rate (the same smoothing as the evaluator's oracle). The sub-model predicts R_x for held-out queries x from:

  in-matrix ("col_*", only if x has an array strain; computed from the view's training rows only)
    col_logit        logit of x's shrunk array-side rate over training queries                (= baselines.query_prop)
    col_negeps_mean, col_eps_q10, col_eps_q50, col_eps_var   (eps clipped to [-1, 0.5])
    col_negsp_mean   mean of -signed log10 P;  col_frac_mild (eps < -0.08), col_frac_sig (signed P < -2),
    col_frac_pos     (eps > 0.08 and P < 0.05)
    col_log_oe       log((k + 1) / (sum_t r_t + 1)): x's positives vs the expectation from the row rates r_t of the
                     training queries observed in x's column (corrects for which queries cover x)
    col_partner_hub  mean R_t over training queries t that hit x (shrunk to the mean R with m = 2 pseudo-partners)
    smf              x's SMF from its array strain
    nb_col / nb_mixed  top-k similarity-weighted mean R_t of training queries most similar to x by array-column
                     Pearson (training queries on the array) / by x's column vs t's ROW ('mixed' kNN similarity)
  external (all x; no GI except the training queries' R_t used as neighbour values)
    gene scalars + ESM PCs (ysl.features.genes.gene_matrix)
    nb_coexp         top-k coexpression(z>0)-weighted mean R_t;  nb_esm: mean R_t of the k nearest by ESM cosine;
    nb_apms          mean R_t over AP-MS partners among training queries (NaN if none), apms_n_partners_train

Cross-fitting (DESIGN §3/§6): training queries' features come from view.without_queries(fold k), K = 10, i.e. with
their own rows (and their fold's rows) removed — the same function `propensity_features` as at test time.
`propensity_predict(view, ctx)` fits on the cross-fitted training features of `view` and predicts `view`'s targets;
stack_lr_prop calls it on every outer fold view (nested cross-fitting) and on the full view.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ysl.models import baselines as B
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView, _ViewSource

R_SMOOTH = 0.5
MIN_SPLIT_TRAIN = 10
EPS_CLIP = B.EPS_CLIP
PROP_DEFAULTS = {"prop_model": "split_enet", "prop_k": 20, "prop_folds": 10, "prop_hgb_iter": 200,
                 "prop_test_view": "matched"}

COL_FEATURES = ["col_logit", "col_negeps_mean", "col_eps_q10", "col_eps_q50", "col_eps_var", "col_negsp_mean",
                "col_frac_mild", "col_frac_sig", "col_frac_pos", "col_log_oe", "col_partner_hub", "smf",
                "nb_col", "nb_col_maxsim", "nb_mixed", "nb_mixed_maxsim"]
EXT_NB_FEATURES = ["nb_coexp", "nb_coexp_maxsim", "nb_esm", "nb_esm_maxsim", "nb_apms", "apms_frac_partners_train"]
FEATURE_SETS = {
    "base": ["col_logit", "smf"],
    "inmatrix": COL_FEATURES,
    "external": "GENE",                         # gene scalars + ESM PCs only (pure cold start)
    "external_nb": "GENE+EXT_NB",               # + external-similarity neighbour propensities
    "all": "GENE+EXT_NB+COL",
}


def _logit(p):
    p = np.clip(np.asarray(p, dtype=np.float64), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


# ----------------------------------------------------------------------------- targets
def row_logit_rate(v: ScreenView) -> np.ndarray:
    """R_t for every training row of v."""
    k = np.nansum(v.Y, axis=1).astype(np.float64)
    n = v.obs.sum(1).astype(np.float64)
    return _logit((k + R_SMOOTH) / (n + 2 * R_SMOOTH))


# ----------------------------------------------------------------------------- feature blocks
def _topk_mean(sim: np.ndarray, vals: np.ndarray, k: int, weighted: bool) -> tuple[np.ndarray, np.ndarray]:
    """Per row of sim [n_t x n_c] (NaN = unavailable): mean of vals over the top-k candidates (weighted by max(sim,0)
    if `weighted`), and the mean top-k similarity. NaN where a row has no available candidate (or zero weight)."""
    S = np.where(np.isnan(sim), -np.inf, sim).astype(np.float64)
    n_av = np.isfinite(S).sum(1)
    kk = min(k, S.shape[1])
    if kk == 0:
        nan = np.full(S.shape[0], np.nan)
        return nan, nan
    idx = np.argsort(-S, axis=1, kind="stable")[:, :kk]
    top = np.take_along_axis(S, idx, axis=1)
    ok = np.isfinite(top)
    w = (np.maximum(top, 0.0) if weighted else np.ones_like(top)) * ok
    ws = w.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (w * vals[idx]).sum(1) / ws
        msim = np.where(ok, top, 0).sum(1) / ok.sum(1)
    bad = (n_av == 0) | ~(ws > 0)
    mean[bad] = np.nan
    msim[n_av == 0] = np.nan
    return mean, msim


def column_features(v: ScreenView, cols: np.ndarray, R: np.ndarray, k: int) -> dict:
    """In-matrix features for targets with array columns `cols` (-1 = off-array -> NaN)."""
    n_t = len(cols)
    out = {c: np.full(n_t, np.nan) for c in COL_FEATURES}
    on = np.flatnonzero(cols >= 0)
    if len(on) == 0:
        return out
    c = cols[on]
    Y, E, SP = v.Y[:, c].astype(np.float64), v.E[:, c].astype(np.float64), v.SP[:, c].astype(np.float64)
    obs = ~np.isnan(Y)
    n = obs.sum(0).astype(np.float64)
    kpos = np.nansum(Y, 0)
    g = B.global_rate(v)
    out["col_logit"][on] = _logit(B.shrunk_rate(kpos, n, g))
    Ec = np.where(obs, np.clip(E, *EPS_CLIP), np.nan)
    SPc = np.where(obs, np.clip(SP, -B.SP_CLIP, B.SP_CLIP), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out["col_negeps_mean"][on] = -np.nanmean(Ec, 0)
            out["col_eps_q10"][on] = np.nanquantile(Ec, 0.10, axis=0)
            out["col_eps_q50"][on] = np.nanquantile(Ec, 0.50, axis=0)
            out["col_eps_var"][on] = np.nanvar(Ec, 0)
            out["col_negsp_mean"][on] = -np.nanmean(SPc, 0)
        nn = np.maximum(n, 1)
        out["col_frac_mild"][on] = (obs & (E < -0.08)).sum(0) / nn
        out["col_frac_sig"][on] = (obs & (SP < -2.0)).sum(0) / nn
        out["col_frac_pos"][on] = (obs & (E > 0.08) & (SP > -np.log10(0.05))).sum(0) / nn
        r = 1.0 / (1.0 + np.exp(-R))
        expected = (obs * r[:, None]).sum(0)
        out["col_log_oe"][on] = np.log((kpos + 1.0) / (expected + 1.0))
        hit = obs & (np.nan_to_num(Y) > 0.5)
        m = 2.0
        out["col_partner_hub"][on] = ((hit * R[:, None]).sum(0) + m * R.mean()) / (hit.sum(0) + m)
    out["smf"][on] = v.array_smf[c]
    for name in ("col_frac_mild", "col_frac_sig", "col_frac_pos", "col_negeps_mean"):
        out[name][on[n == 0]] = np.nan
    # neighbour propensity by in-matrix similarity
    sim_col, cand = B.query_side_similarity(v, c, "col")
    sim_col = np.where(sim_col == 0.0, np.nan, sim_col)  # masked_pearson returns 0 for unusable pairs
    out["nb_col"][on], out["nb_col_maxsim"][on] = _topk_mean(sim_col, R[cand], k, weighted=True)
    sim_mx, cand = B.query_side_similarity(v, c, "mixed")
    sim_mx = np.where(sim_mx == 0.0, np.nan, sim_mx)
    out["nb_mixed"][on], out["nb_mixed_maxsim"][on] = _topk_mean(sim_mx, R[cand], k, weighted=True)
    return out


def external_neighbour_features(v: ScreenView, R: np.ndarray, k: int) -> dict:
    from ysl.features import genes as G
    tc, rc = v.target_orfs_current, v.train_orfs_current
    out = {}
    z = G.similarity("coexp", tc, rc)
    out["nb_coexp"], out["nb_coexp_maxsim"] = _topk_mean(z, R, k, weighted=True)
    s = G.similarity("esm", tc, rc)
    out["nb_esm"], out["nb_esm_maxsim"] = _topk_mean(s, R, k, weighted=False)
    a = G.similarity("apms", tc, rc)
    A = np.nan_to_num(a)
    cnt = A.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["nb_apms"] = np.where(cnt > 0, (A @ R) / np.maximum(cnt, 1), np.nan)
    out["apms_frac_partners_train"] = np.where(np.isnan(a).all(1), np.nan, cnt / max(v.n, 1))
    return out


def propensity_features(v: ScreenView, params: dict) -> pd.DataFrame:
    """THE feature function: one row per target of v (same order as v.target_orfs)."""
    from ysl.features import genes as G
    k = int(params.get("prop_k", PROP_DEFAULTS["prop_k"]))
    R = row_logit_rate(v)
    f = column_features(v, v.target_col, R, k)
    f |= external_neighbour_features(v, R, k)
    df = pd.DataFrame(f, index=pd.Index(v.target_orfs, name="orf"))
    df["on_array"] = (v.target_col >= 0).astype(float)
    gm = G.gene_matrix(v.target_orfs_current, current=True)
    gm.index = df.index
    return pd.concat([df, gm], axis=1)


def feature_columns(df: pd.DataFrame, fset: str) -> list[str]:
    spec = FEATURE_SETS[fset]
    if isinstance(spec, list):
        return list(spec)
    gene = [c for c in df.columns if c not in COL_FEATURES + EXT_NB_FEATURES + ["on_array"]]
    cols = []
    if "GENE" in spec:
        cols += gene
    if "EXT_NB" in spec:
        cols += EXT_NB_FEATURES
    if "COL" in spec:
        cols += COL_FEATURES + ["on_array"]
    return cols


def prop_folds(v: ScreenView, ctx: MethodContext, K: int) -> np.ndarray:
    return (ctx.rng("prop_folds").permutation(v.n) % K).astype(np.int64)


def crossfit_propensity_features(v: ScreenView, ctx: MethodContext) -> pd.DataFrame:
    """Features of v's TRAINING queries, each computed from v.without_queries(its fold) (rows in v.train_orfs order)."""
    K = int(ctx.get("prop_folds", PROP_DEFAULTS["prop_folds"]))
    folds = prop_folds(v, ctx, K)
    parts = []
    for kf in range(K):
        rows = np.flatnonzero(folds == kf)
        fv = v.without_queries(v.train_orfs[rows])
        parts.append(propensity_features(fv, ctx.params))
    return pd.concat(parts).reindex(v.train_orfs)


# ----------------------------------------------------------------------------- regressors
class _Linear:
    """Median imputation (+ missingness indicators for columns with NaN in training) -> standardise -> RidgeCV."""

    def __init__(self, kind: str = "ridge", seed: int = 0):
        self.kind, self.seed = kind, seed

    def _design(self, X: np.ndarray) -> np.ndarray:
        miss = np.isnan(X[:, self.nan_cols]).astype(np.float64)
        Xi = np.where(np.isnan(X), self.med[None, :], X)
        return np.hstack([Xi, miss])

    def fit(self, X: np.ndarray, y: np.ndarray):
        from sklearn.linear_model import ElasticNetCV, RidgeCV
        from sklearn.preprocessing import StandardScaler
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                self.med = np.nan_to_num(np.nanmedian(X, 0))
        self.nan_cols = np.flatnonzero(np.isnan(X).any(0))
        D = self._design(X)
        self.keep = D.std(0) > 1e-12
        self.sc = StandardScaler().fit(D[:, self.keep])
        Z = self.sc.transform(D[:, self.keep])
        if self.kind == "ridge":
            self.m = RidgeCV(alphas=np.logspace(-1, 4, 26)).fit(Z, y)
        else:
            from sklearn.model_selection import KFold
            self.m = ElasticNetCV(l1_ratio=[0.2, 0.5, 0.9], alphas=40, max_iter=5000,
                                  cv=KFold(5, shuffle=True, random_state=self.seed)).fit(Z, y)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.m.predict(self.sc.transform(self._design(X)[:, self.keep]))


class _HGB:
    def __init__(self, seed: int = 0, max_iter: int = 200):
        from sklearn.ensemble import HistGradientBoostingRegressor
        self.m = HistGradientBoostingRegressor(max_iter=max_iter, learning_rate=0.03, max_leaf_nodes=8,
                                               min_samples_leaf=20, l2_regularization=1.0, early_stopping=False,
                                               random_state=seed)

    def fit(self, X, y):
        self.m.fit(X, y)
        return self

    def predict(self, X):
        return self.m.predict(X)


def make_regressor(kind: str, seed: int, params: dict):
    if kind in ("ridge", "enet"):
        return _Linear(kind, seed)
    if kind == "hgb":
        return _HGB(seed, int(params.get("prop_hgb_iter", PROP_DEFAULTS["prop_hgb_iter"])))
    raise ValueError(kind)


def fit_predict_propensity(Xtr: pd.DataFrame, ytr: np.ndarray, Xte: pd.DataFrame, model: str, fset: str,
                           seed: int = 0, params: dict | None = None) -> np.ndarray:
    """Fit a propensity regressor on (Xtr, ytr) and predict Xte. model in
      ridge / enet / hgb                 one model on feature set `fset` for all queries
      split_<ridge|hgb>                  on-array targets: model on `fset` fitted on on-array training queries;
                                         off-array targets: model on 'external_nb' fitted on ALL training queries
      nb_esm / nb_coexp / col_logit ...  a single feature used as the prediction (diagnostic)"""
    params = params or {}
    if model in Xtr.columns:
        return Xte[model].to_numpy(dtype=np.float64)
    if model.startswith("split_"):
        base = model[len("split_"):]
        out = np.empty(len(Xte))
        on_tr = Xtr["on_array"].to_numpy() > 0
        on_te = Xte["on_array"].to_numpy() > 0
        cols = feature_columns(Xtr, fset)
        cold = feature_columns(Xtr, "external_nb")
        if on_te.any() and on_tr.sum() < MIN_SPLIT_TRAIN:      # too few on-array training queries: one model
            m = make_regressor(base, seed, params).fit(Xtr[cols].to_numpy(np.float64), ytr)
            out[on_te] = m.predict(Xte.loc[on_te, cols].to_numpy(np.float64))
        elif on_te.any():
            m = make_regressor(base, seed, params).fit(Xtr.loc[on_tr, cols].to_numpy(np.float64), ytr[on_tr])
            out[on_te] = m.predict(Xte.loc[on_te, cols].to_numpy(np.float64))
        if (~on_te).any():
            m = make_regressor(base, seed, params).fit(Xtr[cold].to_numpy(np.float64), ytr)
            out[~on_te] = m.predict(Xte.loc[~on_te, cold].to_numpy(np.float64))
        return out
    cols = feature_columns(Xtr, fset)
    m = make_regressor(model, seed, params).fit(Xtr[cols].to_numpy(np.float64), ytr)
    return m.predict(Xte[cols].to_numpy(np.float64))


def matched_target_view(v: ScreenView, drop_orfs) -> ScreenView:
    """A view with training rows T \\ drop_orfs and the SAME targets as v (built by the one view constructor from v's
    own training rows; targets are names only). Used so that test features are computed from as many screened
    queries as the cross-fitted training features (matched-view, DESIGN §6) — removes count/extreme-value drift."""
    keep = v.train_orfs[~np.isin(v.train_orfs, np.asarray(drop_orfs, dtype=str))]
    return ScreenView._construct(_ViewSource(v), keep, v.target_orfs)


def propensity_test_features(v: ScreenView, ctx: MethodContext) -> list[pd.DataFrame]:
    """Test-time features of v's targets: one frame from the full view ('full'), or one per matched view
    v minus training fold k, k = 0..K-1 ('matched', default; predictions are averaged over the K frames)."""
    if ctx.get("prop_test_view", PROP_DEFAULTS["prop_test_view"]) == "full":
        return [propensity_features(v, ctx.params)]
    K = int(ctx.get("prop_folds", PROP_DEFAULTS["prop_folds"]))
    folds = prop_folds(v, ctx, K)
    return [propensity_features(matched_target_view(v, v.train_orfs[folds == k]), ctx.params) for k in range(K)]


def predict_frames(Xtr: pd.DataFrame, ytr: np.ndarray, frames: list[pd.DataFrame], model: str, fset: str,
                   seed: int = 0, params: dict | None = None) -> np.ndarray:
    """Fit once, predict every test frame, average."""
    n = len(frames[0])
    pred = fit_predict_propensity(Xtr, ytr, pd.concat(frames), model, fset, seed, params)
    return pred.reshape(len(frames), n).mean(0)


def propensity_predict(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    """Predicted logit query rate for v's targets; the sub-model is trained on v's training queries only, with
    cross-fitted features (nested inside any outer cross-fit, because it reads only `v`)."""
    Xtr = crossfit_propensity_features(v, ctx)
    ytr = row_logit_rate(v)
    frames = propensity_test_features(v, ctx)
    seed = int(ctx.rng("prop_model").integers(2**31 - 1))
    return predict_frames(Xtr, ytr, frames, ctx.get("prop_model", PROP_DEFAULTS["prop_model"]),
                          ctx.get("prop_fset", "all"), seed, ctx.params)


# ----------------------------------------------------------------------------- stacked LR + propensity
PROP_COLS = ("prop_logit", "prop_logit_off")


def _prop_design(prop: np.ndarray, cols: np.ndarray, n_a: int, rows=slice(None)) -> np.ndarray:
    p = np.asarray(prop, dtype=np.float64)[rows]
    off = (np.asarray(cols)[rows] < 0).astype(np.float64)
    P = np.stack([p, p * off], axis=1)                      # [n_rows x 2]
    return np.repeat(P, n_a, axis=0).astype(np.float32)     # row-major like design_matrix's ravel


def stack_prop_features(v: ScreenView, ctx: MethodContext, names) -> dict:
    f = B.pair_features(v, ctx.params, B._knn_needed(names))
    f["prop"] = propensity_predict(v, ctx)
    return f


def crossfit_training_set_prop(v: ScreenView, ctx: MethodContext, names) -> tuple[np.ndarray, np.ndarray]:
    """As baselines.crossfit_training_set, plus the propensity columns from propensity_predict(fold view)."""
    folds = ctx.folds(v)
    knn = B._knn_needed(names)
    loo = B.loo_additive(v, ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)) if ctx.get("additive", "loo") == "loo" \
        else None
    X, ys = [], []
    for k in range(ctx.n_folds):
        rows = np.flatnonzero(folds == k)
        fv = v.without_queries(v.train_orfs[rows])
        f = B.fold_features(fv, ctx.params, knn)
        if loo is not None:
            f = B._loo_override(f, loo, rows)
        prop = propensity_predict(fv, ctx)
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        D = np.hstack([B.design_matrix(f, names), _prop_design(prop, f["cols"], v.n_arrays)])
        X.append(D[keep])
        ys.append(y[keep])
    return np.concatenate(X), np.concatenate(ys).astype(np.int8)


def stack_lr_prop_predict(v: ScreenView, ctx: MethodContext, names=B.LR_FEATURES["stack_lr"]) -> np.ndarray:
    X, y = crossfit_training_set_prop(v, ctx, names)
    sc, clf = B.fit_lr(X, y, ctx)
    del X
    f = stack_prop_features(v, ctx, names)
    n_t, n_a = v.n_targets, v.n_arrays
    out = np.empty((n_t, n_a), dtype=np.float32)
    for i in range(0, n_t, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, n_t))
        D = np.hstack([B.design_matrix(f, names, sl), _prop_design(f["prop"], f["cols"], n_a, sl)])
        out[sl] = clf.decision_function(sc.transform(D)).reshape(-1, n_a)
    ctx.info["lr_coef"] = dict(zip(tuple(names) + PROP_COLS, np.round(clf.coef_[0], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    ctx.info["prop_pred"] = {"mean": float(np.mean(f["prop"])), "sd": float(np.std(f["prop"]))}
    return out


@register("stack_lr_prop", defaults=B._LR_DEFAULTS | B.KNN_DEFAULTS | PROP_DEFAULTS | {"prop_fset": "all"},
          logit_scale=True,
          description="stack_lr + predicted query propensity (logit; + off-array interaction) from the nested "
                      "cross-fitted propensity sub-model (src/ysl/models/propensity.py)")
def stack_lr_prop(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return stack_lr_prop_predict(v, ctx)
