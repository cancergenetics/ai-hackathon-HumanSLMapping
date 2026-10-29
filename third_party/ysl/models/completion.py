"""Within-row ranking methods, rungs 1-2 (DESIGN.md §6): fold-in (B6), two-step Kronecker KRR, tied symmetric low-rank
completion, 2-hop partner similarity (B4b), refined mixed kNN, and the cross-fitted stack `stack_lr_wr`.

Every method reads only the ScreenView. Notation: T = screened (training) queries (rows of the view), T_on = screened
queries that also have an array strain, A = arrays, x = held-out target, c_g = gene g's ARRAY-SIDE column over T
(entries E(t, g), t in T; for a target x this is the only GI information besides the reciprocal, which is part of
c_x: E(a, x) for a in T). Off-array targets (no column) get NaN in feature form; registered single methods fill them
with the method's own mean row over on-array targets (stack_lr_wr handles NaN via an availability flag).

Signals (one per entry of the training rows; NaN where unobserved):
  negeps : -clip(eps, -1, 0.5)        label : Y in {0, 1}        negsp : -clip(signed -log10 P, +-12)
Centring ('row' = subtract each screened query's mean over arrays; 'col' = each array's mean over T; 'both' = row
then col), then missing -> 0.

Methods (all scores for x on the array are functions of c_x and the training rows only):
  B6 fold-in      ridge (dual / kernel form) from a gene's input profile z_g to its query row y_g (over all A),
                  fitted on screened genes; input modes
                    col   : z_g = c_g over T (training genes: T_on)            — test-time distribution exactly
                    mixed : z_t = t's ROW at T_on columns, z_x = c_x at T_on rows (symmetry assumption)
                    sym   : z_t = mean(row, column) at T_on coordinates (all screened genes are training genes)
                  For mixed/sym the input coordinate q equals the target (t, q) itself, so the coordinate of array
                  column a = col(q) is EXCLUDED when predicting column a (exact, per-column Sherman-Morrison downdate).
                  Optional reduced rank: fitted values projected on their top-r right singular vectors (RRR).
  two-step KRR    the fold-in prediction smoothed along arrays with an array kernel K_A = C^T C (C = column-normalised
                  array profiles over T):  P K_A (K_A + l2 I)^-1 = (P C^T)(C C^T + l2 I)^-1 C.  Query kernel = the
                  fold-in kernel (linear or cosine).
  tied low-rank   symmetric gene x gene matrix over A plus off-array screened queries, observed where one gene is a
                  screened query (both orientations; reciprocal pairs averaged), after a two-way additive fit
                  (mu + b_query + b_array). Hard-impute with a rank-r symmetric (indefinite) eigen-approximation and
                  eigenvalue soft-thresholding (randomised subspace iteration, warm-started, seeded). A held-out gene's
                  query factors are tied to its array factors (it is one gene); transductive over the held-out
                  genes' known columns only.
  B4b 2-hop       pp : sum_{t,t'} w(t, x) S(t, t') w(t', a) — x's and a's partners among T (from their columns) are
                       similar screened queries; S = Pearson of screened ROWS over all arrays (diag 0)
                  pg : sum_t w(t, x) sim(t, a) + sum_t w(t, a) sim(t, x) — a resembles x's partners and vice versa
                       ("between-pathway"); sim = 'mixed' Pearson (t's row at T_on vs g's column at T_on)
                  w = positive part of the centred signal, columns L1-normalised (partner distribution).
  kNN refined     mixed kNN (baselines) averaging neighbours' labels (or clipped -eps), array-centred
                  (sum_t w (V(t,a) - mean_a) / (sum_t w obs + m)), k = round(k_frac * n_T).

Cross-fitted stack `stack_lr_wr`: stack_lr's features (baselines.pair_features + exact-LOO additive stats) plus the
within-row scores above; training features for fold k from view.without_queries(fold k) (K = 10), test features from
the full view — the same feature function (`wr_features`) for both.
"""
from __future__ import annotations

import numpy as np

from ysl.models import baselines as B
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView


# ----------------------------------------------------------------------------- signals
def signal(v: ScreenView, kind: str) -> np.ndarray:
    """[n_T x n_A] float64, NaN where unobserved."""
    if kind == "negeps":
        return -np.clip(v.E.astype(np.float64), *B.EPS_CLIP)
    if kind == "label":
        return v.Y.astype(np.float64)
    if kind == "negsp":
        return -np.clip(v.SP.astype(np.float64), -B.SP_CLIP, B.SP_CLIP)
    raise ValueError(kind)


def centre(R: np.ndarray, how: str) -> np.ndarray:
    """Centre a NaN-masked matrix ('none'|'row'|'col'|'both'), then missing -> 0."""
    if how not in ("none", "row", "col", "both"):
        raise ValueError(how)
    R = np.array(R, dtype=np.float64, copy=True)
    if how in ("row", "both"):
        R -= _nanmean(R, 1)[:, None]
    if how in ("col", "both"):
        R -= _nanmean(R, 0)[None, :]
    return np.nan_to_num(R, nan=0.0)


def _nanmean(R: np.ndarray, axis: int) -> np.ndarray:
    n = (~np.isnan(R)).sum(axis)
    return np.where(n > 0, np.nansum(R, axis=axis) / np.maximum(n, 1), 0.0)


def col_mean(R: np.ndarray) -> np.ndarray:
    n = (~np.isnan(R)).sum(0)
    s = np.nansum(R, axis=0)
    return np.where(n > 0, s / np.maximum(n, 1), 0.0)


def _on(v: ScreenView) -> np.ndarray:
    return np.flatnonzero(v.target_col >= 0)


def _full(v: ScreenView, on: np.ndarray, P: np.ndarray) -> np.ndarray:
    S = np.full((v.n_targets, v.n_arrays), np.nan, dtype=np.float32)
    S[on] = P
    return S


# ----------------------------------------------------------------------------- B6 fold-in (+ two-step KRR)
# Recipes selected on S1 dev seed 0 (scripts/41_explore.py; results/m4_completion/explore_seed0.json), then frozen.
FOLDIN_DEFAULTS = {"inp": "sym", "in_sig": "negeps", "in_centre": "row", "norm": False, "target": "label",
                   "lam": 10.0, "rank": 50}
KRR2_DEFAULTS = FOLDIN_DEFAULTS | {"lam2": 0.03}


def foldin_inputs(v: ScreenView, inp: str, in_sig: str, in_centre: str):
    """(Z [n_train_genes x p], Zx [n_on_targets x p], train_rows (indices into v's rows), diag_cols) where
    diag_cols[j] = array column whose target equals input coordinate j for training genes (-1: none)."""
    raw = signal(v, in_sig)
    R = centre(raw, in_centre)
    on = _on(v)
    xc = v.target_col[on]
    t_on = np.flatnonzero(v.train_col >= 0)
    qcols = v.train_col[t_on]
    if inp == "col":
        Z = R[:, qcols].T                                  # c_t over T, for t in T_on   [n_on x n_T]
        Zx = R[:, xc].T
        return Z, Zx, t_on, np.full(Z.shape[1], -1)
    if inp == "mixed":
        Z = R[:, qcols]                                    # t's row at T_on columns     [n_T x n_on]
        Zx = R[t_on][:, xc].T                              # c_x at T_on rows
        return Z, Zx, np.arange(v.n), qcols
    if inp == "sym":
        Rm = np.where(np.isnan(raw), np.nan, R)            # centred, missing kept as NaN for averaging
        rowm = Rm[:, qcols]                                # E(t, q): t's row at T_on columns
        colv = np.full_like(rowm, np.nan)
        colv[t_on] = Rm[t_on][:, v.train_col[t_on]].T      # E(q, t): t's column at T_on rows (t on the array)
        Z = np.nan_to_num(np.nanmean(np.stack([rowm, colv]), axis=0))
        Zx = R[t_on][:, xc].T
        return Z, Zx, np.arange(v.n), qcols
    raise ValueError(inp)


def _normalise_rows(Z: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(Z, axis=1, keepdims=True)
    return Z / np.where(n > 0, n, 1.0)


def foldin_predict(v: ScreenView, params: dict, array_kernel: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fold-in residual predictions for on-array targets: (P [n_on x n_A], on, mu [n_A]); score = mu + P."""
    p = FOLDIN_DEFAULTS | params
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Z, Zx, rows, dcols = foldin_inputs(v, p["inp"], p["in_sig"], p["in_centre"])
    if p["norm"]:
        Z, Zx = _normalise_rows(Z), _normalise_rows(Zx)
    Tsig = signal(v, p["target"])
    mu = col_mean(Tsig)
    Y = np.nan_to_num(Tsig[rows] - mu[None, :])           # array-centred targets of the training genes
    K = Z @ Z.T
    lam = float(p["lam"]) * np.trace(K) / max(len(K), 1)
    Minv = np.linalg.inv(K + lam * np.eye(len(K)))
    alpha = Minv @ Y                                        # [n x n_A]
    Kx = Zx @ Z.T
    P = Kx @ alpha
    excl = np.flatnonzero(dcols >= 0)
    if len(excl):  # exact removal of input coordinate j when predicting array column dcols[j]
        Zj = Z[:, excl]                                     # [n x m]
        a_j = dcols[excl]
        Mz = Minv @ Zj
        d = 1.0 - np.einsum("ij,ij->j", Zj, Mz)
        s = np.einsum("ij,ij->j", Mz, Y[:, a_j])
        al = alpha[:, a_j] + Mz * (s / d)[None, :]
        u = np.einsum("ij,ij->j", Zj, al)
        P[:, a_j] = Kx @ al - Zx[:, excl] * u[None, :]
    if int(p["rank"]) > 0:  # reduced-rank: project on the top-r right singular vectors of the fitted values
        F = K @ alpha
        _, _, Vt = np.linalg.svd(F, full_matrices=False)
        Vr = Vt[: int(p["rank"])]
        P = (P @ Vr.T) @ Vr
    if array_kernel:
        C = centre(signal(v, p["in_sig"]), "row")
        C = C / np.where((n := np.linalg.norm(C, axis=0, keepdims=True)) > 0, n, 1.0)   # columns unit norm
        G = C @ C.T
        l2 = float(p["lam2"]) * np.trace(G) / len(G)
        P = (P @ C.T) @ np.linalg.solve(G + l2 * np.eye(len(G)), C)
    return P, _on(v), mu


# ----------------------------------------------------------------------------- tied symmetric low-rank completion
TIED_DEFAULTS = {"sig": "label", "rank": 10, "tau": 0.35, "n_iter": 30, "n_power": 2}


def two_way_fit(R: np.ndarray, n_pass: int = 5):
    """Additive fit mu + bq[t] + ba[a] on the observed entries (alternating means)."""
    obs = ~np.isnan(R)
    R0 = np.nan_to_num(R)
    mu = R0.sum() / obs.sum()
    bq = np.zeros(R.shape[0])
    ba = np.zeros(R.shape[1])
    no_r, no_c = np.maximum(obs.sum(1), 1), np.maximum(obs.sum(0), 1)
    for _ in range(n_pass):
        bq = ((R0 - mu - ba[None, :]) * obs).sum(1) / no_r
        ba = ((R0 - mu - bq[:, None]) * obs).sum(0) / no_c
    return mu, bq, ba


def symmetric_matrix(v: ScreenView, Res: np.ndarray):
    """Gene x gene symmetric matrix (genes = arrays, then off-array screened queries) and its observed mask."""
    nA = v.n_arrays
    off = np.flatnonzero(v.train_col < 0)
    g = v.train_col.copy()
    g[off] = nA + np.arange(len(off))
    nG = nA + len(off)
    obs = ~np.isnan(Res)
    R0 = np.nan_to_num(Res).astype(np.float32)
    S = np.zeros((nG, nG), dtype=np.float32)
    W = np.zeros((nG, nG), dtype=np.float32)
    S[g[:, None], np.arange(nA)[None, :]] += R0
    W[g[:, None], np.arange(nA)[None, :]] += obs
    S[np.arange(nA)[:, None], g[None, :]] += R0.T
    W[np.arange(nA)[:, None], g[None, :]] += obs.T
    # g is injective, so neither assignment has duplicate indices. For t, t' both screened and on the array the cell
    # (col t, col t') receives E(t, t') from the first and E(t', t) from the second assignment -> averaged by W.
    np.divide(S, W, out=S, where=W > 0)
    np.fill_diagonal(W, 0)
    return S, W > 0, g


def _top_eig_sym(X: np.ndarray, Q: np.ndarray, rank: int, n_power: int):
    """Randomised subspace iteration for a symmetric (indefinite) X: top-|eigenvalue| pairs."""
    for _ in range(n_power):
        Q, _ = np.linalg.qr(X @ Q)
    Bm = Q.T @ X @ Q
    w, U = np.linalg.eigh((Bm + Bm.T) / 2)
    idx = np.argsort(-np.abs(w), kind="stable")[:rank]
    return w[idx], Q @ U[:, idx], Q


def tied_lowrank_predict(v: ScreenView, params: dict, rng: np.random.Generator):
    """Residual low-rank predictions for on-array targets: (P [n_on x n_A], on, ba [n_A], extra)."""
    p = TIED_DEFAULTS | params
    R = signal(v, p["sig"])
    mu, bq, ba = two_way_fit(R)
    Res = R - mu - bq[:, None] - ba[None, :]
    S, M, g = symmetric_matrix(v, Res)
    r, over = int(p["rank"]), 10
    Q = np.linalg.qr(rng.standard_normal((S.shape[0], r + over)).astype(np.float32))[0]
    L = np.zeros_like(S)
    for _ in range(int(p["n_iter"])):
        X = np.where(M, S, L)
        w, U, Q = _top_eig_sym(X, Q, r, int(p["n_power"]))
        tau = float(p["tau"]) * np.abs(w).max()
        w = np.sign(w) * np.maximum(np.abs(w) - tau, 0.0)
        L = ((U * w[None, :]) @ U.T).astype(np.float32)
    on = _on(v)
    P = L[v.target_col[on]][:, : v.n_arrays].astype(np.float64)
    return P, on, ba, {"mu": mu, "bq": bq, "g": g}


# ----------------------------------------------------------------------------- B4b 2-hop partner similarity
TWOHOP_DEFAULTS = {"form": "pp", "sig": "negeps", "centre": "row", "normalise": True}


def _partners(R: np.ndarray) -> np.ndarray:
    W = np.maximum(R, 0.0)
    s = W.sum(0, keepdims=True)
    return W / np.where(s > 0, s, 1.0)


def twohop_predict(v: ScreenView, params: dict) -> tuple[np.ndarray, np.ndarray]:
    p = TWOHOP_DEFAULTS | params
    R = centre(signal(v, p["sig"]), p["centre"])
    W = _partners(R) if p["normalise"] else np.maximum(R, 0.0)     # [n_T x n_A] partner weights of each array
    on = _on(v)
    xc = v.target_col[on]
    if p["form"] == "pp":
        S = B.masked_pearson(v.E.T.astype(np.float64), v.E.T.astype(np.float64)).astype(np.float64)  # [n_T x n_T]
        np.fill_diagonal(S, 0.0)
        P = W[:, xc].T @ S @ W
    elif p["form"] == "pg":
        t_on = np.flatnonzero(v.train_col >= 0)
        sim = B.masked_pearson(v.E[t_on].astype(np.float64),
                               v.E[:, v.train_col[t_on]].T.astype(np.float64)).astype(np.float64)  # [n_A x n_T]
        P = W[:, xc].T @ sim.T + sim[xc] @ W
    elif p["form"] == "direct":
        C = R / np.where((n := np.linalg.norm(R, axis=0, keepdims=True)) > 0, n, 1.0)
        P = C[:, xc].T @ C
    else:
        raise ValueError(p["form"])
    return P, on


# ----------------------------------------------------------------------------- refined mixed kNN
KNNR_DEFAULTS = {"sig": "label", "k_frac": 0.3, "p": 1.0, "m": 0.05, "array_centre": True}


def knn_refined_predict(v: ScreenView, params: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = KNNR_DEFAULTS | params
    on = _on(v)
    sim, cand = B.query_side_similarity(v, v.target_col[on], "mixed")
    k = max(1, int(round(float(p["k_frac"]) * v.n)))
    Wt = B.topk_weights(sim, k, float(p["p"]))
    V = signal(v, p["sig"])[cand]
    mu = col_mean(signal(v, p["sig"]))
    O = (~np.isnan(V)).astype(np.float64)
    V0 = np.nan_to_num(V - mu[None, :]) if p["array_centre"] else np.nan_to_num(V)
    P = (Wt @ V0) / (Wt @ O + float(p["m"]))
    if not p["array_centre"]:
        P = P - mu[None, :]
    return P, on, mu


# ----------------------------------------------------------------------------- one feature function
WR_FEATURES = {  # name -> (kind, params); the within-row scores offered to the stack (residual scale)
    "foldin": ("foldin", {}),
    "foldin_col": ("foldin", {"inp": "col", "rank": 0}),
    "krr2": ("krr2", {}),
    "tied": ("tied", {}),
    "hop_pp": ("twohop", {"form": "pp"}),
    "hop_pg": ("twohop", {"form": "pg"}),
    "knn_lab": ("knnr", {}),
    "knn_eps": ("knnr", {"sig": "negeps", "k_frac": 0.08}),
}


def wr_score(v: ScreenView, kind: str, params: dict, rng_seed: int = 0) -> np.ndarray:
    """Residual within-row score [n_targets x n_A], NaN rows for off-array targets."""
    if kind == "foldin":
        P, on, _ = foldin_predict(v, FOLDIN_DEFAULTS | params)
    elif kind == "krr2":
        P, on, _ = foldin_predict(v, KRR2_DEFAULTS | params, array_kernel=True)
    elif kind == "tied":
        P, on, _, _ = tied_lowrank_predict(v, params, np.random.default_rng(rng_seed))
    elif kind == "twohop":
        P, on = twohop_predict(v, params)
    elif kind == "knnr":
        P, on, _ = knn_refined_predict(v, params)
    else:
        raise ValueError(kind)
    return _full(v, on, P)


def wr_features(v: ScreenView, params: dict, names: tuple[str, ...], rng_seed: int = 0) -> dict:
    """Within-row features of v's targets; params[f'{name}__{key}'] overrides a feature's settings."""
    out = {}
    for nm in names:
        kind, base = WR_FEATURES[nm]
        own = {k.split("__", 1)[1]: val for k, val in params.items() if k.startswith(nm + "__")}
        out[nm] = wr_score(v, kind, base | own, rng_seed)
    return out


# ----------------------------------------------------------------------------- registered single methods
def _fill_off_with_mean_row(S: np.ndarray) -> np.ndarray:
    on = ~np.isnan(S).all(1)
    if not on.any():
        raise ValueError("no on-array targets")
    mrow = S[on].mean(0)
    return np.where(np.isnan(S), mrow[None, :], S).astype(np.float32)


def _single(kind: str, add_mean: bool):
    def fn(v: ScreenView, ctx: MethodContext) -> np.ndarray:
        params = dict(ctx.params)
        if kind in ("foldin", "krr2"):
            P, on, mu = foldin_predict(v, (FOLDIN_DEFAULTS if kind == "foldin" else KRR2_DEFAULTS) | params,
                                       array_kernel=(kind == "krr2"))
            P = P + mu[None, :]
        elif kind == "tied":
            P, on, ba, ex = tied_lowrank_predict(v, params, ctx.rng("tied"))
            # standalone: add the array bias and a query level predicted from x's array bias (linear fit on T_on)
            t_on = np.flatnonzero(v.train_col >= 0)
            beta = np.polyfit(ba[v.train_col[t_on]], ex["bq"][t_on], 1)
            P = ex["mu"] + P + ba[None, :] + np.polyval(beta, ba[v.target_col[on]])[:, None]
        elif kind == "twohop":
            P, on = twohop_predict(v, params)
        elif kind == "knnr":
            P, on, mu = knn_refined_predict(v, params)
            P = P + mu[None, :]
        else:
            raise ValueError(kind)
        return _fill_off_with_mean_row(_full(v, on, P))
    return fn


register("B6_foldin", defaults=dict(FOLDIN_DEFAULTS),
         description="B6 fold-in: kernel ridge from a gene's array-side profile to its query row (off-array -> mean row)")(
    _single("foldin", True))
register("B6_krr2step", defaults=dict(KRR2_DEFAULTS),
         description="two-step Kronecker KRR: fold-in + array-kernel smoothing (array profiles over T)")(
    _single("krr2", True))
register("B6_tied_lowrank", defaults=dict(TIED_DEFAULTS),
         description="tied symmetric low-rank completion (both orientations, two-way biases, hard-impute)")(
    _single("tied", True))
register("B4b_2hop_pp", defaults=TWOHOP_DEFAULTS | {"form": "pp"},
         description="2-hop: x's partners vs a's partners among screened queries (screened-row Pearson)")(
    _single("twohop", False))
register("B4b_2hop_pg", defaults=TWOHOP_DEFAULTS | {"form": "pg"},
         description="2-hop between-pathway: a resembles x's partners and vice versa (mixed Pearson)")(
    _single("twohop", False))
register("B4r_knn_mixed_refined", defaults=dict(KNNR_DEFAULTS),
         description="mixed kNN, array-centred, k = k_frac * n_T (labels; -eps averaging was worse on dev)")(_single("knnr", True))


# ----------------------------------------------------------------------------- stack_lr_wr
STACK_BASE = B.LR_FEATURES["stack_lr"]
WR_STACK_NAMES = ("foldin", "tied", "knn_lab", "hop_pp", "hop_pg")  # overridden by params["wr"]


def _wr_design(fw: dict, names, rows=slice(None)) -> np.ndarray:
    cols = []
    for nm in names:
        X = np.asarray(fw[nm])[rows]
        cols.append(np.where(np.isnan(X), 0.0, X))
    return np.stack([np.asarray(c, dtype=np.float32).ravel() for c in cols], axis=1) if cols else None


def _standardise_wr(fw: dict, scale: dict | None) -> tuple[dict, dict]:
    """Scale each within-row feature by its SD over finite entries (computed per view; scale-free features keep
    the train/test feature distributions comparable across |T| = 0.9 n vs n)."""
    out, sc = {}, {}
    for k, X in fw.items():
        s = float(np.nanstd(X)) if scale is None else scale[k]
        sc[k] = s if s > 0 else 1.0
        out[k] = X / sc[k]
    return out, sc


def stack_wr_predict(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    names = tuple(ctx.get("wr", WR_STACK_NAMES))
    base = STACK_BASE
    knn = B._knn_needed(base)
    folds = ctx.folds(v)
    loo = B.loo_additive(v, ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR))
    seed = int(ctx.rng("tied").integers(2**31))
    Xs, ys = [], []
    for k in range(ctx.n_folds):
        rows = np.flatnonzero(folds == k)
        fv = v.without_queries(v.train_orfs[rows])
        f = B._loo_override(B.fold_features(fv, ctx.params, knn), loo, rows)
        fw, _ = _standardise_wr(wr_features(fv, ctx.params, names, seed), None)
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        Xb = B.design_matrix(f, base)
        Xw = _wr_design(fw, names)
        Xs.append((np.hstack([Xb, Xw]) if Xw is not None else Xb)[keep])
        ys.append(y[keep])
    X, y = np.concatenate(Xs), np.concatenate(ys).astype(np.int8)
    del Xs
    fitted = B.fit_lr(X, y, ctx)
    del X
    f = B.pair_features(v, ctx.params, knn)
    fw, _ = _standardise_wr(wr_features(v, ctx.params, names, seed), None)
    sc, clf = fitted
    n_t, n_a = v.n_targets, v.n_arrays
    S = np.empty((n_t, n_a), dtype=np.float32)
    for i in range(0, n_t, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, n_t))
        Xb = B.design_matrix(f, base, sl)
        Xw = _wr_design(fw, names, sl)
        Xt = np.hstack([Xb, Xw]) if Xw is not None else Xb
        S[sl] = clf.decision_function(sc.transform(Xt)).reshape(-1, n_a)
    ctx.info["lr_coef"] = dict(zip(base + names, np.round(clf.coef_[0], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    return S


register("stack_lr_wr", defaults=B._LR_DEFAULTS | B.KNN_DEFAULTS | {"wr": list(WR_STACK_NAMES)}, logit_scale=True,
         description="stack_lr + within-row completion scores (fold-in, KRR, tied low-rank, 2-hop, kNN-eps); K=10 "
                     "cross-fitted")(stack_wr_predict)
