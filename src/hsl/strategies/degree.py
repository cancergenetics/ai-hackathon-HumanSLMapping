"""Shared query-degree model q̂ (DESIGN v3.1 §6.2 H; review r2 item 14). Used by H and F.

Prior (n < protocol selector.h_switch_n, unfitted): score = pct(1 − SMF) + pct(abundance) for on-array genes and
  2·pct(abundance) for off-array genes (same [0, 2] scale; routing by data availability, not a feature). pct = mid-rank
  empirical CDF: SMF against every non-NaN array SMF of the view (the library), abundance against the pool.
Fitted (n ≥ switch): Gaussian ridge (intercept unpenalised, λ by GCV) on y_s = logit((k_s + 0.5)/(m_s + 1)) of the
  screened rows (k = strong SLs, m = measured cells, visible columns only). Inputs: logit of the shrunk array-side
  rate ĉ_x = (k_x + 100·π̄_S)/(n_x + 100) (column x over the screened rows), array SMF, abundance, AP-MS degree,
  coexpression degree, conservation, paralog flag. Off-array genes (no column, no SMF) are scored by an
  external-only submodel (same response, external inputs only) fitted on all screened rows.
Inputs are standardised on the training rows (SMF: on all non-NaN array SMF values of the view); remaining NaN -> 0
(= the training mean) after standardisation. No missingness indicator enters any model.
Route calibration (knob degree_model.offarray_calibration, default 'quantile'): the two routes have different score
  distributions (e.g. 2·pct(abundance) has twice the variance of pct(1 − SMF) + pct(abundance); a shrunk full model
  vs a less-shrunk external submodel), so un-calibrated routing makes off-array status a de-facto score term (DESIGN
  §1; off-array yeast genes are truly richer in SLs). 'quantile': an off-array gene's percentile among the ON-array
  genes of the scored set under its own route (abundance for the prior, the external submodel for the fitted model)
  is mapped to the same quantile of the on-array genes' main-route scores. 'none' = scores as computed.
"""
from __future__ import annotations

import numpy as np

from hsl.strategies import base


# ============================================================================ GI statistics from the view
def row_counts(view) -> tuple[np.ndarray, np.ndarray]:
    """(k_s, m_s) per screened row over the visible columns."""
    Y = np.asarray(view.Y)[:, np.asarray(view.visible_cols)]
    return (Y == 1).sum(1).astype(np.float64), (Y >= 0).sum(1).astype(np.float64)


def col_counts(view) -> tuple[np.ndarray, np.ndarray, float]:
    """(k_x, n_x) per array column over the screened rows (0 for masked columns) and π̄_S (SL rate of the visible
    screened block)."""
    Y = np.asarray(view.Y)
    vis = np.zeros(Y.shape[1], dtype=bool)
    vis[np.asarray(view.visible_cols)] = True
    k = np.where(vis, (Y == 1).sum(0), 0).astype(np.float64)
    n = np.where(vis, (Y >= 0).sum(0), 0).astype(np.float64)
    pi = float(k.sum() / n.sum()) if n.sum() > 0 else float("nan")
    return k, n, pi


def shrunk_col_rate(view) -> tuple[np.ndarray, float]:
    """ĉ_x for every array column and π̄_S."""
    n0 = float(base.knobs()["spec"]["shrink_n0"])
    k, n, pi = col_counts(view)
    return (k + n0 * pi) / (n + n0), pi


def row_response(view) -> np.ndarray:
    a, b = base.knobs()["spec"]["response_pseudo"]
    k, m = row_counts(view)
    p = (k + float(a)) / (m + float(b))
    return np.log(p / (1 - p))


def screened_rate(view) -> np.ndarray:
    """Observed (smoothed) SL rate of each screened row: q̂_s for screened genes."""
    a, b = base.knobs()["spec"]["response_pseudo"]
    k, m = row_counts(view)
    return (k + float(a)) / (m + float(b))


# ============================================================================ ridge with GCV
class Ridge:
    """Ridge with unpenalised intercept; λ chosen by GCV on a grid λ = n · 10^x."""

    def __init__(self, grid: dict | None = None):
        g = grid or base.knobs()["degree_model"]["ridge_lambda_rel"]
        self.rel = 10.0 ** np.linspace(float(g["log10_min"]), float(g["log10_max"]), int(g["n"]))

    def fit(self, X: np.ndarray, y: np.ndarray) -> "Ridge":
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        n = len(y)
        self.mx, self.my = X.mean(0), float(y.mean())
        Xc, yc = X - self.mx, y - self.my
        if X.shape[1] == 0 or n < 3:
            self.beta, self.lam = np.zeros(X.shape[1]), float("nan")
            return self
        U, d, Vt = np.linalg.svd(Xc, full_matrices=False)
        Uy = U.T @ yc
        r0 = float(yc @ yc - Uy @ Uy)                  # residual outside the column space
        best = (np.inf, None)
        for lam in self.rel * n:
            f = d ** 2 / (d ** 2 + lam)
            rss = r0 + float((((1 - f) * Uy) ** 2).sum())
            df = float(f.sum()) + 1.0
            if n - df <= 0:
                continue
            gcv = n * rss / (n - df) ** 2
            if gcv < best[0]:
                best = (gcv, lam)
        lam = best[1] if best[1] is not None else self.rel[-1] * n
        self.lam = float(lam)
        self.beta = Vt.T @ (d / (d ** 2 + lam) * Uy)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.my + (np.asarray(X, dtype=np.float64) - self.mx) @ self.beta


def _standardise(train: np.ndarray, *others, ref: np.ndarray | None = None):
    R = train if ref is None else ref
    with np.errstate(invalid="ignore"):
        mu = np.nanmean(R, 0) if len(R) else np.zeros(R.shape[1])
        sd = np.nanstd(R, 0) if len(R) else np.ones(R.shape[1])
    mu = np.nan_to_num(mu)
    sd = np.where(np.nan_to_num(sd) > 0, sd, 1.0)
    out = [np.nan_to_num((np.asarray(a, dtype=np.float64) - mu) / sd) for a in (train,) + others]
    return out


# ============================================================================ external inputs
def external_inputs(feats) -> np.ndarray:
    """[n x ≤ 5]: abundance, AP-MS degree, coexpression degree, conservation, paralog flag (NaN kept). A feature
    whose config entry is null is omitted (e.g. the HAP1-minimal scenario drops AP-MS; abundance is required)."""
    f = base.knobs()["degree_model"]["features"]
    cols = [base.abundance(feats)]
    for key in ("apms_degree", "coexp_degree", "conservation"):
        if f.get(key):
            cols.append(feats[f[key]].to_numpy(np.float64))
    if f.get("paralog_count"):
        par = feats[f["paralog_count"]].to_numpy(np.float64)
        cols.append(np.where(np.isfinite(par), (par > 0).astype(np.float64), np.nan))
    return np.column_stack(cols)


# ============================================================================ prior
def prior_scores(view, genes, feats_pool, pool) -> np.ndarray:
    """Unfitted prior for `genes` (subset of `pool`; feats_pool = gene_features(pool))."""
    genes = np.asarray(genes, dtype=str)
    pos = np.searchsorted(pool, genes)
    ab_pool = base.abundance(feats_pool)
    pa = base.midrank_pct(ab_pool[pos], ab_pool)
    pa = np.where(np.isfinite(pa), pa, 0.5)
    col = _cols_of(view, genes)
    smf_all = np.asarray(view.array_smf, dtype=np.float64)
    on = col >= 0
    smf = np.full(len(genes), np.nan)
    smf[on] = smf_all[col[on]]
    ps = 1.0 - base.midrank_pct(smf, smf_all)
    ps = np.where(np.isfinite(ps), ps, 0.5)
    out = np.where(on, ps + pa, 2.0 * pa)
    if _calibrate_routes() and on.any() and (~on).any():
        out[~on] = quantile_map(out[on], ab_pool[pos][on], ab_pool[pos][~on])
    return out


def _calibrate_routes() -> bool:
    mode = str(base.knobs()["degree_model"].get("offarray_calibration", "quantile"))
    if mode not in ("quantile", "none"):
        raise ValueError(f"offarray_calibration {mode!r}")
    return mode == "quantile"


def quantile_map(main_on: np.ndarray, route_on: np.ndarray, route_off: np.ndarray) -> np.ndarray:
    """Off-array scores on the on-array scale: percentile of route_off within route_on (mid-rank; NaN -> median)
    mapped to the same quantile of main_on."""
    p = base.midrank_pct(route_off, route_on)
    p = np.clip(np.where(np.isfinite(p), p, 0.5), 0.0, 1.0)
    return np.quantile(np.asarray(main_on, dtype=np.float64), p)


def _cols_of(view, genes) -> np.ndarray:
    """Array column of each gene (candidates or screened genes of the view; −1 off-array)."""
    m = dict(zip(np.asarray(view.candidates, dtype=str), np.asarray(view.candidate_col)))
    m.update(zip(np.asarray(view.screened_orfs, dtype=str), np.asarray(view.screened_col)))
    return np.array([m[g] for g in np.asarray(genes, dtype=str)], dtype=np.int64)


# ============================================================================ fitted model
def fit_predict(view, genes, feats_pool, pool) -> tuple[np.ndarray, dict]:
    """Predicted logit query-side SL rate for `genes` (⊂ pool) from the screened rows of `view`."""
    genes = np.asarray(genes, dtype=str)
    S = np.asarray(view.screened_orfs, dtype=str)
    y = row_response(view)
    Xe_all = external_inputs(feats_pool)
    Xe_S = Xe_all[np.searchsorted(pool, S)]
    Xe_G = Xe_all[np.searchsorted(pool, genes)]
    c_hat, pi = shrunk_col_rate(view)
    lc = np.log(c_hat / (1 - c_hat))
    smf_all = np.asarray(view.array_smf, dtype=np.float64)
    col_S, col_G = np.asarray(view.screened_col), _cols_of(view, genes)
    on_S, on_G = col_S >= 0, col_G >= 0
    out = np.empty(len(genes))
    info = {"n_train": int(len(S)), "n_train_on": int(on_S.sum()), "pi_S": pi}
    # external-only submodel (all screened rows): off-array genes (and everything if too few on-array rows)
    a, b = _standardise(Xe_S, Xe_G)
    ext = Ridge().fit(a, y)
    out[:] = ext.predict(b)
    info["lam_ext"] = ext.lam
    if on_S.sum() >= int(base.knobs()["degree_model"]["min_train_on_array"]):
        def full(cols, Xe):
            return np.column_stack([lc[cols], smf_all[cols], Xe])
        tr = full(col_S[on_S], Xe_S[on_S])
        te = full(col_G[on_G], Xe_G[on_G])
        a, b = _standardise(tr, te)
        # SMF standardised against the library (all non-NaN array SMF values of the view)
        s_ref = smf_all[np.isfinite(smf_all)]
        mu, sd = (float(s_ref.mean()), float(s_ref.std()) or 1.0) if len(s_ref) else (0.0, 1.0)
        a[:, 1] = np.nan_to_num((tr[:, 1] - mu) / sd)
        b[:, 1] = np.nan_to_num((te[:, 1] - mu) / sd)
        m = Ridge().fit(a, y[on_S])
        ext_pred = out.copy()
        out[on_G] = m.predict(b)
        info["lam_full"] = m.lam
        info["coef_full"] = np.round(m.beta, 4).tolist()
        if _calibrate_routes() and on_G.any() and (~on_G).any():
            out[~on_G] = quantile_map(out[on_G], ext_pred[on_G], ext_pred[~on_G])
    return out, info


def sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.asarray(z, dtype=np.float64)))
