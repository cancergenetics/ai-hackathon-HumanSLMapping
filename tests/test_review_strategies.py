"""Independent review of hsl.strategies (docs/reviews/strategies_code_review.md). Toy worlds only (no real data).

* Brute-force reimplementations (written from DESIGN v3.1 §6.2 / review r2 items 13–15, not from the code): ridge+GCV,
  the degree model q̂, F's objective / embeddings / greedy (fitted and cold), M's GP posterior variance greedy, C's
  facility-location greedy.
* Edge cases: n = 29/30 switch, all candidates off-array, pool exhaustion, batch remainder, static variants ignore
  screened-row content.
* Off-array routing (DESIGN §1, decision D5): in a NULL world genes are exchangeable, so routing a gene to another
  submodel because it lacks a column must not change its selection probability. Known failures are xfail(strict):
  they flip to XPASS (and fail the suite) once fixed, so the marker must then be removed.
"""
from __future__ import annotations

import hsl  # noqa: F401

import numpy as np
import pytest

from hsl import harness, strategies
from hsl.strategies import base, degree, kernels
from strategy_toys import ToyWorld, run_toy

CAP = int(base.protocol_selector()["cluster_cap_per_batch"])
SWITCH = int(base.protocol_selector()["h_switch_n"])


def _logit(p):
    return np.log(p / (1 - p))


def _view(seed, kind, n, **kw):
    w = ToyWorld(seed, kind, **kw)
    S = list(np.random.default_rng(seed).choice(w.queries, n, replace=False)) if n else []
    return w, w.view(S)


def _det_select(name, v, k=10, **kw):
    sel = strategies.make(name, deterministic=True, **kw)
    rng = np.random.default_rng(0)
    sel.reset(v.candidates, rng)
    return sel, sel.select(v, k, rng)


def _bf_greedy(gain_of, n, k, clusters=None):
    """Deterministic greedy: argmax of the exact marginal gain (ties -> lowest index), optional per-cluster cap."""
    picked, counts = [], {}
    for _ in range(k):
        best, bi = -np.inf, None
        for j in range(n):
            if j in picked:
                continue
            if clusters is not None and clusters[j] >= 0 and counts.get(int(clusters[j]), 0) >= CAP:
                continue
            g = gain_of(picked, j)
            tol = 1e-12 * max(1.0, abs(best)) if np.isfinite(best) else 0.0
            if bi is None or g > best + tol:
                best, bi = g, j
        picked.append(bi)
        if clusters is not None:
            counts[int(clusters[bi])] = counts.get(int(clusters[bi]), 0) + 1
    return picked


# ============================================================================ ridge + GCV
def test_ridge_gcv_matches_bruteforce_hat_matrix():
    rng = np.random.default_rng(0)
    n, p = 80, 6
    X = rng.normal(size=(n, p)) * rng.uniform(0.2, 3, p)
    y = X @ rng.normal(size=p) * 0.3 + rng.normal(size=n) + 2.0
    m = degree.Ridge().fit(X, y)
    g = base.knobs()["degree_model"]["ridge_lambda_rel"]
    best = (np.inf, None)
    Xc, yc = X - X.mean(0), y - y.mean()
    for lam in n * 10.0 ** np.linspace(g["log10_min"], g["log10_max"], g["n"]):
        H = np.ones((n, n)) / n + Xc @ np.linalg.solve(Xc.T @ Xc + lam * np.eye(p), Xc.T)   # intercept unpenalised
        r = y - H @ y
        gcv = n * (r @ r) / (n - np.trace(H)) ** 2
        if gcv < best[0]:
            best = (gcv, lam)
    lam = best[1]
    beta = np.linalg.solve(Xc.T @ Xc + lam * np.eye(p), Xc.T @ yc)
    assert np.isclose(m.lam, lam)
    Xt = rng.normal(size=(7, p))
    assert np.allclose(m.predict(Xt), y.mean() + (Xt - X.mean(0)) @ beta)


# ============================================================================ degree model (H, F)
def _std(train, *others, ref=None):
    R = train if ref is None else ref
    mu, sd = np.nanmean(R, 0), np.nanstd(R, 0)
    sd = np.where(sd > 0, sd, 1.0)
    return [np.nan_to_num((a - mu) / sd) for a in (train,) + others]


def _ridge_bf(X, y, Xt):
    g = base.knobs()["degree_model"]["ridge_lambda_rel"]
    n, p = X.shape
    Xc, yc = X - X.mean(0), y - y.mean()
    best = (np.inf, None)
    for lam in n * 10.0 ** np.linspace(g["log10_min"], g["log10_max"], g["n"]):
        H = np.ones((n, n)) / n + Xc @ np.linalg.solve(Xc.T @ Xc + lam * np.eye(p), Xc.T)
        r = y - H @ y
        gcv = n * (r @ r) / (n - np.trace(H)) ** 2
        if gcv < best[0]:
            best = (gcv, lam)
    beta = np.linalg.solve(Xc.T @ Xc + best[1] * np.eye(p), Xc.T @ yc)
    return y.mean() + (Xt - X.mean(0)) @ beta


def test_degree_model_matches_independent_reimplementation():
    w, v = _view(11, "hub", 70)
    f = base.knobs()["degree_model"]["features"]
    feats = w.feats
    cand, S = np.asarray(v.candidates), np.asarray(v.screened_orfs)
    # response logit((k + 0.5)/(m + 1)) over visible columns
    Y = v.Y[:, v.visible_cols]
    y = _logit(((Y == 1).sum(1) + 0.5) / ((Y >= 0).sum(1) + 1.0))
    # shrunk array-side rate
    kx, nx = (v.Y == 1).sum(0).astype(float), (v.Y >= 0).sum(0).astype(float)
    pi = kx.sum() / nx.sum()
    lc = _logit((kx + 100 * pi) / (nx + 100))
    ext_cols = [f["abundance"], f["apms_degree"], f["coexp_degree"], f["conservation"]]

    def ext(g):
        X = feats.loc[g, ext_cols].to_numpy(float)
        return np.column_stack([X, (feats.loc[g, f["paralog_count"]].to_numpy(float) > 0).astype(float)])
    Xs, Xg = ext(S), ext(cand)
    a, b = _std(Xs, Xg)
    out = _ridge_bf(a, y, b)
    on_s, on_g = v.screened_col >= 0, v.candidate_col >= 0
    smf = np.asarray(v.array_smf, float)

    def full(cols, Xe):
        return np.column_stack([lc[cols], smf[cols], Xe])
    tr, te = full(v.screened_col[on_s], Xs[on_s]), full(v.candidate_col[on_g], Xg[on_g])
    a2, b2 = _std(tr, te)
    a2[:, 1], b2[:, 1] = (tr[:, 1] - smf.mean()) / smf.std(), (te[:, 1] - smf.mean()) / smf.std()
    ext_pred = out.copy()
    out[on_g] = _ridge_bf(a2, y[on_s], b2)
    # quantile calibration of the off-array route (D5)
    r = np.sort(ext_pred[on_g])
    pct = (np.searchsorted(r, ext_pred[~on_g], "left") + np.searchsorted(r, ext_pred[~on_g], "right")) / (2 * len(r))
    out[~on_g] = np.quantile(out[on_g], pct)
    pool = base.pool_of(v)
    got, _ = degree.fit_predict(v, cand, v.gene_features(pool), pool)
    assert np.allclose(got, out, atol=1e-9)


# ============================================================================ F
def _P_cols(v):
    pool_cols = set(int(c) for c in np.concatenate([v.candidate_col, v.screened_col]) if c >= 0)
    ao = [c for c in v.visible_cols if int(c) not in pool_cols]
    return np.array(sorted(set(int(c) for c in v.candidate_col if c >= 0) | set(int(c) for c in ao)))


@pytest.mark.parametrize("r", [1, 4])
def test_F_fitted_greedy_matches_bruteforce_logdet(r):
    """Independent p̂, δ, P, M_x; exact log-det marginal gains (slogdet), deterministic greedy with the cluster cap."""
    w, v = _view(21, "hub", 50)
    sel, picks = _det_select("F", v, 10, r=r)
    pool = base.pool_of(v)
    uS, uC = sel.embeddings(v, sel._feats.get(v, pool), pool)
    lo, hi = base.knobs()["spec"]["p_clip"]
    Y = v.Y[:, v.visible_cols]
    q_s = ((Y == 1).sum(1) + 0.5) / ((Y >= 0).sum(1) + 1.0)
    q_c = degree.sigmoid(degree.fit_predict(v, v.candidates, v.gene_features(pool), pool)[0])
    kx, nx = (v.Y == 1).sum(0).astype(float), (v.Y >= 0).sum(0).astype(float)
    pi = kx.sum() / nx.sum()
    c = (kx + 100 * pi) / (nx + 100)
    P = _P_cols(v)
    assert np.array_equal(P, sel.target_cols(v))
    ps = np.clip(np.outer(q_s, c[P] / pi), lo, hi)
    pc = np.clip(np.outer(q_c, c[P] / pi), lo, hi)
    Ws, Wc = ps * (1 - ps), pc * (1 - pc)
    delta = Ws.mean()
    M0 = np.array([delta * np.eye(r) + (uS * Ws[:, x][:, None]).T @ uS for x in range(len(P))])

    def obj(extra):
        M = M0.copy()
        for j in extra:
            M = M + Wc[j][:, None, None] * np.outer(uC[j], uC[j])[None]
        return np.linalg.slogdet(M)[1].sum()
    cache = {}

    def gain(picked, j):
        key = tuple(picked)
        if key not in cache:
            cache[key] = obj(picked)
        return obj(picked + [j]) - cache[key]
    cl = v.clusters(v.candidates)
    bf = _bf_greedy(gain, len(v.candidates), 10, clusters=cl)
    assert picks == [str(v.candidates[i]) for i in bf]


def test_F_embeddings_independent_svd_map_krr_whitening():
    w, v = _view(22, "hub", 60)
    f = strategies.make("F", r=5)
    pool = base.pool_of(v)
    feats = f._feats.get(v, pool)
    uS, uC = f.embeddings(v, feats, pool)
    on_c = v.candidate_col >= 0
    # u = (1, ũ); ũ unit mean square over S ∪ candidates
    assert np.allclose(uS[:, 0], 1) and np.allclose(uC[:, 0], 1)
    assert np.allclose((np.vstack([uS, uC])[:, 1:] ** 2).mean(0), 1.0)
    # screened ũ spans the top-4 left singular space of the row-centred −clip(ε) block
    Ev = np.asarray(v.E, float)[:, v.visible_cols]
    from ysl.models import baselines as B
    R = -np.clip(Ev, *B.EPS_CLIP)
    R = R - np.nanmean(R, 1, keepdims=True)
    R = np.nan_to_num(R)
    Usvd = np.linalg.svd(R, full_matrices=False)[0][:, :4]
    proj = Usvd @ (Usvd.T @ uS[:, 1:])
    assert np.allclose(proj, uS[:, 1:], atol=1e-8)
    # on-array candidates: û = B v with B least squares on screened on-array genes (raw scale)
    col = {int(c): j for j, c in enumerate(v.visible_cols)}
    sp = np.array([col[c] for c in v.screened_col if c >= 0])
    U = uS[:, 1:]                                  # whitened; B is equivariant to the per-column scale
    Vs = R[:, sp].T @ U
    on_s = v.screened_col >= 0
    Bt = np.linalg.lstsq(Vs, U[on_s], rcond=None)[0]
    cp = np.array([col[c] for c in v.candidate_col[on_c]])
    pred = (R[:, cp].T @ U) @ Bt
    # compare directions (b_ridge = 1e-6 relative is negligible)
    assert np.allclose(pred, uC[on_c, 1:], rtol=1e-3, atol=1e-6 * np.abs(pred).max())
    # off-array candidates: KRR from the external kernel into the same space
    K = kernels.external_kernel(v, pool, feats)
    ps_, pc_ = np.searchsorted(pool, v.screened_orfs), np.searchsorted(pool, v.candidates[~on_c])
    Kss = K[np.ix_(ps_, ps_)]
    lam = f.krr_lam * np.trace(Kss) / len(ps_)
    krr = K[np.ix_(pc_, ps_)] @ np.linalg.solve(Kss + lam * np.eye(len(ps_)), U)
    assert np.allclose(krr, uC[~on_c, 1:], atol=1e-8)


def test_F_cold_greedy_matches_bruteforce():
    w, v = _view(23, "hub", 10)
    r = 4
    sel, picks = _det_select("F", v, 10, r=r)
    pool = base.pool_of(v)
    K = np.array(kernels.external_kernel(v, pool, v.gene_features(pool)))
    n = len(K)
    J = np.eye(n) - 1.0 / n
    ev, V = np.linalg.eigh(J @ K @ J)
    o = np.argsort(ev)[::-1][: r - 1]
    E = V[:, o] * np.sqrt(ev[o])
    E = E / np.sqrt((E ** 2).mean(0))                 # sign irrelevant to the log det
    u = np.column_stack([np.ones(n), E])
    uc, us = u[np.searchsorted(pool, v.candidates)], u[np.searchsorted(pool, v.screened_orfs)]
    M0 = np.eye(r) + us.T @ us

    def gain(picked, j):
        M = M0 + uc[picked].T @ uc[picked]
        return np.linalg.slogdet(M + np.outer(uc[j], uc[j]))[1] - np.linalg.slogdet(M)[1]
    bf = _bf_greedy(gain, len(v.candidates), 10, clusters=v.clusters(v.candidates))
    assert picks == [str(v.candidates[i]) for i in bf]


# ============================================================================ M
def test_M_greedy_matches_bruteforce_gp_posterior_variance():
    w, v = _view(31, "hub", 40)
    sel, picks = _det_select("M", v, 10)
    Zs, Zc = kernels.gi_inputs(v)
    lam = sel.lam * float(base.knobs()["kernels"]["w_gi"])
    on = v.candidate_col >= 0

    def var(picked, j):
        Z = np.vstack([Zs] + [Zc[[i]] for i in picked])
        k = Z @ Zc[j]
        return Zc[j] @ Zc[j] - k @ np.linalg.solve(Z @ Z.T + lam * np.eye(len(Z)), k)
    bf = _bf_greedy(lambda p, j: var(p, j) if on[j] else -np.inf, len(v.candidates), 10)
    assert picks == [str(v.candidates[i]) for i in bf]
    assert all(v.candidate_col[list(v.candidates).index(g)] >= 0 for g in picks)


# ============================================================================ C
@pytest.mark.parametrize("normalize", [False, True])
def test_C_greedy_matches_bruteforce_facility_location(normalize):
    w, v = _view(41, "cluster", 30)
    sel, picks = _det_select("C", v, 10, normalize=normalize)
    pool = base.pool_of(v)
    Kx = np.array(kernels.external_kernel(v, pool, v.gene_features(pool)))
    Zs, Zc = kernels.gi_inputs(v)
    Zall = np.zeros((len(pool), Zs.shape[1]))
    Zall[np.searchsorted(pool, v.screened_orfs)] = Zs
    Zall[np.searchsorted(pool, v.candidates)] = Zc
    K = Kx + Zall @ Zall.T                           # mk kernel over the pool
    if normalize:
        d = np.sqrt(np.diag(K))
        K = K / np.outer(d, d)
    ic, is_ = np.searchsorted(pool, v.candidates), np.searchsorted(pool, v.screened_orfs)

    def F(sel_idx):
        cols = list(is_) + [ic[j] for j in sel_idx]
        return np.maximum(K[np.ix_(ic, cols)].max(1), 0).sum() if cols else 0.0
    bf = _bf_greedy(lambda p, j: F(p + [j]) - F(p), len(v.candidates), 10)
    assert picks == [str(v.candidates[i]) for i in bf]


# ============================================================================ edge cases
def test_switch_at_n30():
    for n, ph_h, ph_f in ((SWITCH - 1, "prior", "cold"), (SWITCH, "model", "fitted")):
        w, v = _view(51, "hub", n)
        for name, ph in (("H", ph_h), ("F", ph_f)):
            sel = strategies.make(name)
            rng = np.random.default_rng(0)
            sel.reset(v.candidates, rng)
            sel.select(v, 10, rng)
            assert sel.diagnostics[-1]["phase"] == ph, (name, n, sel.diagnostics[-1])


@pytest.mark.parametrize("name", list(strategies.REGISTRY))
def test_all_candidates_offarray(name):
    w = ToyWorld(52, "hub")
    off = np.array([g for g in w.queries if g not in w.col])
    on = np.array([g for g in w.queries if g in w.col])
    S = list(np.random.default_rng(0).choice(on, 40, replace=False))
    pool = np.union1d(off, S)
    v = w.view(S, pool=pool)
    sel = strategies.make(name)
    rng = np.random.default_rng(0)
    sel.reset(v.candidates, rng)
    picks = sel.select(v, 10, rng)
    assert len(set(picks)) == 10 and np.isin(picks, off).all()
    if name == "M":                               # criterion never picks off-array; batch filled at random
        assert sel.diagnostics[-1].get("n_filled_at_random") == 10


@pytest.mark.parametrize("name", list(strategies.REGISTRY))
@pytest.mark.parametrize("n_cand,k", [(7, 10), (25, 3)])
def test_pool_exhaustion_and_remainder(name, n_cand, k):
    w = ToyWorld(53, "hub")
    rng0 = np.random.default_rng(1)
    S = list(rng0.choice(w.queries, 35, replace=False))
    rest = np.setdiff1d(w.queries, S)
    pool = np.union1d(S, rng0.choice(rest, n_cand, replace=False))
    v = w.view(S, pool=pool)
    sel = strategies.make(name)
    rng = np.random.default_rng(0)
    sel.reset(v.candidates, rng)
    picks = sel.select(v, k, rng)
    assert len(picks) == min(k, n_cand) and len(set(picks)) == len(picks) and np.isin(picks, v.candidates).all()


def _scrambled(v, seed):
    """Copy of a toy view whose screened-row data (E, Y, SP) are replaced by noise; everything else identical."""
    import copy
    sv = copy.copy(v)
    rng = np.random.default_rng(seed)
    sv.Y = (rng.random(v.Y.shape) < 0.05).astype(np.int8)
    sv.E = np.where(sv.Y == 1, -0.3, rng.normal(0, 0.05, v.Y.shape)).astype(np.float32)
    sv.SP = np.where(sv.Y == 1, 5.0, 0.0).astype(np.float32)
    return sv


@pytest.mark.parametrize("name", ["R1", "H-static", "F-static", "C-static"])
def test_static_variants_ignore_screened_row_content(name):
    w = ToyWorld(54, "hub")
    for n in (0, 40, 80):
        S = list(np.random.default_rng(n).choice(w.queries, n, replace=False)) if n else []
        v = w.view(S)
        outs = []
        for view in (v, _scrambled(v, 99)):
            sel = strategies.make(name)
            rng = np.random.default_rng(3)
            sel.reset(w.queries, rng)
            outs.append(sel.select(view, 10, rng))
        assert outs[0] == outs[1], (name, n)


# ============================================================================ off-array routing (DESIGN §1, D5)
def _offarray_share_after_switch(name, n_worlds=12, n_max=70, **kw):
    sh = []
    for s in range(n_worlds):
        w = ToyWorld(7000 + s, "null")
        S = run_toy(w, strategies.make(name, **kw), s, n_max)
        sh += [g not in w.col for g in S[SWITCH:]]
    return float(np.mean(sh)), len(sh)


def _check_base_rate(name, **kw):
    share, n = _offarray_share_after_switch(name, **kw)
    base_rate = 40 / 240
    se = np.sqrt(base_rate * (1 - base_rate) / n)
    assert abs(share - base_rate) < 4 * se, (name, share, base_rate)


@pytest.mark.parametrize("name", ["H", "H-static", "F-static"])
def test_null_world_offarray_share_at_base_rate(name):
    _check_base_rate(name)


@pytest.mark.xfail(strict=True, reason="REVIEW F-1: F's KRR route for off-array candidates is a de-facto off-array "
                                       "score term (0% off-array picks after the switch in a null world)")
def test_null_world_offarray_share_F():
    _check_base_rate("F")


@pytest.mark.xfail(strict=True, reason="REVIEW C-1: C's GI part exists only for on-array candidates, so off-array "
                                       "genes are never picked once S has on-array rows (null world)")
def test_null_world_offarray_share_C():
    _check_base_rate("C")


@pytest.mark.xfail(strict=True, reason="REVIEW R1-1: ridge-imputed SMF puts (nearly) all off-array candidates in one "
                                       "SMF tertile, so the bin round-robin under-samples them")
def test_R1_offarray_spread_over_smf_tertiles():
    fr = []
    for s in range(6):
        w = ToyWorld(7100 + s, "null")
        v = w.view([])
        r1 = strategies.make("R1")
        r1.reset(v.candidates, np.random.default_rng(0))
        b, _ = r1.bins(v)
        off = v.candidate_col < 0
        ns = int(base.knobs()["spec"]["n_abund_bins"])
        fr.append(np.bincount(b[off] // ns, minlength=3).max() / off.sum())
    assert np.mean(fr) < 0.6, fr                    # uniform would be ~1/3


def test_M_round0_may_contain_offarray_genes_spec_note():
    """DESIGN: M round 0 = first R0 batch, while 'off-array never picked by M'. The code follows the former."""
    w = ToyWorld(55, "null")
    seen = any(any(g not in w.col for g in run_toy(w, strategies.make("M"), s, 10)) for s in range(10))
    rnd = any(any(g not in w.col for g in run_toy(w, harness.RandomSelector(), s, 10)) for s in range(10))
    assert seen == rnd


@pytest.mark.xfail(strict=True, reason="REVIEW B-1: the fixed per-batch cluster cap (2 x n_clusters slots) cannot fill "
                                       "large batches (§6.3 batch 110 / one-shot 330, 660 with k = 100 clusters); the "
                                       "remainder is filled at random by Selector.select")
@pytest.mark.parametrize("name", ["H", "F", "H-static", "F-static"])
def test_large_batch_not_filled_at_random(name):
    w = ToyWorld(56, "hub")                         # 10 clusters -> cap admits 20 genes per batch
    v = w.view([])
    sel = strategies.make(name)
    rng = np.random.default_rng(0)
    sel.reset(v.candidates, rng)
    sel.select(v, 60, rng)
    assert sel.diagnostics[-1].get("n_filled_at_random", 0) == 0, sel.diagnostics[-1]
