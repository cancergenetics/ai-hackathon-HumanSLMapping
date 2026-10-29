"""Fold-in family v2 (M4b): multi-kernel fold-in, external-only fold-in (cold start), a neural fold-in, fold-in
ensembles, and the cross-fitted stack `stack_lr_wr2`.

All methods read only the ScreenView (+ label-free external pair features, `ysl.features.pairs`). Notation as in
`ysl.models.completion`: T = screened queries (view rows), T_on = screened queries with an array strain, x = held-out
target, c_x = x's array-side column over T (the only GI information about x; absent for off-array x).

Multi-kernel fold-in (`mk`)
  B6 is kernel ridge over screened genes with the linear GI-profile kernel K_gi = Z Z^T (Z = 'sym' inputs at T_on
  coordinates) from a gene's profile to its array-centred label row. Here the gene kernel is
      K = w_gi K_gi / s_gi + sum_k w_k K_k / s_k            (s = mean diagonal on the training genes)
  with external kernels K_k(g, h) = transformed similarity of genes g, h from the pair-feature store (coexpression,
  mutual rank, ESM cosine, AP-MS neighbourhood cosine / edge, phylogenetic profile, paralog identity, domains). The
  training block K_k(T, T) has unit diagonal for genes covered by the source (0 otherwise) and is projected onto the
  PSD cone; the test block K_k(x, T) is the transformed similarity itself (NaN -> 0). lambda = lam * trace(K) / n
  (lam_mode 'total'; lam = 10 and no external kernel reproduce B6 exactly) or lam * w_gi ('gi'). The GI part of an off-array x is 0, so an off-array x is
  scored by the external kernels alone. The target-equals-input coordinate is excluded exactly as in B6 (per-column
  Sherman-Morrison downdate of the scaled GI part; external kernels contain no GI data).
  w_gi = 0 gives the external-only fold-in (`ext`): no coordinate exclusion needed; scores every target.

Neural fold-in (`nn`; torch)
  encoder: [x's input profile over T_on (masked), availability mask, external gene features] -> MLP -> d-dim h_x;
  decoder: logit(x, a) = h_x . V_a + b_a + c (free array embeddings V, array bias b). Trained on screened genes with
  BCE on observed labels. Per example a random fraction of input coordinates is masked (dropped to 0 with a mask
  flag) and, for ~19% of examples, the whole GI input is dropped (simulated off-array). Target coordinate exclusion is
  EXPLICIT: an output column a = col(q) enters the loss only if the input coordinate q is masked for that example
  (`loss_weights`; tested with autograd). Training input for a screened gene = 'sym' profile (mean of its row and
  column at T_on coordinates), test input for x = c_x at T_on (as in B6).

Within-row feature function for the stack: `wr2_features(view, params, names)` (completion's WR features + those
here), the same function for cross-fit fold views and the full view.
"""
from __future__ import annotations

import warnings

import numpy as np

from ysl.features import pairs as PF
from ysl.models import baselines as B
from ysl.models import completion as CM
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

# ----------------------------------------------------------------------------- external gene kernels
# name -> (pair feature, transform, scale). transforms: 'pos' max(S, 0)/scale clipped to 1; 'lin' S/scale clipped
# to [-1, 1]; 'exp_mr' exp(-(S - 1)/scale) (mutual rank: 1 = best); 'flag' (S > 0).
EXT_KERNELS = {
    "coex": ("coex_z", "pos", 1.0),
    "coexmr": ("coex_mr", "exp_mr", 10.0),
    "coexmr50": ("coex_mr", "exp_mr", 50.0),
    "esm": ("esm_cos", "pos", 0.5),
    "apms": ("apms_cosine", "pos", 1.0),
    "apmse": ("apms_edge", "flag", 1.0),
    "phylo": ("phylo_corr", "pos", 1.0),
    "para": ("para_identity", "pos", 100.0),
    "dom": ("dom_ipr_jaccard", "pos", 1.0),
}


def _transform(S: np.ndarray, tf: str, scale: float) -> np.ndarray:
    S = S.astype(np.float64)
    with np.errstate(invalid="ignore", over="ignore"):
        if tf == "pos":
            K = np.minimum(np.maximum(S, 0.0) / scale, 1.0)
        elif tf == "lin":
            K = np.clip(S / scale, -1.0, 1.0)
        elif tf == "exp_mr":
            K = np.exp(-(S - 1.0) / scale)
        elif tf == "flag":
            K = (S > 0).astype(np.float64)
        else:
            raise ValueError(tf)
    return np.where(np.isnan(S), 0.0, K)


def _psd(K: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh((K + K.T) / 2)
    return (U * np.maximum(w, 0.0)[None, :]) @ U.T


def ext_kernel(v: ScreenView, name: str) -> tuple[np.ndarray, np.ndarray]:
    """(K_TT [n x n] PSD with unit diagonal for covered genes, K_xT [n_targets x n]) for external kernel `name`,
    both divided by the mean diagonal of K_TT. Label-free (external pair features only)."""
    feat, tf, scale = EXT_KERNELS[name]
    ti = PF.index_of(v.train_orfs_current)
    xi = PF.index_of(v.target_orfs_current)
    S_tt = PF.block(feat, ti, ti)
    cov = ~np.isnan(S_tt).all(1)                         # gene covered by the source (some finite similarity)
    K = _transform(S_tt, tf, scale)
    np.fill_diagonal(K, cov.astype(np.float64))
    K = _psd(K)
    Kx = _transform(PF.block(feat, xi, ti), tf, scale)
    s = float(np.trace(K)) / max(len(K), 1)
    s = s if s > 0 else 1.0
    return K / s, Kx / s


# ----------------------------------------------------------------------------- multi-kernel fold-in
MK_DEFAULTS = CM.FOLDIN_DEFAULTS | {"w_gi": 1.0, "kern": {}, "beta": 0.0, "lam_mode": "total"}


def _target(v: ScreenView, target: str, beta: float) -> np.ndarray:
    T = CM.signal(v, target)
    if beta:
        T = T + float(beta) * np.where(np.isnan(T), np.nan, CM.signal(v, "negeps"))
    return T


def mk_predict(v: ScreenView, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Residual predictions for ALL targets (off-array: external kernels only; zero rows if none) and the array mean
    mu: score = mu + P. Requires an input mode whose training genes are all screened genes ('sym' / 'mixed') or 'col'
    (training genes T_on)."""
    p = MK_DEFAULTS | params
    kern = {k: float(w) for k, w in dict(p["kern"] or {}).items() if float(w) != 0.0}
    w_gi = float(p["w_gi"])
    on = np.flatnonzero(v.target_col >= 0)
    Tsig = _target(v, p["target"], p["beta"])
    mu = CM.col_mean(Tsig)
    if w_gi > 0:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            Z, Zx_on, rows, dcols = CM.foldin_inputs(v, p["inp"], p["in_sig"], p["in_centre"])
        if p["norm"]:
            Z, Zx_on = CM._normalise_rows(Z), CM._normalise_rows(Zx_on)
        s = float(np.einsum("ij,ij->", Z, Z)) / max(len(Z), 1)
        c = np.sqrt(w_gi / (s if s > 0 else 1.0))
        Z = Z * c
        Zx = np.zeros((v.n_targets, Z.shape[1]))
        Zx[on] = Zx_on * c
        K = Z @ Z.T
        Kx = Zx @ Z.T
    else:
        rows, dcols = np.arange(v.n), np.full(0, -1)
        K = np.zeros((v.n, v.n))
        Kx = np.zeros((v.n_targets, v.n))
        Z = Zx = None
    for nm, w in kern.items():
        Kt, Kxt = ext_kernel(v, nm)
        K = K + w * Kt[np.ix_(rows, rows)]
        Kx = Kx + w * Kxt[:, rows]
    if not np.any(K):
        return np.zeros((v.n_targets, v.n_arrays)), mu
    Y = np.nan_to_num(Tsig[rows] - mu[None, :])
    # lam_mode 'total': lambda relative to the mean diagonal of the combined kernel (B6 when no external kernel);
    # 'gi': relative to the GI part only (adding external kernels does not raise the ridge on the GI part)
    base = np.trace(K) / max(len(K), 1) if p["lam_mode"] == "total" or w_gi <= 0 else w_gi
    lam = float(p["lam"]) * base
    Minv = np.linalg.inv(K + lam * np.eye(len(K)))
    alpha = Minv @ Y
    P = Kx @ alpha
    excl = np.flatnonzero(dcols >= 0)
    if Z is not None and len(excl):     # exact removal of the GI input coordinate j when predicting column dcols[j]
        Zj = Z[:, excl]
        a_j = dcols[excl]
        Mz = Minv @ Zj
        d = 1.0 - np.einsum("ij,ij->j", Zj, Mz)
        sj = np.einsum("ij,ij->j", Mz, Y[:, a_j])
        al = alpha[:, a_j] + Mz * (sj / d)[None, :]
        u = np.einsum("ij,ij->j", Zj, al)
        P[:, a_j] = Kx @ al - Zx[:, excl] * u[None, :]
    if int(p["rank"]) > 0:
        F = K @ alpha
        _, _, Vt = np.linalg.svd(F, full_matrices=False)
        Vr = Vt[: int(p["rank"])]
        P = (P @ Vr.T) @ Vr
    return P, mu


def kernel_alignment(v: ScreenView, names, params: dict | None = None) -> dict:
    """Centred kernel alignment of each kernel (GI 'gi' + external) with the array-centred label Gram Y Y^T over the
    screened genes (training rows only; diagnostic / weight heuristic)."""
    p = MK_DEFAULTS | (params or {})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Z, _, rows, _ = CM.foldin_inputs(v, p["inp"], p["in_sig"], p["in_centre"])
    Tsig = _target(v, p["target"], p["beta"])
    Y = np.nan_to_num(Tsig[rows] - CM.col_mean(Tsig)[None, :])
    n = len(rows)
    H = np.eye(n) - 1.0 / n
    L = H @ (Y @ Y.T) @ H
    np.fill_diagonal(L, 0.0)                     # off-diagonal alignment (the diagonal is trivially aligned)

    def al(K):
        Kc = H @ K @ H
        np.fill_diagonal(Kc, 0.0)
        return float((Kc * L).sum() / (np.linalg.norm(Kc) * np.linalg.norm(L) + 1e-300))
    out = {"gi": al(Z @ Z.T)}
    for nm in names:
        out[nm] = al(ext_kernel(v, nm)[0][np.ix_(rows, rows)])
    return out


# ----------------------------------------------------------------------------- feature function for the stack
WR2_FEATURES = {  # name -> (kind, params); residual scale, NaN where not defined
    "mk": ("mk", {}),                    # multi-kernel fold-in (params filled in by recipe below)
    "ext": ("ext", {}),                  # external-only fold-in, all targets
    "ext_off": ("ext_off", {}),          # external-only fold-in, off-array targets only (NaN on-array)
    "mk_off": ("mk_off", {}),            # multi-kernel fold-in, off-array targets only
    "nn": ("nn", {}),
}
# Recipes selected on S1 dev seed 0 (scripts/44_explore.py; results/m4b_foldin2/explore_*_seed0.json), then frozen.
MK_RECIPE: dict = {"kern": {"apmse": 0.5, "apms": 0.5, "coex": 0.3, "esm": 0.1}, "lam_mode": "gi", "lam": 10.0,
                   "beta": 1.0, "rank": 100}
EXT_RECIPE: dict = {"kern": {"apms": 1.0, "apmse": 1.0, "coex": 0.5, "esm": 0.3, "para": 0.3, "dom": 0.3},
                    "lam": 5.0, "rank": 0, "beta": 0.0}
NN_RECIPE: dict = {"epochs": 60, "wd": 1e-3, "wd_v": 1e-4, "use_gene": False}


def wr2_score(v: ScreenView, kind: str, params: dict, rng_seed: int = 0) -> np.ndarray:
    on = v.target_col >= 0
    if kind in ("mk", "mk_off"):
        P, _ = mk_predict(v, MK_RECIPE | params)
    elif kind in ("ext", "ext_off"):
        P, _ = mk_predict(v, EXT_RECIPE | {"w_gi": 0.0} | params)
    elif kind == "nn":
        P = nn_residual(v, NN_RECIPE | params, rng_seed)
    else:
        raise ValueError(kind)
    P = np.asarray(P, dtype=np.float32)
    if kind.endswith("_off"):
        P = np.where(on[:, None], np.nan, P).astype(np.float32)
    elif kind == "mk" and not (MK_RECIPE | params).get("kern"):
        P = np.where(on[:, None], P, np.nan).astype(np.float32)
    return P


def wr2_features(v: ScreenView, params: dict, names, rng_seed: int = 0) -> dict:
    """Completion's within-row features (names in completion.WR_FEATURES) plus the ones here; same function for
    train (fold views) and test (full view). params[f'{name}__{key}'] overrides a feature's settings."""
    base = tuple(n for n in names if n in CM.WR_FEATURES)
    out = CM.wr_features(v, params, base, rng_seed) if base else {}
    for nm in names:
        if nm in out:
            continue
        kind, prm = WR2_FEATURES[nm]
        own = {k.split("__", 1)[1]: val for k, val in params.items() if k.startswith(nm + "__")}
        out[nm] = wr2_score(v, kind, prm | own, rng_seed)
    return {n: out[n] for n in names}


# ----------------------------------------------------------------------------- neural fold-in (torch)
NN_DEFAULTS = {"d": 32, "hidden": 0, "drop_in": 0.3, "p_off": 0.19, "epochs": 300, "lr": 3e-3, "wd": 1e-5,
               "wd_v": 1e-5, "aux_eps": 0.0, "use_gene": True, "in_sig": "negeps", "in_centre": "row",
               "n_members": 1, "device": "cpu", "n_threads": 3, "exclude_diag": True}


def nn_inputs(v: ScreenView, p: dict):
    """(Z [n x p] training-gene GI inputs ('sym' at T_on coordinates), Zx [n_targets x p] (c_x at T_on; 0 for
    off-array x), dcols [p] (array column whose target equals input coordinate j), G [n x q], Gx [n_targets x q]
    external gene features (standardised on the training genes; NaN -> 0; availability flags kept))."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Z, Zx_on, rows, dcols = CM.foldin_inputs(v, "sym", p["in_sig"], p["in_centre"])
    assert len(rows) == v.n
    s = float(np.sqrt(np.mean(Z ** 2))) or 1.0
    Zx = np.zeros((v.n_targets, Z.shape[1]))
    Zx[v.target_col >= 0] = Zx_on
    if p["use_gene"]:
        from ysl.features import genes as G
        gt = G.gene_matrix(v.train_orfs_current, current=True).to_numpy(dtype=np.float64)
        gx = G.gene_matrix(v.target_orfs_current, current=True).to_numpy(dtype=np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            m, sd = np.nanmean(gt, 0), np.nanstd(gt, 0)
        m, sd = np.nan_to_num(m), np.where(np.nan_to_num(sd) > 0, sd, 1.0)
        gt, gx = np.nan_to_num((gt - m) / sd), np.nan_to_num((gx - m) / sd)
    else:
        gt, gx = np.zeros((v.n, 0)), np.zeros((v.n_targets, 0))
    return Z / s, Zx / s, dcols, gt, gx


def _nn_model(p_in: int, q: int, n_a: int, p: dict):
    import torch
    from torch import nn

    class FoldinNet(nn.Module):
        def __init__(self):
            super().__init__()
            d, h = int(p["d"]), int(p["hidden"])
            if h > 0:
                self.enc = nn.Sequential(nn.Linear(p_in + q, h), nn.GELU(), nn.Linear(h, d))
            else:
                self.enc = nn.Linear(p_in + q, d)
            self.V = nn.Parameter(torch.randn(n_a, d) * 0.01)
            self.b = nn.Parameter(torch.zeros(n_a))
            self.U = nn.Parameter(torch.zeros(n_a, d))      # auxiliary -eps head
            self.e = nn.Parameter(torch.zeros(n_a))

        def embed(self, z, g):
            return self.enc(torch.cat([z, g], 1))

        def forward(self, z, g):
            h = self.embed(z, g)
            return h @ self.V.T + self.b[None, :], h @ self.U.T + self.e[None, :], h
    return FoldinNet()


def nn_loss(model, tZ, keep, tG, tY, tO, tdc, tE, p: dict):
    """Training loss for one step: inputs masked by `keep` (inverted dropout), BCE on observed cells with the explicit
    target-coordinate exclusion (column dcols[j] of example t counts only if coordinate j is masked for t), optional
    auxiliary -eps regression on the same cells. Returns (loss, weighted per-cell BCE)."""
    import torch
    zin = tZ * keep * (1.0 / (1.0 - float(p["drop_in"])))
    W = tO.clone()
    if p["exclude_diag"]:
        W[:, tdc] = W[:, tdc] * (1.0 - keep)
    logit, epsh, _ = model(zin, tG)
    cell = torch.nn.functional.binary_cross_entropy_with_logits(logit, tY, reduction="none") * W
    loss = cell.sum() / W.sum()
    if float(p["aux_eps"]) > 0:
        loss = loss + float(p["aux_eps"]) * ((epsh - tE) ** 2 * W).sum() / W.sum()
    return loss, cell


def nn_fit(v: ScreenView, p: dict, seed: int):
    """Train one network on the screened genes of v; returns (model, inputs tuple). Deterministic given seed (CPU)."""
    import torch
    p = NN_DEFAULTS | p
    torch.set_num_threads(int(p["n_threads"]))
    torch.manual_seed(int(seed))
    gen = torch.Generator().manual_seed(int(seed) + 1)
    dev = torch.device(p["device"])
    Z, Zx, dcols, gt, gx = nn_inputs(v, p)
    Y = v.Y.astype(np.float64)
    obs = ~np.isnan(Y)
    NE = np.nan_to_num(CM.signal(v, "negeps"))
    rate = (np.nansum(Y, 0) + 0.5) / (obs.sum(0) + 1.0)
    rate = np.clip(rate, 1e-4, 0.5)
    model = _nn_model(Z.shape[1], gt.shape[1], v.n_arrays, p).to(dev)
    with torch.no_grad():
        model.b.copy_(torch.tensor(np.log(rate / (1 - rate)), dtype=torch.float32))
    tZ = torch.tensor(Z, dtype=torch.float32, device=dev)
    tG = torch.tensor(gt, dtype=torch.float32, device=dev)
    tY = torch.tensor(np.nan_to_num(Y), dtype=torch.float32, device=dev)
    tO = torch.tensor(obs, dtype=torch.float32, device=dev)
    tE = torch.tensor(NE, dtype=torch.float32, device=dev)
    tdc = torch.tensor(dcols, dtype=torch.long, device=dev)
    enc_params = list(model.enc.parameters())
    opt = torch.optim.Adam([{"params": enc_params, "weight_decay": float(p["wd"])},
                            {"params": [model.V, model.U], "weight_decay": float(p["wd_v"])},
                            {"params": [model.b, model.e], "weight_decay": 0.0}], lr=float(p["lr"]))
    n, pin = Z.shape
    for _ in range(int(p["epochs"])):
        keep = (torch.rand(n, pin, generator=gen) >= float(p["drop_in"])).float()
        off = (torch.rand(n, 1, generator=gen) < float(p["p_off"])).float()
        keep = (keep * (1 - off)).to(dev)
        loss, _ = nn_loss(model, tZ, keep, tG, tY, tO, tdc, tE, p)
        opt.zero_grad()
        loss.backward()
        opt.step()
    return model, (Z, Zx, dcols, gt, gx)


def nn_predict_logits(v: ScreenView, p: dict, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """(logits [n_targets x n_A], array bias b [n_A]) of one trained member."""
    import torch
    model, (Z, Zx, dcols, gt, gx) = nn_fit(v, p, seed)
    dev = next(model.parameters()).device
    with torch.no_grad():
        logit, _, _ = model(torch.tensor(Zx, dtype=torch.float32, device=dev),
                            torch.tensor(gx, dtype=torch.float32, device=dev))
        return logit.cpu().numpy().astype(np.float64), model.b.detach().cpu().numpy().astype(np.float64)


def _member_seeds(rng: np.random.Generator, n: int) -> list[int]:
    return [int(s) for s in rng.integers(0, 2**31 - 1, size=n)]


def nn_ensemble(v: ScreenView, p: dict, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    p = NN_DEFAULTS | p
    L, Bs = [], []
    for s in _member_seeds(rng, int(p["n_members"])):
        lg, b = nn_predict_logits(v, p, s)
        L.append(lg)
        Bs.append(b)
    return np.mean(L, 0), np.mean(Bs, 0)


def nn_residual(v: ScreenView, params: dict, rng_seed: int) -> np.ndarray:
    """Within-row residual (logit minus array bias) for every target (stack feature)."""
    lg, b = nn_ensemble(v, params, np.random.default_rng(np.random.SeedSequence([int(rng_seed), 7])))
    return lg - b[None, :]


def nn_scores(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    lg, _ = nn_ensemble(v, dict(ctx.params), ctx.rng("nn"))
    return lg.astype(np.float32)


# ----------------------------------------------------------------------------- registered single methods
def _mk_method(ext_only: bool):
    def fn(v: ScreenView, ctx: MethodContext) -> np.ndarray:
        prm = dict(ctx.params) | ({"w_gi": 0.0} if ext_only else {})
        P, mu = mk_predict(v, prm)
        S = (mu[None, :] + P).astype(np.float32)
        if not prm.get("kern"):                     # no external kernel: off-array rows as B6 (mean on-array row)
            on = v.target_col >= 0
            S[~on] = S[on].mean(0)
        return S
    return fn


def nn_method(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return nn_scores(v, ctx)


# ----------------------------------------------------------------------------- stack_lr_wr2
# B6 'foldin' REPLACED by the multi-kernel fold-in 'mk', plus the external-only fold-in for off-array targets.
# (Adding mk next to foldin was identical on dev s0-2: joint delta -0.0000 [-0.0002, +0.0001] -> simpler set ships.)
WR2_STACK_NAMES = ("mk",) + tuple(n for n in CM.WR_STACK_NAMES if n != "foldin") + ("ext_off",)


def stack_wr2_predict(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    """stack_lr_wr with the within-row feature set extended/replaced by fold-in v2 features (`wr2_features`); K = 10
    cross-fit via view.without_queries, the same feature function for training folds and the test view."""
    names = tuple(ctx.get("wr", WR2_STACK_NAMES))
    base = CM.STACK_BASE
    knn = B._knn_needed(base)
    folds = ctx.folds(v)
    loo = B.loo_additive(v, ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR))
    seed = int(ctx.rng("tied").integers(2**31))
    Xs, ys = [], []
    for k in range(ctx.n_folds):
        rows = np.flatnonzero(folds == k)
        fv = v.without_queries(v.train_orfs[rows])
        f = B._loo_override(B.fold_features(fv, ctx.params, knn), loo, rows)
        fw, _ = CM._standardise_wr(wr2_features(fv, ctx.params, names, seed), None)
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        Xs.append(np.hstack([B.design_matrix(f, base), CM._wr_design(fw, names)])[keep])
        ys.append(y[keep])
    X, y = np.concatenate(Xs), np.concatenate(ys).astype(np.int8)
    del Xs
    sc, clf = B.fit_lr(X, y, ctx)
    del X
    f = B.pair_features(v, ctx.params, knn)
    fw, _ = CM._standardise_wr(wr2_features(v, ctx.params, names, seed), None)
    n_t, n_a = v.n_targets, v.n_arrays
    S = np.empty((n_t, n_a), dtype=np.float32)
    for i in range(0, n_t, ctx.chunk_targets):
        sl = slice(i, min(i + ctx.chunk_targets, n_t))
        Xt = np.hstack([B.design_matrix(f, base, sl), CM._wr_design(fw, names, sl)])
        S[sl] = clf.decision_function(sc.transform(Xt)).reshape(-1, n_a)
    ctx.info["lr_coef"] = dict(zip(base + names, np.round(clf.coef_[0], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    return S


register("F2_mk_foldin", defaults=MK_DEFAULTS | MK_RECIPE,
         description="multi-kernel fold-in: GI-profile kernel + external gene kernels (AP-MS, coexpression, ESM); "
                     "label + beta*(-eps) target; off-array x scored by the external kernels")(_mk_method(False))
register("F2_ext_foldin", defaults=MK_DEFAULTS | EXT_RECIPE | {"w_gi": 0.0},
         description="external-only fold-in (cold start): kernel ridge from external gene kernels to the query row")(
    _mk_method(True))
register("F2_nn_foldin", defaults=NN_DEFAULTS | NN_RECIPE,
         description="neural fold-in: masked GI profile (+ gene features) -> embedding; free array embeddings + bias; "
                     "BCE with explicit target-coordinate exclusion and simulated off-array dropout")(nn_method)
register("stack_lr_wr2", defaults=B._LR_DEFAULTS | B.KNN_DEFAULTS | {"wr": list(WR2_STACK_NAMES)}, logit_scale=True,
         description="stack_lr_wr with fold-in v2 within-row features (multi-kernel / external-only fold-in); K=10 "
                     "cross-fitted")(stack_wr2_predict)
