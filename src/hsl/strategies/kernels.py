"""Cached kernel module shared by C, M, F and R1 (DESIGN v3.1 §6.2; review r1 item 2).

* GI part = the predictor's `mk` GI kernel: upstream `completion.foldin_inputs(·, 'sym', 'negeps', 'row')` run on the
  SelectorView's (masked) screened rows — screened gene s: 'sym' profile over the screened on-array coordinates T_on;
  candidate x: its array-side column c_x over T_on (0 for off-array x) — scaled as upstream `mk_predict`
  (c = sqrt(w_gi / mean‖Z_s‖²)). Recomputed every round (depends on S).
* External part = the predictor's tier-1 external kernels (upstream foldin2.EXT_KERNELS transforms/scales and
  foldin2._psd; weights from config/strategies.yaml, default = MK_RECIPE): per kernel the transformed similarity over
  the gene set, unit diagonal for genes covered by the source, PSD projection, divided by its mean diagonal; summed
  with weights. Sources available through the SelectorView: coex = similarity('coexp'); apmse = similarity('apms')
  (edge flag); apms = AP-MS neighbourhood cosine with neighbourhoods restricted to the pool genes (derived); esm =
  cosine of per-dimension standardised ESM-2 PCs (gene_features esm_pc*; the raw ESM cosine exposed by
  similarity('esm') is anisotropic, 0.85–0.99). Screened-set independent -> cached IN MEMORY per (source, gene set,
  weights).
* Kernel PCA of the external kernel (F cold start / F-static) and R1's k-means on the visible-gene universe (dev:
  excludes D_val as well as T), both cached in memory.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from hsl.strategies import base

_CACHE: dict = {}
_MAX_CACHE = 12


def _cached(key, fn):
    if key not in _CACHE:
        if len(_CACHE) >= _MAX_CACHE:
            _CACHE.pop(next(iter(_CACHE)))
        _CACHE[key] = fn()
    return _CACHE[key]


def clear_cache() -> None:
    _CACHE.clear()


# ============================================================================ GI ('sym') inputs
def gi_inputs(view, w_gi: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(Z_S [n_screened x n_on], Z_C [n_candidates x n_on]) scaled GI inputs (rows aligned with view.screened_orfs /
    view.candidates; off-array candidates have zero rows). n_on = 0 when no screened gene is on the array."""
    from ysl.models import completion as CM
    w_gi = float(base.knobs()["kernels"]["w_gi"] if w_gi is None else w_gi)
    n, nc = int(view.n_screened), len(view.candidates)
    scol, ccol = np.asarray(view.screened_col), np.asarray(view.candidate_col)
    if n == 0 or not (scol >= 0).any() or w_gi <= 0:
        return np.zeros((n, 0)), np.zeros((nc, 0))
    duck = SimpleNamespace(E=np.asarray(view.E), train_col=scol, target_col=ccol, n=n)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Z, Zx_on, _, _ = CM.foldin_inputs(duck, "sym", "negeps", "row")
    s = float(np.einsum("ij,ij->", Z, Z)) / max(len(Z), 1)
    c = np.sqrt(w_gi / (s if s > 0 else 1.0))
    Zc = np.zeros((nc, Z.shape[1]))
    Zc[ccol >= 0] = Zx_on
    return Z * c, Zc * c


# ============================================================================ external similarities
def _esm_std_cosine(view, genes, feats) -> np.ndarray:
    cols = [c for c in feats.columns if c.startswith("esm_pc")]
    X = feats[cols].to_numpy(np.float64)
    ok = np.isfinite(X).all(1)
    if ok.sum() < 2:
        return np.full((len(genes), len(genes)), np.nan)
    mu, sd = X[ok].mean(0), X[ok].std(0)
    X = (X - mu) / np.where(sd > 0, sd, 1.0)
    X = np.nan_to_num(X) / np.linalg.norm(np.nan_to_num(X), axis=1, keepdims=True).clip(1e-12)
    S = X @ X.T
    S[~ok, :] = np.nan
    S[:, ~ok] = np.nan
    return S


def _apms_neigh_cosine(view, genes) -> np.ndarray:
    A = np.asarray(view.similarity("apms", genes, genes), dtype=np.float64)
    obs = ~np.isnan(A).all(1)
    A0 = np.nan_to_num(A)
    np.fill_diagonal(A0, 0.0)
    nrm = np.linalg.norm(A0, axis=1)
    S = (A0 @ A0.T) / np.outer(np.where(nrm > 0, nrm, 1.0), np.where(nrm > 0, nrm, 1.0))
    S[~obs, :] = np.nan
    S[:, ~obs] = np.nan
    return S


def raw_similarity(view, kind: str, genes, feats) -> np.ndarray:
    if kind == "coex":
        return np.asarray(view.similarity("coexp", genes, genes), dtype=np.float64)
    if kind == "apmse":
        return np.asarray(view.similarity("apms", genes, genes), dtype=np.float64)
    if kind == "apms":
        return _apms_neigh_cosine(view, genes)
    if kind == "esm":
        return _esm_std_cosine(view, genes, feats)
    raise ValueError(f"external kernel {kind!r} is not derivable from the SelectorView")


def _single_kernel(view, kind, genes, feats) -> np.ndarray:
    from ysl.models import foldin2 as F2
    _, tf, scale = F2.EXT_KERNELS[kind]
    S = raw_similarity(view, kind, genes, feats)
    cov = ~np.isnan(S).all(1)
    K = F2._transform(S, tf, scale)
    np.fill_diagonal(K, cov.astype(np.float64))
    K = F2._psd(K)
    s = float(np.trace(K)) / max(len(K), 1)
    return K / (s if s > 0 else 1.0)


def external_kernel(view, genes, feats, weights: dict | None = None) -> np.ndarray:
    """Weighted sum of external kernels over `genes` (sorted gene set; feats = gene_features(genes)). Cached."""
    genes = np.asarray(genes, dtype=str)
    w = {k: float(v) for k, v in (weights or base.knobs()["kernels"]["external"]).items() if float(v) != 0.0}
    key = ("ext", base.genes_key(view, genes, sorted(w.items())))

    def build():
        K = np.zeros((len(genes), len(genes)))
        for kind, wk in w.items():
            K += wk * _single_kernel(view, kind, genes, feats)
        K.setflags(write=False)
        return K
    return _cached(key, build)


def kernel_pca(K: np.ndarray, dim: int) -> np.ndarray:
    """Top-`dim` kernel-PCA scores of the centred kernel (deterministic LAPACK subset solver; sign fixed so that the
    largest-|loading| entry is positive)."""
    from scipy.linalg import eigh
    n = len(K)
    dim = int(min(dim, n - 1))
    if dim <= 0:
        return np.zeros((n, 0))
    rs = K.mean(1)
    Kc = K - rs[:, None] - rs[None, :] + rs.mean()
    w, V = eigh((Kc + Kc.T) / 2, subset_by_index=[n - dim, n - 1])
    w, V = w[::-1], V[:, ::-1]
    sign = np.sign(V[np.abs(V).argmax(0), np.arange(V.shape[1])])
    sign[sign == 0] = 1.0
    return V * sign * np.sqrt(np.clip(w, 0, None))


def external_kpca(view, genes, feats, dim: int) -> np.ndarray:
    genes = np.asarray(genes, dtype=str)
    w = base.knobs()["kernels"]["external"]
    key = ("kpca", base.genes_key(view, genes, sorted(w.items()), int(dim)))
    return _cached(key, lambda: kernel_pca(external_kernel(view, genes, feats, w), dim))


# ============================================================================ R1 clusters on the visible universe
def visible_universe(view) -> np.ndarray:
    """Candidates ∪ screened ∪ visible (unmasked) array genes: non-test genes, dev: also minus D_val."""
    arr = np.asarray(view.array_orfs, dtype=str)[np.asarray(view.visible_cols)]
    return np.union1d(base.pool_of(view), arr)


def visible_clusters(view, feats_fn) -> tuple[np.ndarray, np.ndarray]:
    """(genes, labels): k-means (k = protocol kmeans_k, seed cluster_seed, n_init 10) on the kernel-PCA embedding
    (protocol cluster_embedding_dim) of the combined external kernel of hsl.features (coexpression z clipped to
    [−1, 1], standardised-ESM cosine, AP-MS edge; mean; unit diagonal) over the visible universe. Cached."""
    from sklearn.cluster import KMeans
    from hsl import config
    pr = base.protocol_selector()
    U = visible_universe(view)
    k, dim, seed = int(pr["kmeans_k"]), int(pr["cluster_embedding_dim"]), config.seed("cluster_seed")
    key = ("kmeans", base.genes_key(view, U, k, dim, seed))

    def build():
        feats = feats_fn(view, U)
        Kc = np.clip(np.nan_to_num(np.asarray(view.similarity("coexp", U, U), dtype=np.float64)), -1, 1)
        Ke = np.nan_to_num(_esm_std_cosine(view, U, feats))
        Ka = np.nan_to_num(np.asarray(view.similarity("apms", U, U), dtype=np.float64))
        K = (Kc + Ke + Ka) / 3.0
        K = (K + K.T) / 2
        np.fill_diagonal(K, 1.0)
        Z = kernel_pca(K, dim)
        lab = KMeans(n_clusters=min(k, len(U)), n_init=10, random_state=seed).fit_predict(Z).astype(np.int64)
        return U, lab
    return _cached(key, build)
