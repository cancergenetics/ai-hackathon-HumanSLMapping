"""Whitelisted per-gene EXTERNAL features for selectors (DESIGN v3 §1, §6.1(e); config/feature_provenance.yaml).

Everything here is label-free and identical for every split (tier-1 stores read through the upstream loaders):
  gene_features(orfs)          ysl.features.genes.gene_matrix minus its `has_*` availability indicators (no
                               missingness indicators reach a selector; missing values stay NaN). Blocks: sequence,
                               abundance (PaxDb), popgen (1002 genomes), AlphaFold, UniProt domain counts, AP-MS node
                               stats (Michaelis 2023), expression level/variability + coexpression hub-ness
                               (COXPRESdb, perturbation studies removed), Compara conservation/paralogs, ESM-2 PCs.
                               Genome coordinates are NOT included.
  similarity(kind, a, b)       'coexp' | 'esm' | 'apms' gene x gene (ysl.features.genes.similarity).
  clusters(orfs)               k-means (k = protocol selector.kmeans_k, seed cluster_seed) on a kernel-PCA embedding
                               of the combined external kernel (coexpression z, standardised-ESM cosine, AP-MS edge),
                               fitted ONCE on all non-test S1 genes (query ∪ array ORFs minus T).
Array SMF is not here: selectors get it from the SelectorView (masked for excluded genes).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import hashlib
from functools import lru_cache

import numpy as np
import pandas as pd

from hsl import config, paths

EXCLUDED_PREFIXES = ("has_",)
KERNELS = ("coexp", "esm", "apms")


def _G():
    from ysl.features import genes as G
    return G


def current_names(orfs) -> np.ndarray:
    return _G().current_names(np.asarray(orfs, dtype=str))


def feature_columns() -> list[str]:
    cols = _G().gene_table_current().columns
    return [c for c in cols if not c.startswith(EXCLUDED_PREFIXES)]


def gene_features(orfs) -> pd.DataFrame:
    """Whitelisted gene scalars for `orfs` (Costanzo names; rows in the given order)."""
    t = _G().gene_matrix(np.asarray(orfs, dtype=str))
    return t[feature_columns()].astype(np.float64)


def similarity(kind: str, orfs_a, orfs_b) -> np.ndarray:
    if kind not in KERNELS:
        raise ValueError(f"kernel {kind!r} not whitelisted ({KERNELS})")
    return _G().similarity(kind, np.asarray(orfs_a, dtype=str), np.asarray(orfs_b, dtype=str), current=False)


def _combined_kernel(genes: np.ndarray) -> np.ndarray:
    """Mean of three label-free kernels over `genes` (NaN -> 0, diagonal 1): coexpression Fisher z clipped to [-1, 1];
    cosine of per-dimension standardised ESM-2 embeddings; AP-MS FDR<5% edge (0/1)."""
    G = _G()
    cur = G.current_names(genes)
    Kc = np.clip(np.nan_to_num(G.similarity("coexp", cur, cur, current=True)), -1, 1)
    e = G._esm_current()
    X = e.reindex(cur).to_numpy(dtype=np.float64)
    mu, sd = np.nanmean(e.to_numpy(np.float64), 0), np.nanstd(e.to_numpy(np.float64), 0)
    X = (X - mu) / np.where(sd > 0, sd, 1)
    X = X / np.linalg.norm(np.nan_to_num(X), axis=1, keepdims=True).clip(1e-12)
    Ke = np.nan_to_num(X @ X.T)
    Ka = np.nan_to_num(G.similarity("apms", cur, cur, current=True))
    K = (Kc + Ke + Ka) / 3.0
    K = (K + K.T) / 2
    np.fill_diagonal(K, 1.0)
    return K


def cluster_universe() -> np.ndarray:
    """All S1 query ∪ array ORFs minus the test queries (names only)."""
    from hsl import views
    b = views.dev_bundle()
    return np.setdiff1d(np.union1d(b.all_query_orfs, b.array_orfs), views.test_queries())


@lru_cache(maxsize=1)
def _cluster_table() -> pd.Series:
    from sklearn.cluster import KMeans
    pr = config.protocol()["selector"]
    k, dim, seed = int(pr["kmeans_k"]), int(pr["cluster_embedding_dim"]), config.seed("cluster_seed")
    genes = cluster_universe()
    key = hashlib.sha256(("|".join(genes) + f"|{k}|{dim}|{seed}").encode()).hexdigest()[:16]
    cache = paths.output_root() / "data" / "cache" / f"clusters_{key}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        return pd.Series(z["label"], index=pd.Index(z["genes"], name="orf"), name="cluster")
    K = _combined_kernel(genes)
    n = len(K)
    J = np.eye(n) - 1.0 / n
    Kc = J @ K @ J
    w, V = np.linalg.eigh(Kc)
    top = np.argsort(w)[::-1][:dim]
    Z = V[:, top] * np.sqrt(np.clip(w[top], 0, None))
    lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(Z).astype(np.int32)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, genes=genes, label=lab)
    return pd.Series(lab, index=pd.Index(genes, name="orf"), name="cluster")


def clusters(orfs) -> np.ndarray:
    """Cluster id per ORF (-1 if not in the non-test universe)."""
    t = _cluster_table()
    return t.reindex(np.asarray(orfs, dtype=str)).fillna(-1).to_numpy(dtype=np.int64)
