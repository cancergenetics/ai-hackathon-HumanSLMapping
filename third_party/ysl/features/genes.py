"""Gene-level EXTERNAL (non-GI) feature matrix and gene-gene similarities for model code (DESIGN.md §7).

Everything here is built from the tier-1 loaders in `ysl.data.features` (no GI data, no GO, no literature columns)
plus ESM-2 embeddings (`load_esm`). Genes are keyed by CURRENT SGD systematic names internally; public functions take
the caller's (Costanzo) ORF names and map them with `to_current()` (e.g. merged ORFs 'YDL038C' -> 'YDL039C').
Missing values stay NaN; every block has a `has_<block>` availability indicator (1.0 / 0.0).

  gene_matrix(orfs)              -> DataFrame [len(orfs) x d] of single-gene scalars (+ ESM PCA components)
  similarity(kind, orfs_a, orfs_b) -> float32 [na x nb]: 'coexp' (mean per-study Fisher z), 'esm' (cosine),
                                    'apms' (1.0 if a Michaelis 2023 FDR<5% edge in either direction, else 0.0);
                                    NaN where a gene is absent from the source.

Features are label-free and identical for every split, so caching them across views cannot leak held-out rows.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from ysl.data import features as F

SCALARS = {  # block -> columns taken from the processed tables (coordinates deliberately excluded)
    "sequence": ["protein_length_aa", "cds_length_nt", "gc_content", "gc3_content", "enc", "cai"],
    "abundance": ["log10_ppm_mueller2020", "log10_ppm_kulak2014"],
    "popgen": ["dnds_yn00_mean", "pangenome_presence_freq", "copy_number_mean", "copy_number_frac_multicopy",
               "lof_freq"],
    "alphafold": ["plddt_mean", "fraction_disordered"],
    "uniprot": ["n_pfam", "n_interpro"],
    "ppi_nodes": ["ppi_bait_observed", "ppi_in_network", "ppi_degree_fdr05", "tabaq_abundance"],
}
N_ESM_PCS = 16
COEXP_TOP = 20          # coexpression hub-ness: mean of a gene's top-20 partner z
COEXP_STRONG = 0.5      # ... and the fraction of measured partners with z > 0.5


# ----------------------------------------------------------------------------- pure helpers (unit-tested)
def paralog_summary(pairs: pd.DataFrame) -> pd.DataFrame:
    """Per gene: number of within-species paralogs, max % identity (gene's own side) to one, number of older
    ('other') paralogs. `pairs` has columns orf_a, orf_b, homology_type, identity_a, identity_b."""
    rows = []
    for side, other, ident in (("orf_a", "orf_b", "identity_a"), ("orf_b", "orf_a", "identity_b")):
        rows.append(pd.DataFrame({"orf": pairs[side].values, "partner": pairs[other].values,
                                  "type": pairs["homology_type"].values, "identity": pairs[ident].values}))
    long = pd.concat(rows, ignore_index=True).drop_duplicates(["orf", "partner", "type"])
    wsp = long[long["type"] == "within_species_paralog"]
    oth = long[long["type"] == "other_paralog"]
    out = pd.DataFrame(index=pd.Index(sorted(set(long["orf"])), name="orf"))
    out["n_paralogs_wsp"] = wsp.groupby("orf").size().reindex(out.index).fillna(0).astype(float)
    out["max_paralog_identity"] = wsp.groupby("orf")["identity"].max().reindex(out.index)  # NaN if none
    out["n_paralogs_other"] = oth.groupby("orf").size().reindex(out.index).fillna(0).astype(float)
    return out


def phylo_summary(profile: pd.DataFrame) -> pd.DataFrame:
    """Per gene: fraction of species with an ortholog (identity > 0) and mean identity over those species."""
    P = profile.to_numpy(dtype=np.float64)
    present = P > 0
    n = present.sum(1)
    with np.errstate(invalid="ignore"):
        mean_id = np.where(n > 0, (P * present).sum(1) / np.maximum(n, 1), np.nan)
    return pd.DataFrame({"phylo_breadth": n / P.shape[1], "phylo_mean_identity": mean_id}, index=profile.index)


def coexp_hubness(z: np.ndarray, top: int = COEXP_TOP, strong: float = COEXP_STRONG) -> tuple[np.ndarray, np.ndarray]:
    """Per gene (row of a symmetric NaN-aware z matrix, diagonal ignored): mean of top-`top` partner z and the
    fraction of measured partners with z > `strong`. NaN for genes with < top measured partners."""
    Z = np.array(z, dtype=np.float32, copy=True)
    np.fill_diagonal(Z, np.nan)
    meas = ~np.isnan(Z)
    nm = meas.sum(1)
    Zf = np.where(meas, Z, -np.inf)
    k = min(top, Z.shape[1] - 1)
    part = -np.partition(-Zf, k - 1, axis=1)[:, :k]
    topmean = np.where(nm >= k, np.where(np.isfinite(part), part, 0).mean(1), np.nan)
    frac = np.where(nm > 0, (np.where(meas, Z, -np.inf) > strong).sum(1) / np.maximum(nm, 1), np.nan)
    return topmean.astype(np.float64), frac.astype(np.float64)


def esm_pca(X: np.ndarray, n_comp: int = N_ESM_PCS) -> np.ndarray:
    """Deterministic PCA scores (SVD of the centred matrix; component signs fixed so the largest |loading| is > 0)."""
    Xc = X.astype(np.float64) - X.mean(0, keepdims=True)
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    sign = np.sign(Vt[np.arange(len(Vt)), np.abs(Vt).argmax(1)])
    return (U[:, :n_comp] * s[:n_comp] * sign[:n_comp])


# ----------------------------------------------------------------------------- cached sources (current names)
@lru_cache(maxsize=1)
def _esm_current() -> pd.DataFrame:
    e = F.load_esm()
    cur = F.to_current(e.index).values
    e = e.set_axis(pd.Index(cur, name="orf"))
    return e[~e.index.duplicated(keep="first")]


@lru_cache(maxsize=1)
def _coexp():
    z, genes, _ = F.load_coexpression("combined")
    return z, pd.Index(genes)


@lru_cache(maxsize=1)
def _apms_edges() -> pd.DataFrame:
    e = F.load_ppi_edges()
    e = e[e["f_fdr_0_05"].astype(bool)]
    return pd.DataFrame({"a": e["bait_orf"].astype(str).values, "b": e["prey_orf"].astype(str).values})


@lru_cache(maxsize=1)
def gene_table_current() -> pd.DataFrame:
    """All external gene scalars keyed by current ORF name (universe = SGD sequence table)."""
    seq = F.load_table("sequence")
    out = pd.DataFrame(index=pd.Index(seq.index.astype(str), name="orf"))
    for block, cols in SCALARS.items():
        t = F.load_table(block)
        out[f"has_{block}"] = out.index.isin(t.index).astype(float)
        out = out.join(t[cols].astype(float), how="left")
    out["log_protein_length"] = np.log(out.pop("protein_length_aa"))
    out["log_cds_length"] = np.log(out.pop("cds_length_nt"))
    out["log1p_ppi_degree"] = np.log1p(out.pop("ppi_degree_fdr05"))
    # expression level / variability per platform (kept, non-perturbation samples)
    for plat in ("rnaseq", "microarray"):
        m, genes, _ = F.load_expression(plat)
        with np.errstate(invalid="ignore"):
            mu, sd = np.nanmean(m, 1), np.nanstd(m, 1)
        df = pd.DataFrame({f"expr_mean_{plat}": mu, f"expr_sd_{plat}": sd}, index=pd.Index(genes.astype(str)))
        df = df[~df.index.duplicated()]
        out[f"has_expr_{plat}"] = out.index.isin(df.index).astype(float)
        out = out.join(df, how="left")
    z, zg = _coexp()
    top, frac = coexp_hubness(z)
    df = pd.DataFrame({"coexp_top20_mean": top, "coexp_frac_strong": frac}, index=pd.Index(zg.astype(str)))
    df = df[~df.index.duplicated()]
    out["has_coexp"] = out.index.isin(df.index[np.isfinite(df["coexp_top20_mean"])]).astype(float)
    out = out.join(df, how="left")
    prof = F.load_phylo_profile()
    out["has_compara"] = out.index.isin(prof.index).astype(float)
    out = out.join(phylo_summary(prof), how="left")
    par = paralog_summary(F.load_paralogs())
    out = out.join(par, how="left")
    has = out["has_compara"] > 0
    for c in ("n_paralogs_wsp", "n_paralogs_other"):   # absent from the pair list but in Compara -> 0 paralogs
        out.loc[has & out[c].isna(), c] = 0.0
    e = _esm_current()
    pcs = esm_pca(e.to_numpy())
    df = pd.DataFrame(pcs, index=e.index, columns=[f"esm_pc{i}" for i in range(pcs.shape[1])])
    out["has_esm"] = out.index.isin(df.index).astype(float)
    out = out.join(df, how="left")
    return out


def current_names(orfs) -> np.ndarray:
    return F.to_current(pd.Index(np.asarray(orfs, dtype=str))).to_numpy(dtype=str)


def gene_matrix(orfs, current: bool = False) -> pd.DataFrame:
    """External gene features for `orfs` (rows in the given order; index = the given names). Genes unknown to every
    source get NaN features and 0 availability flags. `current=True` if names are already current SGD names."""
    orfs = np.asarray(orfs, dtype=str)
    cur = orfs if current else current_names(orfs)
    t = gene_table_current().reindex(cur)
    has = [c for c in t.columns if c.startswith("has_")]
    t[has] = t[has].fillna(0.0)
    t.index = pd.Index(orfs, name="orf")
    return t


def _pos(index: pd.Index, names: np.ndarray) -> np.ndarray:
    p = index.get_indexer(names)
    return p


def similarity(kind: str, orfs_a, orfs_b, current: bool = True) -> np.ndarray:
    """Gene-gene external similarity [na x nb] (float32, NaN where a gene is absent from the source)."""
    a = np.asarray(orfs_a, dtype=str) if current else current_names(orfs_a)
    b = np.asarray(orfs_b, dtype=str) if current else current_names(orfs_b)
    out = np.full((len(a), len(b)), np.nan, dtype=np.float32)
    if kind == "coexp":
        z, zg = _coexp()
        pa, pb = _pos(zg, a), _pos(zg, b)
        ia, ib = np.flatnonzero(pa >= 0), np.flatnonzero(pb >= 0)
        out[np.ix_(ia, ib)] = z[np.ix_(pa[ia], pb[ib])]
    elif kind == "esm":
        e = _esm_current()
        pa, pb = _pos(e.index, a), _pos(e.index, b)
        ia, ib = np.flatnonzero(pa >= 0), np.flatnonzero(pb >= 0)
        X = e.to_numpy(dtype=np.float64)
        Xa, Xb = X[pa[ia]], X[pb[ib]]
        Xa = Xa / np.linalg.norm(Xa, axis=1, keepdims=True)
        Xb = Xb / np.linalg.norm(Xb, axis=1, keepdims=True)
        out[np.ix_(ia, ib)] = (Xa @ Xb.T).astype(np.float32)
    elif kind == "apms":
        e = _apms_edges()
        nodes = set(F.load_table("ppi_nodes").index.astype(str))
        ina = np.array([x in nodes for x in a])
        inb = np.array([x in nodes for x in b])
        out[np.ix_(ina, inb)] = 0.0
        pos_a, pos_b = {}, {}
        for i, x in enumerate(a):
            pos_a.setdefault(x, []).append(i)
        for j, x in enumerate(b):
            pos_b.setdefault(x, []).append(j)
        for u, w in zip(np.r_[e["a"].values, e["b"].values], np.r_[e["b"].values, e["a"].values]):
            if u in pos_a and w in pos_b:
                for i in pos_a[u]:
                    out[i, pos_b[w]] = 1.0
    else:
        raise ValueError(kind)
    return out
