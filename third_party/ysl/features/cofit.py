"""Chemogenomic co-fitness (tier 2, DESIGN.md §7): Hoepfner et al. 2014 HIP/HOP compound-sensitivity scores
(Dryad doi:10.5061/dryad.v5m8v, CC0) -> gene x gene co-fitness similarity (~ DepMap co-essentiality).

No GI data, GO or literature enters anything here: inputs are deletion-pool fitness scores under compounds only.

Fixed recipe (`cofitness`), chosen before any model result:
  1. Per compound (column), robust scaling over genes: z = (s - median) / (1.4826 * MAD), computed on observed genes
     only; compounds with < `min_genes` observed genes or MAD = 0 are dropped. Puts compounds with very different
     potency / dynamic range on one scale, so co-fitness is not dominated by a few strong compounds.
  2. Clip z to [-clip, clip] (clip = 10): keeps strong (real) sensitivities as the dominant signal but limits the
     leverage of single extreme values on a Pearson correlation.
  3. Missing values: genes observed in < `min_frac` (50%) of the kept compounds are dropped (no profile); remaining
     missing entries are set to 0 = the compound median, i.e. "no differential sensitivity". This is conservative:
     it can only shrink |r| towards 0, never create co-fitness, and keeps the similarity a proper Gram matrix (PSD
     before thresholding). Pairwise-complete correlation was not used (noisier for low-overlap pairs, not PSD).
  4. Co-fitness = Pearson correlation of the gene profiles over the kept compounds (diagonal set to NaN).
Merging pools (HOP homozygous = non-essential genes; HIP heterozygous = all genes, the only profile for essential
genes): a gene's profile is its HOP profile if it has one, else its HIP profile (`merge_pools`). Each pool is scaled
separately (step 1) before merging, so profiles share the robust-z scale; compounds are matched by identifier.

Store (`write_store` / `load_store` / `block`): data/processed/features/hoepfner/ with cofit_genes.txt (current
systematic names), cofit_r.npy (float32 [G x G], NaN on the diagonal) and manifest.json (SHA-256 of both).
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
STORE = ROOT / "data" / "processed" / "features" / "hoepfner"
RECIPE = {"min_genes": 100, "clip": 10.0, "min_frac": 0.5, "mad_const": 1.4826, "impute": "zero (compound median)",
          "similarity": "pearson", "pool_merge": "HOP profile if present else HIP"}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------- recipe
def robust_scale_columns(S: np.ndarray, min_genes: int = RECIPE["min_genes"],
                         mad_const: float = RECIPE["mad_const"]) -> tuple[np.ndarray, np.ndarray]:
    """Per-column robust z over observed rows. Returns (Z with NaN kept, kept-column boolean mask)."""
    S = np.asarray(S, dtype=np.float64)
    n_obs = (~np.isnan(S)).sum(0)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(S, 0)
            mad = np.nanmedian(np.abs(S - med[None, :]), 0) * mad_const
    keep = (n_obs >= min_genes) & np.isfinite(mad) & (mad > 0)
    Z = (S[:, keep] - med[keep][None, :]) / mad[keep][None, :]
    return Z, keep


def profiles(Z: np.ndarray, clip: float = RECIPE["clip"], min_frac: float = RECIPE["min_frac"]
             ) -> tuple[np.ndarray, np.ndarray]:
    """Clip, drop genes observed in < min_frac of columns, impute missing as 0. Returns (P [kept genes x C], kept-gene
    mask)."""
    frac = (~np.isnan(Z)).mean(1) if Z.shape[1] else np.zeros(len(Z))
    keep = frac >= min_frac
    P = np.clip(np.nan_to_num(Z[keep], nan=0.0), -clip, clip)
    return P, keep


def pearson_rows(P: np.ndarray) -> np.ndarray:
    """Pearson correlation between rows (0 for constant rows); diagonal NaN. float32."""
    X = P.astype(np.float64)
    X = X - X.mean(1, keepdims=True)
    nrm = np.linalg.norm(X, axis=1)
    ok = nrm > 1e-12
    X[ok] /= nrm[ok, None]
    X[~ok] = 0.0
    R = (X @ X.T).astype(np.float32)
    np.fill_diagonal(R, np.nan)
    return R


def cofitness(S: np.ndarray, **kw) -> tuple[np.ndarray, np.ndarray, dict]:
    """Full recipe on one gene x compound score matrix. Returns (R [kept x kept], kept-gene mask, stats)."""
    Z, kc = robust_scale_columns(S, **{k: kw[k] for k in ("min_genes", "mad_const") if k in kw})
    P, kg = profiles(Z, **{k: kw[k] for k in ("clip", "min_frac") if k in kw})
    return pearson_rows(P), kg, {"n_genes_in": int(len(S)), "n_compounds_in": int(S.shape[1]),
                                 "n_compounds_kept": int(kc.sum()), "n_genes_kept": int(kg.sum()),
                                 "frac_missing_kept": float(np.isnan(Z[kg]).mean()) if kg.any() else float("nan")}


def merge_pools(pools: dict[str, pd.DataFrame], order=("HOP", "HIP"), **kw) -> tuple[pd.DataFrame, dict]:
    """pools: name -> gene x compound DataFrame (index = current ORF, columns = compound ids). Each pool is robust-
    scaled separately; a gene takes its profile from the first pool in `order` where it has one (>= min_frac
    observed). Returns (profile DataFrame [genes x union of kept compounds], per-pool stats incl. source of each
    gene). Missing entries (incl. compounds absent from a gene's pool) -> 0 after merging."""
    clip = kw.get("clip", RECIPE["clip"])
    min_frac = kw.get("min_frac", RECIPE["min_frac"])
    rows, stats, src = [], {}, {}
    for name in order:
        if name not in pools:
            continue
        df = pools[name]
        Z, kc = robust_scale_columns(df.to_numpy(dtype=np.float64), **{k: kw[k] for k in ("min_genes", "mad_const")
                                                                         if k in kw})
        cols = df.columns[kc]
        frac = (~np.isnan(Z)).mean(1) if Z.shape[1] else np.zeros(len(Z))
        ok = (frac >= min_frac) & ~df.index.isin(list(src))
        z = pd.DataFrame(np.clip(Z[ok], -clip, clip), index=df.index[ok], columns=cols)
        rows.append(z)
        for g in z.index:
            src[g] = name
        stats[name] = {"n_genes_in": int(len(df)), "n_compounds_in": int(df.shape[1]),
                       "n_compounds_kept": int(kc.sum()), "n_genes_used": int(ok.sum()),
                       "frac_missing_used": float(np.isnan(Z[ok]).mean()) if ok.any() else float("nan")}
    prof = pd.concat(rows, axis=0, join="outer").fillna(0.0)
    stats["source"] = src
    return prof, stats


# ----------------------------------------------------------------------------------------------- store
def write_store(out_dir: Path, genes: pd.Index, R: np.ndarray, extra: dict | None = None) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if len(set(genes)) != len(genes) or R.shape != (len(genes), len(genes)):
        raise ValueError("genes must be unique and match R")
    if not np.array_equal(R, R.T, equal_nan=True):
        raise ValueError("R not symmetric")
    (out_dir / "cofit_genes.txt").write_text("\n".join(map(str, genes)) + "\n")
    np.save(out_dir / "cofit_r.npy", R.astype(np.float32))
    man = {"genes_sha256": sha256_file(out_dir / "cofit_genes.txt"), "r_sha256": sha256_file(out_dir / "cofit_r.npy"),
           "n_genes": int(len(genes)), "recipe": RECIPE} | (extra or {})
    (out_dir / "manifest.json").write_text(json.dumps(man, indent=1))
    load_store.cache_clear()
    return man


@lru_cache(maxsize=4)
def load_store(out_dir: str = str(STORE)) -> tuple[pd.Index, np.ndarray]:
    """(gene index, R memmap), SHA-256 verified once per process per directory."""
    d = Path(out_dir)
    if not (d / "manifest.json").exists():
        raise FileNotFoundError(f"{d}/manifest.json missing; build with scripts/14_process_hoepfner.py")
    man = json.loads((d / "manifest.json").read_text())
    if sha256_file(d / "cofit_genes.txt") != man["genes_sha256"] or sha256_file(d / "cofit_r.npy") != man["r_sha256"]:
        raise ValueError(f"co-fitness store checksum mismatch in {d}")
    genes = pd.Index([g for g in (d / "cofit_genes.txt").read_text().split("\n") if g], name="orf")
    return genes, np.load(d / "cofit_r.npy", mmap_mode="r")


def block(rows, cols, out_dir: str = str(STORE)) -> np.ndarray:
    """float32 [len(rows) x len(cols)] co-fitness for (current) ORF names; NaN where a gene has no profile or on the
    diagonal (same gene)."""
    genes, R = load_store(str(out_dir))
    ri = genes.get_indexer(pd.Index(np.asarray(rows, dtype=str)))
    ci = genes.get_indexer(pd.Index(np.asarray(cols, dtype=str)))
    out = np.full((len(ri), len(ci)), np.nan, dtype=np.float32)
    r_ok, c_ok = np.flatnonzero(ri >= 0), np.flatnonzero(ci >= 0)
    if len(r_ok) and len(c_ok):
        out[np.ix_(r_ok, c_ok)] = np.asarray(R[ri[r_ok]][:, ci[c_ok]])
    return out
