"""Dense matrices of the Costanzo 2016 NxN data.

HARNESS: `Bundle` / `load_bundle(dataset)` (DESIGN.md §11) — see the section "bundle" below.
LEGACY (pilot record, scripts/20_pilot_baselines.py at commit 403d4f0): `DenseNxN` / `load_dense`, documented here:

    E[i, j], P[i, j]  float32, query strain i x array strain j, NaN = not measured
    Y                 float32 stringent-negative label (ysl.eval.labels), NaN where not measured
Row order = sorted query strain ids (3,608); column order = sorted array strain ids (3,892; one strain per array
ORF in the cellmap release, asserted). Gene-level views use the canonical strain per query ORF
(`costanzo.canonical_query_strains`, recorded in the processed queries parquet as `canonical`).

Index helpers:
    orf_to_col      : array ORF -> column index
    strain_col      : per query strain, column of its ORF on the array (-1 if the ORF has no array strain)
    canon_row       : per query ORF (sorted `query_orfs`), row of its canonical strain
    array_smf       : array-strain single-mutant fitness at 30 C from the strain table (control measurements;
                      constant per array strain in the pair-wise file, checked in the pilot)
Cached as data/processed/costanzo2016_dense_<version>.npz (regenerable; ~120 MB).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from ysl.data import costanzo as C
from ysl.eval import labels as L


def dense_path(version: str = "cellmap"):
    return C.PROCESSED_DIR / f"costanzo2016_dense_{version}.npz"


@dataclass
class DenseNxN:
    E: np.ndarray
    P: np.ndarray
    query_strains: np.ndarray
    query_strain_orf: np.ndarray
    query_is_supp: np.ndarray
    query_canonical: np.ndarray
    array_strains: np.ndarray
    array_orfs: np.ndarray
    array_smf: np.ndarray
    _cache: dict = field(default_factory=dict, repr=False)

    @property
    def n_query_strains(self) -> int:
        return self.E.shape[0]

    @property
    def n_arrays(self) -> int:
        return self.E.shape[1]

    @cached_property
    def observed(self) -> np.ndarray:
        return ~np.isnan(self.E)

    @cached_property
    def Y(self) -> np.ndarray:
        y = L.stringent_negative(self.E, self.P).astype(np.float32)
        y[~self.observed] = np.nan
        return y

    @cached_property
    def orf_to_col(self) -> dict:
        return {o: j for j, o in enumerate(self.array_orfs)}

    @cached_property
    def strain_col(self) -> np.ndarray:
        return np.array([self.orf_to_col.get(o, -1) for o in self.query_strain_orf], dtype=np.int64)

    @cached_property
    def query_orfs(self) -> np.ndarray:
        return np.unique(self.query_strain_orf)

    @cached_property
    def canon_row(self) -> np.ndarray:
        """Row index of the canonical strain for each ORF in `query_orfs`."""
        rows = np.flatnonzero(self.query_canonical)
        orf = self.query_strain_orf[rows]
        order = np.argsort(orf)
        if not np.array_equal(orf[order], self.query_orfs):
            raise ValueError("canonical strains do not cover each query ORF exactly once")
        return rows[order]

    def orf_index(self, orfs) -> np.ndarray:
        """Positions of `orfs` in `query_orfs` (raises if any is unknown)."""
        orfs = np.asarray(orfs, dtype=str)
        pos = np.searchsorted(self.query_orfs, orfs)
        if (pos >= len(self.query_orfs)).any() or not np.array_equal(self.query_orfs[np.minimum(pos, len(self.query_orfs) - 1)], orfs):
            raise ValueError("unknown query ORF")
        return pos

    def strain_rows(self, orfs) -> np.ndarray:
        """All query-strain rows belonging to `orfs` (sorted)."""
        return np.flatnonzero(np.isin(self.query_strain_orf, np.asarray(orfs, dtype=str)))

    def save(self, path) -> None:
        np.savez(path, E=self.E, P=self.P, query_strains=self.query_strains, query_strain_orf=self.query_strain_orf,
                 query_is_supp=self.query_is_supp, query_canonical=self.query_canonical,
                 array_strains=self.array_strains, array_orfs=self.array_orfs, array_smf=self.array_smf)

    @classmethod
    def load_npz(cls, path) -> "DenseNxN":
        z = np.load(path, allow_pickle=False)
        return cls(**{k: z[k] for k in ("E", "P", "query_strains", "query_strain_orf", "query_is_supp",
                                         "query_canonical", "array_strains", "array_orfs", "array_smf")})

    def with_matrices(self, E: np.ndarray, P: np.ndarray) -> "DenseNxN":
        """Copy of the index structure with replaced E/P (used for leakage/poisoning checks)."""
        return DenseNxN(E, P, self.query_strains, self.query_strain_orf, self.query_is_supp, self.query_canonical,
                        self.array_strains, self.array_orfs, self.array_smf)


def build_dense(version: str = "cellmap") -> DenseNxN:
    """Build from the processed parquet with dictionary-encoded strain columns (no pandas joins)."""
    t = pq.read_table(C.processed_path(version), columns=["query_strain", "array_strain", "epsilon", "pvalue"],
                      read_dictionary=["query_strain", "array_strain"])
    codes = {}
    for side in ("query_strain", "array_strain"):
        col = t.column(side).unify_dictionaries().combine_chunks()
        d = np.asarray(col.dictionary.to_pylist(), dtype=str)
        order = np.argsort(d)
        rank = np.empty_like(order)
        rank[order] = np.arange(len(d))
        codes[side] = (rank[col.indices.to_numpy()], d[order])
    qi, qs = codes["query_strain"]
    ai, as_ = codes["array_strain"]
    E = np.full((len(qs), len(as_)), np.nan, dtype=np.float32)
    P = np.full_like(E, np.nan)
    if (np.bincount(qi * len(as_) + ai, minlength=E.size) > 1).any():
        raise ValueError("duplicate (query strain, array strain) pairs")
    E[qi, ai] = t.column("epsilon").to_numpy()
    P[qi, ai] = t.column("pvalue").to_numpy()
    del t

    qidx = pd.read_parquet(C.PROCESSED_DIR / f"costanzo2016_nxn_{version}_queries.parquet").set_index("query_strain")
    aidx = pd.read_parquet(C.PROCESSED_DIR / f"costanzo2016_nxn_{version}_arrays.parquet").set_index("array_strain")
    qidx, aidx = qidx.loc[qs], aidx.loc[as_]
    if aidx["array_orf"].duplicated().any():
        raise ValueError("an array ORF has more than one array strain; column map would be ambiguous")
    strains = pd.read_parquet(C.PROCESSED_DIR / f"costanzo2016_strains_{version}.parquet").set_index("strain")
    smf = strains.loc[as_, "smf_30"].to_numpy(dtype=np.float32)
    if np.isnan(smf).any():
        raise ValueError("missing array SMF in strain table")
    return DenseNxN(E=E, P=P, query_strains=qs, query_strain_orf=qidx["query_orf"].to_numpy(dtype=str),
                    query_is_supp=qidx["query_is_supp"].to_numpy(dtype=bool),
                    query_canonical=qidx["canonical"].to_numpy(dtype=bool),
                    array_strains=as_, array_orfs=aidx["array_orf"].to_numpy(dtype=str), array_smf=smf)


def load_dense(version: str = "cellmap", rebuild: bool = False) -> DenseNxN:
    """LEGACY pilot loader (CellMap). New code uses `load_bundle`."""
    p = dense_path(version)
    if p.exists() and not rebuild:
        return DenseNxN.load_npz(p)
    d = build_dense(version)
    d.save(p)
    return d


# ===================================================================================================== bundle
"""Bundle: the one cached, checksummed dense data object of the harness (DESIGN.md §11).

Strain-level matrices [n_query_strains x n_array_strains] (rows sorted by query strain ID, columns by array strain ID):
    E       float32 epsilon, NaN = not measured
    SP      float32 signed -log10 P = sign(eps) * -log10(max(P, P_FLOOR)) (P is one-sided in the direction of eps)
    Y       int8 stringent-negative label (ysl.eval.labels): 1 / 0, -1 = not measured   (observed = Y >= 0)
    DMF, DMF_sd  float32 double-mutant fitness and its SD
Gene level = the canonical strain row of each query ORF (`canon_row`; S1 has exactly one strain per ORF).
Array metadata: array_strains, array_orfs, array_orfs_current, array_smf (gene SMF from the array strain: the
pair-wise file's constant per-array value, else the strain table's 30 C SMF, else NaN).
Query metadata (strain level; used ONLY by the evaluator / view constructor, never exposed for held-out queries):
    strain_ids, strain_orf, strain_canonical, strain_is_supp, strain_temperature.
`all_query_orfs`: every query ORF of the dataset (names only; sets split sizes). `query_orfs`: ORFs with rows in
this (possibly blinded) bundle. ORF current names (`*_current`) are for external-feature joins only.

Blinding: `load_bundle(dataset)` drops every lockbox query row (DESIGN §4) unless `unblind=True`; the arrays are
read-only. The cache file data/processed/costanzo2016_bundle_<dataset>.npz has a JSON sidecar with its SHA-256 and
the SHA-256 of the raw and processed inputs; loading verifies the checksum and refuses a mismatch.
"""
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

P_FLOOR = 1e-300
BUNDLE_FORMAT = 1
_MATS = ("E", "SP", "Y", "DMF", "DMF_sd")
_STRAIN_META = ("strain_ids", "strain_orf", "strain_canonical", "strain_is_supp", "strain_temperature")
_ARRAY_META = ("array_strains", "array_orfs", "array_orfs_current", "array_smf")
_ORF_META = ("all_query_orfs", "all_query_orfs_current")


def sha256_file(path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def bundle_path(dataset: str = C.DEFAULT_DATASET):
    return C.PROCESSED_DIR / f"costanzo2016_bundle_{C.dataset_of(dataset)}.npz"


def bundle_meta_path(dataset: str = C.DEFAULT_DATASET):
    return bundle_path(dataset).with_suffix(".json")


def signed_neglog10p(eps, p) -> np.ndarray:
    eps = np.asarray(eps, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    return (np.sign(eps) * -np.log10(np.maximum(p, P_FLOOR))).astype(np.float32)


@dataclass
class Bundle:
    dataset: str
    E: np.ndarray
    SP: np.ndarray
    Y: np.ndarray
    DMF: np.ndarray
    DMF_sd: np.ndarray
    strain_ids: np.ndarray
    strain_orf: np.ndarray
    strain_canonical: np.ndarray
    strain_is_supp: np.ndarray
    strain_temperature: np.ndarray
    array_strains: np.ndarray
    array_orfs: np.ndarray
    array_orfs_current: np.ndarray
    array_smf: np.ndarray
    all_query_orfs: np.ndarray
    all_query_orfs_current: np.ndarray
    blinded_lockbox: np.ndarray = field(default_factory=lambda: np.array([], dtype=str))
    unblinded: bool = False
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ shapes / indices
    @property
    def n_arrays(self) -> int:
        return self.E.shape[1]

    @cached_property
    def observed(self) -> np.ndarray:
        return self.Y >= 0

    @cached_property
    def query_orfs(self) -> np.ndarray:
        """Sorted ORFs whose rows are present in this bundle (lockbox ORFs absent when blinded)."""
        return np.unique(self.strain_orf)

    @cached_property
    def canon_row(self) -> np.ndarray:
        rows = np.flatnonzero(self.strain_canonical)
        orf = self.strain_orf[rows]
        order = np.argsort(orf, kind="stable")
        if not np.array_equal(orf[order], self.query_orfs):
            raise ValueError("canonical strains do not cover each query ORF exactly once")
        return rows[order]

    @cached_property
    def orf_to_col(self) -> dict:
        if len(set(self.array_orfs)) != len(self.array_orfs):
            raise ValueError("an array ORF has more than one array strain; ORF -> column map is ambiguous")
        return {o: j for j, o in enumerate(self.array_orfs)}

    @cached_property
    def current_name(self) -> dict:
        """Original query/array ORF name -> current SGD name (for external joins only)."""
        d = dict(zip(self.all_query_orfs, self.all_query_orfs_current))
        d.update(zip(self.array_orfs, self.array_orfs_current))
        return d

    def gene_rows(self, orfs) -> np.ndarray:
        """Canonical-strain row of each ORF in `orfs` (raises if an ORF has no row, e.g. a blinded lockbox ORF)."""
        orfs = np.asarray(orfs, dtype=str)
        pos = np.searchsorted(self.query_orfs, orfs)
        ok = (pos < len(self.query_orfs))
        ok[ok] = self.query_orfs[pos[ok]] == orfs[ok]
        if not ok.all():
            missing = orfs[~ok]
            lb = np.isin(missing, self.blinded_lockbox)
            if lb.any():
                raise PermissionError(f"{int(lb.sum())} requested query ORFs are LOCKBOX ORFs whose rows were removed; "
                                      "pass unblind=True to load_bundle (final evaluation only)")
            raise KeyError(f"unknown query ORFs: {missing[:5].tolist()}")
        return self.canon_row[pos]

    def strain_rows(self, orfs) -> np.ndarray:
        """All strain rows of `orfs` (sorted by strain ID)."""
        self.gene_rows(orfs)  # existence / lockbox check
        return np.flatnonzero(np.isin(self.strain_orf, np.asarray(orfs, dtype=str)))

    # ------------------------------------------------------------------ derived copies
    def _replace(self, rows=None, **mats) -> "Bundle":
        """Copy with selected strain rows and/or replaced matrices (arrays are copied only when changed)."""
        kw = {k: getattr(self, k) for k in _MATS + _STRAIN_META}
        kw.update(mats)
        if rows is not None:
            kw = {k: v[rows] for k, v in kw.items()}
        return Bundle(dataset=self.dataset, **kw, **{k: getattr(self, k) for k in _ARRAY_META + _ORF_META},
                      blinded_lockbox=self.blinded_lockbox, unblinded=self.unblinded, meta=dict(self.meta))

    def blind(self, lockbox) -> "Bundle":
        """Remove every row of the lockbox query ORFs (DESIGN §4)."""
        lockbox = np.sort(np.asarray(lockbox, dtype=str))
        keep = ~np.isin(self.strain_orf, lockbox)
        b = self._replace(rows=keep)
        b.blinded_lockbox, b.unblinded = lockbox, False
        return b

    def with_matrices(self, **mats) -> "Bundle":
        """Copy with replaced strain-level matrices (poisoning / permutation tests). Y must stay consistent."""
        return self._replace(**mats)

    def subset_queries(self, orfs) -> "Bundle":
        """Restrict to the rows of `orfs` and set all_query_orfs = orfs (small real-data bundles for tests)."""
        orfs = np.sort(np.asarray(orfs, dtype=str))
        b = self._replace(rows=np.isin(self.strain_orf, orfs))
        cur = dict(zip(self.all_query_orfs, self.all_query_orfs_current))
        b.all_query_orfs = orfs
        b.all_query_orfs_current = np.array([cur[o] for o in orfs], dtype=str)
        b.blinded_lockbox = self.blinded_lockbox[np.isin(self.blinded_lockbox, orfs)]
        return b

    def freeze(self) -> "Bundle":
        for k in _MATS + _STRAIN_META + _ARRAY_META + _ORF_META:
            getattr(self, k).flags.writeable = False
        return self

    # ------------------------------------------------------------------ io
    def save(self, path) -> None:
        arrays = {k: getattr(self, k) for k in _MATS + _STRAIN_META + _ARRAY_META + _ORF_META}
        tmp = str(path) + ".tmp.npz"
        np.savez(tmp, **arrays)
        import os
        os.replace(tmp, path)

    @classmethod
    def load_npz(cls, path, dataset: str, meta: dict) -> "Bundle":
        z = np.load(path, allow_pickle=False)
        return cls(dataset=dataset, **{k: z[k] for k in _MATS + _STRAIN_META + _ARRAY_META + _ORF_META}, meta=meta)


def build_bundle(dataset: str = C.DEFAULT_DATASET) -> Bundle:
    """Build the strain-level bundle from the processed parquet (scripts/02_process_costanzo.py)."""
    ds = C.dataset_of(dataset)
    src = C.require(C.processed_path(ds))
    cols = ["query_strain", "array_strain", "epsilon", "pvalue", "dmf", "dmf_sd", "array_smf", "temperature"]
    t = pq.read_table(src, columns=cols, read_dictionary=["query_strain", "array_strain"])
    codes = {}
    for side in ("query_strain", "array_strain"):
        col = t.column(side).unify_dictionaries().combine_chunks()
        d = np.asarray(col.dictionary.to_pylist(), dtype=str)
        order = np.argsort(d)
        rank = np.empty_like(order)
        rank[order] = np.arange(len(d))
        codes[side] = (rank[col.indices.to_numpy()], d[order])
    qi, qs = codes["query_strain"]
    ai, as_ = codes["array_strain"]
    if (np.bincount(qi * len(as_) + ai, minlength=len(qs) * len(as_)) > 1).any():
        raise ValueError("duplicate (query strain, array strain) pairs")
    eps = t.column("epsilon").to_numpy()
    p = t.column("pvalue").to_numpy()
    if np.isnan(eps).any() or np.isnan(p).any():
        raise ValueError("NaN epsilon/P in the pair-wise file")
    shape = (len(qs), len(as_))
    E = np.full(shape, np.nan, dtype=np.float32)
    SP, DMF, DSD = E.copy(), E.copy(), E.copy()
    Y = np.full(shape, -1, dtype=np.int8)
    E[qi, ai] = eps
    SP[qi, ai] = signed_neglog10p(eps, p)
    Y[qi, ai] = L.stringent_negative(eps, p).astype(np.int8)  # label from the float64 source values
    DMF[qi, ai] = t.column("dmf").to_numpy()
    DSD[qi, ai] = t.column("dmf_sd").to_numpy()
    # gene SMF of array strains: the pair-wise file's 30 C value (constant per array strain at 30 C, asserted; CellMap
    # 26 C screens carry a 26 C SMF); strain-table smf_30 fallback. S1's pair-wise SMF differs from the strain table
    # (r = 0.992; the S1 eps were computed with it), so the pair-wise value is preferred.
    smf_pw = t.column("array_smf").to_numpy()
    smf = np.full(len(as_), np.nan)
    ok = ~np.isnan(smf_pw) & (t.column("temperature").to_numpy() == 30)
    lo = np.full(len(as_), np.inf)
    hi = np.full(len(as_), -np.inf)
    np.minimum.at(lo, ai[ok], smf_pw[ok])
    np.maximum.at(hi, ai[ok], smf_pw[ok])
    has = np.isfinite(lo)
    if (hi[has] - lo[has] > 1e-9).any():
        raise ValueError("array SMF is not constant per array strain")
    smf[has] = lo[has]
    del t
    strains = pd.read_parquet(C.require(C.strains_path(ds))).set_index("strain")
    st = strains["smf_30"].reindex(as_).to_numpy(dtype=float)
    both = has & ~np.isnan(st)
    smf_source = {"pairwise_30C": int(has.sum()), "strain_table_fallback": int((~has & ~np.isnan(st)).sum()),
                  "missing": int((~has & np.isnan(st)).sum()),
                  "pairwise_vs_strain_table_max_abs_diff": float(np.max(np.abs(smf[both] - st[both]))) if both.any()
                  else None}
    smf[~has] = st[~has]

    qidx = pd.read_parquet(C.require(C.queries_path(ds))).set_index("query_strain").loc[qs]
    aidx = pd.read_parquet(C.require(C.arrays_path(ds))).set_index("array_strain").loc[as_]
    if "query_orf_current" not in qidx or "array_orf_current" not in aidx:
        raise ValueError("processed indices lack current ORF names; re-run scripts/02_process_costanzo.py")
    orfs = qidx.groupby("query_orf")["query_orf_current"].first().sort_index()
    b = Bundle(dataset=ds, E=E, SP=SP, Y=Y, DMF=DMF, DMF_sd=DSD,
               strain_ids=qs, strain_orf=qidx["query_orf"].to_numpy(dtype=str),
               strain_canonical=qidx["canonical"].to_numpy(dtype=bool),
               strain_is_supp=qidx["query_is_supp"].to_numpy(dtype=bool),
               strain_temperature=qidx["temperature"].to_numpy(dtype=np.int16),
               array_strains=as_, array_orfs=aidx["array_orf"].to_numpy(dtype=str),
               array_orfs_current=aidx["array_orf_current"].to_numpy(dtype=str),
               array_smf=smf.astype(np.float32),
               all_query_orfs=orfs.index.to_numpy(dtype=str), all_query_orfs_current=orfs.to_numpy(dtype=str))
    b.canon_row  # validate
    b.orf_to_col
    b.meta = {"format": BUNDLE_FORMAT, "dataset": ds, "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
              "raw_nxn_sha256": sha256_file(C.require(C.nxn_path(ds), "raw NxN file")),
              "processed_parquet_sha256": sha256_file(src),
              "shape": list(shape), "n_observed": int((Y >= 0).sum()), "n_positive": int((Y == 1).sum()),
              "n_query_orfs": int(len(b.all_query_orfs)), "array_smf_source": smf_source, "p_floor": P_FLOOR}
    return b


def write_bundle(b: Bundle) -> dict:
    p = bundle_path(b.dataset)
    b.save(p)
    meta = dict(b.meta, bundle_sha256=sha256_file(p))
    bundle_meta_path(b.dataset).write_text(json.dumps(meta, indent=1))
    return meta


_CACHE: dict = {}


def load_bundle(dataset: str = C.DEFAULT_DATASET, unblind: bool = False, verify: bool = True,
                cache: bool = True) -> Bundle:
    """Load the checksummed bundle. Lockbox query rows are REMOVED unless unblind=True (final evaluation only).
    Fails (never skips/rebuilds silently) if the cache is missing or its checksum does not match."""
    ds = C.dataset_of(dataset)
    key = (ds, bool(unblind))
    if cache and key in _CACHE:
        return _CACHE[key]
    p, mp = C.require(bundle_path(ds), "dense bundle"), C.require(bundle_meta_path(ds), "dense bundle sidecar")
    meta = json.loads(mp.read_text())
    if verify:
        got = sha256_file(p)
        if got != meta["bundle_sha256"]:
            raise ValueError(f"bundle checksum mismatch for {p}: {got} != {meta['bundle_sha256']}; rebuild with "
                             "scripts/02b_build_bundles.py")
    b = Bundle.load_npz(p, ds, meta)
    if unblind:
        b.unblinded = True
    else:
        from ysl.eval.splits import read_lockbox
        b = b.blind(read_lockbox(ds))
    b.freeze()
    if cache:
        _CACHE[key] = b
    return b
