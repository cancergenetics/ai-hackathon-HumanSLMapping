"""Loader for the Costanzo et al. 2016 (Science 353:aaf1420) non-essential x non-essential (NxN) SGA data.

Raw input (see DATA.md, downloaded by scripts/01_download_costanzo.py):
    data/raw/costanzo2016/science_S1_subset_pairwise/SGA_NxN.txt (dataset "s1", PRIMARY: Science Data File S1 NxN)
    data/raw/costanzo2016/cellmap_pairwise/SGA_NxN.txt          (dataset "cellmap", SECONDARY: TheCellMap release)
    data/raw/costanzo2016/*/strain_ids_and_single_mutant_fitness.xlsx

Dataset names (DESIGN.md §1): "s1" (default) and "cellmap". On disk the S1 files keep the historical version stem
"science_subset" (`DATASETS` maps dataset -> version stem); every public loader accepts either name.

The pair-wise file has 11 tab-separated columns:
    Query Strain ID | Query allele name | Array Strain ID | Array allele name | Arraytype/Temp |
    Genetic interaction score (eps) | P-value | Query SMF | Array SMF | Double mutant fitness | DMF std
Strain IDs are "<systematic ORF>_<collection><number>", e.g. YAL002W_sn273, YBL001C_dma94, YBR126C_S26.

Tidy table (one row per tested (query strain, array strain) pair; no aggregation is done here):
    query_strain, query_orf, query_allele, query_gene, query_strain_type, query_is_supp,
    array_strain, array_orf, array_allele, array_gene, array_strain_type, array_is_supp,
    temperature (26/30, the screen temperature; one per query strain), epsilon, pvalue,
    query_smf, array_smf, dmf, dmf_sd
Gene names ("*_gene") are the allele name with any "-supp..." suffix removed, upper-cased; they are the
names supplied by Costanzo et al., not re-mapped to current SGD names. Systematic ORF names are the keys.

Conventions/assumptions (documented in results/data_characterisation/SUMMARY.md):
  * Strain type "sn" = SGA non-essential deletion query; "S" = deletion query strain carrying a characterised
    spontaneous suppressor mutation (allele names "<gene>-suppN"); "dma" = deletion mutant array strain.
    `*_is_supp` is True when the allele name carries a "-supp" suffix (a few "sn" strains do too).
  * Query SMF is missing (NaN) for some query strains in the source; we keep NaN (never impute here).
  * Gene-level views pick ONE canonical query strain per query ORF (see `canonical_query_strains`).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
RAW_DIR = ROOT / "data" / "raw" / "costanzo2016"
PROCESSED_DIR = ROOT / "data" / "processed"

VERSIONS = {
    "cellmap": RAW_DIR / "cellmap_pairwise",
    "science_subset": RAW_DIR / "science_S1_subset_pairwise",
}
DATASETS = {"s1": "science_subset", "cellmap": "cellmap"}  # dataset name -> on-disk version stem
DEFAULT_DATASET = "s1"


def version_of(dataset: str = DEFAULT_DATASET) -> str:
    """Resolve a dataset name ("s1", "cellmap") or a version stem ("science_subset", "cellmap") to the stem."""
    if dataset in DATASETS:
        return DATASETS[dataset]
    if dataset in VERSIONS:
        return dataset
    raise ValueError(f"unknown dataset {dataset!r}; expected one of {sorted(DATASETS)}")


def dataset_of(version: str) -> str:
    v = version_of(version)
    return {stem: name for name, stem in DATASETS.items()}[v]
NXN_FILE = "SGA_NxN.txt"
STRAIN_FILE = "strain_ids_and_single_mutant_fitness.xlsx"

RAW_COLUMNS = [
    "Query Strain ID",
    "Query allele name",
    "Array Strain ID",
    "Array allele name",
    "Arraytype/Temp",
    "Genetic interaction score (ε)",
    "P-value",
    "Query single mutant fitness (SMF)",
    "Array SMF",
    "Double mutant fitness",
    "Double mutant fitness standard deviation",
]
SHORT = ["query_strain", "query_allele", "array_strain", "array_allele", "array_type_temp",
         "epsilon", "pvalue", "query_smf", "array_smf", "dmf", "dmf_sd"]

STRAIN_RE = re.compile(r"^(?P<orf>[^_]+)_(?P<type>[A-Za-z]+)(?P<num>\d+)(?:_(?P<extra>[A-Za-z0-9]+))?$")
# (the Science-subset file contains one id with an extra suffix: YMR186W_sn2174_T30)
SUPP_RE = r"-supp\d*[a-z]?$"

TIDY_COLUMNS = [
    "query_strain", "query_orf", "query_allele", "query_gene", "query_strain_type", "query_is_supp",
    "array_strain", "array_orf", "array_allele", "array_gene", "array_strain_type", "array_is_supp",
    "temperature", "epsilon", "pvalue", "query_smf", "array_smf", "dmf", "dmf_sd",
]


def nxn_path(version: str = DEFAULT_DATASET) -> Path:
    return VERSIONS[version_of(version)] / NXN_FILE


def processed_path(version: str = DEFAULT_DATASET) -> Path:
    return PROCESSED_DIR / f"costanzo2016_nxn_{version_of(version)}.parquet"


def queries_path(version: str = DEFAULT_DATASET) -> Path:
    return PROCESSED_DIR / f"costanzo2016_nxn_{version_of(version)}_queries.parquet"


def arrays_path(version: str = DEFAULT_DATASET) -> Path:
    return PROCESSED_DIR / f"costanzo2016_nxn_{version_of(version)}_arrays.parquet"


def strains_path(version: str = DEFAULT_DATASET) -> Path:
    return PROCESSED_DIR / f"costanzo2016_strains_{version_of(version)}.parquet"


def require(path: Path, what: str = "processed Costanzo data") -> Path:
    """Fail loudly (never skip) when a required data file is missing."""
    if not Path(path).exists():
        raise FileNotFoundError(
            f"{what} missing: {path}. Regenerate with scripts/01_download_costanzo.py, scripts/02_process_costanzo.py "
            f"and scripts/02b_build_bundles.py.")
    return Path(path)


def read_nxn_raw(path: Path | str) -> pd.DataFrame:
    """Read the pair-wise NxN file, validating the header exactly. Returns SHORT-named columns."""
    path = Path(path)
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").rstrip("\r").split("\t")
    if header != RAW_COLUMNS:
        raise ValueError(f"Unexpected header in {path}: {header}")
    dtypes = {c: "string" for c in RAW_COLUMNS[:5]}
    dtypes.update({c: "float64" for c in RAW_COLUMNS[5:]})
    df = pd.read_csv(path, sep="\t", dtype=dtypes, header=0, engine="c", na_values=["", "NaN", "nan", "NA"],
                     keep_default_na=False)
    df.columns = SHORT
    return df


def split_strain_ids(strains: pd.Series) -> pd.DataFrame:
    """Split '<ORF>_<type><num>' into orf / type. Raises if any id does not match."""
    u = pd.Series(strains.unique()).astype(str)
    parts = u.str.extract(STRAIN_RE)
    bad = parts["orf"].isna()
    if bad.any():
        raise ValueError(f"Unparseable strain ids: {u[bad].tolist()[:10]}")
    parts.index = u.values
    return parts.loc[strains.astype(str).values].reset_index(drop=True)


def allele_to_gene(allele: pd.Series) -> pd.Series:
    return allele.astype(str).str.replace(SUPP_RE, "", regex=True).str.upper()


def tidy_nxn(raw: pd.DataFrame) -> pd.DataFrame:
    """Convert SHORT-named raw table to the tidy schema (categoricals for string columns)."""
    out = pd.DataFrame(index=raw.index)
    for side in ("query", "array"):
        strain = raw[f"{side}_strain"].astype(str)
        parts = split_strain_ids(strain)
        out[f"{side}_strain"] = strain.values
        out[f"{side}_orf"] = parts["orf"].values
        out[f"{side}_allele"] = raw[f"{side}_allele"].astype(str).values
        out[f"{side}_gene"] = allele_to_gene(raw[f"{side}_allele"]).values
        out[f"{side}_strain_type"] = parts["type"].values
        out[f"{side}_is_supp"] = raw[f"{side}_allele"].astype(str).str.contains(SUPP_RE, regex=True).values
    att = raw["array_type_temp"].astype(str)
    m = att.str.extract(r"^(?P<arr>[A-Z]+)(?P<t>\d+)$")
    if m.isna().any().any() or set(m["arr"].unique()) != {"DMA"}:
        raise ValueError(f"Unexpected Arraytype/Temp values: {att.unique()[:10]}")
    out["temperature"] = m["t"].astype("int8").values
    for c in ("epsilon", "pvalue", "query_smf", "array_smf", "dmf", "dmf_sd"):
        out[c] = raw[c].astype("float64").values
    for c in out.columns:
        if out[c].dtype == object or str(out[c].dtype) == "string":
            out[c] = out[c].astype("category")
    return out[TIDY_COLUMNS]


def load_nxn(version: str = DEFAULT_DATASET, from_processed: bool = True) -> pd.DataFrame:
    """Tidy NxN table. Uses the processed parquet if present (scripts/02_process_costanzo.py), else parses raw.
    NOTE: contains every query row (held-out queries included). Only `ysl.data`, `ysl.eval` and characterisation
    scripts may call this; models receive a `ScreenView`."""
    p = processed_path(version)
    if from_processed and p.exists():
        return pd.read_parquet(p)
    return tidy_nxn(read_nxn_raw(nxn_path(version)))


def load_strain_table(version: str = DEFAULT_DATASET) -> pd.DataFrame:
    """Strain ids and single-mutant fitness (all collections: sn, S, dma, tsq, tsa, damp)."""
    df = pd.read_excel(VERSIONS[version_of(version)] / STRAIN_FILE, dtype={"Strain ID": str, "Systematic gene name": str,
                                                                "Allele/Gene name": str})
    df.columns = ["strain", "orf", "allele", "smf_26", "smf_26_sd", "smf_30", "smf_30_sd"]
    parts = split_strain_ids(df["strain"])
    if not (parts["orf"].values == df["orf"].values).all():
        raise ValueError("strain id ORF prefix disagrees with 'Systematic gene name'")
    df["strain_type"] = parts["type"].values
    return df


def canonical_query_strains(tidy: pd.DataFrame) -> pd.Series:
    """One query strain per query ORF. Deterministic preference order:
        1. allele is not a suppressor-carrying allele ('-supp' suffix)
        2. strain type 'sn' over 'S'
        3. more arrays tested
        4. lexicographically smallest strain id
    Returns a Series indexed by query_orf -> query_strain.
    """
    q = (tidy.groupby("query_strain", observed=True)
         .agg(query_orf=("query_orf", "first"), is_supp=("query_is_supp", "first"),
              stype=("query_strain_type", "first"), n=("array_strain", "size"))
         .reset_index())
    q["is_sn"] = q["stype"].astype(str) == "sn"
    q = q.sort_values(["query_orf", "is_supp", "is_sn", "n", "query_strain"],
                      ascending=[True, True, False, False, True])
    q["query_orf"] = q["query_orf"].astype(str)
    q["query_strain"] = q["query_strain"].astype(str)
    return q.drop_duplicates("query_orf").set_index("query_orf")["query_strain"]


def query_index(tidy: pd.DataFrame) -> pd.DataFrame:
    """Per-query-strain index: ORF, allele, gene, type, temperature, SMF, #arrays tested, canonical flag."""
    q = (tidy.groupby("query_strain", observed=True)
         .agg(query_orf=("query_orf", "first"), query_allele=("query_allele", "first"),
              query_gene=("query_gene", "first"), query_strain_type=("query_strain_type", "first"),
              query_is_supp=("query_is_supp", "first"), temperature=("temperature", "first"),
              query_smf=("query_smf", "first"), n_arrays=("array_strain", "size"))
         .reset_index())
    for c in ("query_strain", "query_orf", "query_allele", "query_gene", "query_strain_type"):
        q[c] = q[c].astype(str)
    canon = canonical_query_strains(tidy)
    q["canonical"] = q["query_strain"].isin(set(canon.values))
    q["n_strains_for_orf"] = q.groupby("query_orf")["query_strain"].transform("size")
    return q.sort_values("query_strain").reset_index(drop=True)


def add_current_names(df: pd.DataFrame, orf_col: str) -> pd.DataFrame:
    """Add `<orf_col>_current`: the current SGD systematic name (ysl.data.features.to_current; identity when the
    name is not retired). Used ONLY for joins to external features; original names/strain IDs stay the keys for
    splitting and scoring (several retired ORFs map to one current name, so the mapping is many-to-one)."""
    from ysl.data.features import to_current
    m = to_current(sorted(set(df[orf_col].astype(str))))
    df[f"{orf_col}_current"] = df[orf_col].astype(str).map(m).astype(str).values
    return df


def array_index(tidy: pd.DataFrame) -> pd.DataFrame:
    """Per-array-strain index."""
    a = (tidy.groupby("array_strain", observed=True)
         .agg(array_orf=("array_orf", "first"), array_allele=("array_allele", "first"),
              array_gene=("array_gene", "first"), array_strain_type=("array_strain_type", "first"),
              array_is_supp=("array_is_supp", "first"), n_queries=("query_strain", "size"))
         .reset_index())
    for c in ("array_strain", "array_orf", "array_allele", "array_gene", "array_strain_type"):
        a[c] = a[c].astype(str)
    a["n_strains_for_orf"] = a.groupby("array_orf")["array_strain"].transform("size")
    return a.sort_values("array_strain").reset_index(drop=True)


def gene_level(tidy: pd.DataFrame) -> pd.DataFrame:
    """Restrict to canonical query strains -> at most one row per (query_orf, array_orf) when each array ORF
    has one array strain (true for the Costanzo 2016 NxN DMA). Raises if duplicates remain."""
    canon = set(canonical_query_strains(tidy).values)
    g = tidy[tidy["query_strain"].astype(str).isin(canon)].copy()
    if g.duplicated(["query_orf", "array_orf"]).any():
        raise ValueError("duplicate (query_orf, array_orf) after canonical-strain selection")
    for c in g.select_dtypes("category").columns:
        g[c] = g[c].cat.remove_unused_categories()
    return g.reset_index(drop=True)


def read_matrix_cdt(path: Path | str) -> pd.DataFrame:
    """Read Data File S2 clustered .cdt (query rows x array columns, eps rounded to 3 dp) into long form
    with query_strain / array_strain ids reconstructed as '<ORF>_<id>'. Used only for cross-checking."""
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    head = raw.iloc[:6]
    if list(head.iloc[:, 0]) != ["GID", "AID", "ORF", "Name", "Label", "EWEIGHT"]:
        raise ValueError("unexpected cdt layout")
    arr_ids = head.iloc[0, 6:].values
    arr_orfs = head.iloc[2, 6:].values
    body = raw.iloc[6:]
    q_strain = (body.iloc[:, 2] + "_" + body.iloc[:, 1]).values
    vals = body.iloc[:, 6:].replace("", np.nan).astype(float).values
    a_strain = np.array([f"{o}_{i}" for o, i in zip(arr_orfs, arr_ids)])
    qi, ai = np.nonzero(~np.isnan(vals))
    return pd.DataFrame({"query_strain": q_strain[qi], "array_strain": a_strain[ai], "epsilon": vals[qi, ai]})
