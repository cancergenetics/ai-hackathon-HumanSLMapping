"""Sample metadata + genetic-perturbation screening for the COXPRESdb S. cerevisiae expression compendia.

Why: DESIGN.md §7 excludes deletion-mutant expression profiles (Hughes 2000, Kemmeren 2014 "deleteome", and any
sample that is a gene-deletion strain profile) from every expression compendium. COXPRESdb Sce-r (RNA-seq, SRA/ENA)
and Sce-m (Affymetrix, ArrayExpress) aggregate public studies, many of which profile deletion / point / ts / degron /
over-expression mutants next to wild type. Neither Kemmeren 2014 nor Hughes 2000 is in these compendia (two-colour
arrays are not used by COXPRESdb), but ~half of the studies contain mutant samples.

Inputs (data/raw/coxpresdb/, scripts/10_download_features.py):
    <tag>.samples.tsv                     COXPRESdb sample list: study, sample id, study title[, run description]
    sample_metadata/ena/<study>.tsv       ENA read_run report per RNA-seq study (sample titles, strain, ...)
    sample_metadata/arrayexpress/<E-...>.tsv   ArrayExpress SDRF per microarray study
    config/expression_sample_filter.yaml  marker genes treated as background genotype + manual study decisions

Policy (strict, study level): a study is EXCLUDED if its title or the metadata text of ANY of its samples shows a
gene-specific genetic perturbation (deletion/knockout/null, point or ts mutant, degron/anchor-away/CRISPRi depletion,
over-expression of a yeast gene, '::' constructs), after removing standard auxotrophic/background markers
(his3, leu2, ura3, met15, ...). All samples of an excluded study are dropped (including its wild-type samples), so
misclassifying an individual sample cannot let a deletion profile through. Manual overrides (with reasons) live in the
YAML config and are applied after the automatic screen. A per-sample flag is also kept for audit.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import yaml

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
RAW = ROOT / "data" / "raw" / "coxpresdb"
CONFIG = ROOT / "config" / "expression_sample_filter.yaml"
PLATFORMS = {"rnaseq": "Sce-r.c3-0", "microarray": "Sce-m.c4-0"}

SDRF_SKIP = re.compile(r"^(Protocol REF|Term Source REF|Term Accession|Comment \[(ArrayExpress|Derived)|"
                       r"Derived Array|Array Design|Material Type|Label$|Technology Type|Normalization Name|"
                       r"Unit|Performer|Provider|Parameter Value|Scan Name|Comment \[Array Design)", re.I)
ORGANISM = re.compile(r"organism", re.I)
ENA_TEXT_COLS = ["sample_title", "sample_alias", "experiment_title", "strain", "sub_strain", "host_genotype",
                 "sample_description"]

# Word boundaries that treat '_' as a separator (sample names like 'W303_SET2KO_rep1', 'yJC_rpb1-1_RZ_0').
_B, _E = r"(?<![A-Za-z0-9])", r"(?![A-Za-z0-9])"
DEL_RE = re.compile(
    r"(Δ|∆|△|\?(?=[\s_,;)/]|$)|" + _B + r"delta" + _E + r"|delta[-_ ]?[a-z]{3}\d|[a-z]{3}\d+[-_ ]?delta|"
    + _B + r"del(etion|eted|ete|etant)?s?" + _E + r"|knock[- ]?(out|down|ed out)|" + _B + r"KOs?" + _E + r"|"
    r"[a-z0-9]KO" + _E + r"|" + _B + r"null" + _E + r"|::|-/-|" + _B + r"lacking" + _E + r"|loss[- ]of[- ]function|"
    r"disrupt|" + _B + r"(gene|strain|mutant)s? deficien|(?<!glucose )(?<!nitrogen )(?<!carbon )(?<!nutrient )"
    r"(?<!glucose-)(?<!oxygen )deplet|inactivat)", re.I)
MUT_RE = re.compile(
    r"(mutant|mutation|mutated|" + _B + r"alleles?" + _E + r"|temperature[- ]sensitive|" + _B + r"ts" + _E + r"|"
    r"degron|anchor[- ]?away|auxin|" + _B + r"AID" + _E + r"|over-?express|" + _B + r"OE" + _E + r"|CRISPR|"
    r"tet[- ]?off|" + _B + r"tetO|-td" + _E + r"|analog[- ]sensitive|knock[- ]?(in|sideways)|point mutant)", re.I)
# Case-sensitive patterns: yeast mutant-allele notation ('cdc13-1', 'rpb1-1'), amino-acid substitutions ('K1246R',
# 'F221A'), histone residue mutants ('H3K36A', 'H3E50A', 'H4T80A').
POINT_RE = re.compile(
    _B + r"[a-z]{3}\d+-\d+" + _E + r"|" + _B + r"(?!S288C)[ACDEFGHIKLMNPQRSTVWY]\d{1,4}[ACDEFGHIKLMNPQRSTVWY]" + _E
    + r"|H2A|H2B|" + _B + r"H[34][A-Z]?\d{1,3}[A-Z]?")
# NB: the '?' alternative in DEL_RE catches 'Δ' that was mangled to '?' in some titles (e.g. 'lsm1? Saccharomyces').
# Lower-case standard gene names = mutant alleles in yeast nomenclature ('jhd2', 'set2D', 'tps1-').
GENE_TOKEN_RE = re.compile(_B + r"([a-z]{3}\d{1,3})(?=[dDΔ∆△]?(?![a-z0-9]))")


def load_config(path: Path = CONFIG) -> dict:
    return yaml.safe_load(path.read_text())


def coxpresdb_samples(platform: str) -> pd.DataFrame:
    """COXPRESdb sample list -> sample_id, study, study_title, run_description, platform (file order kept)."""
    tag = PLATFORMS[platform]
    df = pd.read_csv(RAW / f"{tag}.samples.tsv", sep="\t", header=None, dtype=str, quoting=3,
                     keep_default_na=False)
    out = pd.DataFrame({"sample_id": df[1], "study": df[0].str.replace(".raw.1", "", regex=False),
                        "study_title": df[2]})
    out["run_description"] = df[3] if df.shape[1] > 3 else ""
    out["platform"] = platform
    return out


def _ena_text(study: str) -> pd.Series:
    p = RAW / "sample_metadata" / "ena" / f"{study}.tsv"
    if not p.exists():
        return pd.Series(dtype=str)
    t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    exp = t["experiment_title"].str.replace(r"^.*?sequencing(;|:)\s*(GSM\d+:\s*)?", "", regex=True)
    txt = t[[c for c in ENA_TEXT_COLS if c != "experiment_title"]].agg(" | ".join, axis=1) + " | " + exp
    return pd.Series(txt.values, index=t["run_accession"].values)


def sample_key(x: str) -> str:
    """Normalise array file names: COXPRESdb ids are R make.names()-style ('+'/'-' -> '.', leading 'X' before a
    digit, '.1' suffix on duplicates), SDRF 'Array Data File' values are raw file names."""
    x = re.sub(r"(\.CEL)\.\d+$", r"\1", str(x), flags=re.I)
    x = re.sub(r"^X(?=\d)", "", x)
    return re.sub(r"[^a-z0-9]", "", x.lower())


def _sdrf_text(study: str) -> tuple[pd.Series, pd.Series]:
    """Text per SDRF row keyed by Array Data File, and by GSM id (fallback)."""
    p = RAW / "sample_metadata" / "arrayexpress" / f"{study}.tsv"
    if not p.exists():
        return pd.Series(dtype=str), pd.Series(dtype=str)
    t = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, header=None)
    header = t.iloc[0].tolist()
    body = t.iloc[1:]
    keep = [i for i, h in enumerate(header) if not SDRF_SKIP.match(h.strip()) and not ORGANISM.search(h)]
    txt = body.iloc[:, keep].agg(" | ".join, axis=1)
    files = [i for i, h in enumerate(header) if h.strip() == "Array Data File"]
    by_file = pd.Series(txt.values, index=body.iloc[:, files[0]].values) if files else pd.Series(dtype=str)
    gsm = txt.str.extract(r"(GSM\d+)")[0]
    if files:
        by_file.index = [sample_key(x) for x in by_file.index]
    by_gsm = pd.Series(txt.values, index=gsm.values)
    return by_file, by_gsm[by_gsm.index.notna()]


def sample_metadata(platform: str) -> pd.DataFrame:
    """COXPRESdb samples joined with their per-sample text (ENA run report / ArrayExpress SDRF)."""
    s = coxpresdb_samples(platform)
    texts = []
    for study, grp in s.groupby("study", sort=False):
        if platform == "rnaseq":
            m = _ena_text(study)
            t = grp["sample_id"].map(m)
        else:
            by_file, by_gsm = _sdrf_text(study)
            by_file = by_file[~by_file.index.duplicated()]
            by_gsm = by_gsm[~by_gsm.index.duplicated()]
            t = grp["sample_id"].map(sample_key).map(by_file)
            t.index = grp.index
            gsm = grp["sample_id"].str.extract(r"(GSM\d+)")[0]
            t = t.fillna(gsm.map(by_gsm))
        texts.append(t)
    s["sample_text"] = pd.concat(texts).reindex(s.index)
    s["has_sample_metadata"] = s["sample_text"].notna()
    s["sample_text"] = s["sample_text"].fillna("")
    return s


def _strip_markers(text: str, markers: list[str]) -> str:
    pat = r"(?<![A-Za-z0-9])(Δ|∆|△|delta[-_ ]?)?(" + "|".join(markers) + r")(?![0-9])(Δ|∆|△|-?delta|[-_:]?\S*)?"
    return re.sub(pat, " ", text, flags=re.I)


def screen_text(text: str, markers: list[str], gene_names: set[str]) -> list[str]:
    """Return the list of perturbation evidence strings found in `text` (empty = no evidence)."""
    t = _strip_markers(text, markers)
    hits = []
    for rx, tag in ((DEL_RE, "del"), (MUT_RE, "mut"), (POINT_RE, "point")):
        for m in rx.finditer(t):
            hits.append(f"{tag}:{m.group(0).strip()}")
    for m in GENE_TOKEN_RE.finditer(t):
        tok = m.group(1)
        if tok in gene_names:
            hits.append(f"gene:{tok}")
    return sorted(set(hits))


def screen(platform: str, gene_names: set[str], config: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample and per-study perturbation screen. Returns (samples, studies)."""
    cfg = config or load_config()
    markers = cfg["background_markers"]
    s = sample_metadata(platform)
    s["sample_hits"] = [";".join(screen_text(x, markers, gene_names)) for x in s["sample_text"]]
    title_hits = {st: ";".join(screen_text(t, markers, gene_names))
                  for st, t in s.drop_duplicates("study")[["study", "study_title"]].itertuples(index=False)}
    s["title_hits"] = s["study"].map(title_hits)
    st = (s.groupby("study", sort=False)
          .agg(study_title=("study_title", "first"), n_samples=("sample_id", "size"),
               n_samples_flagged=("sample_hits", lambda x: int((x != "").sum())),
               n_with_metadata=("has_sample_metadata", "sum"), title_hits=("title_hits", "first"),
               sample_hits_union=("sample_hits", lambda x: ";".join(sorted({h for v in x for h in v.split(";")
                                                                                if h})))))
    st["auto_flag"] = (st["n_samples_flagged"] > 0) | (st["title_hits"] != "")
    overrides = cfg.get("study_overrides", {}) or {}
    st["decision"] = st["auto_flag"].map({True: "exclude", False: "keep"})
    st["reason"] = st["auto_flag"].map({True: "auto: perturbation evidence in title/sample metadata",
                                        False: "auto: no perturbation evidence"})
    for study, ov in overrides.items():
        if study in st.index:
            st.loc[study, "decision"] = ov["decision"]
            st.loc[study, "reason"] = "manual: " + ov["reason"]
    st = st.reset_index()
    s = s.merge(st[["study", "decision", "reason"]], on="study", how="left")
    s["excluded"] = s["decision"] == "exclude"
    s["sample_flagged"] = s["sample_hits"] != ""
    return s, st
