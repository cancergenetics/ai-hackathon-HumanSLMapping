"""External (non-GI) gene-level features for S. cerevisiae, keyed by systematic ORF name (e.g. 'YAL001C').

Tier-1 sources (DESIGN.md §7; download: scripts/10_download_features.py; build: scripts/11_process_features.py;
licences/checksums: DATA.md). Raw inputs live in data/raw/<source>/, processed outputs in data/processed/features/.

Builders (`build_*`) parse raw files into tidy tables; loaders (`load_*`) read the processed outputs. Every gene-level
table has an explicit index named 'orf'. Missing values are NaN (never zero-filled); `gene_feature_table()` adds
`has_<block>` availability indicators. NO GI data, NO GO, NO literature-curated columns enter any table:
  * SGD FASTA headers (curated descriptions) are discarded, only ORF ids + sequences are kept.
  * Michaelis 2023: PubMed / publication-count / APID literature / experimental-system columns, Complex Portal,
    UniProt subcellular location and SGD descriptions are dropped.
  * UniProt: only accession, reviewed flag, ordered-locus names, length, Pfam and InterPro ids (no GO fields).
  * PaxDb: single-experiment datasets only (the 'integrated' set weights datasets by STRING consistency).
  * COXPRESdb: studies with gene-specific genetic perturbations removed (src/ysl/data/expression_meta.py).
"""
from __future__ import annotations

import gzip
import io
import json
import re
import zipfile
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FEAT = PROC / "features"

ORF_RE = re.compile(r"^(Y[A-P][LR]\d{3}[WC](-[A-Z])?|Q\d{4}|R\d{4}[WC]?)$")


# ----------------------------------------------------------------------------------------------- gene universe
def nxn_genes() -> dict[str, pd.Index]:
    """Costanzo 2016 NxN (cellmap) query ORFs, array ORFs and their union (sorted Index objects)."""
    q = pd.read_parquet(PROC / "costanzo2016_nxn_cellmap_queries.parquet", columns=["query_orf"])["query_orf"]
    a = pd.read_parquet(PROC / "costanzo2016_nxn_cellmap_arrays.parquet", columns=["array_orf"])["array_orf"]
    qi, ai = pd.Index(sorted(set(q.astype(str)))), pd.Index(sorted(set(a.astype(str))))
    return {"query": qi, "array": ai, "union": qi.union(ai)}


ORF_LIKE_TYPES = {"ORF", "pseudogene", "blocked reading frame", "transposable element gene"}


@lru_cache(maxsize=1)
def orf_synonyms() -> pd.Series:
    """Retired systematic names -> current SGD systematic name (e.g. merged ORFs: 'YDL038C' -> 'YDL039C').
    Built from the SGD_features alias column (name lookup only). Costanzo 2016 uses some retired names."""
    from ysl.data.sgd import COLS, SGD_FEATURES
    f = pd.read_csv(SGD_FEATURES, sep="\t", header=None, names=COLS, dtype=str, keep_default_na=False)
    f = f[f["feature_type"].isin(ORF_LIKE_TYPES)]
    current = set(f["feature_name"])
    m = {}
    for name, alias in zip(f["feature_name"], f["alias"]):
        for a in alias.split("|"):
            if ORF_RE.match(a) and a not in current:
                m[a] = name
    return pd.Series(m, name="orf_current").sort_index()


def to_current(genes: pd.Index | list[str]) -> pd.Series:
    """Map ORF names to current SGD systematic names (identity when not retired). Index = input names."""
    syn = orf_synonyms()
    g = pd.Index(genes)
    return pd.Series([syn.get(x, x) for x in g], index=g, name="orf_current")


# ----------------------------------------------------------------------------------------------- sequences
def read_fasta_gz(path: Path) -> pd.Series:
    """FASTA.gz -> Series(seq) indexed by the first header token (ORF). Header descriptions are discarded."""
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    s = pd.Series(seqs, name="sequence")
    s.index.name = "orf"
    return s


CODON_TABLE = dict(zip(
    [a + b + c for a, b, c in product("TCAG", repeat=3)],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
SYN = {}
for _c, _aa in CODON_TABLE.items():
    SYN.setdefault(_aa, []).append(_c)


def codon_counts(cds: str) -> Counter:
    cds = cds.upper()
    return Counter(cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3) if set(cds[i:i + 3]) <= set("ACGT"))


def effective_number_of_codons(cnt: Counter) -> float:
    """Wright (1990) ENC from codon counts (stop, Met, Trp excluded). NaN if a degeneracy class is unobserved."""
    fam = {}
    for aa, codons in SYN.items():
        if aa in "*MW":
            continue
        n = sum(cnt[c] for c in codons)
        if n < 2:
            continue
        p = np.array([cnt[c] / n for c in codons])
        F = (n * (p ** 2).sum() - 1) / (n - 1)
        fam.setdefault(len(codons), []).append(F)
    need = {2: 9, 3: 1, 4: 5, 6: 3}
    enc = 2.0
    for k, nk in need.items():
        Fs = [f for f in fam.get(k, []) if f > 0]
        if not Fs:
            if k == 3:  # Ile only: standard fix
                Fs = [np.mean([np.mean(fam.get(2, [np.nan])), np.mean(fam.get(4, [np.nan]))])]
            else:
                return np.nan
        enc += nk / np.mean(Fs)
    return float(min(enc, 61.0))


def cai_weights(cds_by_gene: pd.Series, reference: list[str]) -> dict[str, float]:
    """Sharp & Li (1987) relative adaptiveness w(codon) from a reference gene set (pseudocount 0.5)."""
    tot = Counter()
    for g in reference:
        tot.update(codon_counts(cds_by_gene[g]))
    w = {}
    for aa, codons in SYN.items():
        if aa in "*MW":
            continue
        x = np.array([tot[c] + 0.5 for c in codons])
        for c, v in zip(codons, x / x.max()):
            w[c] = float(v)
    return w


def cai(cnt: Counter, w: dict[str, float]) -> float:
    num, n = 0.0, 0
    for c, k in cnt.items():
        if c in w:
            num += k * np.log(w[c])
            n += k
    return float(np.exp(num / n)) if n else np.nan


def build_sequence_features(cai_reference: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """Protein + CDS sequence scalars. CAI reference = `cai_reference` ORFs (chosen by the caller from an
    independent abundance measurement; see scripts/11_process_features.py). Returns (features, protein seqs)."""
    prot = read_fasta_gz(RAW / "sgd" / "orf_trans_all.fasta.gz")
    cds = read_fasta_gz(RAW / "sgd" / "orf_coding_all.fasta.gz")
    prot = prot.str.rstrip("*")
    w = cai_weights(cds, [g for g in cai_reference if g in cds.index])
    rows = {}
    for g, s in cds.items():
        s = s.upper()
        cnt = codon_counts(s)
        third = s[2::3]
        rows[g] = {
            "cds_length_nt": len(s),
            "gc_content": (s.count("G") + s.count("C")) / max(len(s), 1),
            "gc3_content": (third.count("G") + third.count("C")) / max(len(third), 1),
            "enc": effective_number_of_codons(cnt),
            "cai": cai(cnt, w),
        }
    df = pd.DataFrame.from_dict(rows, orient="index")
    df["protein_length_aa"] = prot.str.len().reindex(df.index)
    df.index.name = "orf"
    return df.sort_index(), prot.sort_index()


# ----------------------------------------------------------------------------------------------- coordinates
def build_coordinates() -> pd.DataFrame:
    from ysl.data.sgd import COLS, SGD_FEATURES
    df = pd.read_csv(SGD_FEATURES, sep="\t", header=None, names=COLS, dtype=str, keep_default_na=False)
    df = df[df["feature_type"].isin(ORF_LIKE_TYPES) & (df["start"] != "")]
    out = pd.DataFrame({
        "orf": df["feature_name"].values,
        "feature_type": df["feature_type"].values,
        "orf_qualifier": df["qualifier"].values,  # Verified / Uncharacterized / Dubious (SGD classification)
        "chromosome": pd.to_numeric(df["chromosome"]).astype(int).values,  # 17 = mitochondrial
        "start": pd.to_numeric(df["start"]).astype(int).values,
        "stop": pd.to_numeric(df["stop"]).astype(int).values,
        "strand": df["strand"].values,
    })
    out["lo"] = out[["start", "stop"]].min(axis=1)
    out["hi"] = out[["start", "stop"]].max(axis=1)
    out["mid"] = (out["lo"] + out["hi"]) / 2
    return out.drop_duplicates("orf").set_index("orf").sort_index()


# ----------------------------------------------------------------------------------------------- abundance
PAXDB_FILES = {
    "mueller2020": "4932-PXD145877_Mueller_Nature_2020_Saccharomyces_cerevisiae.txt",
    "kulak2014": "4932-PXD000270_Kulak_Nat_Methods_2014_S_cerevisiae.txt",
}


def build_abundance() -> pd.DataFrame:
    """log10(ppm) protein abundance from single PaxDb whole-organism MS datasets. NaN = not detected."""
    cols = {}
    for key, fn in PAXDB_FILES.items():
        t = pd.read_csv(RAW / "paxdb" / fn, sep="\t", comment="#", header=None, usecols=[0, 1],
                        names=["id", "ppm"])
        t["orf"] = t["id"].str.replace(r"^4932\.", "", regex=True)
        t = t[t["ppm"] > 0].drop_duplicates("orf").set_index("orf")
        cols[f"log10_ppm_{key}"] = np.log10(t["ppm"])
    df = pd.DataFrame(cols)
    df.index.name = "orf"
    return df.sort_index()


# ----------------------------------------------------------------------------------------------- 1002 genomes
def build_popgen() -> pd.DataFrame:
    """1011-isolate population statistics per reference ORF (Peter et al. 2018)."""
    d = RAW / "y1002"
    dnds = pd.read_csv(d / "gene_dNdS.tab.gz", sep="\t").set_index("Gene")
    dnds = dnds[["YN00_mean", "YN00_median", "NG86_mean", "LWL85_mean", "LPB93_mean"]]
    dnds.columns = ["dnds_yn00_mean", "dnds_yn00_median", "dnds_ng86_mean", "dnds_lwl85_mean", "dnds_lpb93_mean"]

    def ref_cols(m: pd.DataFrame) -> pd.DataFrame:
        orf = m.columns.str.replace(r"^X\d+\.", "", regex=True)
        ncopy = orf.str.extract(r"_NumOfGenes_(\d+)$")[0].astype(float).fillna(1).values
        orf = orf.str.replace(r"_NumOfGenes_\d+$", "", regex=True)
        keep = orf.str.match(ORF_RE.pattern)
        m = m.loc[:, keep]
        m.columns = orf[keep]
        return m, pd.Series(ncopy[keep], index=orf[keep])

    pa, nref = ref_cols(pd.read_csv(d / "genesMatrix_PresenceAbsence.tab.gz", sep="\t", index_col=0))
    cn, _ = ref_cols(pd.read_csv(d / "genesMatrix_CopyNumber.tab.gz", sep="\t", index_col=0))
    pop = pd.DataFrame({
        "pangenome_presence_freq": pa.mean(axis=0),
        "pangenome_n_isolates_called": pa.notna().sum(axis=0),
        "copy_number_mean": cn.mean(axis=0),
        "copy_number_sd": cn.std(axis=0),
        "copy_number_frac_multicopy": (cn > 1).sum(axis=0) / cn.notna().sum(axis=0),
        "pangenome_ref_num_genes": nref,  # '_NumOfGenes_k' suffix: k reference genes collapsed in one ORF cluster
    })
    pop = pop[~pop.index.duplicated()]
    txt = gzip.open(d / "1011LossOfFunction.xls.gz", "rb").read().decode("latin-1")
    lof = pd.read_csv(io.StringIO(txt.replace("\r\n", "\n").replace("\r", "\n")), sep="\t", dtype=str)
    lofm = lof.set_index("NAME").iloc[:, 6:].apply(pd.to_numeric, errors="coerce")
    lof_freq = pd.DataFrame({"lof_freq": lofm.mean(axis=1), "lof_n_isolates": lofm.notna().sum(axis=1)})
    lof_freq = lof_freq[~lof_freq.index.duplicated()]
    df = dnds.join(pop, how="outer").join(lof_freq, how="outer")
    df = df[df.index.str.match(ORF_RE.pattern)]
    df.index.name = "orf"
    return df.sort_index()


# ----------------------------------------------------------------------------------------------- UniProt
def build_uniprot() -> pd.DataFrame:
    """One row per ORF: UniProt accession(s), reviewed, length, Pfam / InterPro ids (';'-joined, sorted)."""
    t = pd.read_csv(RAW / "uniprot" / "UP000002311_domains.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    t.columns = ["accession", "reviewed", "oln", "gene_primary", "length", "pfam", "interpro"]
    rows = []
    for r in t.itertuples(index=False):
        olns = [x for x in re.split(r"[;/\s]+", r.oln) if ORF_RE.match(x)]
        for o in olns:
            rows.append((o, r.accession, r.reviewed == "reviewed", int(r.length), len(olns), r.pfam, r.interpro))
    u = pd.DataFrame(rows, columns=["orf", "uniprot", "reviewed", "uniprot_length", "n_orfs_for_entry", "pfam",
                                    "interpro"])
    u = u.sort_values(["orf", "reviewed", "n_orfs_for_entry"], ascending=[True, False, True])

    def ids(s: pd.Series) -> str:
        return ";".join(sorted({x for v in s for x in v.split(";") if x}))

    g = u.groupby("orf").agg(uniprot=("uniprot", "first"), uniprot_all=("uniprot", ";".join),
                             reviewed=("reviewed", "first"), uniprot_length=("uniprot_length", "first"),
                             pfam=("pfam", ids), interpro=("interpro", ids))
    g["n_pfam"] = g["pfam"].map(lambda s: len(s.split(";")) if s else 0)
    g["n_interpro"] = g["interpro"].map(lambda s: len(s.split(";")) if s else 0)
    return g


# ----------------------------------------------------------------------------------------------- AlphaFold
def build_alphafold(uniprot: pd.DataFrame) -> pd.DataFrame:
    """Canonical AFDB model (entry AF-<acc>-F1) per ORF: mean pLDDT and pLDDT band fractions.
    fraction_disordered = fraction of residues with pLDDT < 50 (AFDB 'very low')."""
    recs = {}
    with gzip.open(RAW / "alphafold" / "afdb_api_UP000002311.jsonl.gz", "rt") as f:
        for line in f:
            r = json.loads(line)
            want = f"AF-{r['accession']}-F1"
            for e in r["entries"]:
                if e.get("entryId") == want:
                    recs[r["accession"]] = {
                        "plddt_mean": e["globalMetricValue"],
                        "plddt_frac_very_low": e["fractionPlddtVeryLow"],
                        "plddt_frac_low": e["fractionPlddtLow"],
                        "plddt_frac_confident": e["fractionPlddtConfident"],
                        "plddt_frac_very_high": e["fractionPlddtVeryHigh"],
                        "afdb_model_length": (e["uniprotEnd"] or 0) - (e["uniprotStart"] or 1) + 1,
                    }
    af = pd.DataFrame.from_dict(recs, orient="index")
    out = uniprot[["uniprot"]].join(af, on="uniprot").drop(columns="uniprot")
    out["fraction_disordered"] = out["plddt_frac_very_low"]
    out = out.dropna(subset=["plddt_mean"])
    out.index.name = "orf"
    return out


# ----------------------------------------------------------------------------------------------- Compara
def build_compara() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ensembl Fungi Compara 116: (paralog pairs, phylogenetic profile).
    Paralogs: one row per unordered S. cerevisiae pair, types within_species_paralog / other_paralog / gene_split.
    Profile: ORF x species, max % identity (query side) of any ortholog in that species; 0 = no ortholog
    called (absence), rows only for ORFs present in the Compara file (others -> NaN at lookup)."""
    usecols = ["gene_stable_id", "identity", "homology_type", "homology_gene_stable_id", "homology_species",
               "homology_identity", "dn", "ds", "is_high_confidence"]
    h = pd.read_csv(RAW / "ensembl_compara" / "Compara.116.protein_default.homologies.tsv.gz", sep="\t",
                    usecols=usecols, na_values=["NULL"],
                    dtype={"gene_stable_id": str, "homology_gene_stable_id": str, "homology_type": str,
                           "homology_species": str})
    para = h[h["homology_type"].isin(["within_species_paralog", "other_paralog", "gene_split"])].copy()
    para = para[para["homology_species"] == "saccharomyces_cerevisiae"]
    a, b = para["gene_stable_id"], para["homology_gene_stable_id"]
    swap = a > b
    para["orf_a"] = np.where(swap, b, a)
    para["orf_b"] = np.where(swap, a, b)
    para["identity_a"] = np.where(swap, para["homology_identity"], para["identity"])
    para["identity_b"] = np.where(swap, para["identity"], para["homology_identity"])
    para = (para.sort_values(["orf_a", "orf_b", "homology_type"])
            .drop_duplicates(["orf_a", "orf_b"])
            [["orf_a", "orf_b", "homology_type", "identity_a", "identity_b", "dn", "ds", "is_high_confidence"]])
    para = para[para["orf_a"] != para["orf_b"]].reset_index(drop=True)
    orth = h[h["homology_type"].str.startswith("ortholog")]
    prof = orth.pivot_table(index="gene_stable_id", columns="homology_species", values="identity", aggfunc="max")
    genes = pd.Index(sorted(set(h["gene_stable_id"])))
    prof = prof.reindex(genes).fillna(0.0).astype("float32")
    prof.index.name = "orf"
    return para, prof


# ----------------------------------------------------------------------------------------------- Michaelis 2023
MICH = RAW / "michaelis2023" / "extracted" / "Interactome_Data"
EDGE_KEEP = ["score_FDR+cor", "score_FDR", "score_cor", "F_FDR_0.05", "F_FDR_0.01", "F_FDR_0.001", "F_FDR_1e-04",
             "R_FDR_0.05", "R_FDR_0.01", "R_FDR_0.001", "R_FDR_1e-04", "cor_val_>=4_[std_dev]",
             "cor_val_>=5_[std_dev]", "Inter-cluster edge"]
EDGE_DROPPED = ["Count of publications", "Pubmed Ids", "Binary (literature evidence APID)", "experimentral System",
                "Highlight Novel In Cytoscape", "Source/Target Gene names"]


def _cys_table(pattern: str) -> pd.DataFrame:
    """Read a typed table from the Cytoscape session (.cys = zip of CSV tables). Used because the decimal
    separators of 'cor_val' and 'TaBaQ' were stripped in the CSV export (e.g. 5.789950827 -> 5789950827)."""
    with zipfile.ZipFile(MICH / "The_Yeast_Interactome.cys") as z:
        name = [x for x in z.namelist() if pattern in x][0]
        raw = z.read(name).decode("utf-8")
    return pd.read_csv(io.StringIO(raw), skiprows=[0, 2, 3, 4], dtype=str)


def build_michaelis() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(edges, nodes). Edges: one row per (bait=source ORF, prey ORF) after expanding ';'-joined protein groups
    (ambiguous paralog peptides) into members (`prey_group_size` > 1). Literature columns dropped.
    Nodes: per ORF network membership, bait observed, GFP library, Markov cluster, TaBaQ abundance."""
    e = pd.read_csv(MICH / "The_Yeast_Interactome_Edges.csv", sep=";", dtype=str)
    cys = _cys_table("SHARED_ATTRS-org.cytoscape.model.CyEdge")[["shared name", "cor_val_[σ]"]]
    e = e.merge(cys, left_on="name", right_on="shared name", how="left", validate="1:1")
    if e["shared name"].isna().any():
        raise ValueError("Cytoscape edge table does not cover all edges")
    out = e[["source", "target"] + EDGE_KEEP].copy()
    out["cor_val_sd"] = pd.to_numeric(e["cor_val_[σ]"])
    for c in ["score_FDR+cor", "score_FDR", "score_cor"]:
        out[c] = out[c].astype(int)
    for c in EDGE_KEEP[3:]:
        out[c] = out[c].map({"TRUE": True, "FALSE": False})
    out.columns = [c.replace("+", "_plus_").replace(">=", "ge").replace("[std_dev]", "sd").replace("-", "_")
                   .replace(".", "_").replace(" ", "_").lower() for c in out.columns]
    out["prey_group_size"] = out["target"].str.count(";") + 1
    out["target"] = out["target"].str.split(";")
    out = out.explode("target").rename(columns={"source": "bait_orf", "target": "prey_orf"})
    out = (out.sort_values("score_fdr_plus_cor", ascending=False)
           .drop_duplicates(["bait_orf", "prey_orf"]).reset_index(drop=True))
    out = out[out["bait_orf"] != out["prey_orf"]].reset_index(drop=True)

    nodes = _cys_table("SHARED_ATTRS-org.cytoscape.model.CyNode")
    nodes = nodes[["shared name", "inGFPlibrary", "Bait detected+enriched (=self loop)", "Cluster ID (Markov)",
                   "TaBaQ (relative abundance)"]].copy()
    nodes.columns = ["name", "in_gfp_library", "bait_self_detected", "markov_cluster", "tabaq_abundance"]
    nodes["name"] = nodes["name"].str.split(";")
    nodes = nodes.explode("name")
    nodes = nodes[nodes["name"].str.match(ORF_RE.pattern)]
    for c in ["in_gfp_library", "bait_self_detected"]:
        nodes[c] = nodes[c].map({"true": True, "false": False})
    nodes["markov_cluster"] = pd.to_numeric(nodes["markov_cluster"])
    nodes["tabaq_abundance"] = pd.to_numeric(nodes["tabaq_abundance"])
    nodes = nodes.sort_values(["name", "bait_self_detected"], ascending=[True, False]).drop_duplicates("name")
    nodes = nodes.set_index("name")
    single = pd.read_csv(MICH / "Single_Nodes.csv", sep=";", dtype=str)
    baits = set(out["bait_orf"]) | set(single["source"])
    allg = pd.Index(sorted(set(nodes.index) | baits | set(out["prey_orf"])), name="orf")
    n = nodes.reindex(allg)
    n["ppi_in_network"] = allg.isin(set(out["bait_orf"]) | set(out["prey_orf"]))
    n["ppi_bait_observed"] = allg.isin(baits)
    n["ppi_degree_fdr05"] = (pd.concat([out.loc[out["f_fdr_0_05"] | out["r_fdr_0_05"], "bait_orf"],
                                        out.loc[out["f_fdr_0_05"] | out["r_fdr_0_05"], "prey_orf"]])
                             .value_counts().reindex(allg).fillna(0).astype(int))
    n["ppi_degree_all"] = pd.concat([out["bait_orf"], out["prey_orf"]]).value_counts().reindex(allg).fillna(0) \
        .astype(int)
    return out, n


# ----------------------------------------------------------------------------------------------- expression
def entrez_to_orf() -> pd.Series:
    gi = pd.read_csv(RAW / "ncbi" / "Saccharomyces_cerevisiae.gene_info.gz", sep="\t", dtype=str,
                     usecols=["GeneID", "LocusTag"])
    gi = gi[gi["LocusTag"].str.match(ORF_RE.pattern)]
    return gi.set_index("GeneID")["LocusTag"]


def read_coxpresdb_matrix(platform: str) -> pd.DataFrame:
    """COXPRESdb ComBat expression (log2 scale) -> genes (ORF) x samples. Unmapped Entrez ids dropped."""
    from ysl.data.expression_meta import PLATFORMS
    tag = PLATFORMS[platform]
    with zipfile.ZipFile(RAW / "coxpresdb" / f"{tag}.expression.combat.txt.zip") as z:
        with z.open(z.namelist()[0]) as f:
            df = pd.read_csv(f, sep="\t", index_col=0, dtype={0: str})
    df.index = df.index.astype(str).map(entrez_to_orf())
    df = df[df.index.notna()]
    df = df[~df.index.duplicated()]
    df.index.name = "orf"
    return df.astype("float32")


def study_fisherz_sums(expr: np.ndarray, studies: np.ndarray, min_samples: int = 6,
                       min_sd: float = 1e-3) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Sum over studies of Fisher-z(Pearson r) computed WITHIN each study (genes x samples input, gene order kept),
    and the number of studies contributing to each pair. A gene that is absent (NaN) or has (near-)zero variance
    in a study contributes nothing for that study. Studies with < min_samples samples are skipped.
    Returns (zsum float64, n_studies uint16, studies used)."""
    g = expr.shape[0]
    zsum = np.zeros((g, g), dtype=np.float64)
    cnt = np.zeros((g, g), dtype=np.uint16)
    used = []
    for st in pd.unique(studies):
        cols = np.where(studies == st)[0]
        if len(cols) < min_samples:
            continue
        x = expr[:, cols].astype(np.float64)
        finite = np.isfinite(x).all(axis=1)  # genes absent from this platform are all-NaN rows
        x[~finite] = 0.0
        x = x - x.mean(axis=1, keepdims=True)
        sd = np.sqrt((x ** 2).sum(axis=1))
        ok = finite & (sd > min_sd * np.sqrt(len(cols)))
        x[ok] /= sd[ok, None]
        x[~ok] = 0.0
        r = x @ x.T
        np.clip(r, -0.999999, 0.999999, out=r)
        z = np.arctanh(r)
        z[~ok, :] = 0.0
        z[:, ~ok] = 0.0
        zsum += z
        cnt += (ok[:, None] & ok[None, :]).astype(np.uint16)
        used.append(st)
    return zsum, cnt, used


def mean_fisherz(zsum: np.ndarray, cnt: np.ndarray) -> np.ndarray:
    """Mean per-study Fisher z (float32); NaN where no study measured both genes."""
    with np.errstate(invalid="ignore", divide="ignore"):
        mz = (zsum / cnt).astype(np.float32)
    mz[cnt == 0] = np.nan
    return mz


def per_study_coexpression(expr: np.ndarray, studies: np.ndarray, min_samples: int = 6,
                           min_sd: float = 1e-3) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Mean over studies of within-study Fisher-z(Pearson r). Returns (mean_z float32, n_studies, studies used)."""
    zsum, cnt, used = study_fisherz_sums(expr, studies, min_samples, min_sd)
    return mean_fisherz(zsum, cnt), cnt, used


# ----------------------------------------------------------------------------------------------- loaders
def _p(name: str) -> Path:
    return FEAT / name


def load_table(name: str) -> pd.DataFrame:
    """Gene-level table by name: sequence, coordinates, abundance, popgen, uniprot, alphafold, ppi_nodes."""
    df = pd.read_parquet(_p(f"{name}.parquet"))
    if df.index.name != "orf":
        raise ValueError(f"{name}: index must be 'orf'")
    return df


def load_protein_sequences() -> pd.Series:
    return pd.read_parquet(_p("protein_sequences.parquet"))["sequence"]


def load_paralogs() -> pd.DataFrame:
    return pd.read_parquet(_p("compara_paralog_pairs.parquet"))


def load_phylo_profile() -> pd.DataFrame:
    """ORF x fungal species max ortholog % identity (0 = no ortholog). Rows only for ORFs in Compara."""
    return pd.read_parquet(_p("compara_phylo_profile.parquet"))


def load_ppi_edges() -> pd.DataFrame:
    return pd.read_parquet(_p("michaelis2023_edges.parquet"))


def load_expression(platform: str = "rnaseq") -> tuple[np.ndarray, pd.Index, pd.DataFrame]:
    """(matrix genes x kept samples float32 log2, gene Index, sample metadata of the kept samples)."""
    m = np.load(_p(f"expression_{platform}.npy"))
    genes = pd.Index(pd.read_parquet(_p(f"expression_{platform}_genes.parquet"))["orf"], name="orf")
    meta = pd.read_parquet(_p("expression_sample_metadata.parquet"))
    meta = meta[(meta["platform"] == platform) & ~meta["excluded"]].reset_index(drop=True)
    if m.shape != (len(genes), len(meta)):
        raise ValueError("expression matrix / metadata shape mismatch")
    return m, genes, meta


def load_coexpression(kind: str = "combined") -> tuple[np.ndarray, pd.Index, np.ndarray]:
    """(mean per-study Fisher-z coexpression float32 [NaN = no study measured both genes], gene Index,
    number of studies per pair uint16). kind in {'combined', 'rnaseq', 'microarray'}."""
    z = np.load(_p(f"coexpression_{kind}_fisherz.npy"))
    n = np.load(_p(f"coexpression_{kind}_nstudies.npy"))
    genes = pd.Index(pd.read_parquet(_p("coexpression_genes.parquet"))["orf"], name="orf")
    if z.shape != (len(genes), len(genes)):
        raise ValueError("coexpression shape mismatch")
    return z, genes, n


SCALAR_BLOCKS = {
    "sequence": ["protein_length_aa", "cds_length_nt", "gc_content", "gc3_content", "enc", "cai"],
    "abundance": ["log10_ppm_mueller2020", "log10_ppm_kulak2014"],
    "popgen": ["dnds_yn00_mean", "pangenome_presence_freq", "copy_number_mean", "copy_number_frac_multicopy",
               "lof_freq"],
    "alphafold": ["plddt_mean", "fraction_disordered"],
    "uniprot": ["n_pfam", "n_interpro"],
    "ppi_nodes": ["ppi_bait_observed", "ppi_in_network", "ppi_degree_fdr05", "tabaq_abundance"],
    "coordinates": ["chromosome", "mid"],
}


def gene_feature_table(genes: pd.Index | None = None) -> pd.DataFrame:
    """Wide table of single-gene scalar features for `genes` (default: NxN union). Missing = NaN; each block gets
    a boolean `has_<block>` indicator (True if the gene has a row in that source)."""
    genes = nxn_genes()["union"] if genes is None else pd.Index(genes, name="orf")
    cur = to_current(genes)
    out = pd.DataFrame({"orf_current": cur.values}, index=genes)
    for block, cols in SCALAR_BLOCKS.items():
        t = load_table(block)
        out[f"has_{block}"] = cur.isin(t.index).values
        vals = t[cols].reindex(cur.values)
        vals.index = genes
        out = out.join(vals, how="left")
    out.index.name = "orf"
    return out


def load_esm(tag: str = "esm2_650m") -> pd.DataFrame:
    """Mean-pooled ESM-2 protein embeddings (scripts/12_esm_embeddings.py): float32 DataFrame, index 'orf'
    (Costanzo NxN ORF names; genes with no SGD protein sequence are absent), columns esm_0..esm_{d-1}."""
    X = np.load(_p(f"{tag}_mean.npy"))
    genes = [g for g in (_p(f"{tag}_genes.txt")).read_text().split("\n") if g]
    return pd.DataFrame(X, index=pd.Index(genes, name="orf"), columns=[f"esm_{i}" for i in range(X.shape[1])])
