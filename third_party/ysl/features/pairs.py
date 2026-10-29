"""External (non-GI) gene x gene PAIR features for the Costanzo 2016 S1 NxN gene universe (DESIGN.md §7).

Gene universe: sorted union of CURRENT SGD names of S1 query ORFs and array ORFs (`genes.txt`; built by
scripts/42_build_pair_features.py). Every feature is a dense symmetric [G x G] matrix over that universe, cached in
data/processed/features/pairs/ with a SHA-256 manifest (`manifest.json`) and packed availability masks (`masks.npz`).
NaN (float) / -1 (int8 flags) = unavailable; a feature is NaN exactly where its availability mask is False.
No GI data, GO, curated complexes or literature-derived resources enter any feature.

Families (name: definition; availability):
  coexpression  coex_z, coex_z_rnaseq, coex_z_microarray: mean within-study Fisher-z Pearson (processed COXPRESdb,
                perturbation studies removed; ysl.data.features.load_coexpression); NaN if no study measured both.
                coex_mr: mutual rank sqrt(rank_i(j) * rank_j(i)) of coex_z among all coexpression genes (lower =
                more coexpressed; Obayashi & Kinoshita 2009).
  AP-MS         Michaelis 2023 (literature columns dropped). apms_edge: best score_fdr_plus_cor over both directions,
                0 = no edge; available if either protein was a bait. apms_shared / apms_jaccard / apms_cosine:
                shared interaction partners and Jaccard / cosine of the undirected neighbourhoods; available if both
                proteins were observed (bait or in network). apms_same_mcl: same Markov cluster of the AP-MS network
                (data-derived, not curated); available if both have a cluster.
  paralogy      Ensembl Compara 116. para_any / para_within / para_other flags (within_species_paralog vs older
                other_paralog/gene_split), para_identity = max of the two % identities (0 if not paralogs);
                available if both genes are in Compara. esm_cos: cosine of ESM-2 650M mean embeddings after
                centring and per-dimension standardisation over all embedded genes (raw cosine is anisotropic).
  genome        SGD coordinates. gen_same_chr, gen_dist_kb (|midpoint difference| in kb; NaN if different
                chromosomes), gen_same_arm (centromere split; mitochondrion = one arm), gen_adj3k / gen_adj10k
                (same chromosome, edge-to-edge gap < 3 kb / < 10 kb; 3 kb = the evaluator's ADJ_BP).
  domains       UniProt Pfam / InterPro ids: dom_pfam_shared, dom_pfam_jaccard, dom_ipr_shared, dom_ipr_jaccard
                (Jaccard 0 when both have no domains); available if both genes have a UniProt entry.
  phylogeny     Compara 334-fungal-species profile. phylo_corr (Pearson of presence, 0 if a profile is constant),
                phylo_jaccard (presence), phylo_ident_corr (Pearson of max % identity); available if both in Compara.
  abundance     PaxDb Mueller 2020 log10 ppm: abund_absdiff |difference|, abund_sum (log of the product).

Loading (no ysl.data import; safe for ysl.models): `genes()`, `index_of(names)`, `load(name)` (memory-mapped stored
array), `block(name, rows, cols)` (float32 with NaN; rows/cols = gene indices, -1 -> NaN), `FEATURES`.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
PAIR_DIR = ROOT / "data" / "processed" / "features" / "pairs"
ADJ_BP = (3000, 10000)


@dataclass(frozen=True)
class PairFeature:
    name: str
    family: str
    dtype: str        # 'f16' | 'f32' | 'i8' (i8: 0/1 flag, -1 = unavailable)
    mask: str         # availability mask name in masks.npz
    direction: int    # +1: larger = expected more interaction-like; -1: smaller (univariate diagnostics only)
    description: str


FEATURES: dict[str, PairFeature] = {f.name: f for f in [
    PairFeature("coex_z", "coexpression", "f16", "coex", 1, "mean within-study Fisher-z, all kept studies"),
    PairFeature("coex_z_rnaseq", "coexpression", "f16", "coex_rnaseq", 1, "mean within-study Fisher-z, RNA-seq"),
    PairFeature("coex_z_microarray", "coexpression", "f16", "coex_microarray", 1,
                "mean within-study Fisher-z, microarray"),
    PairFeature("coex_mr", "coexpression", "f16", "coex", -1, "mutual rank of coex_z (1 = best)"),
    PairFeature("apms_edge", "apms", "f16", "apms_bait", 1, "best AP-MS edge score (0 = none)"),
    PairFeature("apms_shared", "apms", "f16", "apms_obs", 1, "shared AP-MS partners"),
    PairFeature("apms_jaccard", "apms", "f16", "apms_obs", 1, "Jaccard of AP-MS neighbourhoods"),
    PairFeature("apms_cosine", "apms", "f16", "apms_obs", 1, "cosine of AP-MS neighbourhoods"),
    PairFeature("apms_same_mcl", "apms", "i8", "apms_mcl", 1, "same AP-MS Markov cluster"),
    PairFeature("para_any", "paralogy", "i8", "compara", 1, "Compara paralog pair (any type)"),
    PairFeature("para_within", "paralogy", "i8", "compara", 1, "within_species_paralog (recent, WGD-like)"),
    PairFeature("para_other", "paralogy", "i8", "compara", 1, "other_paralog / gene_split (older)"),
    PairFeature("para_identity", "paralogy", "f16", "compara", 1, "max % identity of the paralog pair (0 = none)"),
    PairFeature("esm_cos", "paralogy", "f16", "esm", 1, "cosine of centred, standardised ESM-2 embeddings"),
    PairFeature("gen_same_chr", "genome", "i8", "coords", 1, "same chromosome"),
    PairFeature("gen_dist_kb", "genome", "f16", "coords_same_chr", -1, "|midpoint distance| kb (same chr only)"),
    PairFeature("gen_same_arm", "genome", "i8", "coords", 1, "same chromosome arm"),
    PairFeature("gen_adj3k", "genome", "i8", "coords", 1, "same chr, edge-to-edge gap < 3 kb"),
    PairFeature("gen_adj10k", "genome", "i8", "coords", 1, "same chr, edge-to-edge gap < 10 kb"),
    PairFeature("dom_pfam_shared", "domains", "f16", "uniprot", 1, "shared Pfam ids"),
    PairFeature("dom_pfam_jaccard", "domains", "f16", "uniprot", 1, "Jaccard of Pfam ids"),
    PairFeature("dom_ipr_shared", "domains", "f16", "uniprot", 1, "shared InterPro ids"),
    PairFeature("dom_ipr_jaccard", "domains", "f16", "uniprot", 1, "Jaccard of InterPro ids"),
    PairFeature("phylo_corr", "phylogeny", "f16", "compara_profile", 1, "Pearson of presence profiles"),
    PairFeature("phylo_jaccard", "phylogeny", "f16", "compara_profile", 1, "Jaccard of presence profiles"),
    PairFeature("phylo_ident_corr", "phylogeny", "f16", "compara_profile", 1, "Pearson of % identity profiles"),
    PairFeature("abund_absdiff", "abundance", "f16", "abundance", -1, "|log10 ppm difference| (Mueller 2020)"),
    PairFeature("abund_sum", "abundance", "f16", "abundance", 1, "log10 ppm sum = log of abundance product"),
]}
_NP = {"f16": np.float16, "f32": np.float32, "i8": np.int8}


# ===================================================================================================== helpers
def reindex_square(M: np.ndarray, src: pd.Index, genes: pd.Index, fill=np.nan, dtype=np.float32) -> np.ndarray:
    """[len(src)]^2 matrix over `src` gene order -> [G x G] over `genes` (missing genes -> fill)."""
    pos = src.get_indexer(genes)
    ok = np.flatnonzero(pos >= 0)
    out = np.full((len(genes), len(genes)), fill, dtype=dtype)
    out[np.ix_(ok, ok)] = M[np.ix_(pos[ok], pos[ok])]
    return out


def gene_mask(genes: pd.Index, present) -> np.ndarray:
    """Pair availability: both genes in `present`."""
    v = genes.isin(pd.Index(present))
    return v[:, None] & v[None, :]


def dedup_current(index, to_current) -> tuple[pd.Index, np.ndarray]:
    """Map a source index to current names; keep the first occurrence of each current name.
    Returns (current-name index, positions in the source kept)."""
    cur = to_current(list(index)).to_numpy()
    _, keep = np.unique(cur, return_index=True)
    keep = np.sort(keep)
    return pd.Index(cur[keep]), keep


def mutual_rank(z: np.ndarray) -> np.ndarray:
    """Mutual rank sqrt(r_ij r_ji); r_ij = 1-based rank of j among i's partners by decreasing z (self and NaN last).
    NaN where z is NaN. Ties broken by column order (stable)."""
    n = z.shape[0]
    zz = np.where(np.isnan(z), -np.inf, z).astype(np.float32)
    np.fill_diagonal(zz, -np.inf)
    order = np.argsort(-zz, axis=1, kind="stable")
    r = np.empty((n, n), dtype=np.float32)
    rows = np.arange(n)[:, None]
    r[rows, order] = np.arange(1, n + 1, dtype=np.float32)[None, :]
    mr = np.sqrt(r * r.T)
    mr[np.isnan(z)] = np.nan
    return mr


def row_corr(X: np.ndarray) -> np.ndarray:
    """Pearson correlation between rows of X; 0 where a row is constant."""
    X = X.astype(np.float64)
    X = X - X.mean(1, keepdims=True)
    nrm = np.linalg.norm(X, axis=1)
    ok = nrm > 1e-12
    X[ok] /= nrm[ok, None]
    X[~ok] = 0.0
    return (X @ X.T).astype(np.float32)


def set_overlap(B) -> tuple[np.ndarray, np.ndarray]:
    """Binary gene x item matrix (dense or scipy sparse) -> (shared counts, Jaccard; 0 when both sets empty)."""
    import scipy.sparse as sp
    B = sp.csr_matrix(B, dtype=np.float32)
    S = (B @ B.T).toarray().astype(np.float32)
    n = np.asarray(B.sum(1)).ravel().astype(np.float32)
    union = n[:, None] + n[None, :] - S
    with np.errstate(invalid="ignore", divide="ignore"):
        J = np.where(union > 0, S / union, 0.0).astype(np.float32)
    return S, J


def standardised_cosine(X: np.ndarray) -> np.ndarray:
    """Cosine similarity after centring and per-dimension standardisation of the embedding matrix X [n x d]."""
    X = X.astype(np.float64)
    X = (X - X.mean(0)) / np.where(X.std(0) > 1e-12, X.std(0), 1.0)
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    return (X @ X.T).astype(np.float32)


def genome_pairs(chrom: np.ndarray, lo: np.ndarray, hi: np.ndarray, arm: np.ndarray) -> dict:
    """Genome pair features from per-gene chromosome / lo / hi / arm-id (NaN chrom = unknown)."""
    known = ~np.isnan(chrom)
    av = known[:, None] & known[None, :]
    same = (chrom[:, None] == chrom[None, :]) & av
    mid = (lo + hi) / 2.0
    gap = np.maximum(0.0, np.maximum(lo[None, :] - hi[:, None], lo[:, None] - hi[None, :]))
    out = {"gen_same_chr": same, "gen_same_arm": same & (arm[:, None] == arm[None, :]), "_avail": av,
           "gen_dist_kb": np.where(same, np.abs(mid[:, None] - mid[None, :]) / 1000.0, np.nan).astype(np.float32)}
    for bp in ADJ_BP:
        out[f"gen_adj{bp // 1000}k"] = same & (gap < bp)
    return out


# ===================================================================================================== builders
def build_coexpression(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_coexpression, to_current
    feats, masks = {}, {}
    for kind, name, mname in (("combined", "coex_z", "coex"), ("rnaseq", "coex_z_rnaseq", "coex_rnaseq"),
                              ("microarray", "coex_z_microarray", "coex_microarray")):
        z, src, _ = load_coexpression(kind)
        cur, keep = dedup_current(src, to_current)
        z = z[np.ix_(keep, keep)]
        np.fill_diagonal(z, np.nan)
        feats[name] = reindex_square(z, cur, genes)
        masks[mname] = ~np.isnan(feats[name])
        if kind == "combined":
            feats["coex_mr"] = reindex_square(mutual_rank(z), cur, genes)
        del z
    return feats, masks


def build_apms(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_ppi_edges, load_table, to_current
    e = load_ppi_edges()
    nodes = load_table("ppi_nodes")
    node_cur = to_current(list(nodes.index))
    nodes = nodes.assign(cur=node_cur.values).drop_duplicates("cur").set_index("cur")
    G = len(genes)
    bait = genes.isin(nodes.index[nodes["ppi_bait_observed"].astype(bool)])
    obs = genes.isin(nodes.index[nodes["ppi_bait_observed"].astype(bool) | nodes["ppi_in_network"].astype(bool)])
    b = genes.get_indexer(to_current(list(e["bait_orf"])).values)
    p = genes.get_indexer(to_current(list(e["prey_orf"])).values)
    ok = (b >= 0) & (p >= 0) & (b != p)
    score = np.zeros((G, G), dtype=np.float32)
    sc = e["score_fdr_plus_cor"].to_numpy(dtype=np.float32)[ok]
    np.maximum.at(score, (b[ok], p[ok]), sc)
    score = np.maximum(score, score.T)
    # neighbourhoods use ALL network genes (not only the universe) so degrees are not truncated
    allg = pd.Index(sorted(set(to_current(list(e["bait_orf"])).values) | set(to_current(list(e["prey_orf"])).values)))
    bi, pi = allg.get_indexer(to_current(list(e["bait_orf"])).values), allg.get_indexer(
        to_current(list(e["prey_orf"])).values)
    import scipy.sparse as sp
    keep = bi != pi
    A = sp.coo_matrix((np.ones(keep.sum(), np.float32), (bi[keep], pi[keep])), shape=(len(allg), len(allg))).tocsr()
    A = ((A + A.T) > 0).astype(np.float32).tocsr()
    pos = allg.get_indexer(genes)
    inu = np.flatnonzero(pos >= 0)
    sel = sp.csr_matrix((np.ones(len(inu), np.float32), (inu, pos[inu])), shape=(G, len(allg)))
    rows = sel @ A                                                 # [G x allg] neighbourhood indicator
    shared, jac = set_overlap(rows)
    deg = np.asarray(rows.sum(1)).ravel().astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos = np.where((deg[:, None] * deg[None, :]) > 0, shared / np.sqrt(deg[:, None] * deg[None, :]), 0.0)
    m_bait = bait[:, None] | bait[None, :]
    m_obs = obs[:, None] & obs[None, :]
    mcl = nodes["markov_cluster"].reindex(genes).to_numpy(dtype=np.float64)
    has = ~np.isnan(mcl)
    m_mcl = has[:, None] & has[None, :]
    feats = {"apms_edge": np.where(m_bait, score, np.nan),
             "apms_shared": np.where(m_obs, shared, np.nan),
             "apms_jaccard": np.where(m_obs, jac, np.nan),
             "apms_cosine": np.where(m_obs, cos, np.nan).astype(np.float32),
             "apms_same_mcl": np.where(m_mcl, (mcl[:, None] == mcl[None, :]).astype(np.float32), np.nan)}
    return feats, {"apms_bait": m_bait, "apms_obs": m_obs, "apms_mcl": m_mcl}


def build_paralogy(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_esm, load_paralogs, load_phylo_profile, to_current
    G = len(genes)
    para = load_paralogs()
    compara_genes = set(to_current(list(load_phylo_profile().index)).values)
    m = gene_mask(genes, compara_genes)
    a = genes.get_indexer(to_current(list(para["orf_a"])).values)
    b = genes.get_indexer(to_current(list(para["orf_b"])).values)
    ok = (a >= 0) & (b >= 0) & (a != b)
    a, b = a[ok], b[ok]
    t = para["homology_type"].to_numpy()[ok]
    ident = np.fmax(para["identity_a"].to_numpy(dtype=np.float64), para["identity_b"].to_numpy(dtype=np.float64))[ok]
    flags = {k: np.zeros((G, G), dtype=np.float32) for k in ("para_any", "para_within", "para_other")}
    idm = np.zeros((G, G), dtype=np.float32)
    for k, sel in (("para_any", np.ones(len(a), bool)), ("para_within", t == "within_species_paralog"),
                   ("para_other", t != "within_species_paralog")):
        flags[k][a[sel], b[sel]] = 1
        flags[k][b[sel], a[sel]] = 1
    np.maximum.at(idm, (a, b), np.nan_to_num(ident).astype(np.float32))
    idm = np.maximum(idm, idm.T)
    feats = {k: np.where(m, v, np.nan) for k, v in flags.items()}
    feats["para_identity"] = np.where(m, idm, np.nan)
    esm = load_esm()
    cur, keep = dedup_current(esm.index, to_current)
    C = standardised_cosine(esm.to_numpy()[keep])
    np.fill_diagonal(C, np.nan)
    feats["esm_cos"] = reindex_square(C, cur, genes)
    em = gene_mask(genes, cur)
    np.fill_diagonal(em, False)
    return feats, {"compara": m, "esm": em}


def centromeres() -> pd.Series:
    """Chromosome (int) -> centromere midpoint (bp) from SGD_features (coordinates only)."""
    from ysl.data.sgd import COLS, SGD_FEATURES
    f = pd.read_csv(SGD_FEATURES, sep="\t", header=None, names=COLS, dtype=str, keep_default_na=False)
    f = f[f["feature_type"] == "centromere"]
    mid = (pd.to_numeric(f["start"]) + pd.to_numeric(f["stop"])) / 2
    return pd.Series(mid.values, index=pd.to_numeric(f["chromosome"]).astype(int).values).sort_index()


def build_genome(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_table
    co = load_table("coordinates").reindex(genes)
    chrom = co["chromosome"].to_numpy(dtype=np.float64)
    lo, hi = co["lo"].to_numpy(dtype=np.float64), co["hi"].to_numpy(dtype=np.float64)
    cen = centromeres()
    cmid = pd.Series(chrom).map(cen).to_numpy(dtype=np.float64)        # NaN for mito (17) / unknown
    right = np.where(np.isnan(cmid), 0.0, ((lo + hi) / 2 > cmid).astype(float))
    arm = np.where(np.isnan(chrom), np.nan, chrom * 2 + right)
    g = genome_pairs(chrom, lo, hi, arm)
    av = g.pop("_avail")
    feats = {k: (np.where(av, v.astype(np.float32), np.nan) if v.dtype == bool else v) for k, v in g.items()}
    return feats, {"coords": av, "coords_same_chr": g["gen_same_chr"].copy()}


def build_domains(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_table, to_current
    import scipy.sparse as sp
    u = load_table("uniprot")
    u = u.assign(cur=to_current(list(u.index)).values).drop_duplicates("cur").set_index("cur").reindex(genes)
    m = gene_mask(genes, u.index[u["uniprot"].notna()])
    feats = {}
    for col, tag in (("pfam", "pfam"), ("interpro", "ipr")):
        sets = [s.split(";") if isinstance(s, str) and s else [] for s in u[col]]
        vocab = {d: i for i, d in enumerate(sorted({d for s in sets for d in s}))}
        r = [i for i, s in enumerate(sets) for _ in s]
        c = [vocab[d] for s in sets for d in s]
        B = sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(genes), len(vocab)))
        S, J = set_overlap(B)
        feats[f"dom_{tag}_shared"] = np.where(m, S, np.nan)
        feats[f"dom_{tag}_jaccard"] = np.where(m, J, np.nan)
    return feats, {"uniprot": m}


def build_phylogeny(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_phylo_profile, to_current
    prof = load_phylo_profile()
    cur, keep = dedup_current(prof.index, to_current)
    X = prof.to_numpy(dtype=np.float32)[keep]
    P = (X > 0).astype(np.float32)
    _, J = set_overlap(P)
    feats = {"phylo_corr": reindex_square(row_corr(P), cur, genes),
             "phylo_jaccard": reindex_square(J, cur, genes),
             "phylo_ident_corr": reindex_square(row_corr(X), cur, genes)}
    return feats, {"compara_profile": gene_mask(genes, cur)}


def build_abundance(genes: pd.Index) -> tuple[dict, dict]:
    from ysl.data.features import load_table, to_current
    ab = load_table("abundance")
    ab = ab.assign(cur=to_current(list(ab.index)).values).drop_duplicates("cur").set_index("cur")
    x = ab["log10_ppm_mueller2020"].reindex(genes).to_numpy(dtype=np.float32)
    feats = {"abund_absdiff": np.abs(x[:, None] - x[None, :]), "abund_sum": x[:, None] + x[None, :]}
    has = ~np.isnan(x)
    return feats, {"abundance": has[:, None] & has[None, :]}


BUILDERS = (build_coexpression, build_apms, build_paralogy, build_genome, build_domains, build_phylogeny,
            build_abundance)


# ===================================================================================================== storage
def sha256_file(path: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def encode(name: str, M: np.ndarray) -> np.ndarray:
    spec = FEATURES[name]
    if spec.dtype == "i8":
        out = np.where(np.isnan(M), -1, M).astype(np.int8)
        if not np.isin(out, (-1, 0, 1)).all():
            raise ValueError(f"{name}: flag values outside {{0, 1}}")
        return out
    out = M.astype(_NP[spec.dtype])
    if spec.dtype == "f16" and np.isinf(out[np.isfinite(M)]).any():
        raise ValueError(f"{name}: float16 overflow")
    return out


def decode(name: str, A: np.ndarray) -> np.ndarray:
    if FEATURES[name].dtype == "i8":
        out = A.astype(np.float32)
        out[A < 0] = np.nan
        return out
    return A.astype(np.float32)


def save_features(genes: pd.Index, feats: dict, masks: dict, out_dir: Path = PAIR_DIR, extra: dict | None = None):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "genes.txt").write_text("\n".join(genes) + "\n")
    man = {"n_genes": len(genes), "genes_sha256": sha256_file(out_dir / "genes.txt"), "features": {}, "masks": {}}
    for name, M in feats.items():
        spec = FEATURES[name]
        mk = masks[spec.mask]
        nan = np.isnan(M)
        off = ~np.eye(len(genes), dtype=bool)
        if not np.array_equal(nan[off], ~mk[off]):
            raise ValueError(f"{name}: NaN pattern != availability mask {spec.mask}")
        if not np.array_equal(M, M.T, equal_nan=True):
            raise ValueError(f"{name}: not symmetric")
        p = out_dir / f"{name}.npy"
        np.save(p, encode(name, M))
        man["features"][name] = {"file": p.name, "dtype": spec.dtype, "mask": spec.mask, "family": spec.family,
                                 "direction": spec.direction, "description": spec.description,
                                 "sha256": sha256_file(p)}
    np.savez_compressed(out_dir / "masks.npz", **{k: np.packbits(v.ravel()) for k, v in masks.items()})
    man["masks"] = {"file": "masks.npz", "sha256": sha256_file(out_dir / "masks.npz"), "names": sorted(masks),
                    "packed": "np.packbits over the row-major [G x G] boolean mask"}
    man.update(extra or {})
    (out_dir / "manifest.json").write_text(json.dumps(man, indent=1))
    return man


# ===================================================================================================== loading
@lru_cache(maxsize=4)
def _manifest(out_dir: str = str(PAIR_DIR)) -> dict:
    p = Path(out_dir) / "manifest.json"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing; build with scripts/42_build_pair_features.py")
    return json.loads(p.read_text())


@lru_cache(maxsize=4)
def genes(out_dir: str = str(PAIR_DIR)) -> pd.Index:
    p = Path(out_dir) / "genes.txt"
    if sha256_file(p) != _manifest(out_dir)["genes_sha256"]:
        raise ValueError("pair-feature gene list checksum mismatch")
    return pd.Index([g for g in p.read_text().split("\n") if g], name="orf")


def index_of(names, out_dir: str = str(PAIR_DIR)) -> np.ndarray:
    """Gene-universe index of each (current) ORF name; -1 if absent."""
    return genes(out_dir).get_indexer(pd.Index(np.asarray(names, dtype=str)))


_VERIFIED: set = set()


def load(name: str, out_dir: str = str(PAIR_DIR), verify: bool = True) -> np.ndarray:
    """Stored (encoded) [G x G] array, memory-mapped; SHA-256 verified once per process."""
    meta = _manifest(out_dir)["features"][name]
    p = Path(out_dir) / meta["file"]
    if verify and (out_dir, name) not in _VERIFIED:
        if sha256_file(p) != meta["sha256"]:
            raise ValueError(f"pair feature {name}: checksum mismatch")
        _VERIFIED.add((out_dir, name))
    return np.load(p, mmap_mode="r")


def load_mask(name: str, out_dir: str = str(PAIR_DIR)) -> np.ndarray:
    G = len(genes(out_dir))
    with np.load(Path(out_dir) / "masks.npz") as z:
        return np.unpackbits(z[name], count=G * G).reshape(G, G).astype(bool)


def block(name: str, rows: np.ndarray, cols: np.ndarray, out_dir: str = str(PAIR_DIR)) -> np.ndarray:
    """Decoded float32 [len(rows) x len(cols)] sub-matrix (NaN = unavailable; index -1 -> NaN)."""
    A = load(name, out_dir)
    rows, cols = np.asarray(rows), np.asarray(cols)
    r_ok, c_ok = rows >= 0, cols >= 0
    out = np.full((len(rows), len(cols)), np.nan, dtype=np.float32)
    sub = np.asarray(A[np.maximum(rows, 0)][:, np.maximum(cols, 0)])
    dec = decode(name, sub)
    out[np.ix_(r_ok, c_ok)] = dec[np.ix_(r_ok, c_ok)]
    return out


def feature_names(out_dir: str = str(PAIR_DIR)) -> list[str]:
    return list(_manifest(out_dir)["features"])
