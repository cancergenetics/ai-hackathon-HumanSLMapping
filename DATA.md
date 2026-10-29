# Data register

**No data file is copied into or committed to this repository.** Every dataset below is read **in place, read-only**
from the upstream checkout `YEAST_SLs` (commit `9b5dce9`; location set in `config/paths.yaml`, guarded by
`src/hsl/guard.py`). `/data/` is gitignored; it may only hold files we generate.

Checksums: SHA-256. "Upstream record" = the value in `YEAST_SLs/DATA.md` / the per-source `MANIFEST.json` / the bundle
sidecar. **Verified** = we recomputed it from the file on disk on 2026-10-29 (`scripts/03_data_checksums.py` →
`results/data_checksums/checksums.json`, 461 files hashed; every raw download and every bundle/store file with an
upstream record matched it). Files with no upstream record (intermediate parquet tables) have their first recorded hash
in that JSON. Full provenance of each source (download, filtering, licence notes) is in `YEAST_SLs/DATA.md`; the
licence column is copied from there.

## 1. Genetic interactions (labels and GI inputs) — Costanzo et al. 2016

| item | value |
|---|---|
| Dataset | Costanzo et al. 2016 *Science* 353:aaf1420, **Data File S1** ("representative subset"), pair-wise format, **SGA_NxN.txt** (non-essential × non-essential) + `strain_ids_and_single_mutant_fitness.xlsx`. Upstream ID `costanzo2016-science-S1-subset` (upstream primary). |
| Source URL | https://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/Data%20File%20S1_Raw%20genetic%20interaction%20datasets:%20Pair-wise%20interaction%20format.zip |
| Version | Server Last-Modified 2016-11-28; 431,580,254 bytes; accessed by upstream 2026-10-26 |
| SHA-256 (verified) | zip `2134d990e66f754605e936d1bab70f874e264d1467c9a08534b054570fbe3775`; SGA_NxN.txt `ff6e7e754a05258e0b9e2df965dec7d703408a14b68d4acac1e5f81c67ac3681`; strain xlsx `0806bfb416e2bf343282e9b98e8d6553cd58e229613ce9ee1242beb80bc888b2` |
| Derived files we read | `data/processed/costanzo2016_bundle_s1.npz` `8e89ba3fa3bbea2b2346005456309ec3830176e16ee8ca62df429652e7588a71` (verified against its sidecar on every load by `ysl.data.dense.load_bundle`); sidecar `.json` `5112566e…c3f3`; processed parquet `costanzo2016_nxn_science_subset.parquet` `4b97918dd407da821b9ee512b31e24fcc64d4fd9a421451e140269cf34d9aab8` |
| Licence | Science supplementary material; article © AAAS; no explicit data licence. Academic use with citation; **do not redistribute**. |
| Use | ε, P, labels (P < 0.05 & ε < −0.12), array SMF; split stratification; evaluation (held-out rows only at evaluation time). |

Not used: TheCellMap complete release (`costanzo2016-cellmap-*`, upstream secondary; also on disk, verified), and the
ExN/NxE/ExE/DAmP tables (inside the zips, never extracted upstream).

Upstream S1 lockbox list `YEAST_SLs/config/lockbox_s1.txt` SHA-256 `3493f0cfa6b048c69cd59b7de9152bd5e15a7ea6b95323c18871af4e111b4c73`
(verified against `config/frozen.yaml` by `ysl.eval.splits.read_lockbox`): used only to prefer lockbox genes as test genes.

## 2. External (non-GI) feature stores read by the upstream methods

Which methods read which stores: `stack_lr` none; `stack_lr_wr2` the **pair store** (query-side gene×gene kernels);
`m1_hgb` the pair store (query and target×array blocks) **and the gene-level tables** (propensity sub-model).
`*_t2` methods additionally read the Hoepfner co-fitness store (tier 2; not used by us so far).

| store (under `YEAST_SLs/data/processed/features/`) | key files, SHA-256 (verified unless noted) | built from (upstream DATA.md IDs) |
|---|---|---|
| Pair store `pairs/` (4,426 genes = current names of S1 query ∪ array ORFs) | `genes.txt` `8df297e9…0e694`; `masks.npz` `b92d488b…45a6`; 29 feature `.npy` files each matching `pairs/manifest.json` (manifest `c7f963b0…0889`, no upstream record) | coxpresdb-sce, ncbi-gene-info, michaelis2023-apms, ensembl-fungi-compara116, esm2 (sgd-sequences + esm2-650m-weights), sgd-features, uniprot-UP000002311, paxdb5-single |
| Gene tables | `sequence.parquet` `4e1f0833…b5a`*; `abundance.parquet` `2d5b7290…d031`; `popgen.parquet` `b5d14ec4…579f`*; `alphafold.parquet` `ed632607…5b34`*; `uniprot.parquet` `1d172b37…d899`; `ppi_nodes.parquet` `a310c239…84f`; `expression_{rnaseq,microarray}.npy` `2df23807…2dbf`* / `7b410894…ab13`*; `coexpression_combined_fisherz.npy` `42531b5a…66dd`; `compara_phylo_profile.parquet` `fb41154c…a0f2`; `compara_paralog_pairs.parquet` `03c4a1f1…1188`; `esm2_650m_mean.npy` `3059e074…e234`; `michaelis2023_edges.parquet` `966825d0…2218` | as above + y1002-genomes, alphafold-api-yeast |

\* no upstream record; first hash in `results/data_checksums/checksums.json`. Full 64-character hashes for all files
are in that JSON.

### Raw sources behind the stores (all verified against upstream records)

| upstream ID | source URL | version | SHA-256 (raw file) | licence |
|---|---|---|---|---|
| coxpresdb-sce | https://coxpresdb.jp/download/Sce-r.c3-0/exp/Sce-r.c3-0.expression.combat.txt.zip ; …/Sce-m.c4-0/… | Sce-r.c3-0 (2021-11), Sce-m.c4-0 (2021-06) | r `8310fde5f0d2159a4f83e4ac403edf49d23566c54ae798f2171336e69eb6bd8b`; m `6ff938e906a0dbce1910c399bb0dc7a3f03769bc0fffca243ad211832531cb04` | CC BY |
| coxpresdb-sample-metadata | ENA / ArrayExpress per-study reports | accessed 2026-10-26 | 307 files, manifest digest `d4ff9ed7877598f6bbf6d32a0db1edaa6140077111f9cc1383462818f680a210` | EMBL-EBI terms |
| ncbi-gene-info | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Fungi/Saccharomyces_cerevisiae.gene_info.gz | 2026-10-26 | `92637d185bf5718ad43b192b321a57205eeb49f2772385fb0d4348eb0f72c900` | public domain |
| michaelis2023-apms | https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06739-5/MediaObjects/41586_2023_6739_MOESM4_ESM.rar | Supp. Data 2 (2023-12-06) | `49b0c66fb09bf041615c9d3de5b8645474fb59dc99781f0f4f804e987b369345` | CC BY 4.0 |
| paxdb5-single | https://pax-db.org/downloads/5.0/datasets/4932/ (Mueller 2020, Kulak 2014) | PaxDb 5.0 | Mueller `67378a93c42154d54d1ccdbdb6911bab1bd1a7de630c0336219f61a60f35c6e8`; Kulak `5cf4795921bbeb0a49574fae497fa1a621eef2c9217024dcf04f68da14d97a07` | CC BY 4.0 |
| ensembl-fungi-compara116 | https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-63/tsv/ensembl-compara/homologies/saccharomyces_cerevisiae/Compara.116.protein_default.homologies.tsv.gz | Fungi r63 / Compara 116 | `cb81826edbd3df44e62c505edbf5c9aafbdb0e337553d82b5158f24872314963` | no restrictions |
| uniprot-UP000002311 | https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000002311&… (no GO fields) | UniProt 2026_03 | `2ebe18165191728eea6ae6d7191c30d020dc13d6d0a7e447bbb759ead714d3f4` | CC BY 4.0 |
| y1002-genomes | http://1002genomes.u-strasbg.fr/files/ | 2025-03-31 | dNdS `ed930aeedf57644ffa227fe23b94f9c11de9a708033d44dec94db7a2fcecaa6a`; PresenceAbsence `63b62d6761c378f71360a259add9ad3317befb27cfa45944d0b11433766bf7fc`; CopyNumber `d22fe0a2024921163739553afb612e97e229c99b51a4301021c2433e11b2ee4d`; LoF `067327e1970c48f278dd08abb18f140622458829a6abc320207245c8c4b61c34` | no licence posted; cite; do not redistribute |
| alphafold-api-yeast | https://alphafold.ebi.ac.uk/api/prediction/<accession> | AFDB v6 | `c701712483b32a7184158efa98d95a5bd16cb92dd40a34082095775f3ae5c56d` | CC BY 4.0 |
| sgd-features | https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab | 2026-10-23 | `d7e8d2cc6a392d4e6ba1e5ee170c2879e8b41620d458796abd11bea4d51ebc38` | CC BY 4.0 |
| sgd-sequences | https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz ; …/orf_dna/orf_coding_all.fasta.gz | 2024-06-13 | protein `17e8b47e1ae23178c6000fbc4ab548f102d1b250ef9dff5d811feb3f03dd2c5b`; CDS `a0d329ba5110dcb3661f2fe3ae199f364d222dc5026927fd4451d89df1e7c782` | CC BY 4.0 |
| esm2-650m-weights | https://huggingface.co/facebook/esm2_t33_650M_UR50D | rev `08e4846e537177426273712802403f7ba8261b6c` | model.safetensors `a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0` | MIT |
| hoepfner2014-hiphop (tier 2, **not used so far**) | https://doi.org/10.5061/dryad.v5m8v | Dryad v1/v3 | HOP `99b386a84384eae847657ed41bf222c9550a87ef961f0ab191833c918771ffd7`; HIP `dbc5041defea9c046da0890d5e569f97d5f7afbf50ea0885f539ea8e5980cd24` | CC0 1.0 |

Circularity notes inherited from upstream (to be re-checked in our DESIGN review): coexpression studies with genetic
perturbations were excluded upstream; AP-MS, coexpression, ESM-2, paralogy, domains, phylogenetic profiles, abundance,
population genetics and AlphaFold have human analogues; no GO, no text mining, no STRING, no GI data from other
sources are in these stores.
