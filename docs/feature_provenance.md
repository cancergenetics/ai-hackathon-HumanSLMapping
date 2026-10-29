# Feature provenance (DESIGN v3 §1, §6.1(e))

Machine-checked source: `config/feature_provenance.yaml` (tests/test_provenance.py fails if a selector or predictor
opens a YEAST_SLs file that no store entry for that consumer covers, or if a store is missing from DATA.md). File
lists were recorded with `hsl.guard.record_reads` on 2026-10-29 (stack_lr_wr2 and m1_hgb with `drop=["genome"]` at
n = 30 on dev split 1; SelectorView features and clusters). Full source/version/checksum/licence: DATA.md and
`YEAST_SLs/DATA.md`.

Frozen P* `stack_lr_wr2_lb` (P2, 2026-10-29) opens the stores stack_lr_wr2 opens plus the gene-scalar stores already
cleared for m1_hgb (abundance, alphafold, apms, coexpression, compara, esm2_embeddings, popgen, sequence,
sgd_feature_names, uniprot_domains), recorded with `hsl.guard.record_reads` at n = 30 on dev split 1; `used_by` in
config/feature_provenance.yaml lists it. No new store.

| store | label-free? | GI / text-mining / GO-free? | human analogue | used by predictor | used by selector |
|---|---|---|---|---|---|
| Costanzo 2016 S1 bundle (ε, P, labels, array SMF) | no (it is the label source) | no (GI data) | HAP1 qGI screens; HAP1 single-KO fitness ↔ array SMF | all (screened rows only, via ScreenView) | screened rows (ε, P, labels; T/D_val columns masked) + array SMF (masked for excluded genes); **no DMF** |
| SGD_features.tab (name mapping only) | yes | yes (annotation columns unused) | HGNC/Ensembl IDs | m1_hgb (via loaders) | yes (loaders) |
| COXPRESdb expression + coexpression (perturbation studies removed) | yes | yes | COXPRESdb human / ARCHS4 | stack_lr_wr2 (pair kernel `coex_z`), m1_hgb | yes (expression level/variability, coexpression hub-ness, `coexp` kernel) |
| ESM-2 650M embeddings | yes | yes | same model, human proteome | stack_lr_wr2 (`esm_cos`), m1_hgb | yes (16 PCs, `esm` kernel) |
| PaxDb single-experiment abundance | yes | yes (integrated STRING-weighted set excluded) | PaxDb human | m1_hgb | yes |
| SGD sequences (length, GC, ENC, CAI) | yes | yes | Ensembl human | m1_hgb | yes |
| 1002 Genomes population genetics | yes | yes | gnomAD constraint | m1_hgb | yes |
| AlphaFold (pLDDT, disorder) | yes | yes | AFDB human | m1_hgb | yes |
| UniProt Pfam/InterPro (no GO fields) | yes | yes | UniProt human | stack_lr_wr2 (`dom_ipr_jaccard`), m1_hgb | yes (domain counts) |
| Michaelis 2023 AP-MS (literature columns dropped) | yes | yes | BioPlex / OpenCell | stack_lr_wr2 (`apms_*`), m1_hgb | yes (AP-MS node stats, `apms` kernel); not in the HAP1-minimal scenario |
| Ensembl Compara 116 (phylogenetic profile, paralogs) | yes | yes | Ensembl Compara | stack_lr_wr2 (`para_identity`), m1_hgb | yes (conservation, paralog counts/identity) |
| Pair store, non-genome kernels | yes | yes | as sources | stack_lr_wr2, m1_hgb | no (selectors use hsl.features kernels) |
| Pair store, genome (`gen_*`) | yes | yes | **none** (SGA linkage artefact) | no (m1_hgb `drop=["genome"]`; upstream recipe = sensitivity) | no |
| Hoepfner 2014 co-fitness (tier 2) | yes | yes (phenotype-entangled) | partial | sensitivity only (`*_t2`) | no |
| `stack_lr` | — | — | — | reads no external store (GI + array SMF only) | — |

Selector whitelist (`hsl.features.gene_features`): the upstream gene table (`ysl.features.genes.gene_matrix`) minus
every `has_*` availability indicator; genome coordinates are not in it. Never exposed to selectors: query SMF, DMF/DMF_sd,
n_arrays/n_queries, strain metadata, missingness/off-array indicators, any feature of an excluded (test / D_val) gene.
Inherent caveat: an off-array candidate has NaN array SMF and `candidate_col = -1`; selectors may route such genes
to an external-only submodel (DESIGN §6.2 H) but must not use the indicator as a feature.

Clusters (`hsl.features.clusters`): k-means (k = 100, `cluster_seed`) on a 32-dimensional kernel-PCA embedding of the
mean of three label-free kernels (coexpression z clipped to [−1, 1], cosine of standardised ESM-2 embeddings, AP-MS
FDR < 5 % edge) over all non-test S1 genes; computed once, cached in `data/cache/` (gitignored).
