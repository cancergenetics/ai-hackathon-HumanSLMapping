# LITERATURE — Using the yeast GI network to design human HAP1 SL-mapping screens

Compiled 2026-10-29 by the literature agent for `humanSLmapping`.

**How this was checked.** Primary sources were read as follows:
- PubMed abstracts via NCBI E-utilities.
- Full text of Billmann et al. 2026 (Cell; Europe PMC PMC13281610) and its companion qGI methods preprint (PMC12803112).
- The Billmann 2026 supplementary tables (Data S1, S5, S14, S20) from the Boone-lab supplement site, which I downloaded and parsed myself.
- The Costanzo 2010 Science PDF.
- bioRxiv and arXiv abstracts, plus the GeneDisco and IterPert full texts.

**Markers used below:**
- **[computed]** — a number I derived from a supplementary table. It is not quoted in the paper.
- **[upstream]** — a fact taken from `/Users/colmclaude/claude_experiments/YEAST_SLs/LITERATURE.md` (whose agent read the primary sources). I did not re-check it here.
- **[unverified]** — I could not confirm it from a primary source in this session.

References are numbered [n] and listed at the end.

---

## 0. The twelve most decision-relevant findings

1. **What Billmann et al. 2026 did [1, 2].**
   - **Design:** 222 unique HAP1 query knockout lines and 298 genome-wide screens (including biological replicates for 50 queries), all with the TKOv3 library (~18,000 genes, 4 gRNAs/gene, ~71k gRNAs).
   - **Scoring:** qGI score. The "standard" threshold is |qGI| > 0.3 and FDR < 0.1.
   - **Yield:** 3,934,506 unique gene pairs tested; 88,933 GIs (**47,052 negative**, 41,881 positive); 3,341 "extreme" SLs (qGI < −0.6, FDR < 0.01, plus fitness filters).
   - **Coverage:** the authors say they have assayed "<5% of all possible HAP1 gene pairs".
2. **Billmann's query lines were not chosen for network coverage.** The authors "focused on functionally diverse genes that showed a range of gene expression levels and single mutant fitness defects" [1].
   - Of the 222 non-redundant lines, **185 were bought off the shelf from Horizon Discovery and 37 were made in-house** [computed from Data S1].
   - Many queries come from hypothesis-driven projects (lipid metabolism, mTOR, glycosylation and so on).
   - There was no iterative or model-based selection. The vendor catalogue lists ~7,500 HAP1 KO models [43, vendor claim, unverified count of distinct genes].
3. **HAP1 and yeast labels have similar prevalence and noise. This is the key calibration for translating our simulations.**
   - Prevalence: the mean HAP1 negative-GI density at the standard threshold is **1.20%** of tested pairs [computed, Data S14]. Our yeast strong-SL prevalence is ≈1.1% (NxN S1).
   - HAP1 per-screen accuracy at the standard threshold is **precision 0.61 and recall 0.37** (strict threshold: 0.91 / 0.12), measured against 7 queries screened 4–5 times [1].
   - Yeast single-screen stringent negatives have precision ≈54% and recall ≈39% [3 via upstream].
   - So the yeast "strong SL" label corresponds best to the HAP1 **standard negative GI**. It does not correspond to the HAP1 "extreme SL" (0.14% density), which is closer to yeast ε < −0.35.
4. **Billmann's estimate of the size of the HAP1 network [1].**
   - The candidate universe is ~11,000 expressed genes (log2(TPM+1) > 1), giving ~60.5 M possible pairs.
   - That implies **~1.4 M GIs in total**, ~85,000 extreme SLs and ~45,000 extreme suppressions.
   - Extrapolation was by multiplying density (no saturation model).
   - Of these genes, 1,524 (~14–15%) are HAP1-essential. Their library-side interactions are partly masked.
5. **Yeast queries were also not random.** Costanzo 2010 screened 1,712 queries (334 TS/hypomorphic alleles of essential genes). These "were selected randomly with respect to function; however, preference was given to mutants exhibiting fitness defects" [4].
   - Screening ~30% of the genome as queries gave GI data for ~75% of genes, because nearly all nonessential genes are on the array [4].
   - The 2016 NxN screen covers nearly all nonessential genes (3,283 query strains × 3,827 arrays) [3 via upstream].
6. **Single-mutant fitness, expression and conservation predict GI degree in both species.**
   - Yeast: degree vs SMF r = 0.73 [4]. A feature model predicts degree at r = 0.80 and transfers to *S. pombe* [5].
   - HAP1: hubs are highly and stably expressed, old, constrained genes. "DepMap fitness metrics are predictive of genetic connectivity in the HAP1 network". The authors suggest these features could be used to predict hubs [1].
   - Degree heuristics therefore transfer. They are natural baselines for discovery-oriented query selection.
7. **Individual SL pairs transfer poorly from yeast to human. Network structure transfers better.**
   - Overall: ~4% of HAP1-tested pairs have a yeast ortholog pair, and only ~7% (289/4,189) of HAP1 GIs among conserved pairs are also yeast GIs [1].
   - Negative GIs: in my tabulation of Data S20 (98,275 ortholog pairs tested in both species), P(HAP1 negative | yeast negative) = 3.7% vs 2.0% otherwise, an odds ratio ≈ 1.95 [computed].
   - Earlier work agrees: worm vs yeast GIs are "not highly conserved" [6]; budding vs fission yeast share ~29% of SL interactions [7]; conservation is highest within complexes [8].
8. **No published study has done model-guided, iterative choice of *query genes* for a genome-scale GI map.** I found no SGA/E-MAP or HAP1 paper that did this (absence of evidence, extended searches 2026-10).
   - The closest precedents select *which genes to perturb* so as to predict the profiles of unperturbed genes:
     - IterPert (Perturb-seq) [9]: ~3–5× fewer perturbations than uniform sampling. TypiClust was the best generic baseline, and prior-knowledge kernels were essential in early rounds (cold start).
     - Active learning on target × compound matrices (Murphy lab) [10–13]: 48×48 matrix learned with 29% of experiments; ~60% of hits found after 3% of experiments; stopping rules.
     - BATCHIE (drug combinations) [14]: accurate predictions after exploring ~4% of 1.4 M experiments.
9. **Objective matters.** Strategies that maximise *model accuracy on held-out items* differ from those that maximise *hits discovered*.
   - GeneDisco: random and k-means diversity were relatively good for model improvement but poor for hit discovery; there is "a trade-off between model improvement and hit discovery" [15].
   - Sutherland et al. separate prediction-oriented from search-oriented criteria on low-rank matrices [16].
   - Hacohen et al.: at **low budgets, typical/representative selection beats uncertainty sampling**; the ordering reverses at high budgets [17].
10. **Selecting a query = observing a whole row (≈ column) of a symmetric matrix.** This is column-subset selection / adaptive Nyström sampling.
    - Adaptive column sampling recovers low-rank matrices with fewer observations than uniform sampling [18, 19].
    - Transductive experimental design picks the points that best reconstruct the rest in kernel space [20].
    - These map directly onto our kernel-ridge fold-in predictor.
11. **Cold-start evaluation is mandatory.** Pair-level splits inflate SL/GI prediction badly.
    - Human SL methods: all 12 benchmarked methods drop sharply from pair split (CV1) to one-gene-unseen (CV2) and both-unseen (CV3) splits [21 via upstream].
    - Yeast GI: AUROC ~0.9 → ~0.6 on independent data [22 via upstream].
    - Most published human SL predictors (SL2MF, GRSMF, KG4SL, PiLSL, KR4SL) use GO and/or knowledge graphs, so they are **not admissible** here (§4).
12. **HAP1 has the same "array-side partial profile" structure as yeast.** Every screened query measures nearly every library gene (~17,800 genes), so the profile of an unscreened gene over the screened queries is observed.
    - Exception: essential library genes (1,524), whose interactions are mostly masked. This mirrors yeast, where essential genes are absent from the NxN array.
    - HAP1-specific artefacts mirror yeast linkage and batch effects: chromosome-arm/ploidy biases, a 30-Mbp disomic chr15 region, and "functionally coherent covariation" [2].

---

## 1. Yeast GI mapping

### 1.1 Primary maps

**Tong et al. 2001, Science [23].**
- SGA method: a query mutation crossed into ~4,700 viable deletion mutants.
- 8 queries gave 291 interactions among 204 genes.
- Matters for us: establishes the query × array design, which is the template for HAP1 query × library.

**Tong et al. 2004, Science [24].**
- 132 queries × ~4,700 deletions gave ~1,000 genes and ~4,000 interactions.
- The network has "dense local neighborhoods; therefore, the position of a gene on a partially mapped network is predictive of other genetic interactions".
- This was the earliest explicit statement that a partial screen predicts unscreened interactions.

**Costanzo et al. 2010, Science [4]** (PDF text read).
- 1,712 queries (334 conditional/hypomorphic alleles of essential genes), ~5.4 M pairs, ~170,000 interactions.
- Query choice: "selected randomly with respect to function; however, preference was given to mutants exhibiting fitness defects".
- Coverage: "Despite screening only ~30% of the genome as query genes, we recovered genetic interactions for ~75% of the genome because partial genetic interaction profiles were generated for nearly all nonessential genes".
- Degree vs SMF r = 0.73. Most genes have few interactions and a few are hubs.
- About 2% of arrays have >8× more negative than positive interactions.
- AB/BA precision of stringent negatives is 0.89 (conservative) [upstream].

**Costanzo et al. 2016, Science [3].**
- Scale: >23 M double mutants; ~550,000 negative and ~350,000 positive GIs. Essential genes are densely connected hubs.
- Design [upstream, SOM]:
  - NxN: 3,283 nonessential query strains (natMX) × 3,827-gene deletion array (kanMX).
  - ExN/NxE/ExE: TS alleles (TSA: 786 alleles of 560 essential genes, confirmed in the methods of [1]).
- Screening: each screen run once with 4 replicate colonies. Scoring uses the multiplicative model f_ij = f_i·f_j + ε_ij plus a P-value.
- Thresholds: lenient P < 0.05; intermediate |ε| > 0.08; stringent ε < −0.12 (negative) / ε > 0.16 (positive).
  ε < −0.35 is treated as "most likely SL", giving ~10,000 SL pairs in NxN.
- Reproducibility of a single screen (26 queries × 5 replicates) for stringent negatives: precision ≈54%, recall ≈39%.
- AB/BA reciprocal agreement is "reasonable" but imperfect.
- 319 query strains were removed by a co-batch filter [upstream].
- For us: query rows and array columns of the same gene are different strains. Reciprocity is informative but not a replicate.

**Baryshnikova et al. 2010, Nat Methods [25]** [upstream].
- SGA scoring model and normalisation; batch effects are often as large as GIs.
- Ceiling: the replicate-based limit on correlation with ε is r = 0.67 (as cited in Yu 2016).

**Kuzmin et al. 2018, Science [26].**
- ~200,000 triple mutants. Double-mutant queries were "selected … across a broad spectrum of biological processes, spanning a range of quantitative features of the global digenic interaction network".
- Estimated trigenic network ≈100× larger than digenic.
- Relevance: the only Boone-lab precedent for query sets chosen deliberately to span feature space (stratified/diversity design). It was not iterative.

**Kuzmin et al. 2020, Science [27].**
- Digenic and trigenic interactions for whole-genome-duplication (WGD) paralog pairs.
- About 30% of duplicate pairs interact negatively vs ~3.6% of random pairs [upstream].

**Costanzo et al. 2021, Science [28].**
- 14 conditions × ~30,000 "functionally representative" pairs. Most GIs are unchanged across conditions.
- Relevance: the network is robust to environment, so a single-condition map is a reasonable target (also argued for HAP1 rich vs minimal medium [1]).

**E-MAP strategy: Schuldiner et al. 2005, Cell [29]; Roguev et al. 2008, Science [30]; Ryan et al. 2012, Mol Cell [8].**
- E-MAPs are high-density maps of "logically connected gene subsets" (e.g. 424 early-secretory genes).
- This is the alternative "focused" design: dense within modules, nothing between them.
- Relevance: the opposite extreme to diversity-based query selection.

**Pooled and barcode designs.**
- iSeq / barcode-fusion GI (BFG-GI; Díaz-Mejía et al. 2018, Mol Syst Biol) [31] map GIs in pooled format.
- Compressed Perturb-seq (Yao et al. 2024) decompresses random multi-perturbation pools using sparsity, at an order-of-magnitude lower cost [32].
- I found no compressed-sensing design for genome-scale query × library GI fitness screens [extended search].

### 1.2 Degree distribution, hubs and partial-screen recovery

- **Degree is heavy-tailed and predictable.**
  - Yeast: SMF r = 0.73 [4]. A feature model gives r = 0.80 for unseen genes, transferable to *S. pombe*; predictors are SMF defect, disorder, pleiotropy, PPI degree and low expression variation [5].
  - Essential genes have ~5× more interactions than nonessential genes [1, citing 3].
- **Partial-screen recovery.**
  - The Costanzo 2010 ~30% → ~75% statement [4] is about *genes with any GI data*, not recall of SL pairs.
  - With symmetric handling, screening a fraction f of queries observes 1 − (1 − f)² of unordered pairs among genes that are both query and array (f = 0.05 → 9.75%; f = 0.20 → 36%) [arithmetic].
  - For HAP1 (library ≈ all genes), n queries observe ≈ n × 11,000 of 60.5 M expressed pairs (n = 222 → ~4%; the paper says <5%) [1].
- **Ulitsky et al. 2009 [33]** [upstream]: prioritising pairs with an imputation model recovered 40% of negative GIs at 61% precision by testing 4.2% of pairs. This is a pair-level (not query-level) adaptive design.

---

## 2. Human GI mapping

### 2.1 Billmann, Costanzo et al. 2026 — the target design

Citation: *Cell* 189:3736–3759.e58 (June 2026), doi:10.1016/j.cell.2026.03.044 [1]. Preprint bioRxiv 10.1101/2025.06.30.662193. Data: https://boonelab.ccbr.utoronto.ca/supplement/billmanncostanzo2026/ and Mendeley. All items below were verified in the full text or supplementary tables.

| Item | Value |
|---|---|
| Cell line | HAP1 (Horizon clone C631; near-haploid; male, lost Y). Chosen because it "lacks aneuploidies and is amenable to LOF screens" |
| Query lines | **222 unique query genes** (non-redundant set). 298 screens in the complete dataset; 50 queries have 2–8 biological replicates. Data S1 lists 363 screens including WT controls |
| Query construction | Clonal LOF lines: a 3-gRNA pool (Synthego or TKOv3 guides), limiting dilution and Sanger confirmation, or Horizon catalogue lines. In Data S1, 185/222 are Horizon and 37 in-house [computed]. 209/222 show positive "self-interactions", consistent with complete LOF; 13 suggest partial LOF |
| Essential queries | A few in-house partial alleles of essential genes (e.g. PSMC5 point mutant, POLR2A, TUBB, CCT3, TCP1) [computed from Data S1]. Hypomorphs of essential genes are flagged as future work |
| How queries were chosen | "functionally diverse genes that showed a range of gene expression levels and single mutant fitness defects in HAP1 cells (Figure S4)". No iterative or model-based selection described |
| Library | TKOv3 [34]: ~18,000 protein-coding genes (17,804 in the matrix), 4 gRNAs/gene, ~70k gRNAs. Lentiviral, ≥200-fold coverage (15 M cells per technical replicate × 3), 10–12 population doublings. T0 sequenced at 400×, endpoint at 200× |
| Controls | 39 WT screens (18 rich, 21 minimal medium) for single-mutant fitness. 1,524 HAP1-essential genes (~15% of expressed genes), 81% of them core-essential in DepMap |
| Scoring | qGI = differential LFC (query vs WT) after loess fit, plus corrections for genomically linked/ploidy biases and functionally coherent covariation (SVD) [2]. A 2-gRNA/gene subset "has the potential to perform as well as the full 70,000-gRNA TKOv3 library" [2] — relevant to cost per screen |
| Thresholds | Standard \|qGI\| > 0.3, FDR < 0.1: per-screen **precision 0.61, recall 0.37**. Strict \|qGI\| > 0.6, FDR < 0.01: precision 0.91, recall 0.12. The MCMC (JEDER) gold standard comes from 7 queries × 4–5 replicates. The threshold optimised precision/recall of GO co-annotation |
| Interactions | 3,934,506 unique pairs; 88,933 GIs (47,052 negative, 41,881 positive) |
| Extreme SL | qGI < −0.6, FDR < 0.01, library-gene LFC > −0.5, double-mutant fitness < −1. **3,341 pairs** (0.14% of ~2.44 M expressed pairs) |
| Densities (Data S14) [computed] | Mean negative density per library gene 1.20% (all genes), 1.42% (essential), 2.31% (nonessential with fitness defect). Library negative degree over 223 queries: median 1, 25% zero, 95th percentile 10, 99th 20, max 46 |
| Hubs | Top 5% of genes have ~6× more interactions than average. Hubs are highly and stably expressed, old, mutation-intolerant. >30 gene features are associated with degree (Data S14: HAP1 SMF, DepMap CERES mean/SD, gene age, s_het, pLI, RVIS…). Mitochondrial genes are positive-GI hubs; the authors caution about mitochondrial screen artefacts [35] |
| Paralogs | Expressed paralogs with ≥20% identity are >7-fold enriched for negative GIs; at ≥90% identity, >60-fold |
| PPI | Pairs sharing a PPI are ~3-fold enriched for negative GIs |
| Reciprocity | "agreement among reciprocal gene pairs (query A–library B vs query B–library A)" (Fig. S5F; no number extracted) |
| Profile network | 3,784 genes at PCC > 0.41; 412 modules. Overlap with DepMap co-essentiality networks is low (Jaccard < 0.1 vs ~0.5 DepMap–DepMap). "The DepMap co-essentiality network appears to be approaching saturation", while the HAP1 network is not |
| Conservation | See §3 |
| Size estimate | ~11,000 expressed genes → 60.5 M pairs × 2.26% → ~1.4 M GIs. ~85,000 extreme SLs; ~45,000 suppressions |
| Cost and time per screen | **Not stated** [unverified]. The CHyMErA protocol quotes 15–20 weeks to implement a combinatorial screen [36]. Horizon says in-stock HAP1 KO lines ship in 1–2 weeks and made-to-order lines in ~12 weeks [43, vendor blog, unverified] |

**Companion methods paper: Billmann, Costanzo, …, Moffat, Myers. "Quantitative analysis of genetic interactions in human cells from genome-wide CRISPR-Cas9 screens", bioRxiv 10.1101/2025.06.30.662330 (2025) [2].**
- 363 screens, 39 WT controls, ~125,000 pairs measured in 4–5 replicates.
- Artefacts: whole-chromosome and arm biases; a 30-Mbp chr15 region (disomic in most clones, haploid in eHAP1-derived clones); correcting it in the WDR73 screen kept only 7 of 77 positive GIs and added 11 new negatives.
- "Highly variable" library genes (80, excluded) can flip GI sign between screens.
- Journal version: I could not find a 2023 *Mol Syst Biol* "Billmann et al." qGI paper. The brief's "Billmann 2023 MSB" citation is **[unverified]**. The closest MSB paper is Rahman et al. 2021 (FLEX benchmarking, mitochondrial bias) [35].

### 2.2 Earlier HAP1 and isogenic query-line work

- **Blomen et al. 2015, Science [37].**
  - Gene-trap mutagenesis in haploid HAP1/KBM7: ~2,000 fitness genes.
  - A secretory-pathway SL network "based exclusively on mutations" that "parallels observations made in yeast".
  - Query count and SL numbers [unverified].
- **Hart et al. 2015, Cell [38]:** TKO library screens in 5 lines including HAP1; 1,580 core fitness genes.
- **Hart et al. 2017, G3 [34]:** design of TKOv3; core-essential reference CEG2 (684 genes).
- **Aregger et al. 2020, Nat Metab [39]:**
  - FASN-KO HAP1 × TKOv3 introduced qGI-style scoring.
  - Further queries: LDLR, SREBF1, SREBF2, ACACA (all later included in [1]).
- **Chatr-aryamontri et al. 2025, bioRxiv 10.1101/2025.08.04.665369 [40]:**
  - NALM-6 "gene by genome-wide" design: a single query sgRNA added to cells carrying an inducible genome-wide pool, so no clonal isogenic lines are needed.
  - 26 screens. The authors say it "can be rapidly performed". Data at chemogenix.iric.ca.
  - A cheaper query-line route; query choice criteria not stated.
- **Desjardins et al. 2026, Cell Rep 45:117961 [41]:**
  - Genome-wide CRISPR screens in isogenic lines for 15 cancer drivers, plus DepMap analysis.
  - Validates PELO–SKI-complex dependence in 9p21.3-deleted/MSI-H contexts and FBXW7–PKMYT1, among others.
  - Another query × genome-wide design, chosen for cancer relevance.

### 2.3 Combinatorial (pair-library) CRISPR maps

These are pair-library designs: the experimenter chooses *pairs*, not query lines. They are relevant for validation and paralogs.

| Study | Design | Scale / result |
|---|---|---|
| Horlbeck et al. 2018 Cell [42] | CRISPRi dual-sgRNA, K562 and Jurkat | 222,784 gene pairs (~472 genes [unverified gene count]); GI maps cluster functionally related genes |
| Najm et al. 2018 Nat Biotech [44] | Orthogonal SaCas9/SpCas9 ("Big Papi") | MAPK, apoptosis and TP53 pairs across cell types |
| Shen et al. 2017 Nat Methods; Han et al. 2017 Nat Biotech [45] | Dual-gRNA SpCas9 | Tens of genes / ~21,000 drug-target pairs [details not re-checked] |
| Dede et al. 2020 Genome Biol [46] | enCas12a multiplex, ~400 paralog pairs, 3 lines | 24 SL paralog pairs; 79% seen in ≥2 of 3 lines |
| Thompson et al. 2021 Nat Commun [47] | 1,191 pairs (paralogs + predicted SL) | 105 pairs with fitness loss; 27 across multiple lines |
| Parrish et al. 2021 Cell Rep [48] | pgPEN, >2,000 paralogs, 2 lines | 12% of paralog pairs SL in ≥1 context |
| Ito et al. 2021 Nat Genet [49] | 5,065 paralog pairs, 3,284 genes | DUSP4/DUSP6 digenic dependence |
| Köferle et al. 2022 Cell Rep [50] | CRISPR + DepMap mining | >2,000 candidate paralog dependencies |
| Gonatopoulos-Pournatzis/Aregger CHyMErA [36] | Cas9–Cas12a hybrid guides | Protocol: 15–20 weeks |
| Esmaeili Anvar et al. 2024 Nat Commun [51] | in4mer Cas12a 4-guide arrays; Inzolia library ~49k clones | ~4,000 paralog pairs; 5× smaller library than other GI methods; most paralog SLs are platform-specific |
| Fong et al. 2025 Nat Genet (SCHEMATIC) [52] | Combinatorial KOs across lineages | 1,805 SLs (95% unreported); 34% penetrant across lineages |
| Fielden et al. 2025 Nature (SPIDR) [53] | Dual-guide CRISPRi, 548 DDR genes, RPE1 TP53−/− | ~150,000 gene-level pairs |
| Flister et al. 2025 Cell Rep [54] | Digenic screen of 36,648 paralog pairs; 462 pairs × 49 models | SLs infrequent and variable in penetrance; shared-PPI overlap and essentiality explain strength |

**Lesson for us:** paralog pairs are the best-validated human SL class, but a small fraction of all SLs. In HAP1, ≥20%-identity paralogs are only 7× enriched [1]. A query design should not be judged on paralog recovery alone.

### 2.4 DepMap-based SL inference (human-only context)

- **Lord et al. 2020, eLife [55]:** 220 robust driver-associated GIs; robust SLs are enriched among physically interacting pairs.
- **De Kegel et al. 2021, Cell Syst [56]:** random forest predictor of paralog SL; sequence identity and shared PPIs are the top features (AP ≈ 49% vs 11% chance [upstream]).
- **Other methods:** Benfatto et al. 2021, PARIS [57]; Lee et al. 2018, ISLE [58]; Jerby-Arnon et al. 2014, DAISY [59] [not re-checked].
- **Co-essentiality:** Wainberg et al. 2021 [60]; DepMap/CERES (Tsherniak 2017 [61]) and Chronos (Dempster 2021 [62]).
- **Billmann 2026 "expression dependency" (ED) [1]:** ~218,000 negative and ~247,000 positive ED scores among ~4 M pairs. They complement qGI.
- **Mitochondrial bias:** co-essentiality has a predominant mitochondrial bias reflecting screen dynamics [35].

### 2.5 Numbers for translating a yeast query fraction into HAP1 query lines

| Quantity | Yeast (Costanzo 2016 NxN) | HAP1 (Billmann 2026) |
|---|---|---|
| Genes targetable on array/library | 3,827 nonessential deletions (TSA adds 560 essential genes) | ~17,800 (TKOv3) |
| Expressed / relevant gene universe | ~6,000 genes; NxN query universe ≈ nonessential genes | ~11,000 expressed; ~9,500 expressed nonessential [arithmetic: 11,000 − 1,524] |
| Essential genes | ~1,000–1,100 [unverified exact] | 1,524 (1,231 core) |
| Queries screened | 3,283 strains (≈ all nonessential genes) | 222 unique (≈2.0% of expressed, ≈2.3% of expressed nonessential) |
| Label prevalence | Strong SL (ε < −0.12, P < 0.05) ≈1.1% | Standard negative 1.20%; extreme SL 0.14% |
| Single-screen noise | Precision ≈54%, recall ≈39% (stringent) | Precision 0.61, recall 0.37 (standard) |
| Mean negatives per query | ≈1.1% × 3,827 ≈ 40 strong SLs [arithmetic] | 47,052/222 ≈ 212 standard negatives [arithmetic] |

Translating yeast fractions:
- **By fraction** of a ~9,500–11,000 gene universe: 5% ≈ 475–550 lines, 10% ≈ 950–1,100, 20% ≈ 1,900–2,200.
- **By absolute count:** yeast 5% ≈ 164 queries and 20% ≈ 657 queries. Billmann's 222 lines sit between these.

Which translation is right depends on whether prediction quality scales with the *fraction* of queries or with the *number of queries per functional module*. No published study resolves this. It is a question our simulations can address (see §7).

---

## 3. Yeast → human conservation of GIs/SLs

| Work | Finding | Relevance |
|---|---|---|
| Billmann et al. 2026 [1] | ~4% of HAP1-tested pairs (146,664/3,934,728) have a yeast ortholog pair; ~7% (289/4,189) of HAP1 GIs among conserved pairs also seen in yeast. Negative GIs significantly conserved; positive not. Conserved negatives enriched in vesicle trafficking, mitosis, DNA replication/repair. Module-level negative density correlates (Fig. 7A). "Physiological and evolutionary properties of … hub genes are also conserved". [computed, Data S20] OR ≈ 1.95: P(HAP1 neg \| yeast neg) 3.7% vs 2.0% | **Pair-level transfer is weak, structure-level transfer is good.** Using yeast labels to choose HAP1 queries can only help a little directly. Our yeast simulations inform *strategy* (how many, which kind) rather than *which specific genes* |
| Tischler et al. 2008 Nat Genet [6] | Yeast GIs "are not highly conserved in animals" (*C. elegans*) | Caution against ortholog-based label transfer |
| Dixon et al. 2008 PNAS [7] | *S. pombe* vs *S. cerevisiae* SGA, ~220 genes: ~29% of GIs shared | Upper bound for a short evolutionary distance |
| Roguev et al. 2008 Science [30]; Ryan et al. 2012 Mol Cell [8]; Frost et al. 2012 Cell [63] | Conservation highest within complexes, lower within processes, lowest between processes; functional "repurposing" of some complexes | Within-complex structure (captured by AP-MS/PPI features) is the most transferable signal |
| Ryan et al. 2013 GBE [64] | Complexes are "all or nothing" for essentiality; essentiality flips at complex level between species | Essentiality and fitness priors are module-level properties |
| Koch et al. 2012 Genome Biol [5] | GI degree model (SMF, disorder, PPI degree, expression variation) transfers Sc → Sp | Degree heuristics transfer across species |
| Srivas et al. 2016 Mol Cell [65] | ~169,000 yeast interactions among orthologs of tumour suppressors and drug targets (several genotoxic conditions) used to prioritise thousands of TSG–drug tests in HeLa. Stability across environments and shared function raise the chance of a human hit | Precedent for "yeast-guided" human testing. Gains were modest and conditional |
| Deshpande et al. 2013 Cancer Res [66] | Comparative-genomics prioritisation from yeast GIs found SMARCB1–PSMA4 and ASPSCR1–PSMC2 SLs | Small-scale precedent |
| Jacunski et al. 2015 PLoS CB [67] | SINaTRA: network "connectivity homology" transfers SL models across species (AUROC 0.86 Sc → Sp [upstream]) | Uses curated SL labels; pair-level CV |
| De Kegel & Ryan 2019 [68]; De Kegel 2021 [56] | Paralog buffering explains variable essentiality. 25–35% of yeast paralog pairs are SL vs ~3.5% in human training data [upstream] | Paralog SL rates differ by species. Do not port yeast paralog effect sizes |
| Lehner et al. 2006 Nat Genet [69] | *C. elegans*: chromatin-regulator GI hubs conserved across animals | Hub identity partly conserved |

---

## 4. SL/GI prediction for unscreened queries

The setting here is a *cold-start query with a partial profile*:
- The query row is unseen.
- Its array/library column over screened queries is seen.
- External gene features are available.

This is Park & Marcotte class C2 / Pahikkala setting S2 [70, 71 via upstream].

| Method family | Examples | Inputs / admissible here? | Notes |
|---|---|---|---|
| kNN / collaborative filtering on profiles | Ryan et al. 2010 symmetric kNN for E-MAPs [72]; Ulitsky 2009 GSG features [33] | GI-only, transductive. Admissible (own data) | Strong baseline; uses partial profiles |
| Low-rank / network-guided matrix completion | Žitnik & Zupan 2015 (PPI and GO networks) [73]; Liany 2020 collective MF [74]; inductive MC (Natarajan & Dhillon 2014) [75] | Side information may include GO → use only allowed kernels | Inductive MC handles unseen rows |
| Kernel methods / pairwise KRR | Two-step / Kronecker KRR (Stock 2018) [76]; our upstream fold-in KRR + stacked combiner | Allowed kernels: coexpression, AP-MS, ESM, paralogy, SMF, abundance, chemogenomic co-fitness (tier 2) | Our predictor: AUPRC 0.061 / 0.083 / 0.102 at 10 / 20 / 30% screened, prevalence 0.011 (RESEARCH_TREE) |
| Logistic MF / graph-regularised MF | SL2MF [77] (**GO + PPI**), GRSMF [78] (**GO**) | **Not admissible as published** (GO) | Could be re-run with allowed similarity kernels |
| Knowledge-graph GNNs | KG4SL [79] (KG incl. GO biological processes, compounds, diseases), PiLSL [80] (KG enclosing graphs + omics), KR4SL [upstream] | **Not admissible** (KG / text-derived) | Do well in CV3 benchmarks partly through KG priors |
| Graph autoencoders / multi-view | SLMGAE (multi-view: PPI, GO, …) [21 via upstream] | Not admissible with GO views | Best overall in Feng 2024 under CV1 |
| Feature-based classifiers | ELISL (protein-sequence embeddings + omics) [81]; De Kegel paralog RF [56]; Benstead-Hume Slant (PPI topology) [upstream] | Mostly admissible (sequence, PPI, expression) | ELISL designed for robustness to selection bias |
| Perturbation-response models | GEARS [82] (GO-based gene graph); IterPert/GEARS [9] | GEARS's GO graph not admissible | Shows "predict unseen perturbations" is feasible with priors |

**Benchmarks and pitfalls.**
- **Feng et al. 2024, Nat Commun 15:9058 [21]** [upstream; PubMed shows a 2025 correction]:
  - 12 methods; CV1 (pairs) / CV2 (one gene unseen) / CV3 (both unseen); negatives 1:1 → 1:50.
  - All methods degrade from CV1 to CV3 and as negatives increase; NDCG@10 < 0.01 under CV3 for nearly all methods.
- **Seale et al. 2022 [83]:** SL labels show strong gene selection bias; SL-topology models are the most sensitive to it.
- **Park & Marcotte 2012 [70]:** pair-input CV classes.
- **Pandey et al. 2010 [22]:** yeast GI AUROC 0.84 (CV) → 0.62 (independent).
- **Implication for us:**
  - Evaluate by *query-level* holdout on all tested pairs (no balancing).
  - Report AUPRC lift over prevalence.
  - Include a degree-product baseline: degree drives pooled AUPRC.

---

## 5. Experimental design and active learning for screens

### 5.1 Objectives to keep separate

| Objective | What it rewards | Typical acquisition | Precedent |
|---|---|---|---|
| **A. Model accuracy on held-out queries** (our Q2) | Queries whose profiles reduce prediction error on *other* genes | Representativeness/diversity (k-center core-set, k-means/TypiClust, facility location, DPP/volume sampling, transductive experimental design), expected error/variance reduction, column-subset/leverage sampling | GeneDisco: random and k-means relatively good for model improvement [15]; IterPert: distance maximisation in a prior kernel, TypiClust next best [9]; Hacohen: typical points at low budget [17] |
| **B. SLs discovered by the screens themselves** | Queries with many SLs (hubs) | Exploitation of predicted degree, UCB/Thompson, active search | Kangas 2014: ~60% of hits after 3% of experiments [12]; Garnett 2012 active search [84]; DiscoBAX: diverse high-effect sets [85]; BATCHIE [14] |
| **C. Total SLs known** = screened + confidently predicted | Both | Hybrid / multi-objective | Sutherland 2013 defines criteria for each objective [16] |

These objectives conflict. Hub-seeking (B) over-samples high-SMF and highly expressed genes, as Costanzo 2010 and Billmann 2026 did implicitly. That leaves under-sampled modules poorly predicted.

### 5.2 Batch selection methods (pool-based; batch = 10 queries)

- **Uncertainty / information.**
  - BALD and BatchBALD (Kirsch, van Amersfoort & Gal, NeurIPS 2019) [86]: joint mutual information for batches.
  - Expected model change and expected error reduction (Settles 2009 survey) [87].
- **Diversity / representativeness.**
  - Core-set / k-center greedy (Sener & Savarese, ICLR 2018) [88].
  - BADGE (Ash et al., ICLR 2020) [89].
  - Submodular facility location (Wei, Iyer & Bilmes, ICML 2015) [90].
  - DPPs (Kulesza & Taskar 2012) [91]; batch DPP active learning (Bıyık et al. 2019) [92, not re-checked].
  - TypiClust (Hacohen, Dekel & Weinshall, ICML 2022) [17].
  - LCMD and a benchmark framework for batch active learning in regression (Holzmüller et al., JMLR 2023) [93, not re-checked].
- **Kernel/GP design:**
  - Transductive experimental design picks points that best linearly reconstruct the pool in feature space (Yu, Bi & Tresp, ICML 2006) [20].
  - Mutual-information sensor placement (Krause, Singh & Guestrin, JMLR 2008) [94].
  - V/Σ-optimality on graphs (Ma, Garnett & Schneider, NeurIPS 2013) [95, not re-checked].
  - These are natural for our KRR fold-in, since the KRR posterior variance is computable in closed form.
- **Matrix completion / column selection.**
  - Krishnamurthy & Singh 2013 (NeurIPS): adaptive column sampling. Pick columns whose observed entries are poorly explained by the span of already-selected columns. Noiseless recovery from O(r²n log r) observations, better than uniform sampling [18].
  - Bhargava, Ganti & Nowak 2017 (AISTATS): active PSD completion via adaptive Nyström sampling [19].
  - Nyström/leverage sampling theory: Gittens & Mahoney 2016 (JMLR) [96].
  - Active matrix completion: Chakraborty et al. 2013 (ICDM) [97, not re-checked]; Ruchansky, Crovella & Terzi 2015 (KDD) [98].
  - Active learning *and search* on low-rank matrices (Sutherland, Póczos & Schneider, KDD 2013), including drug–target discovery [16].
  - **Direct mapping to our problem.** Screening query X reveals X's full row. Before screening, X's *column* over the already-screened queries is observed (yeast arrays and HAP1 library both cover nearly all genes).
    The residual of X's partial column, projected onto the span of screened profiles, is therefore a cheap, non-circular score of "how novel X is". This is the adaptive-sampling criterion of [18].

### 5.3 Biological precedents for iterative, model-guided screening

| Work | Setting | Strategy | Result |
|---|---|---|---|
| Huang, Lopez, Hütter, …, Regev. IterPert, bioRxiv 10.1101/2023.12.12.571389 [9] | Perturb-seq K562: 2,058 essential-gene perturbations (205 held out); genome-scale 9,748. Predict responses of *unperturbed* genes with GEARS | Batches of 100 × 5 rounds (genome-scale 300 × 3). Greedy distance maximisation in a fused kernel of priors (other Perturb-seq, optical pooled screens, scRNA coexpression, ESM-2 embeddings, **STRING PPI, BioGPT literature**) + model kernel. Benchmarked against BALD, BatchBALD, Core-Set, BADGE, ACS-FW, LCMD, TypiClust, random | "comparable accuracy at only a third of the number of perturbations profiled as the next best method"; ~1 round ≈ 5 rounds of uniform sampling. Model kernel poorly aligned with truth in early rounds, so priors matter most early. All AL methods beat random here. Authors argue few rounds are practical because each takes ~1 month. **Closest analogue to our Q2.** Its STRING/literature priors are disallowed here |
| Mehrjou et al. 2022 ICLR, GeneDisco [15] | 4 genome-wide CRISPR screens (T-cell, NK) | 9 acquisition functions; batches 16–512; 5 seeds | Random and k-means competitive for test MSE (model improvement), weak for hit discovery; a trade-off between the two |
| Lyle et al. 2023 ICML, DiscoBAX [85] | Gene-perturbation discovery | Bayesian algorithm execution for a diverse high-effect set | Beats AL/BO baselines for diverse discovery |
| Naik et al. 2013 PLoS One; Naik et al. 2016 eLife [10, 11] | Proteins × compounds (imaging); 48 × 48 matrix | Active learning on a matrix model; robot-executed | Learned with 29% of experiments; beat human-like strategies |
| Kangas et al. 2014 BMC Bioinf [12] | 177 assays × 20,000 compounds | Joint target–compound modelling + AL | ~60% of hits after 3% of experiments |
| Temerinac-Ott et al. 2015 BMC Bioinf [13] | Drug–target matrices | Stopping rules from simulated AL traces | Up to 40% savings. Precedent for "when to stop" |
| Tosh, Tansey et al. 2025 Nat Commun, BATCHIE [14] | Drug combinations (206 drugs; 1.4 M experiments) | Information-theoretic Bayesian batch design | Accurate predictions after ~4% of the space; prospective validation |
| Jain et al. 2024 ICLR workshop [99]; Qin et al. 2025 ICML (NAIAD) [100]; Zhu et al. 2025 arXiv [101] | Pair-level combinatorial perturbation (imaging; CRISPR pairs; HIV 356×356 knockdowns) | Active search on prediction error / maximum predicted effect / ensemble uncertainty + exploration | More interactions than random; pair-level, not query-level selection. Zhu 2025 uses the SPOKE KG (disallowed here) |
| Reker & Schneider 2015 [102] | Review: AL in drug discovery | — | Background |

### 5.4 Strategic query selection in GI mapping specifically (what exists)

- **Heuristics, not adaptive design:**
  - Costanzo 2010: random with respect to function, preference for fitness defects [4].
  - Kuzmin 2018: span biological processes and digenic network features [26].
  - Billmann 2026: functionally diverse, range of expression and SMF, availability [1].
  - E-MAPs: dense functional subsets [29].
  - Costanzo 2021: "functionally representative" pairs [28].
- **Pair-level prioritisation:** Ulitsky 2009 [33]; Srivas 2016 (yeast-guided) [65]; Deshpande 2013 [66].
- **Hub prediction for query choice:** Koch 2012 [5]. Billmann 2026 suggests DepMap fitness features "can help" identify hubs [1].
- **Not found:** an iterative, model-guided query-line selection study for SGA or HAP1 GI maps (searched 2026-10). Our Q2 appears to be new.

### 5.5 Known pitfalls of active selection

- **Sampling bias.** Actively acquired training sets are not i.i.d.:
  - Risk estimates on the acquired set are biased (Farquhar, Gal & Rainforth, ICLR 2021) [103].
  - Selection can lock onto wrong hypotheses (Dasgupta 2011) [104].
  - **Always evaluate on the fixed, never-selectable test queries**, not on the shrinking remaining pool. The remaining pool's prevalence and difficulty drift with selection.
- **Model–acquisition coupling.** Sets acquired with one model can transfer poorly to successor models, and AL sometimes underperforms random (Lowell, Lipton & Wallace, EMNLP 2019) [105].
  - Random is a strong baseline at small budgets (Munjal et al., CVPR 2022 [106]; Hacohen 2022 [17]).
  - Re-fit the final predictor on the acquired set and compare against the same predictor on random sets.
- **Cold start.** With few screened queries, model uncertainty is unreliable; IterPert found its model kernel poorly aligned in early rounds [9].
  - Use prior or feature kernels, or diversity, for the first batches. Switch to model-based criteria later.
- **Degree confounding.** Pooled AUPRC rewards predicting hubs. A strategy that screens hubs may lower the test-set prevalence of "easy" positives in unscreened rows, or raise it.
  - Report per-query and pooled AUPRC.
  - Report lift over a degree baseline.
- **Batch effects across rounds.** Real screens drift by batch:
  - SGA co-batch filter removed 319 strains [upstream].
  - HAP1 clonal biases include ploidy and chromosome-arm effects [2].
  - IterPert evaluated batch-effect robustness [9].
  - Simulations should consider adding round-specific noise.
- **Tuning leakage.** Selection hyper-parameters (kernel weights, uncertainty temperature, batch diversity weight) must be tuned only on development-pool queries. Per the project rules, never on test queries.
- **Stopping rules.** Temerinac-Ott et al. [13] estimate accuracy from AL traces. Any stopping rule must use only screened data, e.g. cross-validation among screened queries.

---

## 6. Datasets: yeast source vs human analogue (with circularity flags)

Yeast sources in *italics* are already registered in the upstream DATA.md (versions, checksums, licences there).

| Feature | Yeast source | Human analogue | Licence / notes | Circularity / analogue risk |
|---|---|---|---|---|
| Coexpression | *COXPRESdb Sce-r.c3-0 / Sce-m.c4-0* [107]; SPELL (Hibbs 2007) [108] | COXPRESdb v8 Hsa compendia [107]; ARCHS4 (Lachmann 2018) [109]; GTEx/CCLE | COXPRESdb CC BY; ARCHS4 free [licence unverified] | Compendia include genetic-perturbation samples (deletion/knockdown series). Upstream removed perturbation studies in yeast; do the same for human. Never use Kemmeren 2014 deletion transcriptomes (yeast) or Perturb-seq (human) as "coexpression" |
| Physical interactions (systematic) | *Michaelis 2023 AP-MS: 4,159 pull-downs, 3,927 proteins, 31,004 interactions* [110]; Gavin 2006 / Krogan 2006 AP-MS [111]; Y2H (Yu 2008) [112]; PCA (Tarassov 2008) [113] | BioPlex 3.0: 293T 118,162 interactions / 14,586 proteins; HCT116 5,522 IPs [114]. HuRI: ~53,000 binary PPIs [115]. OpenCell: 1,310 endogenously tagged proteins, AP-MS + imaging [116] | BioPlex, HuRI, OpenCell freely downloadable [licences unverified]; Michaelis CC BY 4.0 | Avoid literature-curated PPI/complex sets (BioGRID physical, IntAct, CORUM, Complex Portal). They are text-curated and correlated with GI literature. Billmann used CORUM only for evaluation |
| Protein language-model embeddings | *ESM-2 650M on SGD proteins* (Lin 2023) | ESM-2 on UniProt human proteome | MIT | Clean. Avoid text-augmented gene embeddings (GenePT, BioGPT, LLM embeddings) |
| Paralogy | *Ensembl Fungi Compara 116*; YGOB ohnologs | Ensembl Compara human paralogs (identity, family size); ohnologs (OHNOLOGS db) | EMBL-EBI open | Clean. HAP1 paralog enrichment: 7× at ≥20% identity, 60× at ≥90% [1] |
| Abundance | *PaxDb 5.0 single datasets (Mueller 2020)* [117]; Ho 2018 unified | PaxDb 5.0 human (cell-line datasets); HAP1 or other cell-line proteomes; HAP1 RNA-seq TPM (Billmann Data S2) | CC BY 4.0 | Do not use PaxDb *integrated* sets (weighted by STRING consistency) |
| Single-mutant fitness | *Costanzo 2016 SMF* (from control screens) | HAP1 SMF from 39 WT TKOv3 screens (Billmann Data S2, rich and minimal) [1]; DepMap Chronos/CERES [61, 62]; Hart 2015 HAP1 [38]; Blomen 2015 [37] | DepMap CC BY 4.0 [unverified]; Billmann data academic | SMF strongly predicts degree in both species. HAP1 SMF comes from WT control screens of the GI study, not from GIs. Acceptable, but say so |
| Co-essentiality (fitness profile across backgrounds) | No exact yeast analogue (closest: chemogenomic co-fitness) | DepMap co-essentiality (Wainberg 2021) [60]; Billmann ED scores [1] | — | **Tier 2 / flag.** Co-essentiality is a genetic-background × gene fitness correlation; its overlap with HAP1 GI profiles is low (Jaccard < 0.1) but non-zero. Billmann ED integrates DepMap with qGI, so it is **GI-derived: banned** |
| Chemogenomic profiles | *Hoepfner 2014 HIP/HOP* (tier 2); Hillenmeyer 2008 (1,144 assays) [118] | Olivieri 2020 RPE1, 31 screens / 27 genotoxins [119]; NALM-6 ChemoGenix [40]; other drug-modifier CRISPR screens | Hoepfner CC0 | Tier 2. Phenotype-entangled with GI profiles. Deletion-collection background mutations create artefactual co-fitness [upstream] |
| Conservation / constraint | *Ensembl phylogenetic profiles; 1002 Yeast Genomes dN/dS and LoF*; AlphaFold pLDDT/disorder | gnomAD constraint (pLI, LOEUF; Karczewski 2020) [120]; gene age (Liebeskind 2016); s_het; AlphaFold | gnomAD open | Clean. Billmann used gene age, s_het, pLI, RVIS as degree features [1] |
| Localisation | GFP collection (Huh 2003) [upstream] | Human Protein Atlas subcellular (Thul 2017) [121]; OpenCell [116] | HPA CC BY-SA 3.0 [unverified version] | HPA locations are image-based (OK). Avoid GO cellular component |
| Orthology (for cross-species analysis only) | Ensembl Compara Sc–Hs | Same; Billmann used 1:1 and N:1 orthologs [1] | — | Using yeast *labels* as features for human prediction is cross-species GI transfer, which the human model should flag. For yeast-internal simulations orthology is irrelevant |
| **Banned** | BioGRID genetic, E-MAPs, DRYGIN, *S. pombe* GIs, STRING, YeastNet, GeneMANIA, GO, SGD phenotypes, text embeddings | BioGRID/SynLethDB (incl. text-mined SLs) [122], STRING, HumanNet, GO, Reactome/KEGG (curated pathways; GI-informed), KG-based (SPOKE, Hetionet) | — | — |

---

## 7. Implications for the design

1. **Calibrate labels across species, not just fractions.**
   - The yeast strong-SL label (ε < −0.12, P < 0.05; ≈1.1% prevalence; single-screen precision ≈54%, recall ≈39%) is a close analogue of the HAP1 standard negative GI (1.2% density; precision 0.61, recall 0.37) [1, 3].
   - State this explicitly so yeast AUPRC curves can be read as HAP1 forecasts.
   - Also report results for a stricter yeast label (ε < −0.35) as an analogue of HAP1 "extreme SL" (0.14%).
2. **Translate budgets two ways and test which is right.**
   - Report results against both the query fraction f and the absolute number of screened queries n.
   - HAP1 equivalents: f × ~9,500–11,000 expressed genes, i.e. 5% ≈ 500, 20% ≈ 2,000 lines, vs n directly (Billmann's 222 ≈ 7% of the yeast NxN query count).
   - Suggested check: subsample the yeast gene universe (e.g. random 50% of genes as both queries and arrays) and see whether AUPRC tracks f or n. This decides whether HAP1 needs ~2,000 or a few hundred lines for "20%-equivalent" performance.
3. **Distinguish the two Q2 objectives and report both.**
   - (A) AUPRC on the fixed test queries; this is the primary endpoint.
   - (B) Number of strong SLs found by the screened queries themselves, and the total "known" (found + predicted at a fixed precision).
   - The literature shows these trade off [15, 16]. A strategy can win one and lose the other.
4. **Baselines at matched budgets.**
   - Random (many seeds).
   - Degree-greedy: highest predicted SL degree from SMF/expression/array-side degree, i.e. Costanzo/Billmann-style hub seeking [4, 5].
   - A "Billmann-like" heuristic: stratify on expression and SMF bins and spread over coexpression/AP-MS clusters [1].
   - Pure diversity: k-center / k-means / TypiClust in the allowed feature kernel [17, 88].
   - Then model-based adaptive strategies (item 5).
5. **Strategies worth trying, in rough order of prior support.** The first two are the most literature-supported for objective A at low budget.
   - Transductive experimental design or facility location in the stacked kernel [20, 90].
   - Adaptive column sampling: pick genes whose *observed array-side partial column* has the largest residual after projection onto screened profiles [18, 19].
   - KRR/GP posterior-variance reduction over the unscreened pool (V-optimal design) [94, 95].
   - DPP/k-DPP batches for diversity [91].
   - Hybrid: diversity for the first rounds (cold start [9, 17]), then uncertainty or expected error reduction.
   - Optional discovery arm: UCB on predicted degree [84, 85].
6. **Batch mechanics.**
   - Batches of 10 with re-fitting each round. Use greedy conditional selection *within* a batch (e.g. update the posterior or residual after each pick) to avoid ten near-duplicates; this is the BatchBALD / k-center rationale [86, 88].
   - Record runtime per round; IterPert notes real rounds take ~1 month [9].
7. **Evaluation hygiene.**
   - Evaluate every strategy on the same fixed 20% test queries, never on the remaining pool [103].
   - Treat AB/BA pairs between a screened and a test gene as observed (as specified) and make sure test-row labels never enter selection.
   - Bootstrap over queries, not pairs.
   - Report pooled and per-query AUPRC, lift over prevalence, and lift over the degree-product baseline (degree dominates pooled metrics [5, upstream]).
8. **Use the array-side partial profile in selection as well as prediction.** It is non-circular for selection, because it consists of measured interactions with *screened* queries only.
   - It is the main advantage yeast SGA and HAP1 genome-wide libraries share [4, 1].
   - In HAP1 it is missing for essential library genes (masked); in yeast, essential genes are not on the NxN array. Model this asymmetry when translating.
9. **Do not over-sell paralogs or ortholog transfer.**
   - Paralogs are strongly enriched but a tiny fraction of pairs [1, 56].
   - Yeast-to-human pair-level conservation is ~2-fold enrichment, ~7% overlap [1, computed].
   - The recommendation should rest on transferable *structure*: degree predictors, module coverage, and the scaling of accuracy with screened queries.
10. **Expected magnitude of gains.**
    - Biology precedents report 3–5× sample savings for model-guided perturbation choice [9] and large gains for hit discovery [12, 14].
    - Generic AL benchmarks often find random competitive at small budgets [15, 105, 106].
    - A realistic prior for "5% chosen ≈ 20% random" is *unlikely* for AUPRC on held-out queries; a more plausible outcome is a 1.3–2× budget saving. This is my judgement, not a literature figure. Pre-register the comparison and report confidence intervals over seeds.
11. **Pitfalls specific to simulation realism.**
    - Single-screen noise: labels are themselves ~50% precise.
    - The reciprocal (AB/BA) ceiling.
    - Batch effects per round: optionally add noise at the round level.
    - Linkage artefacts: chromosomal neighbours in yeast; chromosome-arm and chr15 biases in HAP1 [2].
    - Query availability: Billmann used off-the-shelf Horizon lines for 83% of queries [computed]. A strategy that needs arbitrary genes is constrained in practice (~7,500 catalogue models; essential genes need hypomorphs) [1, 43].
12. **What prediction can realistically contribute (for Q3 framing).**
    - At 20% screened our yeast predictor reaches AUPRC ≈0.08 at prevalence 0.011 (~7× lift). This is useful for *prioritising validation and next queries*, not for replacing screens.
    - The HAP1 network is far from saturated (<5% of pairs assayed [1]). Prediction plus adaptive selection is most valuable for choosing the next query lines and ranking candidate SLs for pair-library validation (combinatorial CRISPR [46–54]).
13. **Feature admissibility audit before use.** Every human-analogue feature must pass the §6 table:
    - No STRING, GO, Reactome/KEGG, curated complexes, or text embeddings.
    - Coexpression without perturbation studies.
    - Co-essentiality and chemogenomics only as tier 2.
    - Never Billmann ED/qGI-derived features.

---

## References

1. Billmann M, Costanzo M, Zhang X, Hassan AZ, …, Moffat J, Myers CL, Boone C, Andrews BJ. Global genetic interaction network of a human cell maps conserved principles and informs functional interpretation of gene co-essentiality profiles. *Cell* 189:3736–3759.e58 (2026). doi:10.1016/j.cell.2026.03.044. PMID 42049019; PMC13281610. Preprint doi:10.1101/2025.06.30.662193. Data: https://boonelab.ccbr.utoronto.ca/supplement/billmanncostanzo2026/
2. Billmann M, Costanzo M, Rahman M, …, Boone C, Moffat J, Myers CL. Quantitative analysis of genetic interactions in human cells from genome-wide CRISPR-Cas9 screens. *bioRxiv* (2025). doi:10.1101/2025.06.30.662330 (PMC12803112). Journal version not found [unverified].
3. Costanzo M, et al. A global genetic interaction network maps a wiring diagram of cellular function. *Science* 353:aaf1420 (2016). doi:10.1126/science.aaf1420
4. Costanzo M, et al. The genetic landscape of a cell. *Science* 327:425–431 (2010). doi:10.1126/science.1180823
5. Koch EN, et al. Conserved rules govern genetic interaction degree across species. *Genome Biol* 13:R57 (2012). doi:10.1186/gb-2012-13-7-r57
6. Tischler J, Lehner B, Fraser AG. Evolutionary plasticity of genetic interaction networks. *Nat Genet* 40:390–391 (2008). doi:10.1038/ng.114
7. Dixon SJ, et al. Significant conservation of synthetic lethal genetic interaction networks between distantly related eukaryotes. *PNAS* 105:16653–16658 (2008). doi:10.1073/pnas.0806261105
8. Ryan CJ, et al. Hierarchical modularity and the evolution of genetic interactomes across species. *Mol Cell* 46:691–704 (2012). doi:10.1016/j.molcel.2012.05.028
9. Huang K, Lopez R, Hütter J-C, Kudo T, Rios A, Regev A. Sequential optimal experimental design of perturbation screens guided by multi-modal priors (IterPert). *bioRxiv* (2023). doi:10.1101/2023.12.12.571389 (RECOMB 2024 venue [unverified])
10. Naik AW, Kangas JD, Langmead CJ, Murphy RF. Efficient modeling and active learning discovery of biological responses. *PLoS One* 8:e83996 (2013). doi:10.1371/journal.pone.0083996
11. Naik AW, Kangas JD, Sullivan DP, Murphy RF. Active machine learning-driven experimentation to determine compound effects on protein patterns. *eLife* 5:e10047 (2016). doi:10.7554/eLife.10047
12. Kangas JD, Naik AW, Murphy RF. Efficient discovery of responses of proteins to compounds using active learning. *BMC Bioinformatics* 15:143 (2014). doi:10.1186/1471-2105-15-143
13. Temerinac-Ott M, Naik AW, Murphy RF. Deciding when to stop: efficient experimentation to learn to predict drug-target interactions. *BMC Bioinformatics* 16:213 (2015). doi:10.1186/s12859-015-0650-9
14. Tosh C, Tec M, …, Tansey W. A Bayesian active learning platform for scalable combination drug screens (BATCHIE). *Nat Commun* 16 (2025). doi:10.1038/s41467-024-55287-7 (author list beyond Tosh/Tansey [unverified])
15. Mehrjou A, Soleymani A, Jesson A, Notin P, Gal Y, Bauer S, Schwab P. GeneDisco: a benchmark for experimental design in drug discovery. *ICLR* (2022). arXiv:2110.11875
16. Sutherland DJ, Póczos B, Schneider J. Active learning and search on low-rank matrices. *KDD* 2013:212–220. doi:10.1145/2487575.2487627 [DOI unverified]
17. Hacohen G, Dekel A, Weinshall D. Active learning on a budget: opposite strategies suit high and low budgets. *ICML* (PMLR 162) (2022). arXiv:2202.02794
18. Krishnamurthy A, Singh A. Low-rank matrix and tensor completion via adaptive sampling. *NeurIPS* 26:836–844 (2013). arXiv:1304.4672
19. Bhargava A, Ganti R, Nowak R. Active positive semidefinite matrix completion: algorithms, theory and applications. *AISTATS* (PMLR 54):1349–1357 (2017)
20. Yu K, Bi J, Tresp V. Active learning via transductive experimental design. *ICML* (2006). doi:10.1145/1143844.1143980 [DOI unverified]
21. Feng Y, et al. Benchmarking machine learning methods for synthetic lethality prediction in cancer. *Nat Commun* 15:9058 (2024). doi:10.1038/s41467-024-52900-7 [upstream]; correction 16:1619 (2025)
22. Pandey G, et al. An integrative multi-network and multi-classifier approach to predict genetic interactions. *PLoS Comput Biol* 6:e1000928 (2010). doi:10.1371/journal.pcbi.1000928 [upstream]
23. Tong AH, et al. Systematic genetic analysis with ordered arrays of yeast deletion mutants. *Science* 294:2364–2368 (2001). doi:10.1126/science.1065810
24. Tong AH, et al. Global mapping of the yeast genetic interaction network. *Science* 303:808–813 (2004). doi:10.1126/science.1091317
25. Baryshnikova A, et al. Quantitative analysis of fitness and genetic interactions in yeast on a genome scale. *Nat Methods* 7:1017–1024 (2010). doi:10.1038/nmeth.1534 [upstream]
26. Kuzmin E, et al. Systematic analysis of complex genetic interactions. *Science* 360:eaao1729 (2018). doi:10.1126/science.aao1729
27. Kuzmin E, et al. Exploring whole-genome duplicate gene retention with complex genetic interaction analysis. *Science* 368:eaaz5667 (2020). doi:10.1126/science.aaz5667
28. Costanzo M, et al. Environmental robustness of the global yeast genetic interaction network. *Science* 372:eabf8424 (2021). doi:10.1126/science.abf8424
29. Schuldiner M, et al. Exploration of the function and organization of the yeast early secretory pathway through an epistatic miniarray profile. *Cell* 123:507–519 (2005). doi:10.1016/j.cell.2005.08.031
30. Roguev A, et al. Conservation and rewiring of functional modules revealed by an epistasis map in fission yeast. *Science* 322:405–410 (2008). doi:10.1126/science.1162609
31. Díaz-Mejía JJ, et al. Mapping DNA damage-dependent genetic interactions in yeast via party mating and barcode fusion genetics. *Mol Syst Biol* 14:e7985 (2018). doi:10.15252/msb.20177985 [volume/DOI not re-checked]
32. Yao D, et al. Scalable genetic screening for regulatory circuits using compressed Perturb-seq. *Nat Biotechnol* 42:1282–1295 (2024). doi:10.1038/s41587-023-01964-9
33. Ulitsky I, Krogan NJ, Shamir R. Towards accurate imputation of quantitative genetic interactions. *Genome Biol* 10:R140 (2009). doi:10.1186/gb-2009-10-12-r140 [upstream]
34. Hart T, et al. Evaluation and design of genome-wide CRISPR/SpCas9 knockout screens. *G3* 7:2719–2727 (2017). doi:10.1534/g3.117.041277
35. Rahman M, et al. A method for benchmarking genetic screens reveals a predominant mitochondrial bias. *Mol Syst Biol* 17:e10013 (2021). doi:10.15252/msb.202010013
36. Aregger M, et al. Application of CHyMErA Cas9-Cas12a combinatorial genome-editing platform for genetic interaction mapping and gene fragment deletion screening. *Nat Protoc* 16:4722–4765 (2021). doi:10.1038/s41596-021-00595-1 (original: Gonatopoulos-Pournatzis T, et al. *Nat Biotechnol* 38:638 (2020) [not re-checked])
37. Blomen VA, et al. Gene essentiality and synthetic lethality in haploid human cells. *Science* 350:1092–1096 (2015). doi:10.1126/science.aac7557
38. Hart T, et al. High-resolution CRISPR screens reveal fitness genes and genotype-specific cancer liabilities. *Cell* 163:1515–1526 (2015). doi:10.1016/j.cell.2015.11.015
39. Aregger M, et al. Systematic mapping of genetic interactions for de novo fatty acid synthesis identifies C12orf49 as a regulator of lipid metabolism. *Nat Metab* 2:499–513 (2020). doi:10.1038/s42255-020-0211-z
40. Chatr-aryamontri A, et al. A platform for the systematic interrogation of genetic interactions in human cells: 26 gene by genome-wide double knockout screens as a proof of concept. *bioRxiv* (2025). doi:10.1101/2025.08.04.665369
41. Desjardins J, et al. Genome-wide CRISPR screens map synthetic lethal interactions across recurrent cancer driver alterations. *Cell Rep* 45:117961 (2026). doi:10.1016/j.celrep.2026.117961
42. Horlbeck MA, et al. Mapping the genetic landscape of human cells. *Cell* 174:953–967 (2018). doi:10.1016/j.cell.2018.06.010
43. Horizon Discovery / Revvity. HAP1 ready-to-go knockout cell lines (vendor page: "over 7,500 knockout models"; delivery 1–2 weeks in stock, ~12 weeks on request). https://horizondiscovery.com/en/engineered-cell-lines/ready-to-go-cell-lines (accessed 2026-10-29) [vendor claim, unverified]
44. Najm FJ, et al. Orthologous CRISPR-Cas9 enzymes for combinatorial genetic screens. *Nat Biotechnol* 36:179–189 (2018). doi:10.1038/nbt.4048
45. Shen JP, et al. Combinatorial CRISPR-Cas9 screens for de novo mapping of genetic interactions. *Nat Methods* 14:573–576 (2017). doi:10.1038/nmeth.4225; Han K, et al. Synergistic drug combinations for cancer identified in a CRISPR screen for pairwise genetic interactions. *Nat Biotechnol* 35:463–474 (2017). doi:10.1038/nbt.3834 [both not re-checked]
46. Dede M, McLaughlin M, Kim E, Hart T. Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. *Genome Biol* 21:262 (2020). doi:10.1186/s13059-020-02173-2
47. Thompson NA, et al. Combinatorial CRISPR screen identifies fitness effects of gene paralogues. *Nat Commun* 12:1302 (2021). doi:10.1038/s41467-021-21478-9
48. Parrish PCR, et al. Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. *Cell Rep* 36:109597 (2021). doi:10.1016/j.celrep.2021.109597
49. Ito T, et al. Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. *Nat Genet* 53:1664–1672 (2021). doi:10.1038/s41588-021-00967-z
50. Köferle A, et al. Interrogation of cancer gene dependencies reveals paralog interactions of autosome and sex chromosome-encoded genes. *Cell Rep* 39:110636 (2022). doi:10.1016/j.celrep.2022.110636
51. Esmaeili Anvar N, et al. Efficient gene knockout and genetic interaction screening using the in4mer CRISPR/Cas12a multiplex knockout platform. *Nat Commun* 15:3577 (2024). doi:10.1038/s41467-024-47795-3
52. Fong SH, et al. A multilineage screen identifies actionable synthetic lethal interactions in human cancers (SCHEMATIC). *Nat Genet* 57:154–164 (2025). doi:10.1038/s41588-024-01971-9
53. Fielden J, et al. Comprehensive interrogation of synthetic lethality in the DNA damage response. *Nature* 640:1093–1102 (2025). doi:10.1038/s41586-025-08815-4
54. Flister MJ, et al. Genetic mapping and predictive modeling of paralog synthetic lethality. *Cell Rep* 44:116512 (2025). doi:10.1016/j.celrep.2025.116512
55. Lord CJ, Quinn N, Ryan CJ. Integrative analysis of large-scale loss-of-function screens identifies robust cancer-associated genetic interactions. *eLife* 9:e58925 (2020). doi:10.7554/eLife.58925
56. De Kegel B, Quinn N, Thompson NA, Adams DJ, Ryan CJ. Comprehensive prediction of robust synthetic lethality between paralog pairs in cancer cell lines. *Cell Syst* 12:1144–1159 (2021). doi:10.1016/j.cels.2021.08.006
57. Benfatto S, et al. Uncovering cancer vulnerabilities by machine learning prediction of synthetic lethality (PARIS). *Mol Cancer* 20:111 (2021). doi:10.1186/s12943-021-01405-8
58. Lee JS, et al. Harnessing synthetic lethality to predict the response to cancer treatment (ISLE). *Nat Commun* 9:2546 (2018). doi:10.1038/s41467-018-04647-1
59. Jerby-Arnon L, et al. Predicting cancer-specific vulnerability via data-driven detection of synthetic lethality (DAISY). *Cell* 158:1199–1209 (2014). doi:10.1016/j.cell.2014.07.027 [not re-checked]
60. Wainberg M, et al. A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nat Genet* 53:638–649 (2021). doi:10.1038/s41588-021-00840-z
61. Tsherniak A, et al. Defining a cancer dependency map. *Cell* 170:564–576 (2017). doi:10.1016/j.cell.2017.06.010 [not re-checked]
62. Dempster JM, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biol* 22:343 (2021). doi:10.1186/s13059-021-02540-7
63. Frost A, et al. Functional repurposing revealed by comparing S. pombe and S. cerevisiae genetic interactions. *Cell* 149:1339–1352 (2012). doi:10.1016/j.cell.2012.04.028
64. Ryan CJ, Krogan NJ, Cunningham P, Cagney G. All or nothing: protein complexes flip essentiality between distantly related eukaryotes. *Genome Biol Evol* 5:1049–1059 (2013). doi:10.1093/gbe/evt074
65. Srivas R, et al. A network of conserved synthetic lethal interactions for exploration of precision cancer therapy. *Mol Cell* 63:514–525 (2016). doi:10.1016/j.molcel.2016.06.022
66. Deshpande R, et al. A comparative genomic approach for identifying synthetic lethal interactions in human cancer. *Cancer Res* 73:6128–6136 (2013). doi:10.1158/0008-5472.CAN-12-3956
67. Jacunski A, Dixon SJ, Tatonetti NP. Connectivity homology enables inter-species network models of synthetic lethality. *PLoS Comput Biol* 11:e1004506 (2015). doi:10.1371/journal.pcbi.1004506
68. De Kegel B, Ryan CJ. Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. *PLoS Genet* 15:e1008466 (2019). doi:10.1371/journal.pgen.1008466 [upstream]
69. Lehner B, et al. Systematic mapping of genetic interactions in Caenorhabditis elegans identifies common modifiers of diverse signaling pathways. *Nat Genet* 38:896–903 (2006). doi:10.1038/ng1844
70. Park Y, Marcotte EM. Flaws in evaluation schemes for pair-input computational predictions. *Nat Methods* 9:1134–1136 (2012). doi:10.1038/nmeth.2259 [upstream]
71. Pahikkala T, et al. Toward more realistic drug–target interaction predictions. *Brief Bioinform* 16:325–337 (2015). doi:10.1093/bib/bbu010 [upstream]
72. Ryan C, Greene D, Cagney G, Cunningham P. Missing value imputation for epistatic MAPs. *BMC Bioinformatics* 11:197 (2010). doi:10.1186/1471-2105-11-197 [upstream]
73. Žitnik M, Zupan B. Data imputation in epistatic MAPs by network-guided matrix completion. *J Comput Biol* 22:595–608 (2015). doi:10.1089/cmb.2014.0158 [upstream]
74. Liany H, Jeyasekharan A, Rajan V. Predicting synthetic lethal interactions using heterogeneous data sources. *Bioinformatics* 36:2209–2216 (2020). doi:10.1093/bioinformatics/btz893
75. Natarajan N, Dhillon IS. Inductive matrix completion for predicting gene–disease associations. *Bioinformatics* 30:i60–i68 (2014). doi:10.1093/bioinformatics/btu269 [upstream]
76. Stock M, et al. A comparative study of pairwise learning methods based on kernel ridge regression. *Neural Comput* 30:2245–2283 (2018). doi:10.1162/neco_a_01096 [upstream]
77. Liu Y, et al. SL²MF: predicting synthetic lethality in human cancers via logistic matrix factorization. *IEEE/ACM TCBB* 17:748–757 (2020). doi:10.1109/TCBB.2019.2909908
78. Huang J, et al. Predicting synthetic lethal interactions in human cancers using graph regularized self-representative matrix factorization (GRSMF). *BMC Bioinformatics* 20:657 (2019). doi:10.1186/s12859-019-3197-3
79. Wang S, et al. KG4SL: knowledge graph neural network for synthetic lethality prediction in human cancers. *Bioinformatics* 37:i418–i425 (2021). doi:10.1093/bioinformatics/btab271
80. Liu X, et al. PiLSL: pairwise interaction learning-based graph neural network for synthetic lethality prediction in human cancers. *Bioinformatics* 38:ii106–ii112 (2022). doi:10.1093/bioinformatics/btac476
81. Tepeli YI, et al. ELISL: early–late integrated synthetic lethality prediction in cancer. *Bioinformatics* 40:btad764 (2024). doi:10.1093/bioinformatics/btad764
82. Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nat Biotechnol* 42:927–935 (2024). doi:10.1038/s41587-023-01905-6
83. Seale C, et al. Overcoming selection bias in synthetic lethality prediction. *Bioinformatics* 38:4360–4368 (2022). doi:10.1093/bioinformatics/btac523 [upstream]
84. Garnett R, Krishnamurthy Y, Xiong X, Schneider J, Mann R. Bayesian optimal active search and surveying. *ICML* (2012). arXiv:1206.6406 [not re-checked]
85. Lyle C, et al. DiscoBAX: discovery of optimal intervention sets in genomic experiment design. *ICML* (PMLR 202) (2023). arXiv:2312.04064
86. Kirsch A, van Amersfoort J, Gal Y. BatchBALD: efficient and diverse batch acquisition for deep Bayesian active learning. *NeurIPS* (2019). arXiv:1906.08158 [not re-checked]
87. Settles B. Active learning literature survey. Computer Sciences Technical Report 1648, University of Wisconsin–Madison (2009) [not re-checked]
88. Sener O, Savarese S. Active learning for convolutional neural networks: a core-set approach. *ICLR* (2018). arXiv:1708.00489 [not re-checked]
89. Ash JT, Zhang C, Krishnamurthy A, Langford J, Agarwal A. Deep batch active learning by diverse, uncertain gradient lower bounds (BADGE). *ICLR* (2020). arXiv:1906.03671 [not re-checked]
90. Wei K, Iyer R, Bilmes J. Submodularity in data subset selection and active learning. *ICML* (2015) [not re-checked]
91. Kulesza A, Taskar B. Determinantal point processes for machine learning. *Found Trends Mach Learn* 5:123–286 (2012). arXiv:1207.6083 [not re-checked]
92. Bıyık E, Wang K, Anari N, Sadigh D. Batch active learning using determinantal point processes. arXiv:1906.07975 (2019) [not re-checked]
93. Holzmüller D, Zaverkin V, Kästner J, Steinwart I. A framework and benchmark for deep batch active learning for regression. *JMLR* 24 (2023). arXiv:2203.09410 [not re-checked]
94. Krause A, Singh A, Guestrin C. Near-optimal sensor placements in Gaussian processes: theory, efficient algorithms and empirical studies. *JMLR* 9:235–284 (2008) [not re-checked]
95. Ma Y, Garnett R, Schneider J. Σ-optimality for active learning on Gaussian random fields. *NeurIPS* (2013) [not re-checked]
96. Gittens A, Mahoney MW. Revisiting the Nyström method for improved large-scale machine learning. *JMLR* 17(117) (2016). arXiv:1303.1849 [not re-checked]
97. Chakraborty S, et al. Active matrix completion. *IEEE ICDM* (2013) [unverified details]
98. Ruchansky N, Crovella M, Terzi E. Matrix completion with queries. *KDD* (2015). doi:10.1145/2783258.2783259 [DOI unverified]
99. Jain M, Denton AK, Whitfield ST, Didolkar AR, Earnshaw B, Hartford J. Active learning to discover pairwise genetic interactions via representation learning. *ICLR 2024 workshop (MLGenX / GEM)*
100. Qin J, Wessels H-H, Fernandez-Granda C, Hao Y. Active learning for efficient discovery of optimal combinatorial perturbations (NAIAD). *ICML* (PMLR 267) (2025)
101. Zhu H, et al. Deep active learning based experimental design to uncover synergistic genetic interactions for host targeted therapeutics. arXiv:2502.01012 (2025)
102. Reker D, Schneider G. Active-learning strategies in computer-assisted drug discovery. *Drug Discov Today* 20:458–465 (2015). doi:10.1016/j.drudis.2014.12.004
103. Farquhar S, Gal Y, Rainforth T. On statistical bias in active learning: how and when to fix it. *ICLR* (2021). arXiv:2101.11665 [not re-checked]
104. Dasgupta S. Two faces of active learning. *Theor Comput Sci* 412:1767–1781 (2011). doi:10.1016/j.tcs.2010.12.054 [DOI unverified]
105. Lowell D, Lipton ZC, Wallace BC. Practical obstacles to deploying active learning. *EMNLP-IJCNLP* (2019). arXiv:1807.04801 [not re-checked]
106. Munjal P, Hayat N, Hayat M, Sourati J, Khan S. Towards robust and reproducible active learning using neural networks. *CVPR* (2022). arXiv:2002.09564 [not re-checked]
107. Obayashi T, et al. COXPRESdb v8: an animal gene coexpression database navigating from a global view to detailed investigations. *Nucleic Acids Res* 51:D80–D87 (2023). doi:10.1093/nar/gkac983 (yeast compendia: see upstream DATA.md)
108. Hibbs MA, et al. Exploring the functional landscape of gene expression: directed search of large microarray compendia (SPELL). *Bioinformatics* 23:2692–2699 (2007). doi:10.1093/bioinformatics/btm403 [not re-checked]
109. Lachmann A, et al. Massive mining of publicly available RNA-seq data from human and mouse (ARCHS4). *Nat Commun* 9:1366 (2018). doi:10.1038/s41467-018-03751-6 [not re-checked]
110. Michaelis AC, et al. The social and structural architecture of the yeast protein interactome. *Nature* 624:192–200 (2023). doi:10.1038/s41586-023-06739-5
111. Gavin AC, et al. Proteome survey reveals modularity of the yeast cell machinery. *Nature* 440:631–636 (2006). doi:10.1038/nature04532; Krogan NJ, et al. Global landscape of protein complexes in the yeast S. cerevisiae. *Nature* 440:637–643 (2006). doi:10.1038/nature04670 [not re-checked]
112. Yu H, et al. High-quality binary protein interaction map of the yeast interactome network. *Science* 322:104–110 (2008). doi:10.1126/science.1158684 [not re-checked]
113. Tarassov K, et al. An in vivo map of the yeast protein interactome. *Science* 320:1465–1470 (2008). doi:10.1126/science.1153878 [not re-checked]
114. Huttlin EL, et al. Dual proteome-scale networks reveal cell-specific remodeling of the human interactome (BioPlex 3.0). *Cell* 184:3022–3040 (2021). doi:10.1016/j.cell.2021.04.011
115. Luck K, et al. A reference map of the human binary protein interactome (HuRI). *Nature* 580:402–408 (2020). doi:10.1038/s41586-020-2188-x
116. Cho NH, et al. OpenCell: endogenous tagging for the cartography of human cellular organization. *Science* 375:eabi6983 (2022). doi:10.1126/science.abi6983
117. Huang Q, et al. PaxDb 5.0: curated protein quantification data suggests adaptive proteome changes in yeasts. *Mol Cell Proteomics* 22:100640 (2023). doi:10.1016/j.mcpro.2023.100640
118. Hillenmeyer ME, et al. The chemical genomic portrait of yeast: uncovering a phenotype for all genes. *Science* 320:362–365 (2008). doi:10.1126/science.1150021 (Hoepfner 2014: see upstream DATA.md)
119. Olivieri M, et al. A genetic map of the response to DNA damage in human cells. *Cell* 182:481–496 (2020). doi:10.1016/j.cell.2020.05.040
120. Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581:434–443 (2020). doi:10.1038/s41586-020-2308-7 [not re-checked]
121. Thul PJ, et al. A subcellular map of the human proteome. *Science* 356:eaal3321 (2017). doi:10.1126/science.aal3321 [not re-checked]
122. Wang J, et al. SynLethDB 2.0: a web-based knowledge graph database on synthetic lethality for novel anticancer drug discovery. *Database* 2022:baac030 (2022). doi:10.1093/database/baac030

*Also used for context:* Replogle JM, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* 185:2559–2575 (2022), doi:10.1016/j.cell.2022.05.013 (human analogue of deletion transcriptomes; GI-adjacent, so not to be used as coexpression). Billmann M, et al. Widespread rewiring of genetic networks upon cancer signaling pathway activation. *Cell Syst* 6:52–64 (2018), doi:10.1016/j.cels.2017.10.015.
