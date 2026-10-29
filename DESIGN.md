# DESIGN — Using the yeast GI network to design human HAP1 synthetic-lethal screens

Version: v1 draft (2026-10-29), before review. Revision history at the end.
Brief: PROMPT.md (fixed). Background: LITERATURE.md. Upstream predictor: YEAST_SLs @ 9b5dce9 (read-only, vendored in
`third_party/ysl`, see docs/recon.md). Data numbers: results/data_characterisation/summary.md.

## 0. Questions and primary endpoints (pre-registered)

| Q | question | primary endpoint |
|---|---|---|
| Q1 | How does SL prediction for unscreened queries scale with the fraction of queries screened (random)? | Pooled AUPRC on E_fixed (§3) vs n screened, per predictor |
| Q2 | Can a smaller, strategically chosen set (iterative, batches of 10) match a larger random set? | (a) AUPRC_strategy(n) − AUPRC_random(n) at matched n; (b) random-equivalent budget n_req(n) and efficiency n_req/n; (c) the brief's example: strategy at n=164 (5%) vs random at n=657 (20%), non-inferiority test |
| Q3 | Implications for HAP1: number and choice of query lines; what prediction contributes | Derived: network coverage curves (§6), predictive lift, profile of selected genes, online stopping estimator accuracy |

Label: strong SL = (ε < −0.12) ∧ (P < 0.05) on Costanzo 2016 Science Data File S1 NxN (gene level, one `_sn` query strain
per ORF, DMA arrays). Prevalence 0.0114. Positive interactions and quantitative ε prediction are out of scope.

## 1. Data and universe
- 3,283 query genes (Q), 3,827 array genes (A); Q∩A = 2,652; array-only 1,175; query-only (off-array) 631.
- Only tier-1 external features with an obvious human analogue (inherited from upstream, all non-GI, no text mining,
  no GO): coexpression (COXPRESdb, genetic-perturbation studies removed) ↔ human COXPRESdb/ARCHS4; AP-MS (Michaelis
  2023) ↔ BioPlex/OpenCell; ESM-2 embeddings ↔ same model on human proteome; paralogs (Ensembl Compara) ↔ same;
  genome coordinates ↔ same; abundance (PaxDb) ↔ PaxDb human; orthology/conservation ↔ same; domains, AlphaFold ↔ same;
  1002-genomes variation ↔ gnomAD constraint; array single-mutant fitness (SMF) ↔ HAP1 single-KO fitness screens.
- Tier 2 (Hoepfner 2014 chemogenomic co-fitness ↔ human CRISPR chemogenomic screens) is a sensitivity analysis only:
  it is a fitness phenotype of the same deletion collection and is entangled with GI-profile similarity.
- Excluded: any GI data from other sources (incl. ExN/ExE/DAmP tables), Costanzo precomputed similarity/degree files,
  STRING/text mining, GO, knockout proteomics or other yeast-only resources.
- Essential genes: S1 NxN contains none (all queries `_sn`, arrays DMA). This matches HAP1, where query lines must be
  viable knockouts; essential genes can only be on the "array" (library) side, which yeast NxN does not model.
  Stated as a limitation; ExN/ExE are not used.

## 2. Fixed held-out test split (done before any development; frozen)
- `scripts/01_make_split.py`, seed 20261029: 656 test query genes (T), 2,627 development pool (D).
- Stratified by number of strong SLs in the query row (zero stratum + 10 quantile bins of nonzero counts); KS p=0.999.
- Within strata, test genes were drawn preferentially from the upstream lockbox (633/656 = 96.5%), because upstream
  hyper-parameters were tuned on its non-lockbox splits. Residual contamination (23 genes) is reported, not hidden.
- T is never screened, never selected, never used for tuning, stopping or strategy design. All strategies are scored on
  the same T. Test rows are scored only in the final runs (§8).

## 3. Observation model, symmetric pairs and evaluation sets
**Screening.** Screening gene g (g ∈ D) reveals its whole query row: ε(g, a) and P(g, a) for all arrays a.
**Observed pairs.** A gene pair {x, y} is observed iff x ∈ S or y ∈ S (S = screened set), regardless of side: if x is
screened, cell (x, y) is known, and the reverse cell (y, x) of an unscreened y is treated as observed too (its
reciprocal measurement ε(x, y) is available to the predictor as a feature, as in upstream).
**Unit of evaluation and symmetric pairs.** The evaluation unit is the directed cell (t, a): query t ∈ T against array
a, labelled by that cell's own measurement. We do NOT merge AB/BA into one gene-pair label, because the two directions
are separate experiments with weak agreement (P(BA SL | AB SL) = 0.145, κ = 0.14); merging would change the label
definition. When both t and a are test genes, (t, a) and (a, t) are two evaluation cells with their own labels.
Self pairs are excluded. Cells with missing ε/P are excluded.
**E_fixed (primary).** All measured cells (t, a) with t ∈ T and a ∈ A \ D (array-only genes and on-array test genes).
No a in E_fixed can ever be screened, so every cell is unobserved under every strategy and budget; E_fixed is identical
for all strategies, budgets and seeds, which removes the confound that selecting hubs would remove hub columns (and
their many SLs) from a budget-dependent evaluation set. Size ≈ 656 × ~1,700 ≈ 1.1 M cells (exact count and
prevalence reported before any scoring).
**E_unobs(S) (secondary).** All measured cells (t, a), t ∈ T, a ∉ S (budget- and strategy-dependent). Reported to show
that conclusions do not hinge on the restriction to non-dev arrays.
**Metric.** Pooled AUPRC (average precision, sklearn definition, ties broken by stable sort) over the evaluation set.
Secondary: per-query mean AUPRC; precision@k per query (k = 10, 50); AUROC; lift over prevalence.

## 4. Predictors
Upstream registered methods, run through `hsl.upstream` (ScreenView from an arbitrary screened set; all
hyper-parameters frozen at upstream values except where §4.1 says otherwise):
- `B1_array_rate` (array degree; GI-only baseline), `stack_lr` (degree + reciprocal ε + kNN; fast),
  `stack_lr_wr2` (multi-kernel fold-in + LR; upstream "simple" method), `m1_hgb` (upstream best tier-1),
  `m1_hgb_t2` (tier-2 sensitivity only).
- Known small-n issue: `m1_hgb` fails at n = 10 with default `prop_folds = 10`; we use `prop_folds = min(10, ⌊n/2⌋)`…
  rule fixed now: `prop_folds = 5` when n < 30, default otherwise. Cross-fit `n_folds = min(10, n)`.
- Memory: n ≥ 2,000 exceeds 6 GB for stack_lr_wr2/m1_hgb; those runs use ≤ 12 GB and run alone.
### 4.1 Predictor development (auto-research)
Allowed only on the inner dev splits (§7). Targets: the low-budget regime (n ≤ 330), where the upstream recipes were
never tuned (they were tuned at 20%), and the effect of non-random (selected) training sets (covariate shift:
selected sets over-represent hubs; candidate fixes: importance weighting by selection propensity, degree-calibrated
offsets). Any new predictor must beat its parent on all inner splits (paired bootstrap CI over validation queries
excludes 0 on ≥ 2/3 inner splits, mean gain ≥ 2%) and pass the leakage poisoning test (§9) before it is used.

## 5. Q1 — scaling with fraction screened (random selection)
- Budgets n ∈ {10, 33, 66, 100, 164, 250, 328, 492, 657, 985, 1313, 1970, 2627} (0.3%–80% of 3,283 queries; 2,627 =
  whole dev pool).
- Random draws from D, seeds `q1_draw_seeds` (10 seeds; nested draws: the n-sets of one seed are nested prefixes of one
  permutation, so curves are monotone in information).
- Predictors: B1, stack_lr, stack_lr_wr2 at all budgets × 10 seeds; m1_hgb at all budgets × 5 seeds.
- Report: AUPRC (mean, SD over draws, 95% CI from a joint test-query bootstrap), per stratum (on-array vs off-array
  test genes; SL-count tertile of the test query), runtime per fit. Fit a saturating curve (e.g. AUPRC(n) = a − b·n^−c)
  descriptively; extrapolation is labelled as such.

## 6. Q2 — strategic iterative selection
**Protocol.** Start S_0 = ∅. Each round chooses a batch of 10 genes from D \ S using only the information an
experimenter would have (§6.1), then reveals their rows. Budget checkpoints n ∈ {33→ use 30, 60, 100, 160, 250, 330,
490, 660} (multiples of 10; random curve uses the same values plus 980, 1310 for interpolation).
**Iterative vs static.** Every adaptive strategy is paired with a static variant that picks all queries at round 0
from external features only, to measure the value of iteration (iteration costs months in HAP1).
### 6.1 What the selector may see (enforced by a `SelectorView`)
Screened rows restricted to columns of non-test arrays; external features and SMF of non-test genes; the candidate list
D \ S. No test gene identifier, no test-gene array column, no unscreened row. Enforced by construction plus a poisoning
test: permuting or poisoning all test-gene data and all unscreened-row data must leave every selection bit-identical.
### 6.2 Candidate strategies (refined on inner dev splits; ≤ 3 frozen for the final)
- **R0 random** (comparator).
- **H hub-first:** pick candidates with the highest predicted query-side SL degree: early rounds from an external
  proxy (low SMF, high abundance/coexpression degree, AP-MS degree — combined by a model fitted only on screened rows
  once ≥ 30 are screened), later from observed array-side SL rate against screened queries (empirical-Bayes shrunk)
  plus the proxy.
- **C coverage/representativeness:** facility-location / k-center over a combined gene kernel (external kernels +
  array-side GI profile), maximising similarity coverage of the unscreened non-test pool D \ S.
- **V variance reduction (optimal design for the fold-in):** treat the kernel-ridge fold-in as a GP; greedily pick genes
  that most reduce the summed posterior variance of the unscreened non-test pool (A-optimal), optionally weighted by
  predicted degree (**VH**, hybrid), since AUPRC is dominated by high-degree rows.
- **U uncertainty:** disagreement among bootstrap/cross-fit predictors on candidate rows.
- **Static variants** of H, C, V (external features only, chosen at round 0).
### 6.3 Comparisons and statistics
- Seeds: random 10 trajectories; each stochastic strategy 5 trajectories (`q2_seeds`); deterministic strategies are run
  once per inner split (dev) and once in the final, with tie-breaking seeded.
- (a) ΔAUPRC(n) = strategy − random at matched n, CI by joint bootstrap over test queries (1,000 reps, same weights for
  all methods) combined with trajectory variance (bootstrap over trajectories within each replicate).
- (b) Random-equivalent budget n_req(n): invert the mean random curve by monotone interpolation in log n; CI from the
  bootstrap. Efficiency = n_req/n.
- (c) Brief's example, non-inferiority: strategy at n=160 vs random at n=660; "matches" iff the lower 95% bound of
  (AUPRC_strat − AUPRC_rand)/AUPRC_rand > −0.05. Holm correction across the ≤ 3 final strategies for (a) at n=160 and
  (c).
- Predictor for Q2: `stack_lr_wr2` (primary, cost) and `m1_hgb` (confirmation at n ∈ {60, 160, 330, 660}); strategy
  rankings must agree for a conclusion to be headline.
## 6.4 Q3 — HAP1 translation
- **Network coverage** (evaluation-only use of the full matrix): fraction of all strong-SL cells and of all SL gene
  pairs (SL in ≥ 1 measured direction) that are directly observed given S, vs n, for random and each strategy.
- **Prediction contribution:** at each n, the number of additional true SLs among the top-k predictions for unscreened
  genes at precision ≥ 0.2 / 0.5 (estimated on T, extrapolated to all unscreened genes), compared to the HAP1
  single-screen precision 0.61 / recall 0.37 (Billmann 2026).
- **Choice of queries:** profile genes chosen by the best strategy (SMF, expression/abundance, AP-MS degree, conservation,
  paralogs) and fit an interpretable model of "selected vs not" restricted to features with human analogues, to give
  HAP1-actionable rules (e.g. prioritise fitness-defective, broadly expressed, complex members).
- **Scaling to HAP1:** two explicit scenarios (no silent assumption): (i) fraction scaling — n_HAP1 = f × N_HAP1 query
  universe (viable KOs among ~9,500–11,000 expressed genes); (ii) count/degree scaling — match the observed network
  coverage curve shape (degree distribution Gini) rather than the fraction. Billmann 2026's 222 queries are placed on
  both scales.
- **Online stopping estimator:** at each round, cross-validated AUPRC on screened queries only (leave-batch-out) —
  what an experimenter can compute. Report its correlation with true test AUPRC and whether a plateau rule on it would
  stop near the true knee. Uses no unscreened information.

## 7. Development protocol (no test information)
- Inner splits of D (3 seeds, `inner_split_seeds`): D_val = 20% of D stratified by SL count; candidates D_cand = rest.
  Inner evaluation set E_fixed_dev = D_val × (A \ D_cand) — the exact analogue of E_fixed.
- All strategy design, hyper-parameters, predictor changes and stopping rules are chosen on inner splits only. Budgets
  in dev are the same absolute n (D_cand has ~2,100 genes, so n ≤ 660 is feasible).
- Freeze (commit hash recorded in RESEARCH_TREE.md) before the final: ≤ 3 strategies, predictor recipes, seeds, budgets,
  endpoints, comparator.

## 8. Final evaluation (one shot)
Run Q1 and Q2 on the test set with D as the candidate pool, seeds `final_*`. Scripts refuse a dirty tree and record
commit, seeds, input checksums, runtimes. Results are reported in full, including non-significant or negative ones.

## 9. Verification
- Toy-data tests with known answers: AUPRC implementation vs sklearn on toy vectors; a planted-hub toy matrix where
  hub-first must beat random and a planted-cluster toy matrix where coverage must beat random (positive controls); a
  null toy (random labels) where no strategy beats random beyond noise (negative control).
- Leakage: poisoning test rows → identical predictions; poisoning test data → identical selections; poisoning
  screened rows → predictions change.
- Sanity: random-curve monotone in n (on average); B1 ≈ upstream B1 at 20%; stack_lr_wr2/m1_hgb at n=657 random
  within the upstream lockbox range (0.076 / 0.083 ± seeds); E_fixed count/prevalence reported.
- Independent rescoring of key numbers by an agent that did not produce them, with no project imports.
- YEAST_SLs read-only: `hsl.guard` audit hook + stamp/`GIT_OPTIONAL_LOCKS=0 git status` at each milestone.

## 10. Compute and runtime reporting
M2 laptop, 8 cores, 24 GB. Jobs ≤ 3 threads, normally ≤ 6 GB; ≤ 2 concurrent. Every result row records fit/predict
and selection wall-clock. Rough cost: Q1 ≈ 10 seeds × 13 budgets × (stack_lr + stack_lr_wr2) ≈ 2 h; m1_hgb 5 × 13 ≈ 3 h;
Q2 dev ≈ 6 strategies × 3 inner splits × (66 selection rounds + 8 evaluations) ≈ several h; final Q2 ≈ 4 × 5–10
trajectories × 8 evaluations.

## 11. Threats and limitations (to be revisited)
Single-direction labels are noisy (AB/BA κ 0.14), capping achievable AUPRC; yeast NxN lacks essential genes; yeast
fractions do not map 1:1 to HAP1 (different genome size, network density, assay); upstream recipe tuned at 20% random;
selection-induced covariate shift; E_fixed excludes dev-pool arrays (checked via E_unobs).

## Revision history
- v1 (2026-10-29): first draft.
