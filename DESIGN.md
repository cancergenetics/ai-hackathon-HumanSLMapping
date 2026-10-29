# DESIGN — Using the yeast GI network to design human HAP1 synthetic-lethal screens

Version: **v3.1** (2026-10-29, frozen design), after review rounds 1 (7 reviewers), 2 (3) and 3 (1); `docs/reviews/design_r*_*.md`;
responses in §14–§15.
Brief: PROMPT.md (fixed). Background: LITERATURE.md. Upstream predictor: YEAST_SLs @ 9b5dce9 (read-only; vendored in
`third_party/ysl`, docs/recon.md). Data numbers: results/data_characterisation/summary.md.

## 0. Questions and pre-registered endpoints

| Q | question | primary endpoint (primary predictor P*, §4) |
|---|---|---|
| Q1 | How does SL prediction for unscreened queries scale with the fraction screened (random)? | Pooled AUPRC on E_fixed (§3) vs n; **n_80, n_90** = budget at which the isotonic, log-n-interpolated mean random curve reaches 80/90% of the measured AUPRC at the full dev pool (n = 2,627), with joint-bootstrap CIs; always shown with absolute AUPRC, lift, precision at 10% recall and the fraction of all strong SLs detected by the screened rows at that n |
| Q2 | Can a smaller strategically chosen set (iterative, batches of 10) match a larger random set? | (c1) **strategy@160 > random@490** (non-inferiority on the budget scale, one-sided); (a) Δ(160) = strategy@160 − random@160 > 0; reported with gap closure **G** and random-equivalent budget n_req(n) |
| Q3 | Implications for HAP1: number/choice of query lines, contribution of prediction | n_80/n_90 translated with the scaling exponent b from the universe-subsampling experiment (§7.3); warm-start gain of the best strategy over R1 (Billmann-like) continuation; SL-discovery frontier; marginal returns per +100 lines; choice-of-query rules on human-available features |

Label (fixed by the brief): strong SL = ε < −0.12 ∧ P < 0.05 on Costanzo 2016 Science Data File S1 NxN (gene level,
one `_sn` query strain per ORF, DMA arrays); prevalence 0.0114. Positive GIs and ε regression are out of scope.
The brief's example "5% vs 20%" is operationalised as budgets **160 (4.9%) vs 660 (20.1%)** of 3,283 query genes.

## 1. Data, features, admissibility
- 3,283 query genes Q, 3,827 array genes A; Q∩A = 2,652; array-only 1,175; query-only (off-array) 631.
- **Feature provenance table** (`docs/feature_provenance.md`, to be written in P0): for each store: label-free?
  GI/text/GO-free? human analogue? used by predictor / selector? Tier 1 (inherited from upstream): coexpression
  (COXPRESdb, genetic-perturbation studies removed ↔ human COXPRESdb/ARCHS4), AP-MS (Michaelis 2023 ↔ BioPlex/OpenCell),
  ESM-2 embeddings (↔ same model, human proteome), paralogs (Ensembl Compara ↔ same), abundance (PaxDb ↔ PaxDb human),
  orthology/conservation, domains, AlphaFold, 1002-genomes variation (↔ gnomAD constraint), array single-mutant fitness
  SMF (↔ HAP1 single-KO fitness).
- **Genome coordinates are excluded from the primary recipes.** Same-chromosome proximity carries an SGA meiotic-linkage
  artefact (dev rows: SL rate 0.26 at 20–50 kb vs 0.0114 inter-chromosomal) with no HAP1 counterpart. `m1_hgb` runs
  with `drop=["genome"]`; the upstream recipe with genome is a reported sensitivity. stack_lr / stack_lr_wr2 do not
  use genome features. Sensitivity: AUPRC on E_fixed minus same-chromosome cells.
- **SMF:** array-strain SMF is a single-mutant phenotype (not a GI), allowed for predictor and selector (human analogue:
  HAP1 fitness). Query-strain SMF of an unscreened gene is excluded everywhere (it only exists once the query strain is
  built/screened). For the 631 off-array genes SMF is missing; selectors may not use a "missing SMF / off-array"
  indicator as a feature or score term (yeast artefact; every HAP1 gene is in the library); routing a candidate to an
  external-only submodel when its data are absent is allowed.
- Tier 2 (Hoepfner 2014 chemogenomic co-fitness) — sensitivity only, never headline.
- Excluded: GI data from any other source (incl. ExN/ExE/DAmP), Costanzo precomputed similarity/degree files,
  STRING/text mining, GO, knockout proteomics and other yeast-only resources.
- Essential genes: S1 NxN has none. Results describe the **nonessential × nonessential submatrix**; HAP1 libraries
  also contain ~1,500 essential genes whose GIs yeast NxN does not model (limitation).
- **Upstream knowledge declared:** 633/656 test genes are the upstream lockbox, whose random-screening AUPRCs
  (m1_hgb 0.061/0.083/0.102 at 10/20/30%) are public in YEAST_SLs. They inform no decision here. Sensitivity: all
  endpoints on T minus the 23 non-lockbox test genes. Absolute AUPRCs of upstream predictors on T are therefore not
  blind estimates; comparisons between strategies with the same predictor are unaffected.

## 2. Fixed held-out test split (frozen, done before development)
`scripts/01_make_split.py`, seed 20261029: 656 test genes T, 2,627 development pool D; stratified by strong-SL count
(zero stratum + 10 quantile bins; KS p = 0.999); 633/656 drawn from the upstream lockbox. T is never screened,
selected or used for tuning/stopping/strategy design. **Dev code cannot load T rows** (`hsl.views.dev_bundle` deletes
them; only `final_bundle` with `HSL_FINAL=1` on a clean frozen commit holds them).

## 3. Observation model, symmetric pairs, evaluation sets, metric
**Screening** gene g ∈ D reveals its whole query row ε(g, ·), P(g, ·). The predictor sees all columns of screened
rows, incl. T (dev: D_val) columns (observed pairs per the brief); only the SelectorView masks them.
**Observed pairs / symmetric handling.** A gene pair {x, y} counts as observed iff x ∈ S or y ∈ S. If x is screened
and y is not, only ε(x, y) (the screened direction) exists; it is available to the predictor as a feature of cell
(y, x) (the "reciprocal" feature), but cell (y, x) is excluded from every evaluation set and is never a training label;
ε(y, x) is never revealed. The evaluation unit is the **directed cell** (t, a), labelled by its own measurement; AB
and BA are not merged (weak agreement: P(BA SL | AB SL) = 0.145, κ = 0.14; merging would change the label). When t
and a are both test genes, (t, a) and (a, t) are two cells. Self pairs and cells with missing ε/P are excluded.
Gene identity (current SGD names; all groups enumerated by script, counts in the manifest): (i) a query pair split
across D/T that is one current gene (e.g. YLR312C-B→YLR313C) → the D member is removed from every candidate pool;
(ii) arrays whose current name equals a dev-pool gene (YFL056C, YHR079C-B, YIR043C, YJL206C-A) are dropped from E_fixed;
(iii) cells whose array is the row gene under current names are dropped (self pairs). All three sets are asserted
empty after application, so E_fixed is invariant.
**E_fixed (primary).** Measured cells (t, a), t ∈ T, a ∈ A \ D (minus the identity exclusions): never observable
under any strategy, so identical for all arms. Pre-registered column strata: **E_AO** = T × array-only (≈1,170
columns, SL rate ≈0.008) and **E_TT** = T × on-array test genes (533 columns; columns are query genes — the closest
analogue of HAP1, where every library gene is a potential query). Counts/prevalence are written to the manifest
(labels only) before any scoring. Note: no E_fixed array is ever screened, so the reciprocal feature is never
available there — E_fixed measures fold-in/propensity prediction.
**E_unobs(S)** (practical, per arm; not paired): (t, a), a ∉ S; reported with prevalence and lift.
**Headline rule:** a Q2 claim is headline only if Δ(160) has the same sign on E_fixed, E_AO and E_TT, with P* and with
m1_hgb; otherwise it is reported as qualified. (E_common was dropped in v3: with R ≥ 20 trajectories it collapses to
E_fixed and depends on R.)
**Metric.** Pooled average precision with **tied scores grouped at one threshold** (sklearn
`average_precision_score` semantics; no stable-sort tie-breaking). Secondary: macro (per-query) AUPRC over test queries
with ≥ 1 positive in the set (zero-positive queries counted and reported; included in pooled); AUROC; precision@k
(k = 10, 50) per query; lift = AUPRC/prevalence. Per-stratum: on-/off-array test genes; SL-count tertile.

## 4. Predictors
- **Primary predictor P\* = `stack_lr_wr2`** (upstream multi-kernel fold-in + LR; ~12–60 s per fit), or its frozen
  successor if one passes the §8 gate before the P2 freeze of P*. All Q1/Q2/Q3 headlines use P*; n_req always inverts the random
  curve of the same predictor.
- Secondary: `m1_hgb` (genome dropped), refit on the same selected sets for 5
  trajectories per arm at {160, 660} (strategies) and {160, 490, 660} (random); baselines `B1_array_rate`, `B2_product` (array rate × query propensity) and `stack_lr`;
  tier-2 `m1_hgb_t2` sensitivity at 660 only.
- **Same recipe in every arm.** Any predictor change (e.g. covariate-shift correction) applies to all arms; the strategy
  is also reported with the unmodified predictor.
- Memory/concurrency: ≤ 4 jobs × 2 threads for n ≤ 660 (≤ 10 GB total), ≤ 2 × 3 threads above; n ≥ 1,970 runs alone
  with ≤ 12 GB.
- Covariate-shift corrections are not part of P*; if tried, it is a post-P4 sensitivity applied identically to all arms.

## 5. Budgets and seeds (single source: `config/budgets.yaml`, `config/seeds.yaml`)
- Random/Q1 grid: n ∈ {10, 30, 60, 100, 120, 160, 210, 220, 330, 410, 490, 660, 820, 980, 1310, 1640, 1970, 2627}
  (220 ≈ Billmann's 222 lines).
- Q2 checkpoints: {30, 60, 100, 160, 220, 330, 490, 660}; trajectories run to 660 (dev: to 660 or pool exhaustion).
- Warm-start checkpoints: {220, 330, 440, 490, 660} (`budgets.yaml: warm_start_checkpoints`).
- Random arm = nested prefixes of seeded permutations of the candidate pool (shared by Q1 and Q2).
- Seeds via `numpy.random.SeedSequence(base).spawn(R)`; trajectory r of every arm uses model seed r (common random
  numbers). Separate RNG streams for selection and evaluation. Dev and final streams disjoint. All named seed lists
  are appended to `config/seeds.yaml` before first use.
- R (trajectories per arm): provisional from the power pilot (§7.1) with R_random = 2·R_strategy (random is shared by
  all strategies); re-estimated from the P4 dev strategy-vs-random runs and fixed at the freeze (may only go up).
- All protocol constants (bins, B, α, margins, thresholds, caps, HAP1 constants) live in frozen `config/protocol.yaml`;
  dev-tuned strategy knobs live in `config/finalists.yaml`. Code never hardcodes them.

## 6. Q2 — strategies
**Protocol.** S_0 = ∅ (cold start) unless stated. Each round picks 10 genes from the candidate pool using only the
SelectorView, then their rows are revealed.
### 6.1 Data access and the SelectorView
- **dev_bundle** = upstream `Bundle.blind(T)` (T rows physically dropped from every row-aligned array); dev views are
  built from it without the upstream lockbox refusal (23 lockbox genes sit in D; tested). Enforced at runtime by
  `hsl.guard`: unless `HSL_FINAL=1`, unblinded bundle loads and direct opens of Costanzo bundle/processed files raise
  outside `hsl.views`.
- **SelectorView** (harness-built; selector code never reads split files; exposes only whitelisted attributes by
  construction): (a) screened rows only; (b) in those rows, columns of T (and, in dev, of D_val) are masked with the
  native "unobserved" code (Y int8 = −1; ε/P = NaN) — masked, not deleted; (c) SMF of T (and D_val) NaN; (d) candidate
  list = pool \ S; (e) whitelisted per-gene features only (tier 1 + array SMF, from `config/feature_provenance.yaml`; a
  test fails if any selector/predictor opens a store absent from it). Not exposed: DMF/DMF_sd, query SMF, n_arrays/
  n_queries, strain metadata, missingness/off-array indicators.
- **Poisoning tests** run on a fixture with a pseudo-test set drawn from D (no unblind exemption for tests): poison
  (a) pseudo-test rows, (b) pseudo-test columns of screened rows, (c) pseudo-test SMF, (d) unscreened pool rows — each
  alone leaves selections bit-identical; adaptive selectors: poisoning a visible column of a screened row changes some
  selection (positive control); static selectors: poisoning any screened row leaves selections identical. Each test has
  a broken twin it must catch.
### 6.2 Strategy set (≤ 3 frozen for the final, chosen on dev)
- **R0 random** (comparator): nested prefixes of seeded permutations.
- **R1 Billmann-like** (status quo, static): round-robin over SMF × abundance tertile bins, spread across k-means
  clusters (k = 100, fixed seed, external-kernel embedding of non-test genes; dev: also excluding D_val); trajectory r
  = seeded permutation of within-bin order.
- **H hub-first:** rank candidates by predicted query-side SL degree q̂. Prior (n < 30, unfitted): rank(1 − SMF) +
  rank(abundance), descending; off-array genes score 2·rank(abundance) (same scale). From n = 30: Gaussian ridge (λ by GCV) on
  logit((k+0.5)/(m+1)) of screened rows; inputs: shrunk array-side rate ĉ_x = (k_x + 100·π̄_S)/(n_x + 100), array SMF,
  abundance, AP-MS degree, coexpression degree, conservation, paralog flag; off-array candidates use an
  external-only submodel. Batch cap ≤ 2 genes per cluster. Dev diagnostic: Spearman(q̂, true degree) on the pool.
- **F Fisher-D (main hybrid):** maximise Σ_{x∈P} log det(δI_r + Σ_{s∈S} p̂_sx(1−p̂_sx) u_s u_sᵀ), P = on-array
  candidates ∪ array-only columns (visible), w_x = 1. u_g = (1, ũ_g): ũ = r−1 truncated-SVD row factors of the visible
  screened block (screened genes: own row factor; on-array candidates: û_g = B·v_g from their column factor v_g with
  an (r−1)×(r−1) map B fitted by least squares on screened on-array genes; off-array candidates and cold start: KRR
  from the external kernel into the same space); ũ whitened to unit mean square; δ = mean p̂(1−p̂);
  p̂_sx = clip(q̂_s·ĉ_x/π̄_S, 1e-4, 0.5). Greedy within a batch via the determinant lemma and Sherman–Morrison; M_x
  recomputed each round. r = 1 (u ≡ 1) reduces exactly to H's ordering (built-in ablation; F uses the same ≤ 2-per-cluster batch cap as H);
  r ∈ {4, 8, 16} tuned on dev. Cold start (n < 30, incl. S = ∅): p̂ ≡ const, so the objective is
  Σ_x log det(I + Σ_s u_s u_sᵀ) with ũ = top r−1 kernel-PCA components of the external kernel (whitened); the SVD space
  is used from n = 30. M_x := δI_r + Σ_{s∈S} p̂_sx(1−p̂_sx) u_s u_sᵀ.
- **C coverage:** facility location over the predictor's `mk` kernel (GI sym kernel over screened on-array genes +
  external kernels), target pool = unscreened non-test candidates; one arm (its static part = external kernel only is
  reported as the round-0 behaviour).
- **M max posterior variance:** pick argmax k(g,g) − k_gS(K_SS+λI)⁻¹k_Sg under the GI-only kernel (adaptive column
  sampling); cheap uncertainty arm. Round 0 = the first batch of the R0 permutation of the same trajectory;
  candidates without a GI column (off-array) are never picked by M (count reported).
- **Static variants** H-static, F-static (prior/external only, chosen at round 0) measure the value of iteration.
- Dropped: U (bootstrap uncertainty; cost and objective mismatch), V (redundant with F on the same geometry).
- Randomised greedy: each pick sampled uniformly from the top-20 by gain (seeded) to give trajectory variance; claims refer to
  this randomised rule; one deterministic (temperature-0) trajectory per finalist is reported alongside.
- Compute rule: no selector fits stack_lr_wr2/m1_hgb per round; allowed per round: kernel/degree modules and at most
  one stack_lr fit. Selection seconds per round logged. Switch points/cut-offs as pool fractions.
### 6.3 Scenarios and sensitivities (best ≤ 2 strategies; descriptive)
- **Warm start:** S_0 = 220 genes chosen by R1 (≈ Billmann's 222), then continue to 660 with R1-continuation (status
  quo), random, and the best strategies. Primary contrast: strategy vs R1-continuation at 440 and 660.
- **Batch size:** 10 (brief), 50, 110 (= 6 rounds to 660) and one-shot (static) at n ∈ {160, 330, 660}; the Q3
  recommendation uses the best strategy feasible in ≤ 6 HAP1 rounds. For batch b, a checkpoint n not divisible by b
  uses the first n − b⌊n/b⌋ genes of the next batch in greedy order.
- **HAP1-feasible pool:** pool minus the bottom array-SMF decile (off-array genes, lacking SMF, are kept); AUPRC lost
  reported.
- **HAP1-minimal features:** selector restricted to fitness + coexpression + protein embeddings (no AP-MS; predictor
  unchanged).
### 6.4 Statistics (confirmatory family fixed now)
- **AP engine:** exact weighted AP (scores sorted once; tie groups kept; per replicate, query weights gathered and
  cumulated at tie-group ends). One replicate array W[rep, arm, r, n] with stratified query weights w^b (11 split strata,
  B = 2,000, same w^b for all arms) drives every inference.
- Estimand Δ(n) = E_traj[AUPRC_s(n)] − E_traj[AUPRC_r(n)] on E_fixed with P*. SE² = V_Q + σ̂²_s/R_s + σ̂²_r/R_r with V_Q = Var
  over replicates of the trajectory-mean contrast; ν = SE⁴/[(σ̂²_s/R_s)²/(R_s−1) + (σ̂²_r/R_r)²/(R_r−1)] (V_Q treated as known;
  ignoring the CRN covariance is conservative); one-sided p = 1 − T_ν(Δ̂/SE). Paired predictor-change contrasts (same
  S) use Var(d_r)/R + V_Q(d̄).
- **Family** (one-sided, Holm, familywise α = 0.025): for each of ≤ 3 frozen strategies, (c1) AUPRC_s(160) >
  AUPRC_r(490) and (a) Δ(160) > 0 → ≤ 6 hypotheses. Everything else is descriptive; curves get Bonferroni bands over
  the checkpoints.
- **Brief's literal example** (strategy@160 vs random@660): estimate + CI of AUPRC_s(160) − AUPRC_r(660); gap closure
  G = (s160 − r160)/(r660 − r160) (joint-bootstrap CI; undefined if the 97.5% interval of r660 − r160 includes 0).
  "Matches 20% random" is tested as the linear contrast H0: s160 − 0.95·r660 ≤ 0 (one-sided 0.025, descriptive) and
  may be stated only if (i) it is rejected, (ii) c1 is rejected, and (iii) the comparison is discriminating: the
  point estimate of G* = 1 − 0.05·r660/(r660 − r160) is ≥ 0.5 (random gap ≥ 2 × margin). G* is reported.
- **n_req(n), n_80, n_90:** joint bootstrap (query weights + trajectories resampled within arm); per replicate PAVA
  on the trajectory-mean random curve over the dense grid, log-n interpolation; value = left endpoint of the reached
  flat segment; censored replicates coded 0/+∞ and reported (one-sided bounds); never extrapolated. Dev uses the dev
  pool maximum in place of 2,627.
- Sensitivity: two-way (query × array) bootstrap (B = 500, per-cell weights; test genes on the array side share their
  query weight) for the ≤ 6 contrasts.

## 7. Q1 / Q3 analyses
### 7.1 Q1 and power pilot
- Random curves for P*, B1, B2_product, stack_lr at all grid budgets × R_random permutations; m1_hgb at
  {30, 160, 490, 660, 1310, 2627} × 5. Report AUPRC (mean, SD, CI), strata, runtime per fit; saturating fit
  a − b·n^−c per permutation, summarised by seed-cluster bootstrap (descriptive only).
- **Power pilot** (dev split 1, random, R = 10, P*): σ_r(n), V_Q,single(n). Provisional V_Q(Δ) := 2·V_Q,single(160)·
  |D_val|/|T| (ρ = 0, conservative); choose R_s with R_r = 2R_s so SE(Δ̂(160)) ≤ g/3.5, g = r210 − r160 (a 1.3× budget
  saving); re-estimated from P4 dev strategy runs; fixed at freeze (may only go up). If V_Q alone exceeds (g/3.5)², Q2
  claims are phrased as bounds.
### 7.2 Discovery axis (co-reported with every Q2 arm)
SLs detected = strong-SL cells in screened rows (cell level, primary); pair level "detected" (SL call in a screened
row) vs "touched" (pair involves a screened gene). Frontier plot: AUPRC (E_fixed) vs SLs detected.
### 7.3 Scaling law (fraction vs count) — decides the HAP1 translation
Dev only: gene sub-universes of 50% and 25% of the non-test genes (same genes removed from query and array sides; 3
draws each; inner fold structure rebuilt inside each; T outside every sub-universe), random curves with P*, ≥ 10
permutations per draw, budget grid scaled so fractions match. External features restricted to the sub-universe and
every derived statistic (degrees, quantile ranks, clusters, kernels' normalisation) recomputed inside it where the
upstream code allows; otherwise the result is labelled a bound. Fit the slope b of log n_80 (and log n_90) on log N
(N = candidate-pool size) by OLS over 9 draws (3 per level: 100% = dev splits 1–3, 50%, 25%; levels {100, 50} only, 6 draws, if the §12
trigger fires; draw-bootstrap CI; N = candidate-pool size, proportional to N_yeast across levels so b is unaffected): b ≈ 0 → count scaling, b ≈ 1 → fraction
scaling. HAP1 translation: n_HAP1 = n_yeast · (N_HAP1/N_yeast)^b, with b = 0 and b = 1 shown as bounds. Report the
positive count of E_fixed_dev at 25% (noise check). Caveat: subsampling mimics fewer genes per module, not more modules.
### 7.4 HAP1 translation and actionability
- **Denominator:** N_yeast = 3,283 constructed query genes ↔ N_HAP1 = viable-KO query universe among expressed
  nonessential genes (range 7,500–9,500 from LITERATURE; protocol.yaml); N_yeast = 4,458 (Q ∪ A) as a sensitivity.
  Costs in lines, screens (× 1.34 replicate factor, Billmann) and sequential rounds; Billmann's 222 lines placed on the
  curves; marginal returns ΔAUPRC and Δ(SLs detected) per +100 lines.
- On-array test genes (E_TT / on-array strata) are the estimate used for HAP1 translation; all-T primary for Q1/Q2.
- **Choice of query lines:** profile genes picked by the best frozen strategy vs random on human-available features;
  interpretable rank rule.
- **Prediction contribution:** precision/recall of P* on E_TT vs the within-yeast replicate reference (reciprocal BA →
  AB on the same test × test cells); score thresholds by leave-batch-out CV on screened rows (not on T); expected true
  SLs among top-k predictions per unscreened gene.
- **Online stopping estimator:** leave-batch-out CV AUPRC on screened rows only; bias/MAE vs true test AUPRC per n,
  within-n correlation across trajectories, stopping-n of a plateau rule (frozen on dev) vs the true knee (smallest n
  with AUPRC ≥ 0.95 × AUPRC(660)).
- Optional, evaluation only, after the REPORT draft and only if compute remains (never feeds back): yeast 222-query
  degree distribution vs Billmann 2026 Data S14; Spearman of our human-feature rank rule vs HAP1 library-side negative
  degree; stricter label ε < −0.35 re-scoring of the final random curve.

## 8. Development protocol (no test information)
- Inner folds: 5 stratified folds of D. **Dev splits k = 1, 2, 3:** D_val = fold k, candidates = D \ fold k (~2,100).
  **Re-split check:** val = folds 4+5, candidates = folds 1–3 (~1,575); not independent (folds 4–5 trained the dev
  models), used only for n ≤ 330; the final is the only independent test.
- E_fixed_dev = D_val × (A \ D_cand \ T) minus identity exclusions — no test-gene columns.
- **Predictor-change gate:** on the shared R0 dev trajectories (R_r,dev = 10 per split, reused for every
  configuration), paired contrasts at n ∈ {30, 160, 330}: mean relative gain ≥ 2% (over budgets and splits), the 95% CI
  at n = 160 excludes 0 on ≥ 2/3 dev splits, same sign on the re-split check.
- **Strategy shortlist:** ≥ 3 trajectories per configuration per split; rank by the mean over the 3 splits of relative
  Δ(160); at most one variant per family (H, F, C/M) advances; ≤ 3 finalists. R1, the static counterparts of finalists, and the warm-start
  and batch-size scenarios for the best ≤ 2 run on T as descriptive arms outside the family; the feasible-pool and
  minimal-features scenarios are dev only.
- **Surrogate** (closed-form fold-in) may be used for tuning only if its within-n Spearman with P* AUPRC across random
  trajectories at n ∈ {30, 160, 330} is ≥ 0.8 and, once ≥ 5 configurations exist, Spearman(surrogate Δ, P* Δ) at
  n = 160 ≥ 0.7.
- **Ledger** (`results/LEDGER.jsonl`, append-only, one line per configuration: id, parent, config hash, split,
  outcome, metrics); the count is reported.
- **Freeze:** `config/finalists.yaml` + `config/frozen.yaml` (hashes of hsl/eval.py, hsl/views.py, hsl/guard.py,
  harness, protocol.yaml, split files, seeds, budgets, config/predictor.yaml, src/hsl/predictors/*.py (added after P2), third_party/ysl/PROVENANCE.md, config/feature_provenance.yaml
  (source; docs/feature_provenance.md rendered from it), DATA.md, yeast_sls_commit.txt, data checksum table).
  Sub-agents never edit frozen files; a failing frozen test is reported, never worked around.

## 9. Phases and gates (recorded in RESEARCH_TREE.md)
- **P0 infrastructure** (eval engine, views/guard, SelectorView, harness, launcher, protocol.yaml, tests, provenance).
- **P1:** dev random curves (P* + baselines; m1_hgb subset) + power pilot + universe subsampling (P* only).
- **P2 predictor development** for n ≤ 330 (n-dependent GI weight, propensity features, small-n degeneracies) on the
  P1 trajectories → **freeze P*** before any strategy work. If P2 changes P*, the power pilot and the
  universe-subsampling curves are re-run with the frozen P* before P3 (the §12 50%-only trigger still applies).
- **P3 ceilings (go/no-go):** oracle hub-first (true query degree; all 3 splits; free) and greedy oracle (closed-form
  fold-in; 200 random candidates per **round of 10**, top-10 by realised E_fixed_dev gain; split 1 only, to n = 160) —
  evaluation-only diagnostics. If both oracles have G < 0.5 at n = 160, de-scope C/M and say so.
- **P4:** strategies on dev (R1, H, H-static, F, F-static, C, M) → shortlist; warm start, batch size, feasible pool,
  minimal features for the best ≤ 2; R fixed → **freeze** → **Final** (Q1 + Q2 on T) → Q3 analyses → REPORT →
  optional items.

## 10. Final evaluation
`scripts/run_final.py` refuses unless `config/finalists.yaml` is committed, the tree is clean, frozen hashes match, and
`results/final/.lock` is absent; writes the lock (git hash); predictions saved before scoring; a crash re-run is a
ledger line. Candidate pool = D (minus identity removals); T evaluated.

## 11. Verification
- Toy tests with known answers: AP literals; all-tied → prevalence; cell-order permutation → identical AP; sklearn
  agreement (unweighted exact; integer-weighted = duplicated cells); stable-sort mutant caught; planted-hub toy (H beats
  random), planted-cluster toy (C beats random), null toy (no strategy beats random beyond noise); n_req on a known
  curve incl. censoring; bootstrap coverage on a toy; Holm on known p-values.
- Leakage: per-channel poisoning (§6.1) with broken twins; prediction poisoning; runtime guard on unblinded loads.
- Data assertions: Y = −1 ⇔ ε or P missing; label rule reproduces Y on D rows; gene counts 3,283/3,827/2,652/1,175/631;
  no E_fixed cell with a dev-pool array; identity-exclusion sets empty after application.
- End-to-end controls at each milestone: shuffled screened labels → AUPRC ≈ prevalence; oracle back door → ≈ 1.
- Sanity (dev only): random curve monotone on average; P* at 660 plausible vs upstream dev numbers; upstream-protocol
  reproduction (cells incl. reciprocal) on dev. T-based upstream comparisons only in the final report, labelled post hoc.
- Independent rescoring from saved predictions (no project imports) by an agent that did not produce them: one dev
  milestone and the final (≥ 3 strategy/n cells + the Q2 verdict). REPORT numbers rendered from `results/claims.json`
  via `{{claim:id}}` placeholders (test: all placeholders resolve).
- YEAST_SLs read-only: `hsl.guard` + `scripts/check_readonly.sh` at each milestone.
- Repo hygiene: a DATA.md row (source, version, checksum, licence) per dataset is asserted at each milestone; raw and
  restricted data paths are git-ignored and a test fails if any is tracked.

## 12. Compute
M2, 8 cores, 24 GB. Concurrency per §4. Rough v3 estimate (round-2 review): 70–90 job-h before cuts. Wall-clock
triggers, applied in order and logged in the tree: if P1 > 6 h wall → universe subsampling to the 50% level only; if
P1–P3 > 14 h → m1_hgb Q1-only budgets {30, 1310} dropped (the §4 Q2 set {160, 490, 660} × 5 and 2627 kept); if P4 >
12 h → stop dev trajectories of strategies eliminated after split 1, drop batch-50; if the projected final (P4 timings
× R) exceeds 10 h, R is set at the pilot minimum *before* the freeze; R is never changed after it. Every result row
records fit/predict and selection wall-clock.

## 13. Threats and limitations
Single-direction labels are noisy (κ 0.14), capping AUPRC; NxN lacks essential genes; yeast→HAP1 scaling tested only by
subsampling; upstream recipes tuned at 20% random; selection-induced covariate shift; E_fixed dominated by array-only
columns (strata reported); T not blind for upstream predictors (declared); randomised-greedy estimand.

## 14. Response to round-1 reviews (A = accepted, P = partly, R = rejected)
| reviewer | item | decision |
|---|---|---|
| ext. PI | 1 n_80/n_90; 2 scaling experiment; 3 gap closure G; 4 oracle ceilings first; 5 R1 + warm start; 6 one primary predictor; 7 discovery axis; 8 reciprocal reference; 9 E_fixed characterisation + agreement rule; 10 batch size; 11 drop U, add M; 12 phase gates; 13 trajectory variance | A (all) |
| ext. PI | 14 availability subset; 15 strict label; 16 per-query rules | P: 14 replaced by HAP1-feasible pool (domain #2); 15 optional secondary; 16 A |
| statistics | 1 tie-aware AP; 2 macro rules; 4 variance decomposition; 5 stratified bootstrap (v3: exact weighted AP, B = 2,000); 8 CRN; 9 power pilot; 11 family/Holm; 12 budget-scale margin (160 vs 490); 13 isotonic censored n_req; 14 no T columns in dev eval; 15 disjoint folds + re-split check (n ≤ 330) + ledger; 16 same recipe; 17 thresholds not on T; 18 stopping metrics | A |
| statistics | 3 E_common (dropped in v3, §15); 6 two-way bootstrap; 7 model-seed variance; 10 m1_hgb on same trajectories; 19 per-seed fit | A |
| data | 1 drop genome; 2 E_AO/E_TT strata + agreement; 3 SelectorView spec + per-channel tests; 4 wording/coverage; 5 identity exclusions; 6 lockbox declaration + 23-gene sensitivity; 7 on-array stratum for HAP1; 8 SMF; 9 provenance table; 10 assertions; 11 seeds | A |
| algorithms | 1 F; 2 shared kernel module; 3 drop U; 4 degree model spec; 5 oracle; 6 surrogate with fidelity check; 7 selector compute rule; 8 B2_product; 9 R1; 10 trajectory variance; 11 batch cap; 12 masking; 13 pool fractions; 14–15 low-n predictor ideas (P2) | A |
| prior work | 1 scaling test; 2 baselines (R1; adaptive column sampling = M; TypiClust ≈ C); 5 n=220; 6 discovery; 7 phased hybrid (F/H + static first); 8 strict label; 9 reframe precision; 11 efficiency currency | A / P (TypiClust folded into C) |
| prior work | 3 like-for-like upstream anchor | P: upstream-protocol reproduction on dev only (cells incl. reciprocal), T comparison post hoc in final |
| prior work | 4, 12 Billmann degree distribution / hub transfer | P: optional, evaluation-only, post-freeze (§7.4) |
| agent failure | 1–5, 7, 8, 10–16, 19, 20 | A |
| agent failure | 18 `.claude` deny rules | R: permission settings are the user's; covered by hsl.guard + check script + agent briefs |
| agent failure | 6 E_fixed vs E_unobs headline rule; 17 lockbox wording; 21 resume block | A |
| HAP1 domain | 1 R1; 2 feasible pool; 3 cost in rounds/screens; 4 scaling test; 5 nonessential scope; 6 feature ablation (HAP1-minimal); 7 strict label (optional); 10 framing; 11 warm start; 12 recommendation shape; 13 E_fixed distributions | A / P |
| HAP1 domain | 8 noise-robustness arm | R: yeast single-direction labels already carry comparable noise (κ 0.14 vs HAP1 precision 0.61); stated as limitation |
| HAP1 domain | 9 HAP1 hub-proxy check | P: optional, evaluation only, post-freeze |


## 15. Response to round-2 reviews
| reviewer | items | decision |
|---|---|---|
| ext. PI | 1 compute + concurrency + triggers; 2 oracle per round; 3 predictor dev before strategies; 4 denominator 3,283; 5 interpolated n_X with absolute context; 6 R1-continuation arm; 7 m1_hgb@490; 8 paired V_Q; 9 noise-aware scaling (→ slope b with bootstrap CI); 10 batch 110; 11 drop V, merge C; 12 trim m1_hgb; 13 deterministic run; 15 marginal returns; 16 optional items last | A |
| ext. PI | 14 E_AO/E_TT distributions | A (P0 manifest, labels/dev rows only) |
| stats/alg | 1 exact weighted AP, B = 2,000; 2 single replicate array + joint bootstrap; 3 Welch formula, Bonferroni bands; 5 linear-contrast "matches" + G* ≥ 0.5; 6 n_req left endpoint, same inversion for n_X; 7 drop E_common; 8 provisional V_Q(Δ), R fixed at freeze; 9 slope b; 10 gate + shortlist spec; 11 re-split restricted to n ≤ 330; 12 within-n surrogate fidelity; 13 F spec; 14 H spec; 15 M redefined; 16 native mask codes, static test; 17 deterministic run | A |
| data/process | 1 Bundle.blind + runtime guard; 2 two-sided identity exclusions; 3 pseudo-test fixture; 4 whitelist by construction, DMF excluded from selectors, provenance YAML; 5 freeze list; 6 protocol.yaml; 7 E_common (→ dropped); 8 re-split naming, dev pool max; 9 sub-universe features; 10 assertions; 11–14 proportionality cuts | A |

## Revision history
- v1 (2026-10-29): first draft.
- v2 (2026-10-29): round-1 revisions (§14).
- v3 (2026-10-29): round-2 revisions (§15).
- v3.1 (2026-10-29): round-3 consistency edits (docs/reviews/design_r3_consistency.md, all accepted). Design frozen.
- **Declared protocol deviations (post-freeze of design, pre-freeze of finalists):**
  - (i) After P2: freeze-hash list extended with config/predictor.yaml and src/hsl/predictors/*.py.
  - (ii) D9: R_s = 20 and R_r = 40, set after P4 dev results existed.
    - The frozen ρ = 0 pilot rule (§7.1) gave no attainable R, because V_Q alone exceeded (g/3.5)².
    - The only finite R the pilot reported, R_s ≈ 35–39, came from a ρ_boot sensitivity. A strict reading of §5
      ('may only go up') would take it as the floor.
    - We deviate because the P4-measured V_Q(Δ) is about 25× the target SE², and the SE is flat in R
      (0.0038 at R_s = 20 vs 0.0037 at 40).
    - A lower R cannot favour strategies over random: it widens CIs and leaves the estimands unchanged.
    - Q2 claims are phrased as bounds (§7.1).
    - Recorded per docs/reviews/p4_review.md §4.
