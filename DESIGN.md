# DESIGN — Using the yeast GI network to design human HAP1 synthetic-lethal screens

Version: **v2** (2026-10-29), after review round 1 (7 reviewers, `docs/reviews/design_r1_*.md`; responses in §14).
Brief: PROMPT.md (fixed). Background: LITERATURE.md. Upstream predictor: YEAST_SLs @ 9b5dce9 (read-only; vendored in
`third_party/ysl`, docs/recon.md). Data numbers: results/data_characterisation/summary.md.

## 0. Questions and pre-registered endpoints

| Q | question | primary endpoint (primary predictor P*, §4) |
|---|---|---|
| Q1 | How does SL prediction for unscreened queries scale with the fraction screened (random)? | Pooled AUPRC on E_fixed (§3) vs n; **n_80, n_90** = smallest budget whose mean AUPRC reaches 80/90% of the measured AUPRC at the full dev pool (n = 2,627) |
| Q2 | Can a smaller strategically chosen set (iterative, batches of 10) match a larger random set? | (c1) **strategy@160 > random@490** (non-inferiority on the budget scale, one-sided); (a) Δ(160) = strategy@160 − random@160 > 0; reported with gap closure **G** and random-equivalent budget n_req(n) |
| Q3 | Implications for HAP1: number/choice of query lines, contribution of prediction | n_80/n_90 under the scaling law chosen by the universe-subsampling experiment (§7.3); warm-start gain over Billmann-like continuation; SL-discovery frontier; choice-of-query rules on human-available features |

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
  indicator (yeast artefact; every HAP1 gene is in the library).
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
**Screening** gene g ∈ D reveals its whole query row ε(g, ·), P(g, ·).
**Observed pairs / symmetric handling.** A gene pair {x, y} counts as observed iff x ∈ S or y ∈ S. If x is screened
and y is not, only ε(x, y) (the screened direction) exists; it is available to the predictor as a feature of cell
(y, x) (the "reciprocal" feature), but cell (y, x) is excluded from every evaluation set and is never a training label;
ε(y, x) is never revealed. The evaluation unit is the **directed cell** (t, a), labelled by its own measurement; AB
and BA are not merged (weak agreement: P(BA SL | AB SL) = 0.145, κ = 0.14; merging would change the label). When t
and a are both test genes, (t, a) and (a, t) are two cells. Self pairs and cells with missing ε/P are excluded.
Gene identity: arrays whose current SGD name equals a dev-pool gene's current name (4 array-only columns:
YFL056C, YHR079C-B, YIR043C, YJL206C-A; list generated by script) are dropped from E_fixed so that E_fixed stays
invariant.
**E_fixed (primary).** Measured cells (t, a), t ∈ T, a ∈ A \ D (minus the identity exclusions): never observable
under any strategy, so identical for all arms. Pre-registered column strata: **E_AO** = T × array-only (≈1,170
columns, SL rate ≈0.008) and **E_TT** = T × on-array test genes (533 columns; columns are query genes — the closest
analogue of HAP1, where every library gene is a potential query). Counts/prevalence are written to the manifest
(labels only) before any scoring. Note: no E_fixed array is ever screened, so the reciprocal feature is never
available there — E_fixed measures fold-in/propensity prediction.
**E_unobs(S)** (practical, per arm; not paired): (t, a), a ∉ S; reported with prevalence and lift.
**E_common(n)** (paired secondary): cells unscreened in every trajectory of the compared arms at n (⊇ E_fixed).
**Headline rule:** a Q2 claim is headline only if Δ has the same sign on E_fixed, E_TT and E_common, with P* and with
m1_hgb; otherwise it is reported as qualified.
**Metric.** Pooled average precision with **tied scores grouped at one threshold** (sklearn
`average_precision_score` semantics; no stable-sort tie-breaking). Secondary: macro (per-query) AUPRC over test queries
with ≥ 1 positive in the set (zero-positive queries counted and reported; included in pooled); AUROC; precision@k
(k = 10, 50) per query; lift = AUPRC/prevalence. Per-stratum: on-/off-array test genes; SL-count tertile.

## 4. Predictors
- **Primary predictor P\* = `stack_lr_wr2`** (upstream multi-kernel fold-in + LR; ~12–60 s per fit), or its frozen
  successor if one passes the §8 gate before the freeze. All Q1/Q2/Q3 headlines use P*; n_req always inverts the random
  curve of the same predictor.
- Secondary: `m1_hgb` (genome dropped; `prop_folds = 5` when n < 30), refit on the same selected sets at checkpoints
  {60, 160, 330, 660}; baselines `B1_array_rate`, `B2_product` (array rate × query propensity) and `stack_lr`;
  tier-2 `m1_hgb_t2` sensitivity at 660 only.
- **Same recipe in every arm.** Any predictor change (e.g. covariate-shift correction) applies to all arms; the strategy
  is also reported with the unmodified predictor.
- Memory: n ≥ 1,970 runs alone with ≤ 12 GB.

## 5. Budgets and seeds (single source: `config/budgets.yaml`, `config/seeds.yaml`)
- Random/Q1 grid: n ∈ {10, 30, 60, 100, 120, 160, 210, 220, 330, 410, 490, 660, 820, 980, 1310, 1970, 2627}
  (220 ≈ Billmann's 222 lines).
- Q2 checkpoints: {30, 60, 100, 160, 220, 330, 490, 660}; trajectories run to 660 (dev: to 660 or pool exhaustion).
- Random arm = nested prefixes of seeded permutations of the candidate pool (shared by Q1 and Q2).
- Seeds via `numpy.random.SeedSequence(base).spawn(R)`; trajectory r of every arm uses model seed r (common random
  numbers). Separate RNG streams for selection and evaluation. Dev and final streams disjoint. All named seed lists
  are appended to `config/seeds.yaml` before first use.
- R (trajectories per arm) fixed from the power pilot (§7.1): expected R_random ≥ 20, R_strategy ≥ 10.

## 6. Q2 — strategies
**Protocol.** S_0 = ∅ (cold start) unless stated. Each round picks 10 genes from the candidate pool using only the
SelectorView, then their rows are revealed.
### 6.1 SelectorView (harness-built; selector code never reads split files)
(a) screened rows only; (b) in those rows, columns of T (and, in dev, of D_val) are NaN-masked in every matrix
(E, P, labels, DMF) — masked, not deleted, so upstream code runs unchanged; (c) SMF of T (and D_val in dev) NaN; (d)
candidate list = pool \ S; (e) feature whitelist (§1, tier 1, array SMF; no query SMF, no missingness/off-array
indicators, no n_arrays/n_queries, no Costanzo strain metadata). Any model inside a selector is fit on this object only.
**Poisoning tests, one per channel:** poison (a) test rows, (b) test columns of screened rows, (c) test SMF, (d)
unscreened pool rows — each alone leaves selections bit-identical; positive control: poisoning a visible column of a
screened row changes some selection; each test has a broken-twin selector it must catch.
### 6.2 Strategy set (≤ 3 frozen for the final, chosen on dev)
- **R0 random** (comparator).
- **R1 Billmann-like** (status quo, static): stratified sampling over SMF × abundance tertile bins, spread across
  external-kernel clusters (k-means on the external kernel embedding), round-robin over bins.
- **H hub-first:** rank candidates by predicted **query-side** SL degree. Degree model: before n = 30 a fixed external
  prior (−SMF rank + abundance rank, pre-stated, unfitted); from n = 30, ridge-logistic model of shrunk query-side SL
  rate fitted on screened rows (inputs: shrunk array-side SL rate over S, array SMF, abundance, AP-MS degree,
  coexpression degree, conservation, paralog flag), refit each round. Within-batch redundancy cap: ≤ 2 genes per
  external-kernel cluster. Diagnostic (dev): Spearman of predicted vs true degree on the unscreened pool per round.
- **F Fisher-D (rank-r logistic design; main hybrid):** maximise Σ_x w_x log det(δI_r + Σ_{s∈S} p̂_sx(1−p̂_sx) u_s u_sᵀ)
  over target pool x (unscreened non-test candidates + visible array-only columns); u_g from truncated SVD of the
  visible screened block (array-side factor; off-array/cold start: kernel PCA of the external kernel); p̂ from the H
  degree model × column rates. Greedy with determinant lemma + Sherman–Morrison; r = 1 reduces to H (built-in
  ablation); r ∈ {4, 8, 16}, w_x ∈ {1, predicted degree} tuned on dev.
- **V GP variance reduction** and **C facility-location coverage** on the predictor's own `mk` kernel (GI sym kernel
  over screened on-array genes + external kernels; one cached kernel module), and **M max-variance / adaptive column
  sampling** (largest residual of the candidate's array-side profile after projection on screened profiles). With
  external kernels only, V/C are label-free (static = adaptive); their adaptivity comes only from the GI part.
- **Static variants** of H, F, C (external features + fixed prior only, chosen at round 0) — value of iteration.
- **Dropped:** U (bootstrap-disagreement uncertainty): objective mismatch at low budget and ~4–5 h per trajectory.
- Deterministic strategies get trajectory variance via seeded randomised greedy (sample each pick from the top-20 by
  gain) and model seeds; each arm has R trajectories.
- Compute rule: no selector may fit stack_lr_wr2 or m1_hgb per round; allowed per round: the kernel/degree modules and
  at most one stack_lr fit. Selection seconds per round are logged.
- Hub cut-offs and switch points are expressed as fractions of the candidate pool (dev pool ~2,100 vs final 2,627).
### 6.3 Scenarios and sensitivities (best ≤ 2 strategies)
- **Warm start:** S_0 = 220 genes chosen by R1 (≈ Billmann's 222), then random vs strategy for the next 440 (to 660).
- **Batch size:** 10 (brief), 50, and one-shot (= static) at n ∈ {160, 330}; report the number of sequential rounds.
- **HAP1-feasible pool:** candidate pool minus the bottom array-SMF decile (genes whose HAP1 KO lines are often
  unavailable/sick); report AUPRC lost.
- **HAP1-minimal features:** selector restricted to fitness + expression/coexpression + protein embeddings (no AP-MS).
### 6.4 Statistics (confirmatory family fixed now)
- Estimand Δ(n) = E_traj[AUPRC_s(n)] − E_traj[AUPRC_r(n)] on E_fixed with P*. Var = V_Q + σ²_s/R_s + σ²_r/R_r: V_Q
  from a **stratified query bootstrap** (resample test queries within the 11 split strata, B = 10,000, same weights
  for all arms, fast binned-AP: ≥ 4,096 score-quantile bins; verified against exact AP, |error| < 1e-4), trajectory
  terms from between-trajectory variance; Welch–Satterthwaite df. Report each component.
- **Family** (one-sided, Holm, familywise α = 0.025): for each of ≤ 3 frozen strategies, (c1) AUPRC_s(160) >
  AUPRC_r(490) and (a) Δ(160) > 0 → ≤ 6 hypotheses. Wald p-values (bootstrap percentile as a check). Everything else is
  descriptive; curves get a simultaneous max-|t| band over checkpoints.
- **Brief's literal example** (strategy@160 vs random@660): reported as estimate + CI and gap closure
  G = (AUPRC_s(160) − AUPRC_r(160)) / (AUPRC_r(660) − AUPRC_r(160)) with CI. "Matches 20% random" may be stated only if
  the CI of AUPRC_s(160) − AUPRC_r(660) has lower bound > −0.05·AUPRC_r(660) **and** the random gap
  AUPRC_r(660) − AUPRC_r(160) exceeds that margin (otherwise the comparison is declared non-discriminating). The budget
  at which random reaches 0.95·AUPRC_r(660) is reported.
- **n_req(n):** per bootstrap replicate, isotonic (PAVA) fit of the trajectory-mean random curve over the dense random
  grid, log-n interpolation, inversion; right/left-censored if outside [AUPRC_r(10), AUPRC_r(2627)] (fraction of
  censored replicates reported; never extrapolated). Efficiency = n_req/n.
- Sensitivity: two-way (query × array) bootstrap (B = 500) for the ≤ 6 contrasts; test genes on the array side share
  their query weight.
- m1_hgb confirmation: same trajectories; "agree" = same sign of Δ(160).

## 7. Q1 / Q3 analyses
### 7.1 Q1 and power pilot
- Random curves for P*, B1, B2_product, stack_lr at all grid budgets × R_random permutations; m1_hgb at
  {30, 60, 160, 330, 660, 1310, 2627} × 5. Report AUPRC (mean, SD, CI), strata, runtime per fit. Saturating fit
  AUPRC(n) = a − b·n^−c per permutation, summarised by seed-cluster bootstrap (descriptive; no extrapolated headline).
- **Power pilot** (dev split 1, random only, R = 10, P*): estimate σ_r(n), V_Q(n), gap g = AUPRC_r(210) − AUPRC_r(160);
  choose R so SE(Δ̂(160)) ≤ g/3.5; if V_Q alone exceeds (g/3.5)² say so and phrase Q2 claims as bounds.
### 7.2 Discovery axis (co-reported with every Q2 arm)
SLs detected = strong-SL cells in screened rows (cell-level primary); pair level: "detected" (SL call in a screened
row) vs "touched" (pair involves a screened gene). Plot AUPRC (E_fixed) vs SLs detected as a frontier.
### 7.3 Scaling law (fraction vs count) — decides the HAP1 translation
Dev only: draw gene sub-universes of 50% and 25% (same genes removed from query and array sides; 3 draws each; inner
split structure rebuilt inside each sub-universe), rerun the random curve with P*. Pre-registered criterion: compute
n_90 per universe; if the coefficient of variation of n_90/N (fraction) across universe sizes is smaller than that of
n_90 (count), adopt fraction scaling, else count scaling; if neither is < 0.25, report both scenarios. Caveat: sub-
sampling mimics fewer genes per module, not more modules, so it bounds rather than proves the HAP1 mapping.
### 7.4 HAP1 translation and actionability
- n_80/n_90 (and strategy equivalents) converted to HAP1 line counts under the chosen law; HAP1 query universe =
  viable KO lines among expressed nonessential genes (~7,500–9,500; LITERATURE), shown with its range. Costs in lines,
  screens (×1.34 replicate factor, Billmann) and sequential rounds. Billmann's 222 lines placed on the curves.
- On-array test genes (E_TT/on-array strata) are the estimate used for HAP1 translation (every HAP1 query is also in
  the library); all-T remains primary for Q1/Q2.
- **Choice of query lines:** profile genes picked by the best frozen strategy vs random on human-available features
  (SMF, abundance, coexpression/AP-MS degree, conservation, paralogs); interpretable rank rule.
- **Prediction contribution:** precision/recall of P* on E_TT against the within-yeast replicate reference (reciprocal
  measurement BA → AB on the same test × test cells); score thresholds chosen by leave-batch-out CV on screened rows (not
  on T); expected true SLs recovered among top-k predictions per unscreened gene.
- **Online stopping estimator:** leave-batch-out CV AUPRC on screened rows only; report bias/MAE vs true test AUPRC per
  n, within-n correlation across trajectories, and stopping-n of a plateau rule (frozen on dev) vs the true knee
  (smallest n with AUPRC ≥ 0.95 × AUPRC(660)).
- Optional, evaluation only, after freeze (never feeds back): degree-distribution shape of 222 random yeast queries vs
  Billmann 2026 Data S14; Spearman between our human-feature rank rule and HAP1 library-side negative degree. HAP1 GI
  data are never a model/selector input.
- Optional secondary: stricter label ε < −0.35 (∧ P < 0.05) for re-scoring the final random curve (analogue of HAP1
  "extreme SL"); the brief's label stays primary.

## 8. Development protocol (no test information)
- Inner folds: 5 stratified folds of D (by SL count). **Development splits k = 1, 2, 3:** D_val = fold k,
  candidates = D \ fold k. **Inner confirmation split:** folds 4+5 as validation (candidates = folds 1–3), evaluated
  once per accepted predictor change and once for the final strategy shortlist.
- E_fixed_dev = D_val × (A \ D_cand \ T) — array-only + D_val arrays; **no test-gene columns**.
- Predictor-change gate (§4.1 of v1): parent vs child paired query-bootstrap CI excludes 0 on ≥ 2/3 dev splits, mean
  gain ≥ 2%, sign holds on the confirmation split; same recipe in all arms. Development focus: n ≤ 330 (upstream was
  tuned at 20%): n-dependent GI kernel weight, propensity features, small-n degeneracies.
- Strategy shortlist: ≥ 3 trajectories per configuration on each dev split; tuning may use the closed-form fold-in
  surrogate only if its Spearman with P* AUPRC across random trajectories × checkpoints is ≥ 0.8.
- **Attempt ledger** (`results/LEDGER.jsonl`, append-only, harness-written): every configuration evaluated, with inner
  results; the count is reported.
- Freeze: `config/finalists.yaml` (≤ 3 strategies + R0 + R1, predictor recipes, budgets, seeds, endpoints) committed;
  hashes in `config/frozen.yaml`.

## 9. Phases and gates (recorded in RESEARCH_TREE.md)
- **P0 infrastructure:** frozen eval module (`hsl/eval.py`: tie-aware AP, E-sets, bootstrap, n_req, tests), dev_bundle,
  SelectorView, budgets/seeds config, run harness (manifest + ledger), launcher (≤ 2 concurrent jobs, timeouts),
  toy/broken-twin tests, read-only check script, feature provenance table.
- **P1:** dev random curves (P* + baselines; m1_hgb subset) + power pilot + universe subsampling.
- **P2 ceilings (go/no-go):** oracle hub-first (true query degree) and greedy oracle (closed-form fold-in, best of 200
  random candidates per pick by realised E_fixed_dev AUPRC) on dev splits — evaluation-only diagnostics, never
  strategies. If the oracles reach G < 0.5 at n = 160, de-scope V/C/M and say so.
- **P3:** static strategies (R1, H-static, F-static, C-static). **P4:** adaptive H, F, V/C/M; warm start; batch size.
- **P5:** predictor development for n ≤ 330 (if a change is accepted, P3–P4 shortlist results are re-run with it).
- **Freeze → Final** (one shot, lock file, clean tree) → Q3 analyses → REPORT.

## 10. Final evaluation
`scripts/run_final.py` refuses unless `config/finalists.yaml` is committed, the tree is clean, frozen hashes match, and
`results/final/.lock` is absent; writes the lock (git hash); predictions are saved before scoring; crash re-runs are
recorded as attempts with reasons. Candidate pool = D; T evaluated. Q1 on T runs after predictor recipes are frozen.

## 11. Verification
- Toy tests with known answers: AP literals; all-tied → prevalence; cell-order permutation → identical AP; sklearn
  agreement on 1e5 cells; stable-sort mutant caught; planted-hub toy (H beats random), planted-cluster toy (C beats
  random), null toy (no strategy beats random beyond noise); n_req inversion on a known curve; bootstrap coverage on a
  toy.
- Leakage: per-channel poisoning (§6.1) with broken twins; prediction poisoning (test rows → identical; screened rows →
  change); dev scripts cannot call the unblinded loader (grep test).
- End-to-end controls at each milestone: shuffled screened labels → AUPRC ≈ prevalence; oracle back door → ≈ 1.
- Sanity on dev: random-curve monotone on average; P* at n = 660 in a plausible range of upstream dev numbers
  (dev only; T-based upstream comparison only in the final report, labelled post hoc).
- Independent rescoring from saved predictions (no project imports) by an agent that did not produce them; all report
  numbers come from generated tables / `results/claims.json`, never hand-typed.
- YEAST_SLs read-only: `hsl.guard` + `scripts/check_readonly.sh` (stamp + `GIT_OPTIONAL_LOCKS=0 git status`) at each
  milestone.

## 12. Compute
M2, 8 cores, 24 GB; jobs ≤ 3 threads, ≤ 2 concurrent, ≤ 6 GB (≤ 12 GB alone for n ≥ 1,970). Estimate (algorithms
review): ≈ 30 job-h ≈ 15 h wall for dev + final; cuts in order: oracle to 1 split; m1_hgb Q1 grid; dev starts for
eliminated strategies. Every result row records fit/predict and selection wall-clock.

## 13. Threats and limitations
Single-direction labels are noisy (κ 0.14), capping AUPRC; NxN lacks essential genes; yeast→HAP1 scaling is an
assumption tested only by subsampling; upstream recipes tuned at 20% random; selection-induced covariate shift; E_fixed
is dominated by array-only columns (strata reported); T is not blind for upstream predictors (declared).

## 14. Response to round-1 reviews (A = accepted, P = partly, R = rejected)
| reviewer | item | decision |
|---|---|---|
| ext. PI | 1 n_80/n_90; 2 scaling experiment; 3 gap closure G; 4 oracle ceilings first; 5 R1 + warm start; 6 one primary predictor; 7 discovery axis; 8 reciprocal reference; 9 E_fixed characterisation + agreement rule; 10 batch size; 11 drop U, add M; 12 phase gates; 13 trajectory variance | A (all) |
| ext. PI | 14 availability subset; 15 strict label; 16 per-query rules | P: 14 replaced by HAP1-feasible pool (domain #2); 15 optional secondary; 16 A |
| statistics | 1 tie-aware AP; 2 macro rules; 4 variance decomposition; 5 stratified B=10k binned; 8 CRN; 9 power pilot; 11 family/Holm; 12 budget-scale margin (160 vs 490); 13 isotonic censored n_req; 14 no T columns in dev eval; 15 disjoint folds + confirmation + ledger; 16 same recipe; 17 thresholds not on T; 18 stopping metrics | A |
| statistics | 3 E_common; 6 two-way bootstrap; 7 model-seed variance; 10 m1_hgb on same trajectories; 19 per-seed fit | A |
| data | 1 drop genome; 2 E_AO/E_TT strata + agreement; 3 SelectorView spec + per-channel tests; 4 wording/coverage; 5 identity exclusions; 6 lockbox declaration + 23-gene sensitivity; 7 on-array stratum for HAP1; 8 SMF; 9 provenance table; 10 assertions; 11 seeds | A |
| algorithms | 1 F; 2 shared kernel module; 3 drop U; 4 degree model spec; 5 oracle; 6 surrogate with fidelity check; 7 selector compute rule; 8 B2_product; 9 R1; 10 trajectory variance; 11 batch cap; 12 masking; 13 pool fractions; 14–15 low-n predictor ideas (P5) | A |
| prior work | 1 scaling test; 2 baselines (R1; adaptive column sampling = M; TypiClust ≈ C); 5 n=220; 6 discovery; 7 phased hybrid (F/H + static first); 8 strict label; 9 reframe precision; 11 efficiency currency | A / P (TypiClust folded into C) |
| prior work | 3 like-for-like upstream anchor | P: upstream-protocol reproduction on dev only (cells incl. reciprocal), T comparison post hoc in final |
| prior work | 4, 12 Billmann degree distribution / hub transfer | P: optional, evaluation-only, post-freeze (§7.4) |
| agent failure | 1–5, 7, 8, 10–16, 19, 20 | A |
| agent failure | 18 `.claude` deny rules | R: permission settings are the user's; covered by hsl.guard + check script + agent briefs |
| agent failure | 6 E_fixed vs E_unobs headline rule; 17 lockbox wording; 21 resume block | A |
| HAP1 domain | 1 R1; 2 feasible pool; 3 cost in rounds/screens; 4 scaling test; 5 nonessential scope; 6 feature ablation (HAP1-minimal); 7 strict label (optional); 10 framing; 11 warm start; 12 recommendation shape; 13 E_fixed distributions | A / P |
| HAP1 domain | 8 noise-robustness arm | R: yeast single-direction labels already carry comparable noise (κ 0.14 vs HAP1 precision 0.61); stated as limitation |
| HAP1 domain | 9 HAP1 hub-proxy check | P: optional, evaluation only, post-freeze |

## Revision history
- v1 (2026-10-29): first draft.
- v2 (2026-10-29): round-1 revisions (§14).
