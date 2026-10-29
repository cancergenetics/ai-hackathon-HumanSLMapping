# Design review r1 — statistics (evaluation sets, inference, power, multiplicity, dev protocol)

Reviewer lens: statistician. Inputs read: PROMPT.md, DESIGN.md v1, results/data_characterisation/summary.md,
config/split/split_manifest.json, docs/recon.md, results/runtime_benchmark. No model run, no test metric computed.
Overall: the architecture (frozen split, E_fixed, paired test-query bootstrap, ≤3 frozen strategies, one-shot final)
is sound. The changes below are about making the inference well defined, adequately powered and hard to game.
Priority: **H** = fix before freeze, **M** = should fix, **L** = nice to have.

## A. Metric definition
1. **(H) Remove "ties broken by stable sort"; use tie-aware AP.** sklearn's `average_precision_score` already treats
   tied scores as one threshold; a stable-sort tie-break makes AP depend on cell order (row-major, i.e. query ORF
   order). This is material: `B1_array_rate` is constant within each array column (656-way ties), and at small n many
   features are 0 (`masked_pearson` → 0, kNN ≈ 0), so large tie blocks are expected. Add a toy test: permuting cell
   order must leave AP bit-identical; AP of an all-tied score vector must equal prevalence.
2. **(H) Fix the per-query (macro) rules now.** Queries with 0 positives in E_fixed (≥ 4 test rows with 0 SLs overall;
   more once arrays are restricted to A\D) are excluded from macro AUPRC and included in pooled; report the count.
   Give macro AUPRC the same inference machinery as pooled (key secondary), because pooled AUPRC rewards between-row
   calibration (query propensity) — which hub-first selection plausibly improves — while macro measures within-row
   ranking. A headline claim must state whether pooled and macro agree in sign at n=160.
3. **(M) E_unobs is not a paired comparison.** E_unobs(S) differs between arms, so ΔAUPRC on it mixes "better
   predictions" with "different, easier/harder cells" (hub-first removes SL-rich columns). Keep E_unobs per arm as
   the practical metric; for a paired secondary use E_common(n) = all (t,a), t ∈ T, a unscreened in every
   trajectory of the compared arms at n (computable from the logged S; E_fixed ⊆ E_common).

## B. Bootstrap and variance (DESIGN §5, §6.3)
4. **(H) Specify the estimand and the variance decomposition explicitly.** Estimand: Δ(n) = E_traj[AUPRC_s(n)] −
   E_traj[AUPRC_r(n)], conditional on D and on the E_fixed array set, for the population of test-like queries.
   Var(Δ̂) = V_Q + σ²_s/R_s + σ²_r/R_r. Compute V_Q by the query bootstrap of the trajectory-averaged Δ (shared
   weights, as planned) and σ²_s, σ²_r as between-trajectory variances on the full E_fixed; CI = Δ̂ ± t_ν·SE with
   Welch–Satterthwaite ν (conservative: the query×trajectory interaction is counted twice). Report V_Q and the
   trajectory terms separately. Replace "bootstrap over trajectories within each replicate": with R = 5 it
   underestimates variance by (R−1)/R and yields ≤ 126 distinct means.
5. **(H) Stratified query bootstrap.** Resample test queries within the 11 split strata (the split was stratified), with
   B = 10,000 for any p-value entering Holm (see 11); 1,000 reps give a resolution of 0.001, too coarse for
   α/6 = 0.0042. Make it cheap: for each score vector, bin scores into ≥ 4,096 quantile bins once, store per-query
   positive/negative counts per bin (656 × 4,096 × 2), and compute each replicate's weighted AP from
   wᵀ·H (tied-bin AP). Verify against exact AP on unweighted data (|error| < 1e-4) in the toy tests.
6. **(M) Array-side dependence: add a two-way (pigeonhole) bootstrap as sensitivity for the headline contrasts.**
   Cells sharing an array are correlated (array SL Gini 0.54; top 10% of arrays hold 41% of SLs) and the query-only
   bootstrap is conditional on the 1,708 E_fixed arrays. For the ≤ 6 confirmatory contrasts: weight cell (t,a) by
   w_t·v_a, where v_a is an independent multinomial weight for array-only genes and v_a = w_a (the same gene's query
   weight) for test genes on the array side, so reciprocal test×test cells (≈ 30% of E_fixed) are resampled as one
   gene. B = 500 exact replicates suffices for a CI. Report whether any conclusion changes.
7. **(M) Model-seed variance for "deterministic" strategies.** Every adaptive strategy whose later rounds use model
   predictions (H after 30, V, VH, U) is stochastic through the cross-fit fold seed; static variants are stochastic
   at evaluation. Run each with ≥ R trajectories (distinct model seeds; seeded tie-breaks), never once.
8. **(M) Common random numbers.** Derive all per-trajectory seeds with `SeedSequence(base).spawn(R)`; trajectory r of
   every arm uses model seed r (paired across arms) — reduces Var(Δ̂) at no cost. Add `q1_draw_seeds`, `q2_seeds`,
   `inner_split_seeds`, `final_*`, `bootstrap_seed` to config/seeds.yaml before first use (none exist yet), with
   dev and final streams disjoint.

## C. Power / number of trajectories (DESIGN §6.3)
9. **(H) Pre-register R from a dev pilot, not a fixed 5/10.** LITERATURE.md's own prior is a 1.3–2× budget saving,
   i.e. "5% ≈ 20%" will likely fail, so the informative endpoint is efficiency n_req/n. Pilot (inner split 1, random
   arm only, R = 10, stack_lr_wr2, n ∈ {120, 160, 210, 330, 490, 660}): estimate σ_r(n), V_Q(n) on E_fixed_dev and the
   random-curve gap g = AUPRC_r(1.3·160) − AUPRC_r(160). Choose R so that SE(Δ̂(160)) ≤ g/(z_{1−α/6} + z_{0.8})
   ≈ g/3.5 (assume σ_s ≈ σ_r), i.e. R ≥ 2σ²_r / ((g/3.5)² − V_Q). If V_Q alone exceeds (g/3.5)², say so: the
   design cannot resolve a 1.3× saving with 656 test queries and the claim must be phrased as a bound.
   Expect R_random ≥ 20 (cheap: nested, ≈ 10 checkpoints × 15–25 s) and R_strategy ≥ 10. Record the calculation.
10. **(M) m1_hgb confirmation on the same trajectories.** Select with stack_lr_wr2 (or stack_lr), then refit m1_hgb only
   at the 4 checkpoints of the already-selected sets (≈ 4 × 2 min per trajectory), so m1_hgb gets all R
   trajectories. "Rankings agree" = same sign of Δ(160) with both predictors and overlapping ranking; define it.

## D. Confirmatory tests, margin and multiplicity
11. **(H) Define the family and the tests.** Family F = {(c) non-inferiority, (a) superiority at n=160} × ≤ 3 frozen
   strategies = ≤ 6 one-sided hypotheses, Holm at familywise α = 0.025 (one-sided). p-values: one-sided Wald from
   item 4 (bootstrap percentile as a check). Report Holm-adjusted p and (1 − 0.05/6) Bonferroni intervals alongside
   the unadjusted 95% CIs. Everything else (other n, strata, E_unobs, iterative-vs-static, Q1) is descriptive;
   cross-budget curves get a simultaneous max-|t| band (bootstrap of the maximum studentised deviation over the 8
   checkpoints), since checkpoints are strongly correlated (nested sets).
12. **(H) The 5% relative AUPRC margin can be passed on a flat curve — anchor it to the random curve.** If
   AUPRC_r saturates, 0.95·AUPRC_r(660) may already be reached by random at, say, 330, so "matches 20% random" would
   be declared for a strategy that is merely as good as 10% random. Primary (c): non-inferiority on the budget scale,
   H0: AUPRC_s(160) ≤ AUPRC_r(n_m) with n_m = 490 (= 0.75·660, frozen now), tested with the same machinery. Keep the
   5% relative rule as a secondary, and report the random-curve budget at which 0.95·AUPRC_r(660) is reached so the
   reader sees how lax it is. Pre-specify: if on inner splits AUPRC_r(660) − AUPRC_r(330) < 2 × 0.05·AUPRC_r(660),
   the relative margin is declared non-discriminating.
13. **(M) n_req inversion.** (i) Fit the random curve per bootstrap replicate as an isotonic (PAVA) fit of the
   trajectory-mean AUPRC over the dense random grid, interpolate in log n, then invert; (ii) add random checkpoints
   at 120, 210, 410, 820 to cut interpolation error in the likely range; (iii) if AUPRC_s(n) > max random (n = 2627
   from Q1 nested draws) or < AUPRC_r(min), report n_req as right/left-censored and give the CI as a one-sided
   bound — never extrapolate the saturating fit for a confirmatory statement; (iv) CI = percentile over replicates of
   the joint (queries + trajectories) bootstrap; report the fraction of censored replicates. Share one set of random
   permutations between Q1 and Q2 (evaluate at the union of budgets) rather than two random arms.

## E. Dev protocol and selection on inner splits (DESIGN §4.1, §7)
14. **(H) E_fixed_dev uses test-gene array columns.** A \ D_cand contains T∩A, so dev evaluation scores cells (v, t),
   t ∈ T — test-gene data used for tuning, and the "poison all test-gene data → identical" check would fail for dev
   numbers. Use E_fixed_dev = D_val × (A \ (D_cand ∪ T)) (array-only + D_val arrays) and extend the poisoning test to
   the dev evaluator.
15. **(H) Guard against forking paths in auto-research.** The 3 inner splits overlap ~80% in D_cand and D_val is
   random each time, so "CI excludes 0 on ≥ 2/3 splits" is weak evidence after many attempts. (i) Make the 3 D_val sets
   disjoint (3 of 5 stratified folds of D) and keep the remaining 2 folds as an inner confirmation split, evaluated
   once per accepted predictor change and once for the final strategy shortlist; (ii) keep an append-only ledger of
   every strategy/predictor/hyper-parameter configuration evaluated, with its inner results, and report the count;
   (iii) shortlist with R ≥ 3 trajectories per configuration (not 1). Final-test inference stays valid regardless
   (fresh test set, family fixed before), but this limits wasted and over-optimistic choices.
16. **(M) Same predictor recipe for every arm.** If a covariate-shift fix (propensity weighting, degree offsets) is
   adopted, the random arm uses the identical recipe (weights ≡ uniform), and the strategy is also reported with the
   unmodified predictor, so a strategy gain is not a predictor-change gain.

## F. Q3 quantities that could be gamed or biased
17. **(M) Precision thresholds must not be chosen on T.** "Top-k at precision ≥ 0.2/0.5 estimated on T" picks the
   threshold on test labels. Choose the score threshold by leave-batch-out CV on screened rows (what an experimenter
   can do) or on inner dev splits; then report realised precision/recall on T. When extrapolating to D \ S under a
   hub strategy, reweight T queries to the degree distribution of D \ S (evaluation-only), since the remaining pool
   is SL-poorer than T.
18. **(M) Online stopping estimator: don't use the across-n correlation.** Both the estimate and the test AUPRC rise
   with n, so their correlation is high by construction. Report (i) bias and MAE of the estimate at each n, (ii) the
   within-n correlation across trajectories, (iii) stopping-n vs a pre-registered true knee (smallest n with test
   AUPRC ≥ 0.95 × the curve's value at 660, per trajectory), with the plateau rule frozen on inner splits.
19. **(L) Q1 curve fit.** Fit a − b·n^−c per seed (nested draws ⇒ within-seed points are correlated) and summarise
   parameters across seeds with a seed-cluster bootstrap; report per-stratum (on/off-array) curves descriptively only.
