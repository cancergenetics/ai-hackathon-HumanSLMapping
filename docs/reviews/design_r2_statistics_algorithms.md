# Design review r2: statistics + algorithms (DESIGN v2 §3, §6, §7.1, §7.3, §8)

Lens: are v2's statistics/strategy specs correct and implementable without guessing? No models run. **H** = fix before
P0 freeze, **M** = should, **L** = optional. Sizes: E_fixed ≈ 1.0M cells, ~10k positives; AUPRC ≈ 0.06–0.10.

## A. Metric and bootstrap engine
1. **(H) Replace quantile-binned AP with exact weighted AP, B = 2,000.** 4,096 quantile bins hold ~250 cells each, so
   the top bin (precision ~0.3–0.5) is scored as one tie block: a downward error of order 1e-3, above the 1e-4
   tolerance, similar to plausible Δ(160), and larger for arms with sharper top rankings. Quantile edges can also split
   a tie group, bringing back r1-1's order-dependence. The confirmatory p-values are Wald, so only a variance is needed
   and B = 2,000 (SE to ~1.6%) suffices. Exact method: sort each vector once and store query index, label and
   tie-group ends; per replicate, gather w[q] → two cumsums → AP = Σ_g (ΔTP_g/TP_tot)·TP_g/(TP_g+FP_g) at group ends,
   ~10–20 ms per 1M cells. Toy test: integer-weighted AP equals sklearn AP on the duplicated cells. Fallback: exact
   groups for the top 2% of cells, tie-respecting bins below, validated on real dev vectors.
2. **(H) One replicate array A[b, arm, r, n]** (AP of trajectory r at budget n under stratified query weights w^b,
   same w^b for all arms) drives everything: V_Q = Var_b of the trajectory-mean contrast, the Welch SE, and a **joint
   bootstrap** for G, "matches", n_req and the bands (each replicate resamples whole trajectories within arm by
   indexing A). v2's n_req replicate is query-only as written and misses trajectory variance. G's CI is the joint
   percentile interval; G is undefined if the 97.5% interval of r660 − r160 includes 0. Cost ≈ (20×17 + 30×8)
   vectors × 2,000 × 15 ms ≈ 5 CPU-h, run in parallel.
3. **(H) Write the Welch formula into §6.4.** SE² = V_Q + a + b, a = σ̂²_s/R_s, b = σ̂²_r/R_r;
   ν = SE⁴/[a²/(R_s−1) + b²/(R_r−1)], V_Q treated as known (ν_Q = ∞); c1 uses σ̂²_r(490); p = 1 − T_ν(Δ̂/SE). Ignoring
   the CRN covariance is conservative; say so. Wald decides; a disagreeing percentile check is reported only.
   **Simpler band:** Bonferroni over the 8 checkpoints instead of max-|t|, which needs its own studentised bootstrap
   for ~10% less width.
4. **No change needed:** Holm over ≤ 6 one-sided hypotheses at 0.025 (valid under any c1/a dependence), the 11-stratum
   bootstrap, CRN seeding, the macro/zero-positive rules, and the B = 500 two-way sensitivity (which needs item 1's
   exact AP, because its weights are per cell).

## B. Q2 effect quantities
5. **(H) "Matches 20% random" is nearly vacuous at its boundary.** With m = 0.05·r660 and "gap > m": if the gap
   equals m, then 0.95·r660 = r160, so "matching" needs only s160 ≥ r160 (G ≥ 0). In general a match implies only
   G ≥ G* = 1 − m/gap. Fix: test H0: s160 − 0.95·r660 ≤ 0 as a linear contrast (no estimated margin treated as a
   constant); call the comparison discriminating only if the point estimate of G* is ≥ 0.5 (gap ≥ 2m, as r1-12 had);
   report G*. Level: one-sided 0.025, unadjusted, descriptive; it may reach the headline only alongside a rejected c1.
6. **(M) n_req details.** Define n_req as the left endpoint of the PAVA flat segment reached. Code censored replicates
   as 0 or +∞, and report a bound that falls in the censored mass as one-sided. In dev, use the dev-pool maximum in
   place of 2,627. **Use the same inversion for Q1's n_80/n_90:** v2's "smallest grid budget" is coarse (e.g. 980 vs
   1,310) and feeds the HAP1 numbers and item 9.
7. **(M) Drop E_common from the headline rule.** The union of screened sets grows with R: with 30 trajectories, ~85%
   of dev genes are screened somewhere at n = 160 and ~98% at n ≥ 490. E_common thus collapses to E_fixed and depends
   on R. Use E_fixed/E_AO/E_TT for sign agreement; keep E_unobs per arm, descriptive.

## C. Power pilot and scaling criterion
8. **(H) A random-only pilot cannot estimate the V_Q that matters,** which is that of the paired strategy − random
   difference (random-vs-random → 0 as R grows). Provisional: V_Q(Δ) := 2·V_Q,single(160)·525/656 (dev fold → |T|;
   ρ = 0, conservative), with R_r = 2·R_s because random is shared by 3 strategies. Then re-estimate V_Q(Δ), σ_s and
   σ_r from the P3–P4 dev strategy-vs-random runs and fix R at the freeze (R may only go up). Record both.
9. **(M) §7.3: replace the CV comparison with a slope.** For sizes 1 : ½ : ¼, CV(n_90) and CV(n_90/N) mirror each
   other; comparing them only asks whether the slope b of log n_90 on log N is above 0.5. Fit b by OLS over the 9
   draws (draw-bootstrap CI) and translate n_HAP1 = n_yeast·(N_HAP1/N_yeast)^b, with b = 0 and b = 1 as bounds: this
   is continuous and carries the uncertainty into Q3. Repeat for n_80, which is more stable (n_90 sits on the flat
   part). N = candidate-pool size; ≥ 10 permutations per draw. Report the 25% E_fixed_dev positive count, since a
   ~130-query D_val may be too noisy.

## D. Dev protocol and gate (§8)
10. **(H) Specify the gate and the shortlist rule.** Gate at n ∈ {30, 160, 330} on the shared R0 dev trajectories:
    mean relative gain ≥ 2% across those budgets, and the 95% CI at n = 160 excludes 0. Predictor changes keep S fixed,
    so the trajectories are paired: use Var(d_r)/R + V_Q(d̄), not the Welch SE. Shortlist: rank by the mean over the
    3 splits of relative Δ(160), at most one variant per family (H, F, V/C/M). Use R_r,dev = 10 per split, reused for
    every configuration. Static counterparts of the finalists run on T as descriptive arms outside the family.
11. **(M) The confirmation split differs structurally.** Its pool (folds 1–3) is ≈ 1,575, so n = 490/660 are 31%/42%
    of it, and folds 4–5 rows trained the models on splits 1–3. Restrict checks to n ≤ 330 (or to equal pool fractions)
    and state the overlap.
12. **(H) The surrogate fidelity check is inflated by the budget trend.** A Spearman across trajectories × checkpoints
    rises with n for both the surrogate and P* (the flaw r1-18 flagged for stopping). Require within-n Spearman across
    trajectories at n ∈ {30, 160, 330}; once ≥ 5 configurations exist, also Spearman ≥ 0.7 between surrogate Δ and
    P* Δ at n = 160.

## E. Strategy specifications
13. **(H) F is not implementable as written.** (i) Use u_g = (1, ũ_g) with r−1 SVD factors in ũ. "r = 1 → H" holds
    only for u ≡ 1, where the gain Σ_x w_x log(1 + p̂_gx/M_x) is monotone in q̂_g (exactly H's ordering); the
    constant carries the D-optimal propensity part. (ii) Use one factor space: screened rows use their real row
    factor; on-array candidates use û_g = B·v_g (column factor v_g, r×r map B fitted by least squares on screened
    on-array genes; reciprocity is weak, ε r = 0.13); off-array candidates and cold start use KRR from the external
    kernel into the same space. v2's "kernel PCA" is a different basis and cannot share the log det. (iii) Whiten ũ to
    unit mean square; δ = mean p̂(1−p̂), not tuned. (iv) p̂_sx = clip(q̂_s·ĉ_x/π̄_S, 1e-4, 0.5). (v) Target pool
    P = on-array candidates ∪ array-only columns, since off-array candidates have no column and carry no information;
    w_x = 1 ("predicted degree" is undefined for array-only columns). (vi) Recompute M_x each round, because the SVD
    basis rotates; Sherman–Morrison only within a batch.
14. **(M) H degree model.** Response: Gaussian ridge on logit((k+0.5)/(m+1)); with ~3,400 cells per row, empirical-
    Bayes shrinkage of the response is negligible. Shrink the array-side-rate input instead:
    ĉ_x = (k_x + 100·π̄_S)/(n_x + 100). λ by GCV. Off-array candidates use a pattern submodel (external inputs only),
    not impute-plus-flag (§1). Prior score = rank(1 − SMF) + rank(abundance), descending; off-array genes use abundance
    only. Clusters for the cap and R1: k-means (k = 100, fixed seed) on **non-test** genes only.
15. **(M) M is degenerate:** |S_on| profiles span R^{|S|}, so the residual → 0. Use k(g,g) − k_gS(K_SS+λI)⁻¹k_Sg (GI only).
16. **(M) SelectorView: mask with the native "unobserved" code** (Y int8 = −1, cannot hold NaN; E/SP/DMF = NaN).
    The positive-control poisoning applies only to adaptive selectors. For static selectors, add: poisoning any screened
    row leaves selections bit-identical.
17. **(L)** Top-20 randomised greedy changes the estimand: say claims refer to it; add 1 deterministic run per finalist.
