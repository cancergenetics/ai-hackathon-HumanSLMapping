# Design review r1 — algorithms & compute (active learning, OED, kernel methods, matrix completion)

Reviewer lens: are the Q2 selection strategies well specified and the right set; what is likely to win for "pooled
AUPRC on E_fixed predicted by fold-in"; implementable algorithms + complexity; compute feasibility (M2, 2 jobs × 3
threads); low-budget predictor fixes. Read: PROMPT.md, DESIGN.md v1, docs/recon.md (+ runtime table),
LITERATURE.md §0/§4/§5/§7, third_party/ysl/models/{foldin2,completion,baselines,registry}.py,
results/data_characterisation/summary.md. No models run, no test metrics computed.

## A. What the predictor needs from the training rows (analysis that drives the proposals)
- `mk` fold-in = KRR from x's array-side column c_x (over screened on-array genes S_on) + external kernels to x's
  query row. Screening g does two things: (i) adds a training pair (sym-profile_g → row_g); (ii) adds an input
  coordinate ε(g, ·) for every unscreened x. A coordinate is informative only if row g has many strong ε — i.e. g is a
  query-side hub. At n = 10–100 most c_x are near-zero unless S contains hubs.
- Pooled AUPRC at 1.1% prevalence is dominated by the row × column propensity structure. Under a rank-1 logistic
  model logit p(s,x) = a_s + b_x, the Fisher information about b_x (test-row propensity) and about the column rate of
  array a from screening s is p(1−p) ≈ p_s·, i.e. ∝ s's query-side degree. Hence **hub-first is exactly the D-optimal
  design for the rank-1 propensity model**; diversity only enters through the rank-r (module) part. Prior: H or a
  hub-weighted diverse design wins at n ≤ 330; pure V/C (variance-only, label-agnostic) are likely to lose because the
  GP variance ignores that rows differ ~10× in SL content (Gini 0.49).
- Caveat that DESIGN v1 misses: query-side vs array-side SL counts of the same gene correlate only Spearman 0.31
  (summary.md). H's "later rounds from observed array-side SL rate" therefore ranks hubs poorly; the degree
  predictor must be a model of *query-side* degree trained on screened rows (array-side rate is one input).
- Second caveat: the propensity map (array-side rate + externals → query-side degree), which m1_hgb learns and which
  drives pooled AUPRC, is a regression whose training design benefits from *spread* over predicted degree, not from
  range-restricted hub-only rows. So pure hub-first may beat random early and then plateau; a hybrid is the likely
  winner by n ≈ 160–330.

## B. Proposed changes (numbered; priority H/M/L; cost = implementation + compute)

1. **Add a principled hybrid "F" (Fisher-D, rank-r logistic design) and freeze it as a primary candidate. [H]**
   - Model: p̂(s,x) from a cheap selector model (prop. 4); embeddings u_g ∈ R^r (r = 8–16) from the truncated SVD
     of the screened block (rows S × visible non-test columns, centred), using each gene's *array-side* factor
     (projected column c_g) as its embedding; off-array candidates: kernel-ridge prediction of u_g from the external
     kernel. Round 0 / static: u from kernel PCA of the external kernel, p̂ from the external degree proxy.
   - Objective over the target pool P = (unscreened dev queries D\S) ∪ (visible array-only columns):
     F(S) = Σ_{x∈P} w_x log det(I_r·δ + Σ_{s∈S} p̂_sx(1−p̂_sx) u_s u_sᵀ). Concave-of-modular in PSD form → monotone
     submodular; greedy batch of 10 with (1−1/e). Gain of g via the determinant lemma: Σ_x w_x log(1 + p̂_gx(1−p̂_gx)
     u_gᵀ M_x⁻¹ u_g); after each pick Sherman–Morrison on every M_x⁻¹.
   - Cost: O(|P|·r²) per candidate → ~2.6k candidates × 3.8k targets × r²=144 ≈ 1.4 G flop ≈ 1 s per pick; ≈ 10–15
     s per round incl. SVD. r = 1 reduces F exactly to H (built-in ablation); w_x = 1 vs predicted degree is the only
     other knob. Replaces the ad-hoc "VH = variance × predicted degree" weighting.
2. **Unify C, V and the Krishnamurthy–Singh column residual in one cached-kernel GP module; keep V and C, add the
   residual as "max-variance". [H]** Kernel = exactly the `mk` kernel (GI 'sym' linear kernel over S_on + MK_RECIPE
   external weights, λ as `lam_mode='gi'`), so selection and prediction share geometry.
   - Precompute once per inner split: external pool×pool kernels for all non-test genes (≈3.8k², 4 kernels, one
     PSD projection each: ~4 × 10 s, 115 MB each float64; cache to disk). GI part per round: Z_P = c over S_on
     (|P| × n_on), O(|P|² n_on).
   - Posterior cov Σ_P|S = K_PP − K_PS(K_SS+λI)⁻¹K_SP: O(|P|² n + n³) ≤ 3.8k²·660 ≈ 10 Gflop ≈ 2–4 s/round. V gain
     of candidate c: ‖Σ_{P,c}‖²/(Σ_cc+λ); after each pick rank-1 update O(|P|²) = 15 M. C = facility location
     Σ_x max_{s∈S} K(x,s) (lazy greedy, O(|P|·|cand|) per pick). Max-variance = argmax Σ_cc (≡ residual of c_g after
     projection on screened profiles under the GI-only kernel = [18]).
   - State explicitly: with external-only kernels V/C are label-free → adaptive and static coincide; their
     adaptivity comes only from the GI part (new coordinates). Say so in §6.2 so the static/adaptive contrast is read
     correctly.
3. **Drop U (bootstrap disagreement) or replace it by max-variance from prop. 2. [H]** Reasons: objective mismatch
   (own-row uncertainty is a hit-discovery/model-change criterion, weak at low budget [17, 105, 106]); cost: B ≥ 5
   refits of stack_lr_wr2 on ~2k targets ≈ 5 × 40–60 s × 66 rounds ≈ 4–5 h per trajectory → >100 job-h over dev
   + final. Max-variance gives the "uncertainty" arm at ~0 cost.
4. **Specify the selector's degree model precisely (shared by H, F, VH). [H]** Response: query-side SL rate of
   screened rows (logit, empirical-Bayes shrunk); inputs: array-side SL rate over S (shrunk; NaN→flag for off-array),
   array SMF, abundance, AP-MS degree, coexpression degree, conservation, paralog flag; model: ridge/LR (not HGB) with
   inner CV on screened rows only; refit each round (≤ 0.1 s). Before n = 30: fixed external prior = −SMF rank +
   abundance rank (pre-stated, no fitting). Report its Spearman with true query degree on D\S per round in dev
   (allowed: dev pool) as a diagnostic. p̂(g,x) for F: rank-1 form σ(â_g + b̂_x) from this model and from column rates.
5. **Add a dev-only oracle-headroom run. [M-H]** On each inner split, greedy forward selection maximising
   AUPRC on E_fixed_dev of the closed-form `F2_mk_foldin` (~2–4 s per fit at n ≤ 330), scoring 200 random candidates
   per pick-batch (batch = top-10 by marginal gain, one refit per candidate). ≈ 200 × 3 s × 16 rounds ≈ 2.7 h per
   split to n = 160 (oracle is evaluation-only, uses D_val labels, never test). Gives the ceiling for Q2(c): if the
   oracle cannot make 160 ≈ 660-random, no honest strategy will, and the result is already interpretable. Also
   profile oracle picks (degree, SMF, module spread) to guide strategy design and the Q3 "choice of queries" rules.
6. **Cheap surrogate for strategy tuning, full model only at checkpoints. [H]** Tune strategy knobs (r, w_x, kernel
   weights, switch round for hybrids) using `F2_mk_foldin` AUPRC on E_fixed_dev; confirm with stack_lr_wr2 at the 8
   checkpoints only. Pre-check surrogate fidelity: Spearman between surrogate and stack_lr_wr2 AUPRC across the 10
   random dev trajectories × checkpoints must be ≥ 0.8, else tune on stack_lr_wr2 at 3 checkpoints (30/160/330).
7. **Rule: no selector may call stack_lr_wr2/m1_hgb per round. [H]** Allowed per-round compute: prop. 2/4 modules
   and at most one `stack_lr` fit (~2–5 s on 2k targets). Record selection seconds per round (§10 already asks).
8. **Add `B2_product` (array rate × query propensity; upstream registered) to Q1 and Q2 as the GI-only
   degree-product baseline. [H, ~0 cost]** B1 cannot rank rows, so it is not the lift reference LITERATURE §4/§7.7 asks
   for. Also add a dev-only diagnostic: AUPRC of oracle-row-degree × array rate on E_fixed_dev, to quantify how much of
   pooled AUPRC is propensity (tells whether F/H should target propensity or modules).
9. **Add a static "Billmann-like" baseline (stratified on SMF × abundance bins, spread over AP-MS/coexpression
   clusters). [M, < 1 h]** It is what HAP1 actually did (LIT §0.2, §7.4); Q3 needs "better than current practice", not
   only "better than random".
10. **Trajectory variance for deterministic strategies. [M]** Run each frozen adaptive strategy from 3 seeded
    starts (first batch = seeded random from the static variant's top-50, or randomised greedy: sample from top-20 by
    gain). Without this, ΔAUPRC CIs for deterministic strategies only reflect test-query bootstrap, while random
    carries trajectory variance — an asymmetric comparison. Cost ×3 on strategy trajectories (see D).
11. **Within-batch redundancy for H. [M]** Top-10 hubs are often one complex (near-identical rows). Specify H's batch
    rule: greedy with a cap of ≤ 2 genes per external-kernel cluster (or use F with r ≥ 8, which handles it).
12. **SelectorView: explicit test-column masking when reusing upstream ScreenView internals. [H]** A ScreenView built
    with targets = D\S still carries test genes' array columns in view.E/Y; SVD (prop. 1) or column rates (prop. 4)
    computed on it would read them. Mask before any selector computation; the poisoning test in §6.1 must poison
    test-gene *columns of screened rows* too (not only test rows).
13. **Pool-size mismatch dev vs final. [L]** D_cand ≈ 2,100 vs D = 2,627: H exhausts the top-10% hubs by n ≈ 210 in
    dev vs 263 in final. Express hub cut-offs/switch points as fractions of the pool, not absolute n.

## C. Low-budget predictor improvements (within §4.1 rules; dev only)
14. **n-dependent GI weight in `mk`. [M]** At n ≤ 100 the GI kernel has ≤ 100 noisy coordinates and fixed external
    weights are under-weighted. One parameter: w_gi(n) = n/(n + n0), n0 ∈ {30, 100, 300}, tuned on inner splits with
    the closed-form F2_mk_foldin (~3 s/fit, 3 splits × 4 budgets × 3 values ≈ 2 min). Also re-check `masked_pearson
    min_n=10` and k = 50 at n ≤ 30 (both degenerate; recon §3).
15. **Query-propensity feature for stack_lr_wr2 at low n. [M]** Upstream `stack_lr_prop` already exists; the
    array-side rate alone is a weak propensity (ρ 0.31). Evaluate stack_lr_prop / wr2+prop at n ∈ {30, 100, 330}; cost
    ≈ stack_lr_wr2 + 10–20 s. Importance weighting for covariate shift: L priority (LR stack is mostly invariant to
    row-selection shift on the features it uses; test only if F/H show a gap between stack_lr_wr2 and m1_hgb).

## D. Compute estimate (job-hours at 3 threads; 2 concurrent ⇒ wall ≈ ½)
| block | assumption | job-h |
|---|---|---|
| Q1 stack_lr + wr2 + B2 | 10 seeds × 13 budgets; wr2 ≈ 580 s/seed (n = 1970/2627 ≈ 130/220 s, run alone ≤ 12 GB) | 1.9 |
| Q1 m1_hgb | 5 seeds × 13 budgets ≈ 2,500 s/seed | 3.5 |
| Q2 dev random | 3 splits × 10 traj × 10 ckpts × ~18 s | 1.5 |
| Q2 dev strategies | 3 splits × ~8 strategies × 3 starts × (8 ckpts × 18 s + 66 rounds × ≤ 5 s) ≈ 8 min | 9.6 |
| Q2 dev oracle (prop. 5) | 3 × 2.7 h | 8.1 |
| Q2 final | (3 frozen × 3 starts + 3 static + 10 random) × 8 ckpts wr2 + m1_hgb at 4 ckpts | 5.5 |
| total | | ≈ 30 job-h ≈ 15 h wall |
Cuts if needed (in order): oracle to 1 split (−5.4); m1_hgb Q1 to budgets {33, 164, 328, 657, 1313} × 3 seeds
(−2.5); dev starts 3 → 1 for strategies eliminated after split 1 (−4). Infeasible variant to avoid: full-stack
per-round selectors (U, or model-based degree from stack_lr_wr2 on all candidates): ≈ 1–5 h per trajectory.

## E. Expected outcome (stated as a prior, to be tested)
F ≥ H > random > V ≈ C at n ≤ 160 for pooled AUPRC; V/C competitive on per-query AUPRC; static H/F close to adaptive
(most signal is in the external degree proxy), which would argue for few HAP1 rounds. "5% chosen ≈ 20% random" is
unlikely; 1.5–2× budget efficiency is plausible. The oracle (prop. 5) bounds this.
