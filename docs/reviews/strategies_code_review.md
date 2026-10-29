# Code review: strategy selectors (R1, H, H-static, F, F-static, C, M) and oracles

Reviewer: independent (did not write the code). Date 2026-10-29. Scope: `src/hsl/strategies/*`, `src/hsl/oracles.py`,
`config/strategies.yaml`, `tests/test_strategies.py`, `tests/strategy_toys.py`, against DESIGN v3.1 §6.1–§6.3, review
r2 items 13–16 and decision D5. New tests: `tests/test_review_strategies.py` (toy data only; 42 pass, 7 strict xfail
that mark the open issues below and will turn into failures once those issues are fixed). Real-data diagnostics used one dev view (`dev1`, one seeded random S at
n ∈ {0, 20, 60, 160}). They are first-pick gain rankings and bin counts, not strategy-vs-random evaluations.

## Verdict: **safe with fixes**

The numerics are right. F's log-det gains and Sherman–Morrison downdates, the cold-start log det, M's GP posterior
variance with within-batch pseudo-observations, C's facility-location greedy, the degree-model ridge with GCV, and the
shrunk ĉ / logit response all match independent brute-force reimplementations exactly, including the deterministic
greedy with the cluster cap. Information hygiene is good. However, **F and C, as implemented, almost never pick
off-array genes once screened rows exist. R1 under-samples them by about half.** In every case the cause is routing:
the gene lacks a column, not anything about the gene itself. D5 fixed this failure for H, and DESIGN §1 forbids it.
Fix F-1, C-1 and R1-1 before any P4 dev comparison involving F, C or R1 (R1 is also the warm-start S0). H, H-static,
F-static and M can be relied on as specified, with the caveats below.

## Issues

| id | sev | where | issue |
|---|---|---|---|
| F-1 | **High** | `fisher.py:128-144` | KRR route for off-array candidates acts as an off-array score term |
| C-1 | **High** | `coverage.py:36-41` | GI kernel part exists only for on-array candidates, so off-array genes are never picked |
| R1-1 | Medium | `r1.py:39-56` | ridge-imputed SMF puts off-array genes in one SMF tertile, so they are under-sampled |
| B-1 | Medium | `base.py:113-128, 141-147` | fixed per-batch cluster cap cannot fill large batches; the rest is filled at random |
| M-1 | Medium | `maxvar.py:33-44`, `strategies.yaml:44` | λ = 10 makes the posterior variance ≈ the prior variance, so M ≈ array-side column-norm ranking |
| H-1 | Low | `degree.py:191-213` | residual off-array enrichment of fitted H (30% of top-50 vs base 20%) |
| L-1 | Low | `harness.py:341,358` | selector diagnostics (phase, `n_filled_at_random`, λ's) are not persisted in trajectory logs |
| L-2 | Low | RESEARCH_TREE D4 vs `kernels.py:156-184` | D4 says R1 clusters include D_val; code (correctly, per DESIGN) excludes it |
| L-3 | Low | `maxvar.py:70-76`, `base.py:141-147` | M's round 0 (R0 batch) and its exhaustion fill can include off-array genes |
| L-4 | Low | `degree.py:133-150` | H prior differs from DESIGN text (D5 quantile map replaces 2·rank(abundance)); DESIGN not amended |
| L-5 | Low | `kernels.py:65-89` | kernel approximations (AP-MS neighbourhood over pool, ESM-16, no pair stores) |
| L-6 | Info | `fisher.py:172`, `degree.py:185` | π̄_S = 0 → NaN gains → whole batch filled at random (cannot happen on real data at n ≥ 30) |

### F-1 (High): F's off-array route is a de-facto score term
*Mechanism.* On-array candidates get ũ = B·v_g from their column factor. Off-array candidates get ũ from a KRR on the
external kernel. The KRR route systematically yields smaller u_gᵀM_x⁻¹u_g (shrinkage and direction), so the gain
Σ_x log(1 + w_gx u_gᵀM_x⁻¹u_g) is lower for every off-array gene.
*Evidence.*
- Null toy world (genes exchangeable): off-array share of F's picks after the n = 30 switch is 0.000 (r = 8) and
  ≈0.02 (r = 4), against a base rate of 0.167. H gives 0.15–0.16. See `test_null_world_offarray_share_F` (xfail).
- dev1, n = 60 and 160: off-array share of F's top-200 first-pick gains is 0.000 and 0.010 (base 0.195). F with
  r = 1 (≡ H) gives 0.225.
- **Route-swap control** (dev1, n = 160): sending 300 random *on-array* candidates through the KRR route drops their
  share of the top-200 from 0.155 (base) to 0.010. Median gain: B-route 24.8, KRR-route on-array 17.3, off-array 18.3.
  The route causes the effect, not the genes.
- Rescaling the KRR-route ũ to the on-array per-component mean square does not help (still 0.0/0.01).

*Consequence.* F's selections are effectively on-array only. In yeast that discards the off-array genes, which are
SL-richer (D5: 52 vs 37). This is a yeast-specific artefact (every HAP1 gene is in the library), and it biases F vs
random/H in the Q2 contrasts.
*Fix (needs a logged decision, like D5).* Preferred: use **one route for all candidates**, i.e. KRR from the external
kernel into the SVD space (on-array column information still enters through q̂ and ĉ). Alternative: use the
expected outer product E[uuᵀ] = ûûᵀ + Σ_post for both routes, so each route's uncertainty is restored. A third
option is D5-style quantile mapping of off-array first-pick gains onto the on-array distribution under a common
route. Acceptance: `test_null_world_offarray_share_F` must XPASS (then remove the marker), and the dev route-swap
share must be within noise of base.

### C-1 (High): C never picks off-array genes once S has on-array rows
*Mechanism.* K = K_ext + Z Zᵀ, where the GI rows Z are zero for off-array candidates. On-array candidates therefore
have a larger self-coverage K(c,c) and extra similarity mass to every other on-array candidate, which inflates
Σ_x max(0, K(x,c) − cov_x).
*Evidence.*
- Null toy, per-round off-array share: 0.125 in round 0 (external only), then 0.000 in every later round
  (`test_null_world_offarray_share_C`, xfail).
- dev1: off-array share of the top-50 gains is 0.00 at n = 20, 60 and 160 (base 0.195).
- The existing knob `C.normalize: true` over-corrects: 0.68 (n = 60) and 0.86 (n = 160) of the top-50.

*Consequence.* This is the same §1 problem as F-1. DESIGN states C's "target pool = unscreened candidates", but in
effect only on-array candidates are ever selected.
*Fix (logged decision).*
- Option (a): stratify by route. Give each batch off-array slots in proportion to their share of the candidates
  (systematic or seeded-Bernoulli), and run the greedy on the joint objective within each stratum. This is
  routing-only and neutral.
- Option (b): impute off-array GI features by KRR from on-array candidates' Z_c and match their norm distribution
  to the on-array norms.

Verify with the null test and the dev share.

### R1-1 (Medium): off-array genes collapse into one SMF tertile
*Mechanism.* Off-array SMF is replaced by an external-feature ridge prediction. That prediction regresses to the
mean, and with the library's skewed SMF distribution the mean lies in a single tertile.
*Evidence.*
- dev1, n = 0: 347 of 410 off-array candidates fall in the lowest-SMF tertile (bins 0–2 hold 319/326/398 genes vs
  134–228 elsewhere).
- Because bins are visited round-robin, R1 picks 10.5% off-array in its first 220 and 13.2% in its first 660
  (base 19.5%).
- Null toy: all 40 off-array genes fall in the middle tertile (`test_R1_offarray_spread_over_smf_tertiles`, xfail).

*Consequence.* This biases the status-quo comparator, the warm-start S0 (220 R1 genes) and the primary warm-start
contrast (strategy vs R1-continuation).
*Fix.* Map off-array imputed SMF onto the on-array distribution, as D5 does for H. Apply the same ridge to on-array
candidates, take each off-array gene's percentile among those predictions, and map it to that quantile of the
on-array candidates' actual SMF. Alternatively, assign tertiles by rank within the off-array set.

### B-1 (Medium): cluster cap in large batches
`capped_update` allows ≤ 2 picks per cluster per batch. With k = 100 clusters at most 200 genes per batch can be
picked by gain. For §6.3's one-shot at n = 330 or 660 (and, near exhaustion, any batch), `Selector.select` silently
fills the rest **uniformly at random**. On a toy with 10 clusters, a 60-gene batch gets 40 random fills
(`test_large_batch_not_filled_at_random`, xfail ×4).
*Fix.* Scale the cap with batch size (2·⌈b/10⌉), or apply it per consecutive block of 10 in greedy order. When all
clusters are full, reset the counts rather than fill at random. Record which rule was used.

### M-1 (Medium, interpretive): M is effectively "largest GI column norm"
λ = lam·w_gi = 10, while the K_SS eigenvalues on dev1 at n = 60 have median 0.23 (max 13). The conditioning on S is
negligible: Spearman(posterior variance, prior k(g,g) = ‖c_g‖²) = 0.998 over on-array candidates. M therefore ranks
candidates by the magnitude of their array-side column over the screened rows, which is a measure of array-side GI
hubness, not uncertainty. This matches the spec, which gives no λ (the code uses the upstream default). Either add λ
to the dev grid ({0.1, 1, 10}) or describe M accurately in the report.

### H-1 (Low): residual off-array enrichment in fitted H
With D5 calibration, the off-array share of fitted H's top-50 is 0.30 and 0.32 on dev1 (n = 60, 160), and 0.23 in
the top-200, against a base of 0.195. The prior is unbiased: 0.168 ± 0.006 over 150 null worlds, and null-world
trajectories are at base rate. This is consistent with D5's intent: the external features themselves rank those genes
high. But the external submodel is trained on all screened rows, including off-array rows that carry the yeast SL
excess, so features correlated with off-array status can still carry the artefact. Monitor the per-round off-array
share in P4. Optional ablation: fit the external submodel on on-array rows only.

### L-1 … L-6
- **L-1.** `Trajectory.run` logs only picks and timings. Selector `diagnostics` are dropped, including M's off-array
  count (which DESIGN requires to be reported), `n_filled_at_random`, phase, and the GCV λ's. Fix: append
  `getattr(sel, "diagnostics", [])[-1]` to each round record.
- **L-2.** D4 says the R1 and H-cap clusters include D_val. Actually the H/F cap uses `features.clusters` (non-test
  genes, incl. D_val), while R1 uses `kernels.visible_clusters` (excludes D_val, per DESIGN §6.2). Correct D4.
- **L-3.** "M round 0 = R0 first batch" and "off-array never picked by M" conflict: round 0 can contain off-array
  genes, and when no on-array candidates remain the batch is filled at random (`test_all_candidates_offarray[M]`).
  Document this.
- **L-4.** The DESIGN §6.2 H prior text ("off-array: 2·rank(abundance)") is superseded by D5. Add it to the
  deviations list.
- **L-5.** None of the three kernel approximations leaks labels or is tied to off-array status.
  - The AP-MS neighbourhood is computed over the pool only, so it depends on the pool. On dev1 it correlates 0.88
    (non-zero pairs) with neighbourhoods over the visible universe. 566 of 1,539 AP-MS-covered pool genes have no pool
    partner (39% of on-array vs 27% of off-array genes). The final pool (2,627) is larger than dev1's (2,101), so the
    geometry that C and F (cold start/KRR) are tuned on in dev shifts in the final. Recommendation: build
    neighbourhoods over the visible universe.
  - ESM-16 standardised cosine weights all 16 PCs equally. It has weight 0.1 in the selector kernel but 1/3 in R1's
    clustering.
  - With no pair-store kernels, C's kernel is not literally the predictor's `mk` kernel. This is an interpretive
    caveat.
- **L-6.** If no screened cell is a strong SL, π̄_S = 0 and F's p̂ is NaN. Only possible on toys. A guard
  (fall back to the cold start) would be cheap.

## Checks that passed
- **Spec fidelity.**
  - H: prior with D5, h_switch_n = 30, response logit((k+0.5)/(m+1)), ĉ with n0 = 100 (entered as a logit), GCV
    ridge with unpenalised intercept, external-only submodel, quantile calibration, cap 2.
  - F: u = (1, ũ), top r−1 left singular vectors of the row-centred −clip(ε) visible block, LS map B, KRR,
    unit-mean-square whitening over S ∪ C, δ = mean p̂(1−p̂), p̂ clip [1e-4, 0.5], P = on-array candidates ∪
    array-only visible columns, determinant lemma + Sherman–Morrison within the batch, M_x rebuilt each round,
    r = 1 ≡ H (existing test). Cold start: kernel-PCA log det.
  - C: facility location, target = candidates, round 0 on the external kernel.
  - M: GI-only posterior variance, round 0 = R0 batch from the same selection stream.
  - R1: tertile bins × cluster round-robin with seeded within-bin order.
- **Brute force** (new tests): ridge GCV against the explicit hat matrix; full degree model reimplemented; F fitted
  greedy (r = 1, 4) against slogdet, including the cap; F embeddings against an independent SVD, lstsq B, KRR and the
  whitening; F cold greedy against an independent dense kernel PCA; M greedy against an exact joint-GP solve each
  step; C greedy against the explicit facility-location objective (normalize on and off). Mutation checks confirmed
  these tests catch broken updates and embeddings. (An overall scale on ĉ is, correctly, order-invariant.)
- **Hygiene.**
  - Selectors read only SelectorView attributes (candidates, candidate_col, screened_orfs, screened_col,
    n_screened, visible_cols, array_smf, array_orfs, E, Y, gene_features, similarity, clusters). None reads SP, DMF,
    query SMF, split files or bundles. Missing values are mean-imputed after standardisation; no indicator is
    created.
  - Static variants (R1, H-static, F-static, C-static) give identical picks when E/Y/SP are scrambled (new test),
    in addition to the existing TrapRows test.
  - The oracle separation holds (AST check plus fresh-process `sys.modules` check). Note that `hsl.strategies`
    imports `hsl.harness` → `hsl.views`, so isolation from bundle loaders is by convention and tests, not physical.
- **Edge cases.** n = 0; the 29/30 switch (phase diagnostics); all candidates off-array; pool exhaustion (7 left,
  k = 10); remainder batch (k = 3) all return valid batches. No exact ties among the top-20 gains on dev1 (ties
  would break by index, i.e. alphabetical ORF order).
