# Design review r3: consistency, brief coverage, implementability (fresh external reviewer)
Read: PROMPT.md, DESIGN v3, r2 reviews (skimmed), config/{budgets,protocol,seeds}.yaml. No model run, no label read.
Verdict: **needs edits (listed)**. All are text/config patches. No new analyses. Brief coverage is complete apart from item 16.

## A. Contradictions (must fix before freeze)
1. **R can go down after the freeze.** §5 says "fixed at the freeze (may only go up)", but §12 says "if the final > 10 h → R at the pilot minimum".
   §12 replace with: "if the projected final (P4 timings × R) exceeds 10 h, R is set at the pilot minimum *before* the freeze; R is never changed after it."
2. **The warm-start contrast budget is not on any checkpoint list.** §6.3 contrasts "at 440 and 660", but §5 Q2 checkpoints have no 440.
   §5 add: "- Warm-start checkpoints: {220, 330, 440, 490, 660} (`budgets.yaml: warm_start_checkpoints`)."
3. **P\* can change after the P1 runs that use it.** P1 runs the power pilot and universe subsampling "with P\*", then P2 may replace P\*.
   §9 P2 append: "If P2 changes P\*, the power pilot and the universe-subsampling curves are re-run with the frozen P\* before P3 (the §12 50%-only trigger still applies)."
4. **m1_hgb budget sets disagree in three places.** §4 gives strategies {160, 660} and random {160, 490, 660} × 5. §7.1 gives {30, 160, 490, 660, 1310, 2627} × 5.
   §12 cuts this to "{160, 660, 2627} × 3", which drops 490 and the 5 trajectories that §4 needs.
   `budgets.yaml` has stale `m1_hgb_q1_grid: [30, 60, 160, 330, 660, 1310, 2627]` and `m1_hgb_q2_checkpoints: [60, 160, 330, 660]`.
   §12 replace with: "→ m1_hgb Q1-only budgets {30, 1310} dropped (the §4 Q2 set {160, 490, 660} × 5 and 2627 kept)".
   Rewrite the two budgets.yaml keys to match §4/§7.1.
5. **The off-array indicator ban (§1) conflicts with off-array routing.** §1 says "selectors may not use a 'missing SMF / off-array' indicator", but H, F and the HAP1-feasible pool all route off-array genes separately.
   §1 replace with: "...may not use such an indicator as a feature or score term; routing a candidate to an external-only submodel when its data are absent is allowed."
   §6.2 H prior: replace "off-array genes use abundance only" with "off-array genes score 2·rank(abundance) (same scale)".
6. **"PROVENANCE.md" in the §8 freeze list does not exist.** §1 names `docs/feature_provenance.md` and §6.1 names `config/feature_provenance.yaml`.
   Replace "PROVENANCE.md" with "config/feature_provenance.yaml (source; docs/feature_provenance.md rendered from it), DATA.md".
7. **Stale terms in §14 and config.** §14 statistics: "5 stratified B=10k binned" → "5 stratified bootstrap (v3: exact
   weighted AP, B = 2,000)"; "15 disjoint folds + confirmation + ledger" → "...+ re-split check (n ≤ 330) + ledger";
   "3 E_common" → "3 E_common (dropped in v3, §15)"; algorithms "(P5)" → "(P2)". protocol.yaml: `confirm_val_folds`
   comment "confirmation split" → "re-split check (n ≤ 330 only)"; `match_margin_frac` comment → "linear contrast H0:
   s160 − 0.95·r660 ≤ 0 (§6.4)"; add `g_star_min: 0.5`, gate 2% and surrogate 0.8/0.7 (§5 says they live there).
8. **§4 has a moot clause.** "`prop_folds = 5` when n < 30" can never trigger, because m1_hgb never runs below n = 30. Delete it.
   Replace "if one passes the §8 gate before the freeze" with "...before the P2 freeze of P\*", because "the freeze" otherwise means the P4 freeze.

## B. Blocking ambiguities for an implementer
9. **F at cold start is undefined.** π̄_S and ĉ_x do not exist at S = ∅. q̂ is only a rank for n < 30. "KRR into the same space" has no SVD space yet.
   §6.2 F append: "For n < 30 (incl. S = ∅): p̂ ≡ const, so with δ = mean p̂(1−p̂) the objective is Σ_x log det(I + Σ_s u_s u_sᵀ). ũ = top r−1 kernel-PCA components of the external kernel (whitened). The SVD space is used from n = 30. M_x := δI_r + Σ_{s∈S} p̂_sx(1−p̂_sx)u_s u_sᵀ."
   Also state whether F uses H's ≤ 2-per-cluster cap. "Reduces exactly to H" holds only if it does.
10. **M is undefined at round 0 and for off-array candidates**, which have no GI column. Append to M:
    "Round 0 = the first batch of the R0 permutation of the same trajectory. Candidates without a GI column (off-array) are never picked by M (reported)."
11. **Variance sources unspecified.** Randomised greedy: "sampled from the top-20" → "sampled uniformly from the top-20".
    R1 is static with no randomness: add "R1 trajectory r = seeded permutation of within-bin order".
12. **Batch sizes 50/110 do not hit the checkpoints.** n = 160 and 330 are not multiples of 50 or 110. §6.3 append:
    "For batch b, a checkpoint n not divisible by b uses the first n − b⌊n/b⌋ genes of the next batch in greedy order."
13. **It is not stated which arms run on T.** Q3 needs warm-start and batch-size results, but §9 runs §6.3 only in P4 (dev). §8 lets only "static counterparts" run on T.
    §8 replace "Static counterparts of finalists run on T" with: "R1, the static counterparts of finalists, and the warm-start and batch-size scenarios for the best ≤ 2 run on T as descriptive arms outside the family; the feasible-pool and minimal-features scenarios are dev only."
14. **§7.3 "OLS over the 9 draws" is ambiguous** (50% and 25% × 3 = 6). Replace with "over 9 draws (3 per level: 100% =
    dev splits 1–3, 50%, 25%; levels {100, 50} only, 6 draws, if the §12 trigger fires)".
15. **Smaller clarifications.** §3, after "reveals its whole query row": add "The predictor sees all columns of screened
    rows, incl. T (dev: D_val) columns (observed pairs per the brief); only the SelectorView masks them." §6.3 feasible pool:
    add "(off-array genes, lacking SMF, are kept)"; HAP1-minimal: add "(predictor unchanged)". §9 P3: "If the oracles reach
    G < 0.5" → "If both oracles have G < 0.5". Clusters "on non-test genes" → "...(dev: also excluding D_val)".

## C. Brief coverage
Covered: 20% stratified split fixed first (§2, 656/3,283); observed-regardless-of-side + symmetric rule (§3); AUPRC on
held-out queries (§3); Q1 (§7.1); Q2 batches of 10, matched budgets vs random, c1/a, G (§6); Q3 (§7.3–7.4); runtime
(§7.1, §12); no unscreened-query info in selection/tuning/stopping (§6.1, §7.4, §8; once item 1 is fixed); no GI/text/GO,
human analogues only (§1); seeds, freeze, ledger, independent rescoring.
16. **Gap: the brief's "never commit restricted data" and DATA.md upkeep are not in DESIGN.**
    §11 append: "- Repo hygiene: DATA.md row (source, version, checksum, license) per dataset is asserted at each milestone; raw/restricted data paths are git-ignored and a test fails if any is tracked."

## D. Optional (no change needed if declined)
- Symbol collisions: a, b in SE² vs slope b vs saturating-fit b (write SE² = V_Q + σ̂²_s/R_s + σ̂²_r/R_r); replicate
  array "A[b, …]" vs gene set A (rename W[b, arm, r, n]).
- §7.4 uses N_yeast = 3,283 while §7.3 uses "N = candidate-pool size". Both are fine because pool ∝ universe across levels, so b is unchanged. Say so in one clause.
