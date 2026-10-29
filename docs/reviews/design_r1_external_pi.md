# Design review r1 — external PI

Reviewed: DESIGN.md v1 (2026-10-29), against PROMPT.md, LITERATURE.md §7, docs/recon.md, data characterisation.
Overall: the evaluation core (frozen T, E_fixed, SelectorView + poisoning tests, dev-only tuning) is sound. My concerns
are (i) the Q3 numbers the consortium will quote have no pre-specified definition and no experiment that decides the
fraction-vs-count translation, (ii) the brief's "5% vs 20%" test can pass or fail for trivial reasons, (iii) the plan
has no early go/no-go evidence on how much Q2 can possibly gain, and (iv) the scenario does not match HAP1's real
starting point. Proposed changes, highest priority first.

## High priority

1. **Pre-register what "enough query lines" means (§6.4, §0 Q3).** Add to §0 a Q3 primary endpoint:
   n_X = smallest budget whose AUPRC (E_fixed, primary predictor) reaches X% of AUPRC at the full dev pool
   (n = 2,627, measured, not extrapolated), X ∈ {80, 90}; plus the marginal-return curve ΔAUPRC per +10 queries and
   Δ(observed strong SLs) per +10 queries. Report n_X for random and for the best frozen strategy, with bootstrap CIs.
   *Why:* without this, "estimated query numbers" will be chosen post hoc from a smooth curve; the saturating-fit
   extrapolation in §5 must not be the source of the headline number.

2. **Add an experiment that decides fraction- vs count-scaling (§6.4 "Scaling to HAP1").** Universe-subsampling on
   dev only: draw random gene sub-universes of 50% and 25% (same genes removed from query and array sides; 3 draws
   each), re-run the random Q1 curve with the primary predictor on inner splits, and test whether curves collapse on
   f = n/N or on n (pre-register the collapse criterion, e.g. smaller between-universe spread of n_90). Use the
   result to pick (not just list) scenario (i) or (ii); state that subsampling mimics fewer genes per module, not more
   modules, so it bounds rather than proves the HAP1 mapping.
   *Why:* LITERATURE §2.5/§7.2 shows this choice moves the HAP1 recommendation from ~200–600 to ~2,000 lines. It is
   the single most consequential number for a consortium, and the design currently offers two scenarios with no data
   to choose between them. Cheap (stack_lr_wr2, random only).

3. **Guard the brief's 5%-vs-20% example against a vacuous pass (§6.3 c).** Add a gap-closure statistic:
   G = (AUPRC_strat(160) − AUPRC_rand(160)) / (AUPRC_rand(660) − AUPRC_rand(160)), with CI, and require that the
   random gap AUPRC_rand(660) − AUPRC_rand(160) exceeds the 5% non-inferiority margin (lower CI bound) before
   reporting "matches". State the margin's justification (e.g. ≈ seed-to-seed SD of random at 660).
   *Why:* if the random curve is near-flat between 5% and 20%, non-inferiority passes for any strategy and the
   headline "5% chosen = 20% random" would be true but meaningless; G makes the result interpretable either way.

4. **Run cheap ceilings first, as a go/no-go for Q2 investment (new §6.0; sequencing).** On inner splits only, and
   labelled evaluation-only diagnostics (never part of any strategy): (a) oracle hub-first = rank D_cand by true
   query-row SL count; (b) greedy oracle with stack_lr choosing each pick from 100 random candidates by realised
   E_fixed_dev AUPRC. If even the oracles reach G < 0.5 at n = 160, de-scope V/VH/U and say so.
   *Why:* bounds what any strategy can achieve, tells us early whether the answer to Q2 is "no in principle" or
   "possible", and protects against over-investing in sophisticated selectors (literature prior: 1.3–2× savings).

5. **Add a status-quo comparator and a warm-start scenario (§6.2, §6 protocol).**
   (a) **R1 "Billmann-like"**: stratified sampling over SMF × expression-level bins, spread across coexpression/AP-MS
   clusters (LITERATURE §2.1, §7.4); static. (b) **Warm start**: S_0 = 160 genes chosen by R1, then compare random vs
   the best ≤ 2 strategies for the next 160–500 queries.
   *Why:* a HAP1 consortium is not at S_0 = ∅; it holds ~222 non-random lines. The actionable question is "how should
   the next lines be chosen, and is it better than continuing as we do?" Random-from-scratch alone does not answer it.

6. **Pre-specify ONE primary predictor for the Q1/Q2/Q3 headlines (§4, §5, §6.3).** State now: primary =
   stack_lr_wr2 (or the frozen successor from §4.1 if it passes the gate before freeze); others are secondary.
   n_req inversion (§6.3 b) must use the random curve of the same predictor as the strategy.
   *Why:* four predictors × two evaluation sets give eight Q1 curves; picking the headline afterwards is a
   garden-of-forking-paths risk.

## Medium priority

7. **Make the discovery objective a co-reported Q2 axis, not only a Q3 add-on (§6.4 coverage → §6.3).** For every
   strategy and budget, report (A) AUPRC on E_fixed and (B) strong SLs directly observed by S (cells and unordered
   pairs), and plot A vs B as a frontier. *Why:* LITERATURE §7.3/§5.1 — hub-first will win B and may lose A;
   consortia value both, and a single-axis headline would mislead.

8. **Prediction-vs-screen in the same units (§6.4 "Prediction contribution").** Drop the direct comparison with
   HAP1 precision 0.61/recall 0.37 (that is vs a replicate consensus; our labels are single-screen). Replace with a
   within-yeast reference: precision/recall of the reciprocal measurement (BA SL predicting AB label) on the
   test×test on-array cells of E_fixed, alongside the predictor's PR curve on the same cells.
   *Why:* tells the reader whether prediction at n screened is worth more or less than one extra noisy measurement —
   the realistic framing of "what computation contributes".

9. **Characterise E_fixed and pre-register the E_fixed/E_unobs disagreement rule (§3).** Before scoring, report
   E_fixed's array-degree distribution, SMF and share of each test row's SLs vs all arrays (array-only genes may be
   a biased class). Rule: a Q2 strategy conclusion is headline only if ΔAUPRC has the same sign on E_unobs.
   *Why:* E_fixed drops ~55% of arrays; a consortium will predict over the whole library.

10. **Batch-size sensitivity for the best strategy (§6).** Keep batch 10 as the brief's protocol; add batches of 50
    and one-shot (already = static) for the best strategy at n ∈ {160, 330}. *Why:* 66 rounds of 10 is not feasible in
    HAP1 (rounds take months); the consortium needs to know how much of the gain survives realistic round counts.

11. **Trim the strategy set (§6.2).** Drop U (uncertainty) unless the ceilings in #4 justify it — weakest support at low
    budget (Hacohen; IterPert cold start) and costly. Add LITERATURE's adaptive column sampling (largest residual of the
    candidate's array-side profile after projection on screened profiles) as a cheap alternative to V.
    *Why:* fewer, better-motivated arms; the freeze already caps at 3.

12. **Phase order with decision gates (§7/§10).** Write the sequence explicitly: P1 dev random curves (all budgets,
    primary predictor) + #2 subsampling; P2 ceilings (#4); P3 static heuristics (H, R1, C); P4 adaptive (V/VH or
    column sampling) + warm start; P5 predictor development §4.1 limited to the n ≤ 330 regime; freeze; final.
    Each gate records in RESEARCH_TREE.md whether to continue. *Why:* most decision-relevant results first; avoids
    re-running Q2 when the predictor changes late (predictor changes after P4 require re-running P3–P4 on dev).

13. **Deterministic strategies need trajectory variation (§6.3).** Run deterministic strategies with 5 seeded
    perturbations (random first batch of 10 or 80% candidate subsampling per round) rather than once. *Why:* a single
    trajectory confounds strategy with one lucky/unlucky early batch; CI from test-query bootstrap alone understates
    uncertainty relative to random's 10 trajectories.

## Low priority

14. **Query-availability sensitivity (§6.4).** Re-run the best strategy restricted to a random 70% "available" subset
    of D (proxy for catalogue lines / constructible KOs). *Why:* LITERATURE §7.11; tells whether gains rely on
    specific hard-to-get genes.
15. **Stricter-label sensitivity (§0).** Optional secondary for the final random curve only: ε < −0.35, P < 0.05 as the
    analogue of HAP1 "extreme SL". Must be labelled secondary; the brief's label stays primary. Skip if compute-bound.
16. **Per-query metrics (§3).** State handling of zero-SL test queries (excluded from per-query AUPRC, counted in
    pooled) and list seeds explicitly in config/seeds.yaml referenced by name in DESIGN.

## No change needed
- §1 data/feature admissibility (tier 1/tier 2 split and exclusions are correct and within scope).
- §2 test split (stratified, frozen, lockbox-preferred, contamination reported).
- §3 observation model and AB/BA handling (directed cells, no merge) — clearly stated as the brief requires.
- §6.1 SelectorView and poisoning tests; §8 one-shot final; §9 verification plan (add oracle/G checks to toy tests).
- Scope: no drift beyond the brief; keep positive GIs, ε regression, ortholog pair transfer out.
