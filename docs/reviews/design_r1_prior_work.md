# Design review r1 - prior work lens (DESIGN.md v1 vs LITERATURE.md, docs/recon.md)

Verdict: the design is well grounded (fixed E_fixed per Farquhar/Lowell, cold-start awareness, query-level holdout,
prevalence-matched label). No contradiction with the literature found. The gaps are (a) missing baselines that
reproduce what yeast/HAP1 groups actually did, (b) the fraction-vs-count translation is stated but never tested,
(c) the upstream sanity anchor is not like-for-like, (d) several cheap anchors that make results comparable to
Billmann 2026 and to IterPert/BATCHIE-style claims. Numbers checked: HAP1 precision 0.61/recall 0.37, 222 queries,
1.20% negative density, 11,000 expressed / 1,524 essential (~9,500 nonessential), 47,052 negatives, <5% coverage
(222 x 11,000 / 60.5M = 4.0%; 1-(1-222/9,500)^2 = 4.6%) are all consistent between DESIGN and LITERATURE and
arithmetically reproducible. Yeast figures (prevalence 0.0114, 39.3 SLs/query, 5%=164, 20%=657) agree with
results/data_characterisation/summary.md.

## Proposed changes
1. **[HIGH] Test fraction-vs-count scaling empirically (LITERATURE 7.2, absent from DESIGN 6.4).** Q3 currently
   lists two scaling "scenarios" but no evidence for choosing between them; this decides whether HAP1 needs ~200
   or ~2,000 lines for "20%-equivalent". Add a universe-subsampling experiment: restrict the gene universe to a random
   50% (and 25%) of genes as queries AND arrays (screened from D, evaluated on T intersect universe), run stack_lr /
   stack_lr_wr2 random curves, and plot AUPRC against both n and f. If curves collapse on n, HAP1 scales by count;
   if on f, by fraction; report the intermediate case honestly. Cheap (stack_lr, 5 seeds). Pre-register the collapse
   criterion (e.g. which axis gives smaller between-universe spread of AUPRC at matched x).
2. **[HIGH] Add the baselines that correspond to real practice (LITERATURE 7.4; DESIGN 6.2 has only R0/H/C/V/U).**
   - **P (Costanzo-2010-like):** random with preference for low SMF (fitness-defective genes).
   - **S (Billmann-like):** stratified random over expression x SMF bins, spread over coexpression/AP-MS clusters.
   - **TypiClust / k-means in the feature kernel** as the named diversity baseline (best generic AL baseline in
     IterPert; strong in GeneDisco model-improvement); C should be run in at least one of these forms, not only
     facility location.
   - **ACS (adaptive column sampling, Krishnamurthy-Singh 2013; Bhargava 2017):** residual of the candidate's
     observed array-side partial column after projection on screened profiles. It is the best-supported
     non-circular criterion for objective A in the literature and is not in the list (V is related but uses
     the KRR posterior variance on feature kernels, not the GI column). Leverage-score / Nystrom sampling is its
     static analogue.
   Without P and S the report cannot say whether Q2 gains are over "what labs do" rather than over uniform random.
   Keep the frozen-final cap (<=3) but allow these cheap heuristics as additional reference rows (no tuning).
3. **[HIGH] Make the upstream sanity anchor like-for-like (DESIGN 9 vs recon.md 5).** Upstream's 0.083 at 20% scores
   every observed cell of target rows, including ~14% of cells whose array is screened and whose reciprocal epsilon
   is a feature. E_fixed and E_unobs exclude those by construction, so the check "random n=657 within upstream
   range (0.076/0.083)" will not be comparable. Add one explicit reproduction row on the upstream evaluator
   (all observed cells, dev-side only, not test) and report both numbers plus the drop attributable to excluding
   reciprocal-known cells. Also state which of E_fixed / E_unobs / upstream-style is the headline for
   cross-study comparison (recommend E_unobs as the brief-faithful one, E_fixed as the budget-invariant primary).
4. **[HIGH] Calibrate simulation realism against Billmann's own distribution (cheap, uses HAP1 Data S14 and
   dev-pool yeast only).** Subsample 222 random dev queries in yeast and compare array-side degree distribution with
   HAP1 (median 1, 25% zero, p95 10, p99 20, max 46 over 223 queries; mean ~2.7). Compute the HAP1 Gini the DESIGN
   promises for scenario (ii) from Data S14 rather than leaving it implicit. If shapes agree, state this as the
   justification for reading yeast AUPRC curves as HAP1 forecasts; if not, use the count/degree-matched scale.
5. **[MED] Anchor Billmann's 222 explicitly in budgets.** Add n=220 to Q1 and Q2 checkpoints (222 = 6.8% of 3,283 on
   the fraction-of-yeast scale; 2.3% of ~9,500 on the HAP1-fraction scale; pair coverage 4.6% matches yeast f~2.3%,
   i.e. n~75). Report all three anchors on the Q1/Q2 curves; the current grid (100, 164, 250) brackets but never
   lands on 222. Also align Q1 and Q2 grids (Q1 uses 164/657, Q2 160/660): the matched-budget random comparator must
   be a measured random point at the same n, not an interpolation. Use 160/660 (multiples of 10) in both.
6. **[MED] Report both Q2 objectives (LITERATURE 5.1, 7.3; GeneDisco, Sutherland 2013).** DESIGN's Q2 primary
   endpoint is objective A only. Hub-first (H) will lose on A and win on B (hits found directly, "total known" =
   observed + predicted at fixed precision); the literature warns these trade off. Promote "total strong SLs known vs
   n (screened + predicted at precision >= 0.5)" to a co-primary descriptive endpoint of Q2, not only Q3, so strategy
   recommendations state which objective they serve. Add a degree-product baseline (predicted query degree from SMF/
   expression x array rate; Koch 2012) because B1 uses array degree only and degree dominates pooled AUPRC; report lift
   over it and the n=0-feature-only SMF baseline value in Q3.
7. **[MED] Follow cold-start guidance in the strategy set (Hacohen 2022; IterPert; LITERATURE 5.5).** At n<~100 typical/
   diverse selection beats uncertainty/model-based; U (disagreement) is least supported at low budget and should be
   dropped or run only after n>=100. Add an explicit phased hybrid: C (feature kernel) until n~100, then V/H/ACS,
   with switch point tuned on inner splits only. H's switch at 30 screened rows is likely too early; tune it.
   Specify greedy within-batch conditioning (update residual/posterior after each of the 10 picks) for C, H, U to avoid
   duplicates (BatchBALD/k-center rationale).
8. **[MED] Add a stricter-label sensitivity (LITERATURE 7.1).** Score the same predictions against yeast epsilon < -0.35
   (~HAP1 "extreme SL", 0.14%); the brief fixes the primary label, this is a secondary analogue. Also state that the
   yeast strong-SL label corresponds to the HAP1 *standard negative* (1.20%, |qGI|>0.3, FDR<0.1) not the extreme SL,
   and that single-screen precision 0.54 (yeast) vs 0.61 (HAP1) sets a shared label-noise ceiling.
9. **[MED] Fix the apples-to-oranges comparison in DESIGN 6.4.** "Predicted precision at top-k vs HAP1 single-screen
   precision 0.61/recall 0.37": 0.61/0.37 is replicate reproducibility of a *screen*, while measured precision of a
   prediction against single-screen labels is itself capped near label precision (~0.54). Express the result as
   "k at which predictions reach single-replicate-equivalent precision" and note the ceiling; do not claim predictions
   are "as good as a screen". Cite Temerinac-Ott 2015 and Kangas 2014 as precedents for the online stopping rule.
10. **[MED] Anchor the "choice of queries" rules to HAP1 data and practice (Q3).** Compare selected-gene feature profile
    with Billmann Data S14 degree correlates (HAP1 SMF, DepMap CERES mean/SD, gene age, s_het/pLI) and Koch 2012
    (SMF, disorder, PPI degree, expression variation): report rank agreement of feature importances. Add practical
    constraints from LITERATURE 11: 83% of Billmann queries were Horizon catalogue lines (~7,500 vendor models),
    essential genes need hypomorphs, and a ~1-month round time (IterPert) implies few rounds. Hence add a batch-size/
    rounds sensitivity at one or two budgets (batch 10, 50, 100, one-shot static) so the HAP1 recommendation can be
    stated as "rounds needed", directly comparable to IterPert's "1 round ~ 5 rounds of uniform".
11. **[LOW] Report efficiency in the literature's currency.** IterPert ("a third of the perturbations"), BATCHIE (~4% of
    space), Kangas (60% of hits at 3%): present Q2 as n_req/n fold-savings (already planned) and as the fraction of
    final random-20% AUPRC reached at 5%. Also temper expectation in the pre-registration: LITERATURE 7.10 judges
    1.3-2x savings plausible and 4x unlikely; IterPert's 3-5x relied on STRING/BioGPT priors that are banned here. State
    this prior before running so a null result at 5% vs 20% is read as expected, not as failure.
12. **[LOW] Optional external-validity check (evaluation only, no tuning).** Using Billmann Data S20/S14, test whether
    yeast query-side degree (or predicted degree) of orthologous genes rank-correlates with HAP1 library-side degree
    (hub conservation is a reported result). Supports transferring "choice" rules (change 10) while Billmann's low
    pair-level transfer (OR~1.95, ~7% overlap) argues against ortholog-label features; keep orthology out of
    predictors.

## Things that prior work suggests will NOT work / should be down-weighted
- Uncertainty sampling and model-kernel selection at n<100 (cold start; IterPert, Hacohen).
- Hub-first as an accuracy strategy (GeneDisco; Sutherland): expect it to be good for B, weak for A.
- Importance weighting for selection bias (4.1) may add variance; keep only if it passes the paired bootstrap gate.
- Pair-level conservation (yeast SL -> HAP1 SL) as a feature: ~7% overlap; do not claim pair transfer.
