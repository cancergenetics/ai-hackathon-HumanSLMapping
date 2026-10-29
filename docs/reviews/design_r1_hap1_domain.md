# Design review r1 — HAP1 domain / actionability lens

Reviewed: PROMPT.md, DESIGN.md v1, LITERATURE.md, data_characterisation/summary.md. Web check: Horizon states a HAP1 KO
bank of ">7,500 cell lines" (lines, not distinct genes; vendor claim) and does not sell highly essential or multi-copy
genes. All items stay within Q1-Q3; none alters label, split, metric or primary endpoints.

Verdict: yeast machinery (§3-§9) is sound. Q3 (§6.4) describes the yeast->HAP1 translation but does not yet compute the
HAP1 realities (who can be a query, round time, feature availability). Changes below make Q3 quantitative.

## Priority 1 (change before dev starts)

1. **Add a "Billmann-like" static comparator B0 to Q2.** Round-0 selection stratified on SMF and expression bins,
   spread over coexpression/AP-MS clusters (LITERATURE §6 lists it; DESIGN §6.2 dropped it). Why: the consortium's
   real counterfactual is the "functionally diverse, range of expression and fitness" heuristic behind the 222 lines,
   not uniform random. Comparator only, not one of the <=3 frozen strategies.

2. **HAP1-feasible candidate pool as a pre-registered sensitivity.** Define D_feas = D minus the sickest yeast queries
   (bottom SMF decile, ~ the 10-15% essential/near-essential share in HAP1: 1,524/~11,000); run frozen strategies and
   random on it. Why: hub-seeking by low SMF selects exactly the genes whose HAP1 KOs are unavailable (not sold),
   slow-growing (noisy qGI, poor QC) or only hypomorphs (13/222 Billmann lines were partial LOF). Report the AUPRC lost;
   if H/VH win only via the sickest 10%, the recommendation must say so. Add query-SMF tertile as a reporting stratum.

3. **Report cost in HAP1 units: lines, screens, sequential rounds.** Billmann: 298 screens/222 lines = 1.34 (replicates).
   Add a batch-size sweep at fixed n (batch 10/30/60) and a "3-round" schedule (static batch of ~n/2 from features, then
   2-3 adaptive batches). Why: batch 10 is a fine simulation unit but not what a consortium runs. A round = line
   procurement/making (Horizon in-stock ~1-2 wk, made-to-order ~12 wk; vendor blog, unverified) + clonal QC
   (Sanger, ploidy/chr15) + TKOv3 screen at >=200x, 10-12 doublings (~3-4 wk) + sequencing + qGI scoring (loess/SVD need a
   screen panel): plausibly 2-3 months, so 16 serial rounds for n=164 (>2 yr) unless several labs run in parallel.
   Decision output: "rounds to reach 90% of batch-10 performance". Label times as assumptions. Cheap (same code).

4. **Make Q3 scaling testable.** Add the universe-subsample experiment from LITERATURE §6 (absent from DESIGN §6.4):
   subsample genes (50%, 25%) on both sides, run the random Q1 curve at fixed n and at fixed fraction, report whether
   AUPRC tracks n or f. Caveat to write: HAP1 has ~3x more queryable and ~4.7x more library genes than yeast NxN, so
   this only approximates HAP1. Fix the denominator now: HAP1 query universe = expressed AND viable AND available
   (Horizon >7,500 lines; expressed nonessential ~9,500); quote both, i.e. 5% ~ 375-475 and 20% ~ 1,500-1,900 lines
   (not 500/2,000 off ~11,000 expressed genes, which includes essentials that cannot be queries).

5. **Restate what the HAP1 forecast covers.** Yeast NxN has no essential array genes; ~14% of HAP1 library genes
   (1,524) are essential, with interactions largely masked (library negative density 1.42%). State in §6.4/§11 that
   all AUPRC/coverage numbers apply to the nonessential-library x nonessential-query submatrix (~86% of library genes);
   essential-library columns are an unmodelled part of the HAP1 network. Compute coverage only on that submatrix.

## Priority 2 (Q3 outputs / sensitivities)

6. **Feature-availability ablations for the selectors (not just a list).**
   - SMF: continuous for every yeast gene; HAP1 fitness is ~0 for most nonessentials (mild-defect LFC is noisy). Add a
     "HAP1-like SMF" arm (quantised essential-like/defect/neutral + noise) and an "SMF removed" arm for static H/V/VH.
   - AP-MS: Michaelis covers ~3,900 yeast proteins; BioPlex is 293T (14.6k proteins), HCT116 5.5k IPs, OpenCell 1,310
     baits; none is HAP1. Add arms masking AP-MS to ~60% bait coverage (BioPlex-like), ~10% (OpenCell-like), and removed.
   - Coexpression/ESM/paralogs/gnomAD/PaxDb are fine; use HAP1 RNA-seq (Billmann Data S2) alongside COXPRESdb/ARCHS4.
   Output: selector gain vs random under full / HAP1-like / minimal features = which datasets labs must obtain.

7. **Stricter-label sensitivity (ε < -0.35, P < 0.05)** for Q1 curves and one Q2 comparison: analogue of HAP1 "extreme
   SL" (0.14%); recommended in LITERATURE §6, omitted in DESIGN. Primary label unchanged. In §11 state that absolute AUPRC
   does not transfer (yeast precision ~0.54/recall ~0.39 vs HAP1 0.61/0.37, plus HAP1 clonal/arm/chr15/mito
   artefacts); the transferable results are relative (efficiency n_req/n, fraction of random-20% performance).

8. **Noise-robustness arm for screened rows.** HAP1 calls at the standard threshold carry ~40% false positives. Degrade
   screened (never test) rows (ε noise or P resampling calibrated to HAP1 precision/recall) and check strategy rankings
   hold; hub-count-driven selectors are the likeliest to be fragile.

9. **Retrospective HAP1 check of the hub-proxy rule (evaluation only, after freeze).** Express the static selector as a
   species-neutral rank rule (expression x fitness defect x constraint x ...). Post-freeze, compute Spearman between the
   human-feature proxy and HAP1 library-side negative degree (Billmann Data S14, ~17.8k genes) and top-decile hub
   enrichment: the only direct evidence the yeast-derived rule ports. Data S14 is GI data: it must never enter selection,
   training or tuning; PI decides if this evaluation-only use is acceptable, else flag transfer as untested.

10. **Q3 prediction-contribution framing.** A HAP1 screen has precision 0.61; a predicted list at 0.2-0.5 does not
    substitute for a screen. Frame prediction as (a) prioritising validation and (b) imputing the ~95% unscreened pairs;
    report direct yield (true SLs per screened line) beside extra SLs recovered at precision >= 0.2/0.5 for the same n,
    the lift over B1 (so it is not just "hubs have more SLs"), and results as a fraction of the noise ceiling (AB/BA
    agreement, κ 0.14) to avoid over-selling.

## Priority 3 (presentation / cheap)

11. **Warm-start scenario.** HAP1 starts from 222 non-random lines, not S_0 = ∅. One Q2 scenario: S_0 = B0 with n_0 ~ 200,
    then adaptive batches; report marginal gain per extra round over B0 alone (cheap given item 1).

12. **Fix the shape of the recommendation now:** (a) n and rounds under both scalings with CIs; (b) eligibility filter
    (expressed, not core-essential, line available/derivable, no extreme fitness defect); (c) ranked selection rule
    (item 9-checked); (d) stopping rule (§6.4 estimator) with its correlation to true AUPRC; (e) explicit "cannot be
    answered from yeast": essential-library columns and pair-level transfer (only ~7% of HAP1 GIs among conserved pairs
    are yeast GIs; OR ~1.95), so advice concerns strategy, not specific pairs.

13. **Minor:** E_fixed holds array-only genes (1,175) plus test genes; report their SMF/degree distributions vs all genes.
    Always give yeast counts (164/657) next to HAP1-scaled counts.

## Not needed
Split, label, observation model, symmetric-pair handling, AUPRC definition and leakage tests are appropriate for HAP1.
