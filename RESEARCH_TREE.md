# RESEARCH TREE — Designing human HAP1 SL-mapping screens from the yeast GI network

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-30 evening, milestone M7: project complete)

The one-shot final on T (656 held-out queries) ran 4,019/4,019 units and was independently rescored with 0 mismatches.

Q1 (random screening):
- P\* AUPRC rises from 0.031 (30 queries) to 0.042 (160), 0.063 (660) and 0.092 (whole 2,627 pool), with no plateau.
- n_80 / n_90 = 1,143 / 1,761 queries (44% / 67% of the pool).

Q2 (strategies vs random):
- Adaptive selection in batches of 10 beats random at a matched budget of 160 [C]: H +31%, C +31%, F +26%; Holm p ≤ 0.0014.
- It was not shown to beat 490 random (c1), and "5% matches 20%" is not met. G = 0.55–0.65.
- The random-equivalent budget is 2.4–2.7× [D].
- The gain is pooled-AUPRC, concentrated on SL-rich targets: macro-AUPRC +5–10%, AUROC unchanged.
- A Billmann-like R1 panel is indistinguishable from random.
- A warm start (+27–40%) and batch 110 both keep most of the gain.

Q3 (HAP1):
- With b ≈ 1, reaching 80% of whole-pool performance needs ≈ 3,300 HAP1 lines at random vs ≈ 1,000 adaptively (4 rounds of ≈ 300).
- Billmann's 222 lines ≈ 36–43% of whole-pool performance.
- A static three-feature rule fixed on dev (low fitness, conserved, stable expression) recovered 91–108% of H's gain on T (post-final descriptive, 3 runs).
- Prediction at 20% adaptive reaches the precision of a reverse-direction remeasurement but not of a screen.
- Online CV on adaptively screened lines is badly biased, and the frozen stopping rule did not beat a fixed budget.

Deliverables are complete:
- REPORT.md (final, reviewed; 6 figures);
- results/final, results/q3_final and the independent rescorer;
- reproducible code; DATA.md; all reviews in docs/reviews/.

Remaining ideas are listed in REPORT.md §8. The most valuable is a pre-registered confirmation of the static three-feature rule against H, and an H-vs-C comparison. No compute is running.

**Pre-stated readings** (written before scoring, 2026-10-30 00:45; outcome in brackets):
- (1) If (a) is rejected for no finalist: no confirmed advantage. [Not triggered: (a) rejected for all 3.]
- (2) (a) rejected, c1 not: "beats random at matched budget", with G and n_req as estimates. [APPLIED.]
- (3) "Matches 20%" only under the 3-condition rule; never "approaches". [Not met; honoured.]
- (4) E_unobs and E_TT reported prominently. [Done: REPORT §4.3.]

## Decision log
- 2026-10-30 D11: FREEZE and FINAL START.
  - Freeze: finalists.yaml frozen (bc449a7); frozen.yaml written and committed (9694404; 81 files).
  - First `start` failed in preflight before any lock: hashing the upstream Costanzo raw files tripped the label gate.
    This path is final-only and untested by the dry run.
  - `start` was re-run with HSL_FINAL=1, which opens the gate for hashing only; start reads no labels. Preflight passed:
    81 frozen hashes, 461 upstream checksums, YEAST_SLs clean.
  - Lock 00:22:47 at HEAD 9694404: 4,019 units, 394 jobs.
  - Until scoring finishes, no edits to frozen files; only notes (tree, lessons) are committed.
- 2026-10-29 D10 (pre-freeze, before any T contact). Final cost projected at 14.5–28 h wall, over §12's 10 h. §12's
  'R at the pilot minimum' is undefined here, because the pilot gave no finite R. The closest rule in its spirit:
  - Random arm at R = 10 for n > 660 and R = 5 at the pool max. R_r = 40 stays at n ≤ 660 (all Q2 contrasts).
  - Memory-aware concurrency: ≤ 6 jobs, a free-memory gate, and ≤ 2 jobs at n ≥ 1970.
  - All arms kept.
  - budgets.yaml m1_hgb Q2 key aligned to DESIGN §4 {160, 490, 660} (stale per r3 review).
  - The final-build agent's 19 ambiguity resolutions are accepted (commit 579bac5 report): extended freeze list, pool-max
    R = 5, R0 warm-start continuation, no C-static, reference predictor R = 10, ignored run ledger, n_80 anchored on
    the isotonic fit, Holm-adjusted c1 in the matches rule, and others.
  - Seeds final_bootstrap_seed and final_dryrun_{traj_base,bootstrap}_seed appended.
  - The final starts only after batch6 is committed (a clean tree is required throughout the run) and after the
    independent review of the pipeline code.
- 2026-10-29 D9 (pre-freeze, dev only):
  - Proposed finalists H, C, F_r16 (§8 rule), pending independent P4 review.
  - R_s = 20, R_r = 40. Reason: measured dev V_Q(Δ) ≈ 0.9–3e-5 dominates SE², while σ_s²/R_s + σ_r²/R_r < 3% of SE²
    at these R. The frozen ρ = 0 pilot rule gave no attainable R (V_Q alone > (g/3.5)²), so it set no floor, and Q2
    claims are phrased as bounds per §5. The pilot's ρ_boot = 0.993 figure (R_s ≈ 39) was a sensitivity, not the
    provisional R.
  - Deterministic runs of C and F are not off-array-neutral (argmax under route calibration). They stay descriptive
    with that caveat; strategy code unchanged.
  - M is not a finalist (C wins the C/M slot by the rule; M also never picks off-array genes).
- 2026-10-29 D8: P* frozen after P2 (stack_lr_wr2_lb). Old P* (stack_lr_wr2) results at n ≤ 1640 are kept as the
  reference-predictor sensitivity; every downstream curve (pilot, universe, n ≥ 1970, P3, P4) uses the frozen P*. The
  25% universe level is un-deferred (P2 freed the CPU). New-P* reruns start only after the independent review passes.
- 2026-10-29 D1: Test split drawn preferentially from the upstream lockbox (upstream tuned on non-lockbox rows).
- 2026-10-29 D7: User freed disk (~500 GB+). Disk policy relaxed going forward: prediction vectors may be saved for all
  predictors (useful for baseline bootstraps and independent rescoring); not retroactive (no regeneration of deleted
  baseline vectors). The 10 GB free-disk guard stays. D0's decision to read upstream data in place is kept (it also
  guarantees provenance).
- 2026-10-29 D6: Scheduling: P2 (critical path to P* freeze) has priority; P1 holds n ≥ 820 and universe jobs until P2
  finishes. All jobs export BLAS thread caps (load reached 38–120 on 8 cores from oversubscription). The §12 P1 > 6 h
  trigger (50%-only universe) is deferred, not applied: the overrun came from PI scheduling contention, not P1 cost;
  the 25% level runs later if compute allows (decided before any universe result).
- 2026-10-29 D5: H/F degree model: predictions of the off-array (external-only) submodel are quantile-calibrated to the
  on-array model's distribution (default on). Without it H picked 40–73% off-array genes (base 19.5%): off-array
  genes have more SLs in yeast (52 vs 37, a yeast strain-collection artefact), so the submodel split acted as the
  forbidden off-array score term (DESIGN §1).
- 2026-10-29 Seeds appended by the strategy agent (see config/seeds.yaml tail; dated there).
- 2026-10-29 D4 (corrected after strategy review L-2): H/F batch-cap clusters (`features.clusters`) are computed once on
  all non-test genes incl. D_val (label-free external features only; no label leakage); R1 uses `visible_clusters`,
  which excludes D_val as DESIGN §6.2 says. H prior text in DESIGN §6.2 is superseded by D5.
- 2026-10-29 D3: Rejected `.claude` permission deny rules proposed by a reviewer (user's settings); guard + check script instead.
- 2026-10-29 D2: Evaluation unit = directed cell (t,a) with its own label; AB/BA not merged (κ 0.14, merging would
  change the label). Primary set E_fixed is strategy-invariant so hub-selection cannot shrink the evaluation set.
- 2026-10-29 D0: Do not copy upstream data (11 GB; 38 GB free disk). Read upstream processed data in place; any code
  that could write into YEAST_SLs is vendored/redirected. A stamp file (`.yeast_readonly_stamp`) plus upstream
  `git status` are checked at each milestone to prove YEAST_SLs was not modified.

## Tree
- [x] M0 Setup: repo, PROMPT.md, tree, LESSONS, DATA.md; uv env; vendored ysl (third_party/ysl, path patches only); write guard hsl.guard
- [x] M1 Literature + reconnaissance
  - [x] LITERATURE.md (606 lines, 122 refs). Key: Billmann 2026 = 222 HAP1 query KOs (298 screens), TKOv3, qGI
        |qGI|>0.3 FDR<0.1, 47k negative GIs, <5% of pairs covered; queries chosen for functional diversity, not
        iteratively. HAP1 negative density 1.20% ≈ yeast strong-SL prevalence 1.1% → our label ≈ HAP1 standard negative.
        No prior iterative model-guided query selection for GI maps (closest: IterPert, BATCHIE, GeneDisco, adaptive
        column sampling). Yeast 5%/20% ≈ 164/657 queries; HAP1 equivalents depend on scaling assumption (count vs fraction).
  - [x] Upstream recon (docs/recon.md): arbitrary splits via hsl.upstream.make_view/run_method. S1: 3,283 queries,
        3,827 arrays, 2,652 both, 1,175 array-only, 631 query-only; prevalence 0.0114; top 10% queries hold 34% of SLs;
        AB/BA agreement weak (P(BA SL|AB SL)=0.145, κ 0.14). Runtime (s) n=30/650/1300: stack_lr 1.6/4.6/8.9,
        stack_lr_wr2 14/25/59, m1_hgb 89/139/234; m1_hgb needs prop_folds=5 at n=10; n=2600 >6 GB.
- [x] M2 Fixed test split: 656 test / 2,627 dev pool (config/split/), stratified (KS p=0.999), 96.5% from upstream lockbox
- [x] M3 DESIGN.md + review rounds (3 rounds, 11 reviews)
  - [x] DESIGN v1: primary eval set E_fixed = test rows × arrays never screenable (A \ D), strategy-invariant
  - [x] Review round 1 (7 reviewers) → DESIGN v2: tie-aware AP; dev_bundle w/o T rows; E_fixed_dev w/o T columns;
        P* = stack_lr_wr2; confirmatory family (strategy@160 > random@490; Δ(160) > 0; Holm); gap closure G;
        oracle ceilings as go/no-go; R1 Billmann-like + warm start (S0=220); F (Fisher-D) hybrid, U dropped;
        universe-subsampling decides fraction vs count scaling; genome features dropped (SGA linkage artefact).
  - [x] Review round 2 (external PI, statistics+algorithms, data+process) → v3: exact weighted AP (B=2,000) with one
        replicate array; E_common dropped; "matches" = linear contrast s160 − 0.95·r660 with G* ≥ 0.5 discriminating
        rule; predictor development moved before strategies (freeze P* first); F/H/M fully specified; V dropped;
        warm start vs R1-continuation; batch 110 (= 6 HAP1 rounds); scaling slope b; Bundle.blind + runtime guard.
  - [x] Review round 3 (consistency) → v3.1, all 16 edits accepted. **Design frozen** (commit after 5270a0e).
- [x] P0 infrastructure (8ed970f): eval engine (tie-grouped AP, exact WeightedAP bootstrap, Welch, PAVA n_req, Holm),
      views (Bundle.blind dev bundle + runtime label gate), SelectorView, harness/Trajectory, launcher; 97+5 tests pass.
      E_fixed (final, labels only): 999,102 cells, prevalence 0.00939 (E_AO 0.00806, E_TT 0.01232); 15 identity-excluded
      arrays (rule as written; DESIGN text said 4 — accepted, conservative). E_fixed_dev ≈ 745k cells, prev 0.009–0.010.
      Dev smoke (split 1, random): stack_lr_wr2 AUPRC 0.029 @30, 0.036 @100.
  - [x] Independent P0 code review: safe with fixes; no leakage, AP/stats correct. Fixed M1 (censored CI), L2, L3, disk
        guard (M2). L1 (dev k-means clusters include D_val; label-free external features) accepted as deviation D4.
- [~] P1 dev random curves + power pilot + upstream-protocol sanity + universe subsampling
  - [x] Interim n ≤ 660 (19e7a53; 3 splits × R=10): P* AUPRC (E_fixed_dev, prev ≈ 0.0095) 0.029 @30, 0.041 @160,
        0.049 @330, 0.060 @660 (SD of split means ≈ 0.002–0.005); stack_lr 0.030/0.039/0.043/0.048; B1 0.042 @660.
        P* ≈ stack_lr at n ≤ 100 (−1.2% @30) → target for P2. m1_hgb (dev1) 2–9% above P*. Random gap r660 − r160
        ≈ 0.019 ≫ 5% margin → the "5% vs 20%" comparison is discriminating.
  - [x] Power pilot (dev1, R=10, P*): r160 0.0389 (σ_r 0.0021), r210 0.0411, g = 0.0021. Frozen provisional rule
        (ρ=0) gives V_Q(Δ) ≫ (g/3.5)² → "bounds only" provisionally; shared-weight correlation (ρ_boot 0.994) implies
        R_s ≈ 21 suffices. Per DESIGN, V_Q(Δ) is measured directly in P4 and R fixed at freeze.
  - [x] Upstream-protocol sanity: harness reproduces upstream dev split 0 on its T-free subset (P* 0.0774 vs recorded
        0.0775; no T rows touched — upstream dev splits exclude the lockbox). dev1 is a hard query fold (upstream's
        own predictions restricted to dev1 rows: 0.0673; dev2 0.0864; dev3 0.0750) → absolute dev AUPRC varies ±10%
        across folds; paired contrasts unaffected. E_fixed_dev @660 (dev1): P* 0.0544, stack_lr 0.0450, B1 0.0417.
- [x] Strategy selectors (b71a4df, 767563c): R1, H, H-static, F, F-static, C, M (+ diagnostic C-static), oracles in
      src/hsl/oracles.py (dev-only, import-guarded). Toy controls pass (planted hub/cluster, null; F r=1 ≡ H order);
      per-channel poisoning + broken twins pass; selection ≤ 1 s/round (F 0.6–0.9 s). 208 fast tests pass.
      Kernel approximations in selectors: AP-MS neighbourhood over pool genes, ESM via 16 standardised components,
      pair-store kernels unavailable to selectors (SelectorView unchanged).
  - [x] Independent strategy code review (docs/reviews/strategies_code_review.md): safe with fixes — math verified by
        brute force, hygiene clean; F and C nearly never picked off-array genes (inverse of D5), R1 SMF imputation,
        cluster cap at large batches, M λ too large. Fixes (ad45deb): route-stratified calibration (off-array genes
        ranked among themselves, mapped onto the on-array score distribution) for H/F/C; R1 random SMF tertile for
        off-array; cap max(2, ⌈b/k⌉); M λ = 0.3·mean diag. All strategies within ±5 pp of the 19.5% off-array base on
        dev1; 267 tests pass. PI checked the one edited reviewer test: reference updated to the new spec, not weakened.
- [x] P2 predictor development (n ≤ 330; gate per DESIGN §8) → freeze P*
  - [!] P* frozen = stack_lr_wr2_lb (stack_lr_wr2 + array-side dense ε degree + EB array-degree prior + w_gi=n/(n+300)):
    mean rel. gain +4.45% at n∈{30,160,330} (3/3 splits CI>0 at 160; re-split +4.1%; +1.1% at 660); 10 configs (cap),
    4 passed. Oracle decomposition: low-n loss is degree estimation (oracle query degree +92% at n=30, inadmissible).
    Commit c1154c7; results/p2/summary.md. Review ACCEPT WITH FIXES (fixes applied): read the gain as +4.0–4.5%, the best of
    10 adaptive dev configs with no independent confirmation (+0.001–0.003 AUPRC absolute); leakage, gate and
    reproduction verified independently.
  - [~] Re-run power pilot + universe subsampling (50% and 25%) + random curves with frozen P* (DESIGN §9)
    - [x] Random arm, frozen P*: dev1–3 × R=10, n ≤ 660 (12 budgets), 0 failures (results/p1/curves.csv).
    - [x] Power pilot (frozen P*, dev1): r160 0.0416, r210 0.0433, r490 0.0509, r660 0.0553; g = 0.0017 [0.0009, 0.0025].
      Under the frozen ρ = 0 rule, Q2 claims are bounds only (V_Q(Δ) 4.5e-6 ≫ (g/3.5)²); the shared-weight sensitivity
      gives R_s ≈ 35–39. V_Q(Δ) is measured directly in P4 and R is fixed at the freeze.
    - [x] Universe subsampling with frozen P* (50%, 25%) and n 820–1640 (batch5, commit 3418dfd, 0 failures).
      Frozen P* vs pre-P2 P*: +4.9% @30, +5.7% @160, +1.1% @660, +0.5% @1310 (AUPRC 0.0722 at 1310).
      Provisional slope (50/25% only): b(n_80) = 1.08 [0.72, 1.58], b(n_90) = 1.13 [0.79, 1.54]. n_80/N ≈ 0.35–0.37,
      so fraction scaling; reported as a bound because external stores are not restricted to the sub-universe.
    - [x] n ≥ 1970 + pool max (batch6), commit 1431962. P* at the pool max = 0.0800 (stack_lr 0.0559, B1 0.0432).
      - [!] n_80 / n_90 (dev, P*): 810 / 1279 (0.39 / 0.61 of the pool). The baselines saturate early (stack_lr n_80 ≈ 384–506, B1 ≈ 110).
      - [!] Universe slope over 9 draws: b(n_80) = 1.06 [0.85, 1.28], b(n_90) = 1.06 [0.89, 1.24]. That is fraction
        scaling (n_80/N ≈ 0.35–0.38), reported as a bound.
- [x] P3 oracle ceilings → GO (C, M kept). Commit fb18bb0, results/p3/summary.md.
  - [!] Greedy oracle (dev1, in-sample): s160 0.0698 vs r160 0.0416, G = 2.05 [1.55, 2.75]. Upper bound only, because it
    selects on the scored cells; a held-out-half version is queued.
  - [!] Hub-first oracle (true query degree): pooled G = 0.49 [0.27, 0.72]; split-dependent (dev1 0.02, dev2 0.69,
    dev3 0.64). Its top picks are 60% off-array (base rate 19.5%, D5 yeast artefact), so on-array and E_TT strata are requested.
  - [x] §11 end-to-end controls: shuffled labels → 1.05× prevalence; back door → 1.000. Seed p3_e2e_control_seed appended.
- [x] P4 strategies on dev (commit 01ec9ad; results/p4/summary.md, shortlist_proposal.md)
  - [!] Adaptive selection gives large gains at n = 160 (mean over 3 splits): H +32.7% (G 0.84 [0.63, 1.05]),
    C +29.9% (G 0.77), F_r16 +28.0% (G 0.72), M +28.5%. Static versions gain much less: H-static +13.4%, F-static −6.5%.
    R1 (Billmann-like diverse panel) ≈ random (+0.7%). Same sign on E_AO, E_VV and on-array rows; held at n = 660 (+23–28%).
  - [!] s160 ≈ r490 on dev (c1 is borderline); 'matches 20%' (vs 0.95·r660) is not met on any split.
  - [x] Held-out greedy oracle G = 0.61 [0.21, 0.98] (attainable ceiling); the in-sample G = 2.05 is not attainable.
  - [!] H beats the true-degree hub oracle on dev1 and dev2 despite degree Spearman 0.02–0.21. Mechanism to be checked in review.
  - [x] Scenarios (H, C):
    - warm start (220 R1): +25% vs R1-continuation at 440/660;
    - batch 110: −5% at 160 (still +26% vs random); batch 50: −1 to −4%;
    - feasible pool: −6 to −9%;
    - minimal features: ±2%.
  - [x] Independent P4 review: ACCEPT WITH FIXES (fixes applied, commit 1938d80).
    - No leakage: selections are bit-identical under 4 poisonings, and the positive controls change them.
    - The statistics recompute exactly.
    - [!] Mechanism: adaptive arms screen moderately SL-rich genes (1.5–1.65× random's degree). AUPRC tracks the SLs
      revealed, and the gain comes through each on-array target's own column. The true-degree oracle overshoots to
      extreme hubs. Both oracles are relabelled as references, not ceilings.
    - The R deviation is declared in DESIGN.
    - The gain holds with stack_lr_wr2 and with m1_hgb on stored selections.
  - [x] run_final.py build + dry run on dev (commit 579bac5): 179/179 units; refusals and crash-resume tested; 17 tests.
  - [x] Independent review of the final pipeline: READY WITH FIXES. All fixes applied (43af996, e7d9a69); 321 tests pass.
    The projection is 10.2–13.6 h at an 18 GB gate.
  - [x] Freeze (9694404)
  - [x] Final on T (one shot): 4,019/4,019 units, 0 failures (00:22–12:29), scored at 13:39; frozen check ok
    (commit 8c1818a).
    - [!] Confirmatory: (a) Δ(160) > 0 is REJECTED for H, C and F_r16 (Holm p 0.0003 / 0.0001 / 0.0014).
      - s160 = 0.0555 / 0.0555 / 0.0534 vs r160 = 0.0423, i.e. +31% / +31% / +26%.
      - c1 (s160 > r490 = 0.0576) is NOT rejected: Δ −0.002 / −0.002 / −0.004.
      - Per pre-stated reading (2): the claim is "beats random at matched budget"; c1 is not met.
    - [!] G = 0.65 [0.41, 0.90] (H), 0.65 [0.43, 0.87] (C), 0.55 [0.29, 0.80] (F), with G* = 0.85.
      'Matches 20%' is NOT allowed (s160 − 0.95·r660 < 0). The headline sign rule holds for all 3 (P* and m1_hgb;
      all, AO, TT).
    - [!] n_req(160) (descriptive): 434 [317, 583] (H), 433 (C), 384 (F).
      - That is 2.4–2.7× the random budget (ratio CIs ≈ 1.7–3.6×), about 13% of the query set.
      - It is NOT '3×' or '≈ r490' (c1 not shown).
      - On E_TT, s160 − r490 < 0 with CIs excluding 0 for C and F.
    - [!] R1 (Billmann-like) ≈ random at every n. H-static gives +12–22%; F-static ≈ random (+1%).
    - [!] Warm start (220 R1, then adaptive): +27–35% vs R1-continuation at 440 and +38–40% at 660.
    - [!] Batch size: H at batch 110 gives +33% vs random at 160 (batch 10: +31%), batch 50 +42%. But C at batch 110
      loses 6–7% vs batch 10 at 160 and 330 (CIs exclude 0). Read: 'modest or no loss', not 'no loss'.
    - [!] Q1 (T): P* AUPRC is 0.0917 at the full pool (2,627). n_80 = 1,143 [970, 1,311] (0.44 of the pool);
      n_90 = 1,761 [1,571, 1,899] (0.67). The baselines saturate early (B1 n_80 = 119).
    - Sensitivities agree:
      - lockbox: H +32%;
      - no same-chromosome: H +33%;
      - two-way bootstrap: CI excludes 0;
      - on-array rows: H +40%;
      - off-array rows: +2–6% (yeast artefact rows).
    - [!] E_unobs (the brief's literal set; unpaired; prevalence differs by arm):
      - n = 160: H 0.0521 vs R0 0.0457.
      - n = 660: H = R0 = F = 0.0695, while C is +12%.
      - Strategies observe 1.63× (160) and 1.43× (660) random's SLs, so the unobserved remainder is sparser.
      - Do not claim 'better prediction of unobserved pairs' from lift alone.
    - [!] Where the gain is (n = 160):
      - pooled AUPRC +26–31%, but per-query macro-AUPRC only +5–10%;
      - AUROC unchanged;
      - the gain is concentrated in the top SL-count tertile (+22–24%); F tertile 1 is −10%, n.s.
    - Dev → T: no winner's curse on rel Δ(160) (dev 30.2% vs T 29.5%). G fell 0.84 → 0.65 (H) because T's random
      curve is steeper. The finalist ranking is not tested (H ≈ C > F at 160, F > H > C at 660).
    - Do not quote the saturating-fit parameters (c at its bound). Deterministic arms have R = 1, so their CIs omit
      trajectory variance.
    - Review: docs/reviews/final_results_review.md (claims list).
  - [x] Q3 on the final (b3be29b; results/q3_final):
    - [!] HAP1 lines for 80% of whole-pool performance: random ≈ 3,321 [2,570, 4,374] vs H ≈ 999 (37 rounds) or H b110 ≈ 1,037
      (4 rounds). The b = 0 floor is 1,211. The denominator lever is ±27%.
    - [!] Static three-feature rule (dev-fixed coefficients at bde1622; T evaluation code written post-scoring): 91–108% of H's
      gain on T (dev 72–81%); fitness-only 77–97%. Post-final descriptive, 3 runs. Confirmation is needed.
    - [!] Prediction: on-array precision at reference recall is 0.224 (H@660) vs 0.178 for the reciprocal measurement; top-10 has
      ≈ 1 true SL per gene (10–14× random).
    - [-] Online stopping: LBO-CV is inflated up to +143% under H. The frozen rule's |log2 err| is 0.23 vs 0.12 for a fixed
      n = 490 (hindsight). Negative result.
  - [x] REPORT.md final (940f64f); reviews: docs/reviews/report_review.md (READY WITH FIXES, all applied).
    - [x] Independent rescoring (no hsl import): 0 mismatches.
      - 12,738 rows match exactly, including all 4,019 unit APs and 6,299 strata APs.
      - 154 Monte-Carlo rows are within tolerance.
      - All 6 Holm decisions are identical.
  - [x] Independent rescorer (analysis/rescore, no hsl import; commit f72da0d), 44 tests:
    - The dry run matches exactly: 179 unit APs, 378 strata APs and all point contrasts.
    - CIs and p-values agree within Monte-Carlo error, and all decisions agree.
    - A negative control with 4 planted errors flagged all 4.
    - Two DESIGN readings are reported both ways (G-defined interval; n_X anchor).
- [x] M4 Q1: performance vs. fraction screened (random)
- [x] M5 Q2: strategic iterative selection (batches of 10) vs. random at matched budgets
- [x] M6 Q3: translation to HAP1 (query numbers, choice of query lines, role of prediction)
- [x] M7 REPORT.md, figures, final verification
