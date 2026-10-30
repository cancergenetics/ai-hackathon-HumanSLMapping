# RESEARCH TREE — Designing human HAP1 SL-mapping screens from the yeast GI network

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-29, M0 setup)

Project started. PROMPT.md saved verbatim. Upstream method (YEAST_SLs, commit 9b5dce9, read-only) is a kernel-ridge
fold-in + stacked combiner that predicts stringent negative GIs (ε < −0.12, P < 0.05) for unscreened yeast queries
from the Costanzo 2016 NxN matrix (lockbox AUPRC 0.083 at 20% screened, prevalence 0.011; 10/20/30% → 0.061/0.083/0.102).
It already reads Costanzo S1/CellMap bundles and tier-1 external features (coexpression, AP-MS, ESM-2, paralogs, etc.).

Next: literature review (LITERATURE.md) and technical reconnaissance of the upstream code (how to call it with an
arbitrary screened set, runtime vs. number screened, data characterisation), in parallel; then the fixed test split,
then DESIGN.md and up to three review rounds.

## Decision log
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
  - [~] Freeze (frozen.yaml) → Final on T (one shot)
- [ ] M4 Q1: performance vs. fraction screened (random)
- [ ] M5 Q2: strategic iterative selection (batches of 10) vs. random at matched budgets
- [ ] M6 Q3: translation to HAP1 (query numbers, choice of query lines, role of prediction)
- [ ] M7 REPORT.md, figures, final verification
