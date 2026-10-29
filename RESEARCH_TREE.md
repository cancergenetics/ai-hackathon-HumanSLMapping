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
- 2026-10-29 D1: Test split drawn preferentially from the upstream lockbox (upstream tuned on non-lockbox rows).
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
  - [~] Independent P0 code review
- [~] P1 dev random curves + power pilot + upstream-protocol sanity + universe subsampling
- [~] Strategy selectors implementation (R1, H, F, C, M, static variants) + toy controls — no dev evaluation before P* freeze
- [ ] M4 Q1: performance vs. fraction screened (random)
- [ ] M5 Q2: strategic iterative selection (batches of 10) vs. random at matched budgets
- [ ] M6 Q3: translation to HAP1 (query numbers, choice of query lines, role of prediction)
- [ ] M7 REPORT.md, figures, final verification
