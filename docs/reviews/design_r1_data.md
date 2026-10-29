# Design review r1 — data handling, labels, leakage (reviewer: data scientist)

Scope read: PROMPT.md, DESIGN.md v1, docs/recon.md, DATA.md, results/data_characterisation/summary.md,
src/hsl/upstream.py, third_party/ysl/{models/view.py, models/baselines.py, models/combiner.py, models/foldin2.py,
features/genes.py, features/pairs.py, data/expression_meta.py, eval/labels.py, eval/splits.py, runner.py}, upstream
YEAST_SLs RESEARCH_TREE/DESIGN (read-only). Extra numbers below come from a throw-away script over DEV-POOL ROWS ONLY
(no test row was read, no model run); the PI should re-derive them in a committed script before relying on them.

## What is sound (no change needed)
- **ScreenView implements the observation model correctly.** It physically holds only the screened rows (E, SP, Y,
  DMF, DMF_sd); targets carry only name + array-column index. For a test target t the predictor sees ε(s, t) for
  screened s via the array-side column `E[:, target_col]` and, for a cell (t, a) with a screened, the reciprocal
  ε(a, t) (`baselines.reciprocal`) — both are "pairs involving a screened gene" and are allowed. No unscreened row,
  no query SMF, no per-target observed mask enters. Models touch only the view plus label-free external stores
  (no `load_bundle` in `models/`; `guarded_call` blocks the Costanzo files while a method runs). Cross-fit fold
  views use the same constructor. `make_view`'s final-kind Split only relaxes the lockbox refusal, nothing else.
- **E_fixed is correctly defined:** A \ D = 1,175 array-only + 533 on-array test genes; none can ever be screened,
  so no E_fixed cell ever has a reciprocal or is observed under any S. Directed cells with their own labels is the
  right choice given κ = 0.14 and matches upstream; self pairs never observed.
- Test split: lockbox preference is legitimate (upstream lockbox = seeded random draw); stratification uses test
  labels only to build the split, which the brief mandates; no stratum id reaches any feature.

## Proposed changes

1. **[HIGH] Remove the genome-coordinate pair features from the primary m1_hgb recipe (SGA linkage artefact).**
   m1_hgb uses `gen_same_chr`, `gen_dist_kb` (combiner PAIRS, group "genome"). Dev rows: SL prevalence by
   query–array distance on the same chromosome is 0.26 (20–50 kb, n = 19 observed), 0.040 (50–100 kb), 0.019
   (100–200 kb), 0.0148 (> 200 kb) vs 0.0114 on different chromosomes; < 20 kb is 100% masked and same-chr cells are
   ~70% missing. This is meiotic linkage in SGA, which has no HAP1 counterpart (the coordinate *source* has an
   analogue; the *signal* does not). Change: primary m1_hgb runs with `params={"drop": ["genome"]}` (supported, no
   code change); upstream recipe with genome kept as a reported sensitivity. stack_lr/stack_lr_wr2 do not use
   genome features (`foldin2.EXT_KERNELS`), so no change there. Add one evaluation sensitivity: AUPRC on E_fixed
   minus same-chromosome cells. Update §1 ("genome coordinates ↔ same") accordingly.

2. **[HIGH] Pre-register column strata of E_fixed; it is dominated by a non-representative array population.**
   69% of E_fixed columns are array-only genes, whose SL rate in dev rows is 0.0084 (median per-array rate 0.0046)
   vs 0.0126 for on-array test-gene columns and 0.0129 for dev-gene columns. Expected E_fixed prevalence ≈ 0.0097,
   not 0.0114. Change: report E_fixed split into E_AO = T × array-only and E_TT = T × on-array T (533 cols; columns
   are query genes, the closest analogue of HAP1 where every library gene is a potential query), with E_unobs(S) as
   already planned. Headline Q2 claims require the same sign on E_fixed, E_TT and E_unobs (pre-registered rule).

3. **[HIGH] Specify the SelectorView exactly and test each channel separately.** §6.1 is right in intent but not
   operational. Define it as a ScreenView-shaped object built by the harness (selector code never reads
   `config/split/test_queries.txt`): (a) only screened rows; (b) columns of T (and, in dev, of D_val) NaN-masked in
   E/SP/Y/DMF/DMF_sd of those rows — mask, do not delete, so upstream code runs unchanged; (c) `array_smf` of T NaN;
   (d) targets = D \ S only. Any predictor used inside a selector (H after n ≥ 30, U, VH) is fit on this object and
   never reuses a fit made for evaluation; evaluation fits write outputs only after a trajectory finishes. Use
   separate seeded RNG streams for selection and evaluation so adding a checkpoint cannot change a trajectory.
   Poisoning tests: poison (a) test rows, (b) test columns of screened rows, (c) test SMF, (d) unscreened dev rows —
   each alone must leave selections bit-identical; positive control: poisoning a visible column of a screened row
   must change some selection.

4. **[MEDIUM] Tighten the symmetric-pair wording and the coverage definition (§3, §6.4).** "The reverse cell (y, x)
   ... is treated as observed too" is ambiguous. State: ε(y, x) is never revealed; (y, x) is excluded from every
   evaluation set and is never a training label; only ε(x, y) (screened direction) is available, as a feature.
   For Q3 coverage, an SL pair is *detected* only if an SL call exists in a screened row; a pair SL only in the
   unscreened direction is "touched but not detected" (with κ = 0.14 this is most of the SL-in-either-direction
   pairs: 5,168 both vs 60,372 one way). Report cell-level coverage (primary) and pair-level detected/touched; do not
   use "SL in ≥ 1 measured direction" as the numerator.

5. **[MEDIUM] Resolve gene identity by current name where it affects the observation model.** Five dev query ORFs
   map (current SGD name) to an array strain with a different original ORF: YFL057C→YFL056C, YHR079C-A→YHR079C-B,
   YIR044C→YIR043C, YJL205C→YJL206C-A (these four arrays are array-only, i.e. in E_fixed), YLR312C-B→YLR313C (itself
   a query). Screening one makes its E_fixed column a "pair involving a screened gene", breaking E_fixed invariance;
   upstream matches by original ORF so the predictor never links them. Fix (pick one, state it): drop those 4 array
   columns from E_fixed / E_TT, or remove the 5 genes from the candidate pool. Also have the E_fixed build script
   report how many E_fixed columns share a current name with another array (27 such array groups exist).

6. **[MEDIUM] Declare upstream-lockbox knowledge and the 23 contaminated test genes.** T ≈ upstream lockbox
   (633/656), whose random-screening AUPRCs (m1_hgb 0.061/0.083/0.102 at 10/20/30%) are already in YEAST_SLs. State
   that these numbers inform no decision (predictor choice, budgets, strategies); the §9 sanity check against them
   runs only in the final, after freeze. Add a pre-registered sensitivity: all endpoints on T minus the 23
   non-lockbox test genes (upstream hyper-parameters were tuned with them as dev rows).

7. **[MEDIUM] Off-array test genes are a yeast artefact for HAP1.** 123 test genes have no array column (no GI
   information at all, no SMF); in HAP1 every query gene is in the genome-wide library and has a parental fitness
   value. Pre-register the on-array test stratum (533) as the estimate used for Q3/HAP1 translation (all-T remains
   primary for Q1/Q2). Keep query-strain SMF excluded (it comes from the unscreened query's own strain — ambiguous
   under the brief) and say so.

8. **[LOW] SMF: acceptable; document it.** Array SMF is a single-mutant phenotype of the array strain from the same
   SGA data release, not a GI; analogue = HAP1 parental-line fitness of every library gene (measured as the control
   arm of each screen, and in DepMap-type screens). Costanzo array SMF pools colony sizes across screens, so test
   screens contribute negligibly; ε is defined relative to SMF, so SMF→SL is partly the multiplicative model's
   scale — HAP1 scoring has the same property. Write this in §1 rather than leave it implicit.

9. **[LOW] Add a feature-provenance table to §1** (store → label-free? GI/text/GO-free? human analogue? used by
   predictor / selector?). Points to record: COXPRESdb strict study-level perturbation filter (good); UniProt only
   Pfam/InterPro counts/Jaccard, no GO fields, InterPro assignment is sequence-based; Michaelis 2023 edge FDR — state
   whether thresholds were calibrated on curated complexes (no GI, acceptable either way); `ppi_bait_observed` is a
   technical coverage flag; Hoepfner (tier 2) never enters any primary predictor or selector.

10. **[LOW] Label/missingness assertions and E_fixed counting.** Assert Y == −1 ⇔ ε or P missing and that
    `ysl.eval.labels.stringent_negative(E, P)` reproduces Y on dev rows. Build the E_fixed column set from structure
    (A \ D, T, item 5) now; report its dev-row-estimated prevalence now and the test-row count/prevalence only in
    the final run (keeps "test rows read only in §8" literal). State that linkage-masked cells (mostly same-chr)
    are excluded from every evaluation set.

11. **[LOW] Seeds and small-n edge cases.** `config/seeds.yaml` lacks `q1_draw_seeds`, `q2_seeds`,
    `inner_split_seeds`, `final_*` referenced in §5–§8: append them before first use. `run_method` defaults
    `split_seed=0` (part of the method RNG): record what is passed. `run_method` bypasses `predict.predict`'s refusal
    of zero-SL screened rows; 16 dev genes have zero SLs and can be selected — add a test that all predictors run
    with zero-SL screened rows at n = 10, and never drop such rows silently.
