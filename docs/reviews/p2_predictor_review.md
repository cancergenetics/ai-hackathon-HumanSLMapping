# Review: P2 frozen predictor P* = `stack_lr_wr2_lb` (commit c1154c7)

Reviewer: independent; I did not write this code. Date: 2026-10-29. Scope: `config/predictor.yaml`, `src/hsl/predictors/*`,
`scripts/30_p2_predictor.py`, `tests/test_p2_predictors.py`, `results/p2/*`, LEDGER `p2_*`, checked against DESIGN v3.1
§1, §3, §8 and §9. Dev data only: no test-query rows were read. New tests are in `tests/test_review_p2_predictor.py`:
7 pass and 1 strict xfail (P-1). I also re-ran `tests/test_p2_predictors.py` and `tests/test_provenance.py`; all pass
except the hygiene test (H-1).

## Verdict: **ACCEPT WITH FIXES**

The improvement over the parent is real and leak-free. I reproduced it exactly from saved predictions, and the frozen
recipe reproduces bit for bit. The fixes are about bookkeeping and wording: provenance, repo hygiene and the freeze
list. The quoted "+4.45%" also needs to be framed as a best-of-10 adaptive dev estimate; it is not a confirmed effect
size. None of the fixes changes P*.

## 1. Leakage / admissibility: PASS (one bookkeeping issue)

- **Structural.** `ScreenView` holds no target rows (targets by name only). `dense_arr` reads only screened rows'
  `v.E`/`v.Y`; `aprior` reads `v.Y`, `v.array_smf`, `gene_matrix(array genes)`. No query SMF, D_val/T rows or
  per-target mask; `target_col` feeds only the `q_*` columns, which P* drops (`lowbudget.py:249`).
- **aprior inputs.** These are the `genes.SCALARS` blocks (sequence, abundance, popgen, alphafold, uniprot counts,
  AP-MS node stats), plus expression mean/sd, coexpression hub-ness, Compara phylo/paralogs, 16 ESM PCs and `has_*`
  flags. Each comes from a store marked `label_free` and `gi_text_go_free` with a human analogue, and each is already
  cleared for m1_hgb. There is no GO, no text mining, no genome data and no GI-derived column.
  `test_pstar_opens_only_admissible_stores` records every file P* opens and finds none outside those stores.
- **P-1 (M, `config/feature_provenance.yaml`, `docs/feature_provenance.md`).** `used_by` does not list
  `stack_lr_wr2_lb` for any store. Before P2, stack_lr_wr2 read only the pair store; P* now reads about 10 gene-scalar
  stores. `test_provenance.py` audits only `stack_lr` and `stack_lr_wr2`, so this gap is invisible. Failure scenario:
  a later edit to the `SCALARS` list could add an inadmissible column to P*, and no test would catch it. Fix: add the
  consumer to the yaml and the rendered doc, and add `stack_lr_wr2_lb` to the parametrised provenance test. My strict
  xfail `test_pstar_listed_as_consumer_in_provenance` will then start to XPASS; remove its marker at that point.
- **dense_arr LOQO is exact** (`lowbudget.py:166-177`: row t removed from column and global-prior sums; existing
  recomputation test passes). Training uses n−1 rows, test n rows, as upstream `loo_additive`.
- **aprior own label: excluded.** For training rows in fold k, the GLM is fit on `fv = v.without_queries(fold k)`
  (`lowbudget.py:274`). New test `test_fold_view_array_prior_ignores_fold_rows`: garbling the fold rows leaves that
  column bit-identical, and garbling a non-fold screened row changes it. The only train/test asymmetry is that training
  rows' priors are fit on about 0.9n rows and test priors on n rows. That affects fitting and calibration only; it
  cannot bias held-out evaluation, because E_fixed_dev cells never enter any GLM or LR.
- **Poisoning: PASS.** The existing n = 60 test passes. My n = 30 test also passes: garbling D_val rows or unscreened
  pool rows gives bit-identical predictions, and garbling one screened row changes them. At cell level
  (`test_one_visible_cell_changes_array_side_features`), flipping one screened cell changes that column's `arr_mild`
  more than 100× any other column's (the others move only through the global prior) and also changes the GLM prior.
- L: `_memo` keys (`lowbudget.py:82-95`) hash only ORF lists and `_ADESIGN` only array names (not data/SMF); a cache
  spanning two bundles (e.g. a poisoning test inside `feature_cache()`) would silently reuse stale features. No
  current path does this. Fix: add the bundle sha (and SMF bytes) to the keys.

## 2. Gate correctness: PASS (with caveats on interpretation)

- Recomputed from saved `.npz` predictions (`E.average_precision`, fresh `p1.query_weights(split,'p2_gate')`,
  `E.WeightedAP`, `E.paired`; missing parent predictions from P1's stack_lr_wr2 files). AUPRCs match runs.csv to
  < 1e-12; bootstrap replicates match `boot/*.npy`. P* vs parent, dev1 n = 160: Δ = +0.00266 [+0.00200, +0.00332],
  +6.84%, 10/10 trajectories positive (= `contrasts.csv`). Re-split (8 perms with predictions): +3.9/+5.5/+2.9% at
  n = 30/160/330; from the R = 10 boot data +4.14% (reported +4.14%). Every verdict recomputed from boot + runs.csv
  matches `verdicts.csv`. Parent AUPRC = P1 stack_lr_wr2 in all 120 cells (shared trajectories genuine).
- Verdict logic (`30_p2_predictor.py:314-325`) matches DESIGN §8; no metric, label or test changed. The n = 660
  "≤ 1% loss" rule and the 10-config cap are **not in DESIGN §8** (stricter, harmless), but line 45 cites "DESIGN §8"
  for the cap and −0.01 is hard-coded (line 325) instead of read from `protocol.yaml` (L).
- Futility stopping (prop, dense_prop at R = 2–3) is legitimate: it can only remove passes; both count toward the 10.
- **G-1 (M, interpretation).** The process was adaptive and the selection is a maximum:
  - Configurations 7–10 were designed after interim looks (`contrasts_interim_R2/R1.csv`) that already showed every
    split, *including the re-split* and n = 660, on the same R0 trajectories, so the re-split was not a blind check.
  - The rule "largest mean gain among passing" is not in DESIGN; it was chosen at P2.
  - The gi300 component failed the gate on its own (+0.39%; CI > 0 on 1/3 splits). Its increment over
    `dense_arr_aprior` is +0.45% on average, with CI > 0 at n = 160 only on dev2 and the re-split, and roughly 0 on
    dev1 and dev3. It enters P* by riding on passing components.

  Effect: the *existence* of a gain is solid. Every dense_arr variant passes, the CIs at n = 160 are about 8 SE from 0,
  and aprior adds a consistent 1–2% over dense_arr. The *magnitude* +4.45% is optimistic by perhaps 0.3–0.5 pp; the
  less-selected `dense_arr_aprior` gives +4.00%. Fix (wording only): report "+4.0–4.5% (best of 10 adaptive dev
  configurations; no independent confirmation)" and state that gi300 is retained because it is non-negative
  everywhere, not because it passed on its own.

## 3. Reproducibility: PASS

- Fresh `harness.run_predictor(RunSpec('dev1', perm0[:30]), 'stack_lr_wr2_lb', yaml params, model_seed from
  p1.random_traj)` produces E_fixed predictions **bit-identical** to the stored gate run (AUPRC 0.0308956; 22.8 s,
  2 threads). `dense_arr_aprior` (run at an earlier git head, 7dc95f9) is also bit-identical when re-run on the
  current code, so I found no code drift.
- The seeds come from `config/seeds.yaml` (`dev_traj_base_seed`, `bootstrap_seed` with crc32 spawn 'p2_gate'). The
  Poisson GLM and spline are deterministic. `config/predictor.yaml` plus `lowbudget.py` at c1154c7 is enough to
  rebuild P*, and the existing yaml-vs-registry test passes.
- **F-1 (M, freeze).** P* constants live in code (`MILD_EPS`, `NEGEPS_PRIOR_M`, `APRIOR_ALPHA`, `APRIOR_KNOTS`) and
  depend on `ysl.features.genes.SCALARS`. The DESIGN §8 freeze-hash list does not include
  `src/hsl/predictors/*.py` or `config/predictor.yaml`. Fix: add both to `config/frozen.yaml` when it is written.
- L: every P2 run has `git_dirty: True` (4 heads); the bit-identical re-runs above mitigate this.

## 4. Robustness: PASS (one pre-existing upstream limitation)

- n = 10 random with off-array targets: finite scores; `array_prior` finite on an all-zero-SL view.
- n = 10 all-zero-SL screened set: **P* raises ZeroDivisionError, and so does the parent** stack_lr_wr2 (and
  stack_lr), in upstream `baselines.shrunk_rate` because the global rate is 0. This is not a P2 regression (pinned by
  `test_pstar_n10_zero_sl_fails_like_parent`). It is practically unreachable for random sets: only 12/2,100 pool genes
  have zero SLs. A selector that favours low-degree genes, or a 9-row fold view at n = 10, could hit it (L; a guard
  would need a wrapper outside `third_party`).
- n ≥ 1,970: no size-specific paths; `dense_loo` holds ~8 float64 n × 3,827 arrays (80 MB each at n = 2,627), so
  ~0.7–1 GB extra peak on the parent's > 6 GB: within the 12 GB budget, but measure before the P1 re-run (L).
- Runtime claim (21–42 s, about the parent's) plausible: my n = 30 fit took 22.8 s.

## 5. Other reasons to temper "+4.45%"

- Gain is on on-array target rows (+5–7% at n = 160); off-array rows gain less, sometimes negative (dev2 n = 30 −1.9%,
  dev3 n = 660 −1.2%). On-array is the HAP1-relevant number, but report it. Absolute gain is small (+0.001–0.0027
  AUPRC, e.g. 0.039 → 0.042 at n = 160); +1.1% at n = 660, not significant on dev3 at 660 (p = 0.14).
- It was validated only on R0 random trajectories. dense_arr and aprior estimate array degree from whichever rows are
  screened; under hub-first selection those rows are degree-enriched. DESIGN applies the same P* to every arm, so Q2
  remains fair, but the P2 gain is not evidence for strategy arms (descriptive note).

## Required fixes (none change P* or its numbers)

1. **H-1 (M):** c1154c7 tracks 1,418 `results/p2/boot/*.npy` files. `test_no_data_or_predictions_tracked` now
   **fails** at HEAD. Either untrack them and add `results/p2/boot/` to `.gitignore`, or make an explicit exception
   in the test and in DESIGN.
2. **P-1 (M):** add `stack_lr_wr2_lb` to `used_by` in the provenance yaml and doc, and to the parametrised provenance
   test.
3. **F-1 (M):** add `src/hsl/predictors/*.py` and `config/predictor.yaml` to the freeze-hash list.
4. **G-1 (M, wording):** in summary.md, predictor.yaml and RESEARCH_TREE, state that this is a best-of-10 adaptive
   selection with re-split and n = 660 visible at the interim looks, give the +4.0% (`dense_arr_aprior`) alongside it,
   and document the ad hoc rules (n = 660 ≤ 1% loss, 10-config cap, max-gain selection), moving −0.01 into
   `protocol.yaml`.
5. Optional (L): data-aware cache keys; peak-RSS measurement at n = 2,627; a zero-global-rate guard in a wrapper.
