# P0 code review (independent)

Reviewer: independent agent; I did not write the code under review. Scope: `src/hsl/{eval,views,guard,harness,features}.py`,
`config/{protocol,budgets,seeds,feature_provenance}.yaml`, `tests/`, against DESIGN v3.1 §3, §6.1, §6.4, §8, §11.
Date: 2026-10-29. I used dev data only. No test-query rows or labels were read, and `HSL_FINAL` was never set.

## Verdict: **safe with fixes**

I found no leakage path and no error in the AP metric or the confirmatory statistics. Two medium issues need
fixing before P1 work that depends on them:
- **M1** must be fixed before any n_req, n_80 or n_90 confidence interval is computed.
- **M2** must be handled before the P1 launch.

The low issues are robustness and spec-conformance gaps. None of them biases the P0 results produced so far.

## What I checked independently (tests/test_review_p0.py, 24 tests; scratch scripts not committed)

| area | check (my own reference code) | result |
|---|---|---|
| pooled AP | brute-force distinct-threshold AP and sklearn, 300 random heavy-tie cases, unit and fractional weights | equal to 1e-12 |
| adversarial AP | all tied → prevalence; all positive; single cell; zero positives → NaN; empty; ±0.0 form one tie; NaN scores and −1 labels **refused** (missing cells are never counted as negatives) | pass |
| WeightedAP | equals AP on duplicated cells (integer query weights incl. 0), fractional weights, single query, all ties, positives with weight 0 → NaN | pass |
| bootstrap | stratum totals exact; multinomial variance; ReplicateEngine uses the **same w^b for every arm and budget**; chunked = per-replicate; missing (r, n) stay NaN | pass |
| Welch / paired | SE² = V_Q + σ²_s/R_s + σ²_r/R_r, ν formula, one-sided p with margin, paired Var(d)/R + V_Q (DESIGN §6.4) recomputed by hand | pass |
| Holm | step-down reference, 500 random families incl. ties | pass |
| PAVA + n_req | brute-force fine-grid inversion of an isotonic, log-n-interpolated curve, 400 random non-monotone curves; left endpoint of flat segments; left/right censoring 0/+∞; unsorted grid refused | pass (see note N1) |
| p@k, Bonferroni | expectation over random tie-breaks; t quantile 1 − α/2K | pass |
| dev_bundle | every row-aligned array (E, SP, Y, DMF, DMF_sd, 5 strain arrays) has 2,627 rows; no T strain; blinded_lockbox = T | pass |
| label gate | blocks np.load (incl. mmap), np.memmap, os.open, Path.read_bytes, pq.read_table, pd.read_parquet, pq.ParquetFile, upstream.load_s1, and label_access from another module | pass (see I1) |
| E_fixed_dev (dev1–3, resplit) | no T columns and no D_cand columns, by name and current name; rows are D_val only; labels are each cell's own; **independent recount of the cell set matches exactly** | pass |
| final E_fixed columns | no dev-pool array by name or current name (checked on names only) | pass |
| SelectorView (dev1, n = 60) | masked columns = cols_matching(T ∪ D_val); E/SP NaN and Y −1 there; visible columns unchanged; SMF NaN for excluded genes; public attribute set is exactly the whitelist; features/clusters refuse D_val and T genes; screen_view has DMF all-NaN | pass |
| prediction poisoning, stack_lr, n = 60, dev1 | full score matrix bit-identical when D_val rows are garbled (E, SP, Y, **DMF, DMF_sd and the missingness pattern**), and when unscreened pool rows are garbled; changes when 5 screened rows are garbled (positive control); harness scores = full matrix at E_fixed cells (row and column alignment) | pass |
| trajectories | adding checkpoints plus stack_lr evaluations leaves randomised-greedy selections unchanged; random selector gives the same nested prefixes for batch sizes 10 and 50 and for max_n 100 and 300; top-20 picks are uniform (χ²) | pass |

Fast suite: `.venv/bin/python -m pytest -m "not slow"`: **117 passed, 5 skipped (HSL_FINAL-only), 4 xfailed**. The 4
xfails are my strict pins of issues M1, L1, L2 and L3. Each will turn into a failure once its fix lands, and the
xfail marker must then be removed.

## Issues

### M1 (medium): `percentile_ci` returns NaN instead of +∞ for right-censored n_req replicates
- `src/hsl/eval.py:417-423`. `np.quantile` uses linear interpolation. When both neighbours of the quantile
  position are +inf, it computes inf − inf = NaN.
- Failing example: `E.percentile_ci(np.r_[np.full(950, 100.), np.full(50, np.inf)])` returns `(100.0, nan)`. The
  correct result is `(100.0, inf)`. This happens whenever at least about 2.5% of the replicates are right-censored,
  which is the typical case for n_req of a strong strategy or n_90 near the pool maximum.
- Why it matters: DESIGN §6.4 asks for censored replicates to be "reported (one-sided bounds)". A NaN upper bound
  reads as "undefined" when the truth is "unbounded, i.e. > n_max". The same applies to G when its denominator is 0.
- Fix: use an order-statistic quantile, `np.quantile(x, [a, 1 - a], method="inverted_cdf")` (or `"lower"`/`"higher"`),
  or map ±inf to finite sentinels beyond the data range before interpolating and map them back afterwards. Also
  return the censored fractions. Pinned by `test_percentile_ci_keeps_censored_bounds`.

### M2 (medium, operational): saved predictions are about 3 MB per run against 38 GiB free disk
- `src/hsl/harness.py:222-243`. Each run saves float32 scores for every E_fixed_dev cell (dev1: 748k cells =
  2.99 MB, see `results/predictions/dev1/*.npz`). Final E_fixed will be larger.
- Rough P1 count: 4 predictors × ~17 grid budgets × 10 trajectories × 3 splits ≈ 2,000 runs ≈ 6 GB. Universe
  subsampling, P2 re-runs on the shared trajectories, and P4 add several times that. `df` shows 38 GiB free (96%
  used).
- Fix (PI decision): save predictions only where §10/§11 need them (the final run and the designated
  independent-rescoring dev milestone). Elsewhere keep `save=False` and store metrics plus replicate APs. Add a
  free-space check to the launcher. `np.savez_compressed` gains little on float32 logits.

### L1 (low, spec deviation): k-means cluster universe includes D_val in dev
- `src/hsl/features.py:77-81`. `cluster_universe()` = (Q ∪ A) \ T for every split. DESIGN §6.2 (R1, and the H/F
  cluster cap) says "dev: also excluding D_val".
- The kernel is label-free, so this does not leak labels. It does make dev clusters transductive in D_val, unlike
  the final, where T is excluded. Example: `np.isin(dev_split(1).targets, cluster_universe()).sum() == 526`
  (526 of 526 D_val genes).
- Fix: `cluster_universe(split)` = non-excluded genes, with the cache key per split digest. Pinned by
  `test_cluster_universe_excludes_dval_in_dev`.

### L2 (low): the last batch overshoots `max_n`
- `src/hsl/harness.py:339-343`. `k = min(b_size, len(cand))` is not clipped to `max_n − len(S)`. Example: batch 50
  with max_n 655 selects 700 genes.
- Checkpoints still use prefixes, so metrics are unaffected. However, one extra selection round runs (it may include
  a stack_lr fit), and `log["selected"]` is longer than max_n. Any discovery-axis count (§7.2) that uses
  `selected` would then be overstated.
- Fix: `k = min(b_size, len(cand), max_n - len(S))`, or keep the round but store `selected[:max_n]`. Pinned by
  `test_trajectory_does_not_overshoot_max_n`.

### L3 (low): checkpoints above the pool size are dropped silently
- `src/hsl/harness.py:326-327`. `cps` keeps only `c <= max_n = min(max_n, len(pool))`. A dev Q1 run with grid point
  2627 (dev pool 2,101) silently records no point at the dev pool maximum. DESIGN §6.4 uses that point as the
  n_80/n_90 reference: "Dev uses the dev pool maximum in place of 2,627".
- Fix: map any checkpoint `> len(pool)` to `len(pool)` and log the mapping, or raise. Pinned by
  `test_pool_max_checkpoint_not_silently_dropped`.

### L4 (low): E-set cache keyed by `id(bundle)`
- `src/hsl/harness.py:89`. Once a bundle is garbage-collected, a new bundle (for example a poisoned or shuffled copy
  in a test or control) can reuse its id and receive the old bundle's labels and cells.
- No reuse happened in my 200-try probe (the bundles are large), and production uses the cached singleton
  `dev_bundle()`, so the risk is low. The failure would be silent, though.
- Fix: key on `(b.meta["bundle_sha256"], b.meta.get("hsl_view"), split.digest(), name, …)` plus a content digest
  for modified bundles, or use a `weakref.WeakKeyDictionary`.

### L5 (low): no assertion that score-matrix rows align with `split.targets`
- `src/hsl/harness.py:185`. `ESet.scores` indexes `Smat[rows, cols]`. It assumes `view.target_orfs` (np.unique,
  sorted) equals `split.targets` (sorted in `_make_split`). That holds today: my check confirmed alignment on dev1.
  A `Split` built by hand with unsorted targets (as the test fixture does) would silently misalign.
- Fix: `assert np.array_equal(view.target_orfs, split.targets) and np.array_equal(view.array_orfs, b.array_orfs)`.

### L6 (low): non-atomic shared cell-index write
- `src/hsl/harness.py:235-236`. `if not cp.exists(): np.savez(cp, …)`. With 4 concurrent jobs on the same split,
  two writers can race on `_cells/<digest>.npz`. A reader may also see a partial file.
- Fix: write to a temporary file and `os.replace` it into place. Apply the same to the per-run file.

### L7 (low): inconsistent NaN handling for missing trajectories
- `src/hsl/eval.py:278-284` vs `:322`. `welch`/`paired` return NaN if any trajectory estimate is NaN, and `R`
  counts NaN trajectories. `joint_resample` uses `nanmean` and silently drops them.
- Fix: decide on one policy, preferably to raise on missing cells in confirmatory contrasts, and report `R_eff`.

### L8 (low): process-global thread limits change after the first predictor call
- `hsl.upstream.run_method` → `ysl.runner.set_thread_limits` calls `threadpool_limits(n)` without a context
  manager (`third_party/ysl/runner.py:170`). After the first checkpoint, OpenMP is limited to 3 threads for the rest
  of the process (observed: `[('openmp', 3)]`).
- Selection results could in principle differ by floating-point rounding between runs with and without
  checkpoints. My test saw no change.
- Fix: set the thread limits once at `Trajectory.run()` start, before any selection.

### Info / notes
- **I1** `guard.label_access` (`src/hsl/guard.py:274-276`) checks the caller's `__name__`. `exec(code, {"__name__":
  "hsl.views"})` opens the gate. This is an honour-system limit, and the docstring and grep test should say so.
- **I2** `views.efixed_label_summary()` (`src/hsl/views.py:393-405`) reads T rows without `HSL_FINAL`. It returns
  aggregate label counts only, as DESIGN §3 allows. The grep test whitelists only `scripts/04_eset_manifest.py`.
  Consider requiring an explicit flag so it cannot be called again after P0.
- **I3** `welch`/`paired` hardcode 0.975 for `ci95` (`eval.py:291,311`) instead of reading `stats.ci_level`. This
  conflicts with the "code never hardcodes" rule.
- **I4** DESIGN §3 gives YLR312C-B→YLR313C as an example of a D/T split pair. Both genes are in D, and the D/T
  shared current-name set is empty, so rule (i) is a no-op for the final (`final pool_removed = []`). Correct
  behaviour, but the DESIGN example is wrong.
- **I5** `randomised_greedy` breaks ties at the top-k boundary by lower index (`harness.py:276`). For exactly tied
  gains it is not uniform over all tied candidates. This is irrelevant for continuous gains; document it for
  integer gains.
- **N1** Not a defect: when two adjacent PAVA blocks are equal up to 1 ulp (an unpooled value 0.12 next to a pooled
  mean 0.12000000000000001), "left endpoint" depends on float rounding. This only arises when the target equals a
  fitted value exactly.
- **Upstream `YEAST_SLs/results/`** is not under the label gate. I listed its file structure without reading
  values: it holds aggregate lockbox metrics (declared in DESIGN §1). Per-cell upstream predictions under
  `data/processed/predictions/` are gated.

## Verified correct (no action)
- Tie-grouped AP; WeightedAP sparse segment algebra; shared-W replicate array; Welch, paired and Holm formulas;
  isotonic n_req.
- dev_bundle blinding; the guard on all audited read channels; E_fixed_dev / E_fixed column rules incl. identity
  rules (ii) and (iii); SelectorView masking with native codes, SMF masking and the whitelist.
- Predictions for D_val that do not depend on D_val or unscreened rows (stack_lr, n = 60).
- Selection and evaluation streams that are independent (the model seed is derived from the evaluation child of the
  trajectory seed, so it is the same across arms and gives CRN); random selector as nested prefixes; batch remainder
  taken in greedy order.
- Manifest content (git, seeds, params, checksums, runtime, RSS, threads) and one ledger line per run, including
  failures.
- No dtype overflow (int64 counts, float64 cumulative sums, int32 indices well within range).
