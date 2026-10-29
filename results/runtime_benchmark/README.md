# Runtime benchmark — summary (runtime and memory only; NO accuracy computed or looked at)

Raw: `runtime.json`/`runtime.md` (main grid, 2026-10-29 08:19–08:45) and `runtime_m1_hgb_propfolds5.{json,md}`
(small-n workaround). Script: `scripts/02_runtime_benchmark.py`. Screened = nested random draws from the dev pool
(`runtime_benchmark_draw_seed: 1`); targets = 656 test queries (names only); 3,827 arrays; 3 BLAS/OMP threads;
Apple M2 (8 cores, 24 GB); one fresh process per cell; watchdog cap 6 GB RSS. fit+predict excludes bundle load (~1 s)
and view construction (< 1 s).

fit+predict wall time in s (peak RSS GB):

| n screened | stack_lr | stack_lr_wr2 | m1_hgb |
|---|---|---|---|
| 10 | 1.6 (0.6) | 12.2 (1.2) | **fails** (IndexError); with `prop_folds=5`: 86.1 (2.2) |
| 30 | 1.6 (0.7) | 14.0 (1.2) | 89.3 (2.1); with `prop_folds=5`: 84.6 (2.0) |
| 100 | 1.8 (0.8) | 13.1 (1.8) | 117.9 (2.0) |
| 300 | 2.8 (0.9) | 15.8 (1.8) | 111.9 (2.4) |
| 650 | 4.6 (1.5) | 25.1 (2.4) | 139.4 (2.9) |
| 1300 | 8.9 (2.4) | 58.8 (3.5) | 234.4 (4.9) |
| 2600 | 22.5 (4.1) | **killed** at 6.0 GB after 285 s | **killed** at 6.2 GB after 333 s |

Notes
- **m1_hgb at n = 10 fails** inside upstream code: each outer cross-fit fold view has 9 screened queries, the nested
  propensity sub-model uses `prop_folds = 10` → one empty inner fold → `ysl.features.genes.similarity("apms", [], …)`
  builds `np.ix_` from an empty float array → IndexError. Parameter-only workaround (no code patch): `prop_folds=5`
  (must be ≤ the smallest outer-fold view size, ⌊0.9 n⌋). Any n < ~11 is affected; n = 30 runs with the defaults.
- stack_lr and stack_lr_wr2 run at n = 10 (NumPy "Degrees of freedom <= 0" warnings at n ≤ 100 are harmless).
- **Memory, not time, limits n = 2600** for stack_lr_wr2 and m1_hgb under the 6 GB cap (the cross-fitted training
  design matrix over ~9 M observed screened cells). Options if needed: raise the cap for those cells, or subsample
  training cells (method change → must be justified/validated). Upstream ran ≤ 5.2 GB at n ≈ 657–1,078.
- Time is roughly flat below n ≈ 300 (fixed costs: 656 targets × 3,827 arrays prediction grid, feature-store reads,
  10 matched test views for m1_hgb) and grows ~linearly above. For iterative screening in batches of 10 from 10 to
  ~650 queries (~64 refits per trajectory): stack_lr ≈ 3 min, stack_lr_wr2 ≈ 15–20 min, m1_hgb ≈ 2 h per trajectory
  (at 3 threads, 656 targets; more targets — e.g. scoring the unscreened dev pool for selection — add proportionally
  to the predict part).
