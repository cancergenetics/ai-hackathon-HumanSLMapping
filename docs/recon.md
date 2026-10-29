# Upstream API reconnaissance (YEAST_SLs @ 9b5dce9)

Scope: how to run the upstream methods on our own splits, what they need, how the S1 data are represented, and what
the code could write. Numbers: `results/data_characterisation/summary.json`, `results/runtime_benchmark/runtime.json`.

## 1. Using the upstream package from our code

- `third_party/ysl/` = vendored copy of `YEAST_SLs/src/ysl` (PROVENANCE.md: per-file SHA-256; only path roots
  patched). Data and upstream config (seeds, lockbox, frozen hashes, recipes) are read **in place** from the root in
  `config/paths.yaml` (`yeast_sls_root: ../YEAST_SLs`).
- `import hsl` first in every entry point: installs `hsl.guard` (audit-hook write guard on YEAST_SLs) and exports
  `YSL_UPSTREAM_ROOT`/`YSL_OUTPUT_ROOT` for the vendored package. `ysl._paths` also imports `hsl`, so any `ysl`
  import is guarded.
- Env: `uv sync` (Python 3.12; `uv.lock` seeded from upstream's lock — identical pins for 3.12).

## 2. Running ANY registered method on an ARBITRARY screened set and target set

Upstream entry points and why we do not use them as-is:

| entry point | problem for us |
|---|---|
| `ysl.runner.run(...)` | fixed upstream dev/final splits only; writes results/ledger/predictions; calls `git` and `pytest` with cwd = upstream root (blocked by our guard) |
| `ysl.predict.costanzo_view(ds, screened, targets)` | loads the **blinded** bundle and **refuses lockbox ORFs** as screened or target (our test set is 96% lockbox) |
| `ysl.predict.predict(view, ...)` | fine, but restricts `model_seed`/`split_seed` to upstream `config/seeds.yaml` |

**Clean path (implemented in `src/hsl/upstream.py`):**

```python
import hsl                                     # guard + paths (must be first)
from hsl import upstream as U
b = U.load_s1()                                # ysl.data.dense.load_bundle("s1", unblind=True); SHA-verified, ~1 s, ~0.3 GB
view = U.make_view(b, screened, targets)       # Split("final", seed, train=screened, test=targets, lockbox, "s1", N)
                                               #   -> ScreenView.from_split(b, split)  (same constructor the runner uses)
S, info = U.run_method(view, "stack_lr_wr2", model_seed=17, n_threads=3)   # float32 [n_targets x 3827]
```

- `from_split` allows lockbox ORFs only for kind `"final"` with an unblinded bundle → that is why `make_view` builds
  a final-kind Split. Equivalent lower level: `ScreenView._construct(_BundleSource(b), screened, targets)`.
- `run_method` = what `predict.predict` does (registered defaults ← params, `MethodContext`, `set_thread_limits`,
  `runner.guarded_call`), minus the upstream seed whitelist. Rows of `S` follow `view.target_orfs` (sorted);
  columns follow `b.array_strains` (sorted strain IDs; `b.array_orfs` gives the ORFs).
- **Targets are arbitrary**: any query ORFs not in the screened set (test genes, unscreened dev genes for active
  selection, or both). Only target NAMES and their array column index enter the view; target rows never do.
- **Our own evaluation**: labels of targets from `b.Y[b.gene_rows(targets)]` (int8: 1/0/−1 = unobserved), outside the
  method call. Upstream's `ysl.eval.evaluator.Evaluator(b, split)` also accepts our final-kind split (pooled AUPRC
  over every observed cell of the target rows, macro AUPRC, ORF-cluster bootstrap, strata on/off-array & reciprocal).
- **User-matrix mode** (`ysl.predict.user_view(matrix_tsv, targets)`) also works on arbitrary screens without the
  bundle (long table query_orf/array_orf/epsilon/pvalue); arrays then = arrays of the given rows, sorted by ORF unless
  `array_strain` given (column order changes the `tied` feature's random start: r = 0.99995 vs bundle mode).

## 3. Method requirements

Methods are `fit_predict(view, ctx) -> [n_targets x n_arrays]` in `ysl.models.registry` (32 production methods;
`registry.production_methods()`). The three of interest:

| method | inputs (beyond screened rows) | model | external stores |
|---|---|---|---|
| `stack_lr` | array rate (LOO), query propensity = x's array-side rate, on-array flag, reciprocal (avail, −ε, label, −signed log P), 3 kNN scores (query-side col/mixed, array-side), array SMF | LR on cross-fitted cells | none (GI only) |
| `stack_lr_wr2` | stack_lr + within-row scores: `mk` multi-kernel fold-in (GI-profile kernel + AP-MS, AP-MS-edge, coexpression, ESM kernels; rank 100), `tied` low-rank (rank 10), refined kNN, 2-hop pp/pg, `ext_off` external-only fold-in (off-array x) | LR | pair store (query×target kernels) |
| `m1_hgb` | stack_lr_wr2 features + predicted query propensity (nested cross-fitted sub-model) + 15 external pair features (target×array) | HistGradientBoosting, K-fold additive stats, test features averaged over 10 matched views | pair store + gene-level tables |

`*_t2` variants add the Hoepfner 2014 co-fitness store (tier 2). External-feature universe = 4,426 current ORF
names (S1 query ∪ array); targets outside it are refused by `predict` (not by `run_method`) — all S1 queries are inside.

**Cross-fitting.** `ctx.n_folds = 10`: training-query folds `rng("folds").permutation(n) % 10`
(rng = SeedSequence(model_seed, split_seed, crc32(tag))). Training cells of fold k get features from
`view.without_queries(fold k)`, built by the same constructor as the test view. m1_hgb's propensity sub-model is
nested (`prop_folds = 10` inside each outer fold). Cost ≈ (K+1) feature builds (stack_lr_wr2) and ≈ K×K for m1_hgb.

**Small screened sets (matters for iterative screening from 10 queries).**
- No explicit minimum n in the methods. `predict.predict` refuses screened rows with zero SLs (not an issue: 10
  random S1 queries carry ~400 SLs). With n < K = 10 some folds are empty → `without_queries([])` → untested; keep
  n ≥ 10 or pass `n_folds=min(10, n)`.
- Benchmark (runtime only, §6): **all three methods run and return finite scores at n = 10, 30, 100**.
  NumPy "Degrees of freedom <= 0" warnings at n ≤ 100 (variance over ≤ 1 value in a fold view) are harmless.
- Degeneracies to expect at small n (from the code; accuracy not measured here): `masked_pearson(min_n=10)` sets
  query-side similarities over < 10 shared screened coordinates to 0, so the query-side kNN features are ~0 when
  fold views hold ≤ 9 on-array screened queries; `k = 50` neighbours and the `mk` rank 100 are clipped to n;
  propensity sub-model falls back to one pooled model when < 10 on-array training queries (`MIN_SPLIT_TRAIN`); the
  HGB `min_samples_leaf = 1000` with ~30k training cells at n = 10 gives very shallow trees. **All hyper-parameters
  were tuned upstream at n ≈ 657 (20%) and are untested at 10–300.**

## 4. What the S1 bundle holds (`data/processed/costanzo2016_bundle_s1.npz`, built from Data File S1 SGA_NxN.txt)

- Strain level: E (ε), SP (sign(ε)·−log10 P), Y (label int8: 1/0/−1), DMF, DMF_sd, each [3,283 query strains ×
  3,827 array strains]; strain metadata (ID, ORF, canonical flag, suppressor flag, temperature); array metadata
  (strain, ORF, current ORF name, array SMF); `all_query_orfs`.
- **S1: exactly one strain per query ORF** (all `_sn`, 30 °C, no suppressor alleles) → strain = gene. Upstream's
  gene-level rule (`canonical_query_strains`: non-suppressor > `sn` > most arrays > smallest ID) only matters for
  TheCellMap (3,608 strains / 3,597 ORFs). Arrays: 3,827 DMA strains, one per ORF.
- Genes: 3,283 query ORFs, 3,827 array ORFs, 2,652 both (80.8% of queries on-array), 1,175 array-only, 631
  query-only (off-array). 89.9% of the grid observed (missingness = linkage mask etc.); prevalence 1.142%.
- SMF: array SMF for all 3,827 arrays (pair-wise file's 30 °C value, 4 from the strain table); query SMF is not
  exposed (upstream: x's SMF only via its array strain).
- **Essential genes / ExN / NxE / ExE / DAmP: not present** in the bundle or on disk unzipped. The raw zips do
  contain them (S1 zip: SGA_ExN.txt, SGA_ExE.txt, SGA_DAmP.txt, hsp90_data.txt; CellMap zip: SGA_ExN_NxE.txt,
  SGA_ExE.txt, SGA_DAmP.txt); upstream extracted only NxN + strain table. Using them would require reading inside
  the zip (read-only is possible with `zipfile`), but none of the upstream methods/stores cover TS alleles.

## 5. Reciprocal (AB/BA) measurements

- Stored as two independent cells: (A as query, B as array) = `Y[row(A), col(B)]`, (B as query, A as array) =
  `Y[row(B), col(A)]`. 89.5% of measured unordered pairs among the 2,652 on-array queries are measured both ways
  (52.8% of all observed cells belong to such a pair); self pairs are never observed.
- Agreement is weak: P(BA SL | AB SL) = 0.145, Jaccard 0.079, κ 0.136, ε Pearson 0.13 (5,168 pairs SL both ways,
  60,372 one way). So a reciprocal cell is a noisy second measurement, not a copy.
- In the view: a target x's **array-side column** `view.E[:, view.target_col[i]]` gives ε(t, x) for every screened t
  (the only GI information about x); for a target cell (x, a) with a's ORF screened, the reciprocal is
  `view.Y[row(a), col(x)]` → `stack_lr`'s `recip_*` features (≈ 14% of test cells at 20% screened upstream).
- Upstream evaluation scores **every observed cell of a held-out query row**, including cells whose reciprocal (a
  screened, x on array) is known and used as a feature, and cells whose array is another held-out gene. It does not
  merge AB/BA. → Our DESIGN must state how symmetric pairs are handled (PROMPT: "pairs involving a screened gene are
  observed regardless of side"): e.g. exclude from AUPRC the test cells (x, a) with a screened, or evaluate at
  unordered-pair level with a stated AB/BA combination rule.

## 6. Runtime (fit+predict only; no accuracy; 656 test-query targets; 3 threads; Apple M2, 8 cores)

See `results/runtime_benchmark/runtime.md` for the full table (filled by `scripts/02_runtime_benchmark.py`).
RUNTIME_TABLE_PLACEHOLDER

Scaling is with n_screened × n_arrays (training cells) plus n_targets × n_arrays at predict time; m1_hgb also
re-builds test features over 10 matched views.

## 7. What the code writes, and where

- Method `fit_predict`: **writes nothing**. Feature stores are opened read-only (`np.load(mmap_mode="r")`, parquet
  reads); the bundle is opened read-only. Inside `guarded_call` the Costanzo raw/processed/bundle/prediction files
  cannot even be opened (upstream data guard, still active in the vendored copy).
- Upstream writers that would target YEAST_SLs — all refused by `hsl.guard`, none needed by us:
  `dense.write_bundle`/`Bundle.save`, `features.pairs.save_features`, `features.cofit` store builder,
  `eval.splits.write_lockbox`; the upstream `scripts/*` (not vendored).
- `ysl.runner`: results/ledger/LOCK/predictions → redirected to `results/ysl_runner/` and `data/ysl_predictions/`
  in humanSLmapping (patch); its `git`/`pytest` subprocesses run with cwd = YEAST_SLs → blocked by the guard. Do not
  use `runner.run`/`run_final_batch`.
- `ysl.predict.write_table/write_report`: user-given paths only.
- Python bytecode: vendored sources → `third_party/ysl/**/__pycache__` in our tree (gitignored); nothing is
  imported from YEAST_SLs. Upstream `tests/`, `uv sync`, scripts: never run.
