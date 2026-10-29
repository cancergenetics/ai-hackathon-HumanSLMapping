# Independent review: final evaluation pipeline (pre-lock)

Reviewer: independent agent, 2026-10-29. Scope: commits 579bac5, abbe89d and 7e7803b, plus the working-tree edits in
flight at about 23:35 (memory gate, `calibrate`, per-predictor scaling). I read `src/hsl/final.py`,
`scripts/{run_final,score_final,make_frozen,launch}.py`, `config/finalists.yaml`, the reused `hsl.{harness,views,eval,
q2dev}` code and `results/final_dryrun/`. T was not touched, HSL_FINAL was never set, nothing was written to
`results/final/`.

New test file: `tests/test_review_final.py` (6 pass, 1 xfail documenting the ci defect; about 8 s). The existing
`tests/test_run_final.py` passes (17/17).

**Verdict: READY WITH FIXES.** No leakage and no statistical error in the confirmatory path. The fixes are about
running a 15–35 h one-shot job safely in a multi-agent repo (§4), one CI defect and pre-registration bookkeeping.
Make F1–F6 before `make_frozen.py`; every code fix touches a file on the freeze list.

## 1. Test isolation — PASS
- **Selection never sees T.** Every arm (R0, R1, finalists, static, det, warm start, batch) goes through
  `harness.Trajectory` → `SelectorView.build(b, S, final_split)` (`views.py:438-482`):
  - screened rows only, with S ⊂ pool asserted;
  - T columns and T SMF masked, including identity-matched current names (`cols_matching`);
  - visible features exclude T genes;
  - pool = D: `final_split` gives 2,627 genes and 0 rule-(i) removals (YLR312C-B and YLR313C are both in D);
    checked label-free by `test_final_split_pool_excludes_T_and_efixed_columns_exclude_D`.
- **Warm start:** S0 = first 220 genes of R1 trajectory r, all from the pool.
- **Selection precedes labels in every process.** It runs once per (arm, r) at the start of its job, before any
  `run_unit`, so no E_fixed/T label array is in memory yet (`_ESETS` empty), and it is saved before evaluation.
- **No hidden data access.** Strategies, predictors and features open no bundle, split or label file. The predictor
  view holds training rows only (`ScreenView._construct` + `_BundleSource.rows(train)`); it may see T columns of
  screened rows (allowed, DESIGN §3).
- **E_fixed on T** = T × (A \ D) − rule (ii) − rule (iii) (`views.e_fixed`, frozen, P0-reviewed). `score()` asserts
  one `eset_digest` for all rows (strategy-invariant). 633/656 T genes are upstream-lockbox (the 23-gene sensitivity).
- **Labels reach only the evaluation path:** `_labels` job → `run/labels.npz`, plus the inline AP in `run_unit`.
  Optional hardening: selections could run on `dev_bundle()` (no T rows at all), making isolation physical.

## 2. Pre-registration fidelity
- **PASS — confirmatory family** (`score_final.py:135-156`): per finalist (c1) s160 > r490 and (a) s160 > r160;
  Holm over 6 at α = 0.025, one-sided; Welch SE² = V_Q + σ̂²_s/R_s + σ̂²_r/R_r with the §6.4 ν (V_Q known)
  (`eval.welch`). Recomputed from scratch: match to 1e-9.
- **PASS — G, G*, "matches."** G: joint bootstrap (shared W^b, trajectories resampled within arm). "Defined" = 2.5%
  quantile of (r660 − r160) > 0. G* uses the point means. The matches claim needs all three: s − 0.95·r660 at
  p ≤ 0.025, c1 rejected after Holm (D10), G* ≥ 0.5.
- **PASS with defect — n_req, n_80, n_90:** PAVA on the dense random grid, log-n interpolation, left endpoint,
  censoring coded 0/∞, joint bootstrap; n_80 anchored on the PAVA value at the pool maximum (D10). CI defect: F2.
- **PASS — Bonferroni bands:** per arm over its checkpoints, two-sided 95%.
- **PASS (L) — headline rule** (`score_final.py:338-347`): Δ(160) same sign on E_fixed, E_AO, E_TT with P* and
  m1_hgb, as DESIGN §3 says (E_unobs is not in it). The code also requires the shared sign to be positive, so a
  uniformly negative result is labelled "qualified". Acceptable; state it.
- **PASS — sensitivities:** lockbox-only rows (633); no same-chromosome cells (current systematic names; unresolved
  counted in the manifest); two-way bootstrap (B = 500, w_t·v_a, array columns that are T genes reuse that gene's
  query weight) on the 6 family contrasts.
- **PASS — finalist knobs:** H defaults, C `normalize: false` (= default), F r = 16; quantile calibration, top-k 20,
  cluster cap 2. These match `shortlist_proposal.md` §1, the P4 runs (`q2dev.parse`) and the strategy code/yaml
  (unchanged since ad45deb); `check_global_knobs` asserts them at run time. F-static at r = 16 ran on dev only in the
  dry run (descriptive).
- **PASS — seeds:**
  - trajectory r of every arm = `SeedSequence(final_traj_base_seed = 920297).spawn(.)[r]`; selection = child 0,
    model seed = child 1, shared across arms (CRN, §5); disjoint from dev (410293) and from the dry run;
  - bootstrap: named children of the new `final_bootstrap_seed`.
- **PASS — predictor:** pstar = `stack_lr_wr2_lb {extra: [dense_arr, aprior], gi_n0: 300}`, ref = `stack_lr_wr2 {}`,
  both read from `config/predictor.yaml` (confirmed in dry-run rows); m1_hgb `drop: [genome]`.
- **ISSUE M (F6) — D10 not declared in DESIGN.** D10 sets the random arm to R = 10 at n > 660 and R = 5 at the pool
  maximum (§7.1 says R_random at all budgets), and replaces the §4/§12 concurrency rule ("n ≥ 1970 alone") with the
  memory gate. Both are only in RESEARCH_TREE; the revision history lists (i) and (ii) only. The projection is still
  12.5–35 h, so §12's 10 h trigger still fires as written. Also log whether the §12 trims (m1_hgb Q1 {30, 1310},
  batch 50) were evaluated.

## 3. Statistical correctness — PASS (recomputed independently)
Recomputed from the saved dry-run predictions and `labels.npz` with `hsl.eval` primitives plus numpy/sklearn:

| quantity | recomputed | matches |
|---|---|---|
| unit AP | = sklearn `average_precision_score` (1e-9) | `auprc_inline` |
| Welch H (a), from-scratch SE²/ν/p | Δ 0.00705, SE 0.00401, ν 209.4, p 0.0401 | `confirmatory.csv` |
| G_H | 0.458 [0.049, 0.834]; G* 0.821 | `summary.json` |
| n_80 (pstar) | 359.3 | `q1_nx.csv` |

- **Bootstrap:** one `BootEngine.W` (`E.bootstrap_weights` on `labels.npz:query_strata`) drives every arm and
  stratum; the strata are T's 11 split strata (0–10), checked in a test.
- **ISSUE M (F2) — `q2dev.ci` drops +∞ replicates** (`q2dev.py:439-441`, finite values only). It feeds the n_req CI
  (`q2dev.py:509`) and the G CI (`q2dev.py:460`). Dry run: H@490 is 97% right-censored yet its CI is [624, 974];
  H@660 shows [—, —]. Fix: `E.percentile_ci` on non-NaN values (keeps ±∞, P0 fix M1). Q1 n_80/n_90 already use it.
- **L — pool-max resampling:** `joint_resample` draws 40 indices but only 5 trajectories exist at the pool maximum;
  in ~0.5% of replicates ((35/40)^40) none is drawn, giving NaN that is silently dropped (`score_final.py:265-272`).
  Resample within the available trajectories per column, or report the count.
- **L — saturating fit** pools 18-point (r < 10) and 12-point (n ≤ 660) trajectories (`score_final.py:274-298`), so
  the asymptote mixes extrapolations. Fit full-grid trajectories only, or report separately (descriptive).
- **L — secondary E_fixed metrics missing** from the score step: macro AUPRC, AUROC, p@10/50 and the SL-count-tertile
  stratum (DESIGN §3). Recomputable from saved predictions + `labels.npz`; add now or label post-hoc-from-saved.

## 4. Robustness of the one-shot run
- **ISSUE H (F1) — a dirty tree or moved HEAD makes every later job refuse.** Each job requires a clean tree
  (untracked files count) and HEAD == locked commit (`guard_final`, `final.py:475-490`; `views.final_bundle`,
  `views.py:135`). Over 15–35 h, any edit to a tracked file, any new non-ignored file or any commit makes every later
  job exit 2, and the launcher drains the queue as failures. Triggers include:
  - a RESEARCH_TREE commit; `results/LEDGER.jsonl` appends by other jobs;
  - `run_final.py cost` / `calibrate` (write tracked `results/final_dryrun/*`); `score --partial` (F3);
  - launcher stdout redirected into `logs/` (only `*.log` is ignored);
  - my `tests/test_review_final.py` until it is committed.

  Results cannot be corrupted (jobs refuse before loading data), but time is lost, and after a commit `resume`
  refuses until HEAD is moved back. **Prevent:**
  1. Before `start`, commit or ignore all outstanding paths (`logs/`, `results/p1/*`, `results/LEDGER.jsonl`,
     batch6).
  2. Freeze all repo writes and commits by every agent (incl. PI tree edits, `cost`, `calibrate`) until `score` ends.
  3. Send launcher output to `results/final/run/`.
  4. Document recovery: `git checkout <locked sha>` (detached), clean the tree, `resume`.
  5. Optional: stop launch at the first exit-2 refusal.
  6. Optional: run from a separate clone; its venv must point at the clone's `src` (`.venv/_hsl.pth` pins the main
     tree), and `data/cache/clusters_*.npz` must be copied in.
- **ISSUE M (F3) — `score --partial` in final mode** (`score_final.py:621-693`) writes tracked outputs to
  `results/final/` and appends to tracked `results/LEDGER.jsonl`: it dirties the tree mid-run and is an unplanned
  look at T results. Refuse it in final mode while units are incomplete (unless an explicit abandon flag), or write
  partial outputs only under `run/`.
- **ISSUE M (F4) — score does not verify frozen code.** `score()` checks lock and plan hash but not `check_frozen()`,
  and records no git HEAD for the scoring code (`score_final.py`, `q2dev.py`, `eval.py` are all frozen). Call
  `check_frozen`, store `views.git_state()` in `summary.json`, refuse on mismatch in final mode.
- **ISSUE M (F5) — no protection against two launchers.** `resume --launch` while the first launcher is alive runs the
  same jobs twice: duplicate units, concurrent `np.savez` to the same `run_id` file (risk of a corrupt prediction
  file), and an `ensure_selection` race (`final.py:527-561`, no file lock). Add a launcher PID lock in `run/` checked
  by `resume`, and an O_EXCL claim per selection file.
- **L — truncated JSONL line:** `p1.read_rows` (`p1.py:332-338`) raises on a partial trailing line (process killed
  mid-append), breaking `done_uids`, `progress`, `resume` and `score`. Skip and log a malformed last line.
- **L — pool-max memory/time unmeasured:** m1_hgb at n = 2627 has never run; RSS hint ~8.3 GB vs kill cap ~12.4 GB,
  and a unit that hits the cap fails on every resume. Run `calibrate` at the dev pool maximum (~2,100) before the
  freeze. Projected ~1,500 s is under the 7,200 s timeout.
- **L — memory gate** now counts inactive pages (which include compressible anonymous memory). RSS caps still
  protect; watch for swap at max_big = 2.
- **PASS:** disk (~35 GB needed: ~4 MB E_fixed per unit, ~8 MB more with E_unobs; 572 GB free, 10 GB guard);
  predictions saved before scoring; a unit is done only when row and file exist; selections written atomically and
  never recomputed (dry-run crash test); score needs no HSL_FINAL and refuses until all units are done; the lock is
  O_EXCL (if `start` dies between lock and `plan.json`, recovery is manual — document it).

## 5. Other (L)
- Copy `run/LEDGER_final.jsonl` (crash/resume history, §10; ignored) to `results/final/` at score time.
- Record the sha256 of `data/cache/clusters_*.npz` in the lock: it feeds the H/F/C cluster cap but is unhashed.
- Add `check_readonly.sh` and the upstream checksum check to the pre-start checklist.
- `features.cluster_universe()` loads `dev_bundle()` in final selection jobs (harmless; memory and time).
- Confirm the §7.4 Q3 analyses need only saved artefacts (predictions, `labels.npz` with both directions of T×T
  cells, selections), so T never reopens after the lock.

## Required before the lock
- **F1** (H, procedural): repo freeze — clean tree, no agent writes or commits during the run, launcher output under
  `run/`, written recovery steps.
- **F2** (M): n_req and G CIs via `E.percentile_ci`, keeping ±∞.
- **F3** (M): refuse `score --partial` in final mode, or confine its output to `run/`.
- **F4** (M): `check_frozen` plus recorded git HEAD in `score`.
- **F5** (M): launcher PID lock and per-selection O_EXCL.
- **F6** (M): declare D10 (random R by n, memory gate) in the DESIGN revision history; log the §12 trigger outcomes.

The L items are recommended; none changes a confirmatory number.
