# Design review r1 — AI-agent failure modes (DESIGN.md v1)

Lens: where PI + sub-agents can silently go wrong. MUST = before dev code runs; SHOULD = before the first headline
number; COULD = cheap extra. Read: PROMPT, DESIGN, LESSONS, tree, recon.md, hsl/guard.py, seeds, split manifest, upstream
LESSONS + r1 failure-mode review (generic items there transfer; only design-specific ones below). v1 is already good on
leakage-by-construction; gaps follow.
## A. Test-set contact and metric drift
1. **Dev code must not hold test rows (MUST).** `U.load_s1(unblind=True)` returns all 3,283 rows incl. T; SelectorView
   and poisoning tests protect only code that uses them. Add `hsl.views.dev_bundle(b)` (T rows deleted, not masked)
   as the only data entry for dev scripts; `final_bundle` needs `HSL_FINAL=1` + a clean frozen commit. Grep test: dev
   scripts never call `load_s1`.
2. **Inner dev evaluation must also exclude T columns (MUST).** §7 `D_val x (A \ D_cand)` includes on-array T genes as
   columns; (d,t) labels are the reciprocal of test cells (t,d), so tuning on them touches T labels weakly. Use
   `D_val x (A \ D_cand \ T)`.
3. **No test scoring before freeze (MUST).** §9 "n=657 random within upstream lockbox range" implies scoring T in dev.
   Run that check on inner splits during dev; on T only in the final report, labelled post-hoc. Counts/prevalence of
   E_fixed (labels, not scores) may go in the manifest.
4. **Frozen evaluation module + hashes (MUST).** Before dev write PI-owned `hsl/eval.py` (AUPRC, E_fixed/E_unobs, bootstrap,
   n_req, non-inferiority). `config/frozen.yaml` holds sha256 of it, `config/split/*`, `seeds.yaml`, `budgets.yaml`,
   `third_party/ysl/PROVENANCE.md`, upstream commit file; `tests/test_frozen.py` fails on mismatch. Sub-agents own
   `hsl/strategies|predictors/*`, `scripts/exp_*`, agent tests; a failing PI test is reported, never edited around.
5. **AUPRC tie rule is self-contradictory (MUST).** §3: "sklearn definition, ties broken by stable sort"; sklearn groups
   tied scores at one threshold, a stable sort does not. At small n scores tie heavily (B1 is constant down a column;
   kNN features ~0 at n<=100), so a stable sort makes AP depend on cell order (optimistic if positives sort first). Use
   tie-grouped AP. Tests: all-equal -> AP = prevalence; shuffled cell order -> identical AP; hand-computed literals;
   agreement with sklearn on 1e5 cells; a stable-sort/ascending mutant must be caught.
6. **E_fixed drops the reciprocal feature and is not the brief's literal set (SHOULD).** No array in E_fixed is ever
   screened, so reciprocal eps(a,t) (a major upstream feature) is never available: the headline is "fold-in only".
   The brief's "pairs involving a screened gene are observed" is E_unobs(S). Both are pre-registered, but the risk is
   choosing the favourable one later. Rule now: a Q2 claim is headline only if sign and non-inferiority verdict agree
   on both sets; report E_unobs prevalence per strategy (hub columns removed) and compare by lift/ratio-to-random.
7. **One-shot final with lock (MUST).** `run_final.py` refuses unless `config/finalists.yaml` (<=3 strategies,
   predictors, budgets, comparator, endpoints) is committed, tree clean, `results/final/.lock` absent. Lock holds
   the git hash; crash re-runs record `attempt: n` + reason; predictions saved before scoring. Q1-on-T may run
   once predictor recipes are frozen, not before (it informs nothing, but avoids "peeking" suspicion).
## B. Budgets, seeds, comparator
8. **One budget/seed config (MUST).** v1 is inconsistent: §0 uses 164 vs 657, §6.3 uses 160 vs 660; Q1 grid {33,66,...}
   vs Q2 {30,60,...}. Create `config/budgets.yaml` (Q1 grid a superset of Q2's, decide which of 160/164 and 657/660
   is the headline "5%/20%"); random comparator computed on the union. Pre-add to append-only `seeds.yaml`:
   `q1_draw_seeds, q2_seeds, inner_split_seeds, final_*, bootstrap_seed, model_seed`; derive tie-break, pool-subsample,
   fit, bootstrap RNGs with `SeedSequence(seed, key)`; same `model_seed` for paired strategy-vs-random fits. Flag the
   `prop_folds` rule (recipe changes at n=30) on plots.
9. **Deterministic strategies have no trajectory variance (SHOULD).** One run of H/C/V gives a CI over test queries only,
   not over candidate-pool luck. Give every strategy (incl. static) 5 trajectories via seeded random round-0 batch +
   random 80% candidate-pool subsample; final CI combines both sources (as §6.3 does for stochastic ones).
10. **n_req inversion is unstable (MUST).** Flat/noisy random curve or strategy above the random value at 2,627 makes
    n_req/n explode or be censored. Report n_req with bootstrap CI, censor at ">= 2,627", suppress efficiency when the
    random-slope CI includes 0. Claims rest on ΔAUPRC and the non-inferiority test (c); n_req is descriptive. State
    that the 5% relative margin is fixed (not tuned); report the lower bound itself, not only the verdict.
## C. Number provenance and sub-agent reporting
11. **Harness-written manifest + ledger for every run, dev included (MUST).** `hsl.run_experiment` is the only writer
    of `metrics.json` and `manifest.json` (git head, dirty flag + sha256 of src/scripts/config, resolved config hash,
    data/split/frozen hashes, all seeds, threads, library versions, selection vs fit wall-clock separately, pytest
    pass/fail/SKIP counts, attempt) and appends to `results/LEDGER.jsonl` incl. failed/abandoned runs and parent id.
    Report states total variants tried (multiplicity behind "<=3 finalists"). Results append-only (`_r2`). Skips
    fail under `HSL_STRICT=1`. Headline/final refuses a dirty tree.
12. **No hand-typed numbers (MUST).** `scripts/make_tables.py` writes tables, figures and `results/claims.json` (id,
    value, CI, manifest id); REPORT numbers come from generated blocks or cite a claim id; a test greps REPORT.md for
    numbers absent from claims.json. A [!] tree node needs `evidence:` (result path + key) and `verified-by:`.
13. **Independent check on predictions, not summaries (MUST).** `.gitignore` drops `results/**/predictions/`, so save
    compact prediction arrays locally with sha256 in the manifest. A different agent, importing nothing from `hsl`,
    rescores AUPRC and ΔAUPRC from labels + saved predictions for >= 3 random (strategy, n) cells and the Q2(c)
    verdict; a mismatch blocks the report. Agent summaries paste manifest metric lines; PI spot-checks 1 in N.
14. **Vacuous-test guard (MUST).** Each leakage/metric/selection test ships a broken twin it must catch (selector that
    peeks at T columns or unscreened rows must fail the poisoning test; full-matrix strategy must "win" planted toys).
    Otherwise "selections bit-identical" can pass because the selector ignores its input (e.g. nanmean swallowing
    poison). End-to-end controls (shuffled screened labels -> ~prevalence; oracle back door -> ~1) at every milestone,
    numbers in the manifest.
## D. Feature/selector hygiene specific to this design
15. **Selector feature whitelist (MUST).** Selector may use only features with a named human analogue (§1) and
    array-keyed SMF; never `query_smf`, `dmf`, `n_arrays/n_queries`, or the observed-cell mask of unscreened rows
    (query_smf of an unscreened candidate exists only after screening). "Off-array / missing SMF" (631 genes) is a
    yeast-only artefact a selector could learn and HAP1 cannot use: ban such indicators in the selector, report them
    if used in the predictor.
16. **Stale/cached features (SHOULD).** Upstream keeps `D._CACHE` (bundle) and `lru_cache` on the cofit store. Rule: no
    on-disk cache of anything depending on the screened set; per-process caches for static data only; manifest records
    sha256 of every external feature file loaded; feature coverage report (match rate to 3,283 query / 3,827 array
    ORFs, fail below threshold) + univariate signal per feature, so a join bug is not reported as "feature useless".
17. **Lockbox reuse (SHOULD).** 96.5% of T was upstream's lockbox, on which upstream's method ranking was confirmed;
    predictor choice is not blind. Keep comparisons within-predictor and word absolute AUPRC as not an unbiased
    estimate of a new method.
## E. Process, compute, context, read-only upstream
18. **Guard covers Python only (SHOULD).** `hsl.guard` does not see sub-agent shell commands (`git status`, redirects,
    `uv run` with cwd in YEAST_SLs). Add `.claude/settings.json` deny rules (Write/Edit on `YEAST_SLs/**`; Bash `cd
    *YEAST_SLs*`), a brief line "no cwd or write paths in YEAST_SLs; read by absolute path", and
    `scripts/check_readonly.sh` (stamp `find -newer` + `GIT_OPTIONAL_LOCKS=0 git status`) run after every agent batch and
    before every commit. No chmod of upstream (itself a modification).
19. **Brief template + file ownership (SHOULD).** Template: goal; owned files; frozen files; output <=15 lines (result
    path, manifest id, metric lines); forbidden actions. One agent per strategy/predictor file; shared modules
    PI-owned; changes to anything a result depends on go in `*_v2`. Never `git commit -a` with agents' uncommitted files.
20. **Compute guard (SHOULD).** One launcher (`scripts/launch.py`) enforces per-job timeout, <=2 concurrent jobs and the
    RAM cap (m1_hgb at n>=2,000 needs 12 GB) instead of relying on instructions; each new method/strategy gets a
    1-split pilot with cost estimate; CPU-hour budget line in the tree. Resumable jobs with status files. Profile greedy
    V/VH (2,100 candidates x 66 rounds) early; it may need a low-rank kernel.
21. **Resume protocol and claim tiers (COULD).** Tree gets an "open agents / unmerged work" block and a resume
    checklist (PROMPT, tree summary + decisions, LEDGER tail, LESSONS, `check_stale.py`, `git status`). Tag claims
    `confirmed` (T, pre-registered), `exploratory` (dev) or `translational` (assumption-driven, e.g. yeast AUPRC vs
    HAP1 precision 0.61/recall 0.37, fraction vs count scaling); never mix tiers in one sentence, and never write
    "5% matches 20%" without margin, CI, predictor and evaluation set.
Top three: 1/2/4 (dev/test separation, frozen eval hashes); 5 (tie-correct AUPRC + mutants); 11/12/13 (manifest, generated tables, independent rescoring).
