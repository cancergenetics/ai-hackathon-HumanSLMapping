# Design review r2 — leakage controls, provenance, process (data + AI-agent failure modes)

Read: PROMPT, DESIGN v2, r1 data + failure-mode reviews, recon.md, src/hsl/{upstream,guard}.py, ysl/models/view.py,
ysl/data/dense.py, config/. No model run, no label read. v2 closes r1 on paper; below are implementability gaps on
the real upstream objects, a few silent assumptions, and cuts. P = priority.
## No change needed
SelectorView as a masking `_BundleSource` feeding the unchanged `ScreenView._construct` is implementable; per-channel
poisoning with positive control, E_fixed_dev without T columns, E_AO/E_TT strata, lockbox declaration, tie-grouped AP,
final lock + frozen hashes + independent rescoring are adequate and proportionate (trim per 11-14).
## Proposed changes
1. **[HIGH] `dev_bundle` = `Bundle.blind(T)`, enforced at runtime, not by grep.** `blind` already drops T rows from every
   row-aligned array (E, SP, Y, DMF, DMF_sd, strain_*) and makes `gene_rows(T)` raise. Snag: it sets `unblinded=False`,
   so `ScreenView.from_split` raises for the 23 lockbox genes that sit in D (T holds only 633/656); build dev views via
   `_construct(_BundleSource(dev_b), ...)` and test on those 23. A grep for `load_s1` misses aliases, `importlib`,
   `pd.read_parquet`, `np.load(bundle.npz)`, direct `ysl.data.dense.load_bundle(unblind=True)`; upstream's `data_guard`
   blocks those files only while a method runs. Extend `hsl.guard`: unless `HSL_FINAL=1`, `load_bundle(unblind=True)`
   raises and opening Costanzo raw/processed/bundle files raises unless a thread-local flag set inside `hsl.views` is on.
2. **[HIGH] Identity exclusions are one-sided.** §3 checks only array current name = dev gene. The same collision exists
   query-query (r1 item 5: YLR312C-B -> YLR313C, itself a query): if the two strains fall on opposite sides of the split,
   screening the dev strain reveals the test gene's row, and cell (t, a) with a = the same gene is a self pair that the
   ORF-level rule misses. P0 script enumerates every current-name group over query and array ORFs (incl. the 27
   array-array groups), writes counts to the manifest and applies pre-registered rules: D/T-crossing query pair ->
   remove the D member from the pool; array = dev query -> drop from E_fixed (existing); array = row's own gene -> drop
   cell. Assert all three sets empty afterwards.
3. **[MEDIUM] Poisoning tests must not need T.** With item 1 the dev bundle has no T rows to poison. Run the six-channel
   tests on a fixture (`subset_queries` of real D genes or a toy bundle with a designated pseudo-test set); say so in
   §6.1/§11 so nobody adds an unblind exemption for tests.
4. **[MEDIUM] Whitelist by construction; resolve the DMF ambiguity.** The whitelist is prose and the provenance table is
   "to be written". Make `SelectorView` expose only whitelisted attributes (no strain metadata, obs-mask helpers,
   off-array indicator) and keep the provenance table as checked YAML (a store opened by any selector/predictor but absent
   from it fails a test). §6.1(b) masks DMF while r1 FM-15 banned `dmf`: state whether DMF/DMF_sd of screened rows are
   selector inputs (observed-pair data; HAP1 analogue = double-KO fitness).
5. **[MEDIUM] Freeze list is unspecified.** §10 says "frozen hashes" only. Include leakage/measurement code, not just
   `eval.py`: `hsl/views.py`, `hsl/guard.py`, `hsl/run_experiment.py`, `config/protocol.yaml` (6), split files, seeds,
   budgets, PROVENANCE.md, `yeast_sls_commit.txt`, data-checksum table. Sub-agents cannot edit them; a failing frozen
   test is reported, never worked around.
6. **[MEDIUM] Magic numbers into one frozen `config/protocol.yaml`.** Literals now scattered: H switch n = 30 (a count,
   while §6.2 says cut-offs are pool fractions), top-20 randomised greedy, cluster cap 2, r grid, 4,096 bins, |err| 1e-4,
   B = 10,000/500, margins 5%/2%, Spearman 0.8, CV 0.25, G 0.5, knee 0.95, alpha 0.025, replicate 1.34, HAP1 universe
   7,500-9,500, SMF decile, thread/RAM caps. Protocol constants frozen there; dev-tuned ones (r, w, switch) in
   `finalists.yaml`. Scripts read, never hardcode.
7. **[MEDIUM] `E_common(n)` is R-dependent.** "Unscreened in every trajectory of the compared arms" with R >= 20
   permutations covers 1 - 0.939^20 ~ 72% of D at n = 160 (E_common keeps ~28% of D columns) and ~all at n >= 490
   (E_common ~ E_fixed). The headline sign-agreement rule then changes with R. Define it per matched trajectory pair
   (same index r, same model seed) and average, or drop it.
8. **[MEDIUM] "Confirmation split" is not independent.** Folds 4+5 are candidates (training rows) in dev splits 1-3, so
   their labels shaped every tuned choice. Call it a "re-split check"; the final lock is the only independent test.
   Also dev pools are ~2,102 genes, not 2,627: n_80/n_90 in dev must reference the largest candidate count per split.
9. **[MEDIUM] Universe subsampling (§7.3) leaves features on the full universe.** Removing genes from the GI matrix only
   keeps quantile ranks, k-means/PCA clusters, degrees and precomputed query x target pair kernels computed on all 4,426
   genes. Either mask external lookups to the sub-universe and recompute every derived statistic inside it, or label the
   result a lower bound on information. State the budget grid per universe (fractions must match) and that T is outside
   every sub-universe.
10. **[LOW] r1 #10/#11 assertions never reached v2; pre-split artefacts.** P0 tests: Y == -1 <=> E or P NaN; `stringent_negative(E, P)`
    reproduces Y on D rows; gene counts 3,283/3,827/2,652/1,175/631; E_fixed built from structure with assert (no cell
    with a screened array); all predictors run with zero-SL screened rows (16 dev genes) at n = 10; manifest records `split_seed`/`n_folds` passed to
    `run_method`. Label `results/data_characterisation/` (aggregate only) "all-row, pre-split"; derive tuning numbers from `dev_bundle`.
## Proportionality — cut or simplify
11. **Manifest:** drop per-run sha256 of src/scripts/config (git head + dirty flag + diff hash suffice) and per-run
    hashing of external feature files; verify `results/data_checksums/checksums.json` once per session and at final and
    record its hash.
12. **Ledger:** one line per configuration (not per fit/trajectory): id, parent, config hash, split, outcome. Its job is
    the multiplicity count and abandoned variants.
13. **claims.json:** keep it, but replace "grep REPORT.md for numbers absent from claims.json" (false positives on
    years, n, parameters) with `make_report.py` rendering `{{claim:id}}` placeholders; test only that all placeholders
    resolve.
14. **Final lock / twins / rescoring:** keep lock + clean tree + hashes + predictions saved before scoring; drop separate
    "attempt n + reason" bookkeeping (a ledger line suffices). Broken twins only for leakage, tie-rule, selection and
    n_req tests (~6). Independent rescoring: final run + one dev milestone, >= 3 (strategy, n) cells + Q2 verdict.
