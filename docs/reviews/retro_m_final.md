# Milestone retrospective at final start (fresh reviewer, 2026-10-30)

Read: PROMPT, tree, LESSONS, DESIGN §7/§11-13, 16 reviews (skimmed), git log, dry-run summary. No code run; results/final/ not read (only the job ledger tail seen).

## 1. Stop doing (6)
1. **Freezing pre-registered rules that were never evaluated on pilot output.** The §7.1 power rule (rho=0) gave no attainable R, so R_s=20/R_r=40
   was set after P4 results existed (D9). The §12 "R at the pilot minimum" action was undefined, so the final was run over budget (D10). Both are
   now declared deviations. Before any future freeze, evaluate every rule and trigger on the pilot numbers and confirm it returns a value.
2. **Developing a predictor while downstream arms run on the old one.** P* changed (D8) after P1 had run pre-P2 curves, so the pilot, universe, random
   and n>=1970 arms were redone (commits ea31501, 1431962). The gain was +4.0-4.5% (+0.001-0.003 AUPRC) as best of 10 adaptive configs with no independent
   confirmation (P2 review). Hold P1 until P* is frozen (the design already ordered it this way), and cap predictor tuning at a prespecified gain.
3. **Quoting oracle "ceilings" and best-of-k dev gains as headlines.** The P3 go/no-go used the in-sample greedy G=2.05. The held-out value is 0.61 [0.21, 0.98].
   "H beats the true-degree oracle" turned out to be overshoot, not a ceiling (P4 review). "+32.7%" is the best of about 7 strategies x knobs on 3 dev
   splits (per-split range +15 to +52%). Report these as dev-selected and expect shrinkage on T (the P4 review says so).
4. **Prose resource rules.** Thread and core budgeting took 4 lesson rewrites (D6, 6058ee9, 505e25e, lesson 6): load 120, then 40, then 44% idle.
   Enforce the budget in launch.py (memory gate now exists) and stop restating it in the tree and in lessons.
5. **Letting reviews find final-only and procedural defects late.** The first final `start` failed on an untested path (D11). The pipeline review found F1-F6
   (clean-tree coupling, `score --partial` peeking, no `check_frozen` in score, no launcher lock, D10 undeclared in DESIGN). Two strategy "High" bugs (F-1, C-1: off-array
   share inverted) were found only at review although DESIGN §1 already named the failure. Require a build-time checklist (final-only path, selector off-array share, partial-score refusal) before requesting review.
6. **Tree bookkeeping by append only.** The top "Status summary" still reads "2026-10-29, M0 setup" although the brief requires a two-paragraph status at each
   milestone. P1 is still `[~]` although batch6 is done (1431962). Decisions D10 and D11 live only in the tree, not in DESIGN. Update the summary and status marks in the same commit as each decision.

## 2. What worked well (6)
1. **Leakage prevented by construction, then attacked.** Bundle.blind, the runtime label gate and SelectorView are in place. The P4 poison replay gave bit-identical selections under 4 poisonings,
   while the positive controls changed 42-62 of 160 picks. Shuffled-label and back-door controls gave 1.05x prevalence and 1.000. The broken-twin tests catch vacuous passes.
2. **Independent review at every milestone, with real catches and accepted fixes.** Examples: 11 design reviews (exact weighted AP, G*, tie-aware AP, predictor-before-strategies);
   P0 (censored CI); strategies (F-1/C-1, off-array calibration D5); P2 (1,418 tracked .npy, gain reframed); P4 (mechanism corrected, oracles relabelled); final pipeline (F1-F6).
   The PI also checked the one reviewer test that was edited.
3. **Dry run plus frozen lock for the one-shot.** The dev dry run had 179/179 units, crash/resume and refusal tests. It exercised `score`. The freeze covers 81 hashed files
   and 461 upstream checksums, the lock is O_EXCL, and a seeds ledger is kept. The dry-run `frozen check: NOT OK` shows the check path is live.
4. **Honest handling of negative and qualifying results.** R1 (Billmann-like) is about random (+0.7%), F-static is -6.5%, "matches 20%" is not met on any dev split, c1 is borderline,
   the dev1 fold is hard (0.0673), and Q2 is phrased as bounds when the pilot rule fails. Deviations are written with the direction-of-bias argument (DESIGN §15 (ii)-(iv)).
5. **Mechanistic diagnostics before belief.** The harness reproduces upstream (0.0774 vs 0.0775). The gain is traced to SL-rich (1.5-1.65x degree) screened genes acting
   through target-side columns (+35-44% on-array vs +5-12% off-array). B1_array_rate is worse than random on the same selections. The yeast off-array artefact (52 vs 37 SLs) is caught.
6. **Paired design.** Common random numbers, shared bootstrap weights, V_Q measured directly, universe-subsampling to decide fraction vs count scaling
   (b = 1.06 [0.85, 1.28]), and scenarios (warm start +25%, batch 110 -5%, feasible pool -6 to -9%, minimal features +-2%) give Q3 concrete inputs.

## 3. LESSONS.md pruning
Redundant, stale or too specific: #1 (disk part stale since 38 GB to 500 GB; merge with #3), #2 (one-off gitignore), #3, #4 and #7 (one family), #5 (very specific),
#6 (three rewrites in one line), #8 (one incident). Proposed list of 7:
1. YEAST_SLs read-only: read in place behind `hsl.guard`, upstream git only with `GIT_OPTIONAL_LOCKS=0`, run `check_readonly.sh` per milestone. (Evidence: stamp trip; 2026-10-29.)
2. Repo hygiene with agents active: anchored ignore patterns (`/data/`), stage explicit paths, anchored tree edits verified by grep. (Evidence: lost vendored data; lost two tree entries.)
3. Policy changes (disk, threads) state retroactivity and who may delete what. (Evidence: P1 bulk delete flagged by the safety monitor.)
4. Compute: one global core/RAM budget enforced by the launcher; size by measured CPU/RSS, not thread caps; one PID per background wait. (Evidence: load 120, 44% idle, two PIDs.)
5. Dry runs must execute every final-only path (deep preflight, label gate, scoring, resume). (Evidence: first final `start`; F1-F6.)
6. Run the pre-registered rule/trigger on pilot output before freezing it. (Evidence: D9, D10.)
7. Dev-selected winners and oracles are labelled as such; hold an unused dev split for any predictor change. (Evidence: P2 +4.45% unconfirmed; in-sample G 2.05 vs 0.61.)
Missing and worth adding: #6 and #7 above; "a decision lives in DESIGN and the tree summary, not only the decision log" (D10 and D11 absent from DESIGN; top summary stale).

## 4. Risks to the remaining steps (scoring, rescoring, Q3, REPORT)
Brief deliverables checklist:
- Written summary with figures: **not started** (REPORT.md and M4-M7 are absent). Only dev and dry-run figures exist.
- Recommended screening strategy and estimated query numbers: **not yet computed** on T (Q3 inputs are dev n_80/n_90 = 810/1279 of 2,627, and b = 1.06).
- Reproducible code: present (321 tests, scripts 00-41, run_final, score_final).
- DATA.md with source, version, checksum, licence: present.
- Independent rescoring and claim-rendered numbers: **no scripts exist**. DESIGN §11 promises a no-`hsl`-import rescorer by an agent that did not produce the results
  (>= 3 strategy/n cells plus the Q2 verdict) and `results/claims.json` with `{{claim:id}}` placeholders. `grep` finds no `claims` code. The only recompute
  (`scripts/review/p4_recompute_stats.py`) imports `hsl.eval`, `p1` and `q2dev`, so it is not independent. Write and test the rescorer on the dry run now, while the final runs.
- Q2 claim risk:
  - The dry run shows the confirmatory SE is about 0.004 with Holm at alpha 0.025 over 6 tests, and (a) did not reject for H, C or F (Holm p 0.05-0.16).
    The dev per-split gains are heterogeneous. Pre-state in REPORT that a null on T is possible and what it means (bounds, not "no effect").
  - "5% matches 20% random" is not met on dev. Do not let the report soften it to "approaches". Quote G, G* and r490-equivalence as written.
  - The gain runs through target-side columns of screened rows on E_fixed. The brief's literal set (pairs involving a screened gene count as observed) is E_unobs,
    which is not in the headline sign rule (DESIGN §3, contrary to r1 failure-modes item 6). Report E_unobs and E_TT prominently.
- Q3 overstatement risk:
  - The fraction-scaling slope (b about 1) comes from random sub-universes that mechanically thin genes per module, with external stores not restricted. The tree says "bound",
    and the REPORT must keep b = 0 and b = 1 as both shown.
  - n_80/n_90 are relative to the saturating pool-max AUPRC (0.080, absolute precision low), and the stronger baselines saturate early (stack_lr n_80 about 384-506).
  - HAP1 differences are only partly modelled: essential-library columns are excluded, the noise-robustness arm (HAP1 about 40% false positives) was declined in DESIGN §14, and the ε<-0.35 label check is optional.
  - HAP1 lines (7,500-9,500 denominator) times 0.39-0.61 gives thousands of lines, against 222 in Billmann. Explain this gap explicitly rather than presenting it as a recommendation.
- Undeclared or tree-only deviations to record in DESIGN before scoring:
  - D11 (`HSL_FINAL=1` opened the label gate for hashing at `start`); D4/L1 (clusters include D_val); 15 vs 4 identity-excluded arrays;
  - headline rule requiring a positive shared sign (the pipeline review said "state it");
  - the §12 "pilot minimum" wording (D10 is declared in (iii), the rule text is still wrong); the post-lock HEAD move (0ef629d over lock 9694404) is covered only by (iv).
  - Also confirm that F2 (infinite CI replicates) and the pool-max resampling NaN (about 0.5% of replicates, dropped silently) were fixed or logged in 43af996.
- Process risks during the run:
  - Any edit to a frozen file (config, eval, score_final, q2dev) between now and scoring breaks `check_frozen`. This retro file is a docs note and is not on the list.
  - m1_hgb at n=2,627 had never run (RSS hint 8.3 GB vs 12.4 GB kill cap), and a unit that hits the cap fails on every resume. Watch `LEDGER_final.jsonl` for it.
  - `resume` after HEAD has moved was untested under the (iv) rule. Recovery steps were required by F1 and should be confirmed written.
  - R_r = 5 at the pool max and R=10 at n>660 give wide CIs for Q1 n_80/n_90, so censored values must be reported as bounds.
