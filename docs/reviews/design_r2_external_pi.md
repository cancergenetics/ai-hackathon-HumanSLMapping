# Design review r2 — external PI

Reviewed: DESIGN.md v2 against PROMPT.md, round-1 reviews + §14, LITERATURE §7, data characterisation, runtime benchmark.
Verdict: the round-1 concerns are resolved in substance (n_X endpoints, scaling experiment, G, oracle gate, R1 + warm
start, one primary predictor, discovery axis, E_common/headline rule). Both rejections in §14 are justified: `.claude`
deny rules are the user's settings, and yeast single-direction labels already carry HAP1-like noise (precision ≈ 0.54
vs 0.61). The remaining problems: (i) §12 still quotes the round-1 estimate, which predates v2's extra work, and one
ambiguous spec could make the oracle 10× more expensive; (ii) three numbers or comparators do not line up; (iii) the
Q3 numbers can still be misread by a consortium.

## High priority
1. **Recompute §12 for v2 and raise throughput.** The "30 job-h" figure predates R_random 10→20, R_strategy 3→10,
   m1_hgb on every final trajectory, universe subsampling, warm start/batch/feasible/minimal scenarios, the
   confirmation split and a full P3–P4 re-run after P5. My rough estimate from the runtime benchmark is 70–90 job-h,
   or 35–45 h wall at 2 concurrent jobs, which is beyond 1–2 days. Fixes: (a) allow ≤ 4 concurrent jobs × 2 threads
   for jobs with n ≤ 660 (stack_lr_wr2 peaks at 2.4 GB RSS; 4 jobs ≈ 10 GB of 24 GB) and keep ≤ 2 for larger n;
   (b) apply items 2, 3, 11 and 12; (c) write the table into §12 with a wall-clock trigger for each cut (e.g. if dev
   P1–P4 > 20 h wall, apply the next cut). *Why:* without this the run either overruns or gets cut ad hoc late.
2. **Fix the oracle spec (§9 P2).** v2 says "best of 200 random candidates **per pick**". At about 3 s per fold-in fit
   that is 200 × 3 s × 160 picks ≈ 27 h per split. The algorithms review costed 200 candidates **per round** (top-10
   by realised gain, 16 rounds ≈ 2.7 h). Write "per round of 10". Run the greedy oracle on dev split 1 only and the
   oracle hub-first (free) on all 3 splits.
3. **Move predictor development (P5) to right after P1, and freeze P\* before P3.** The §8 gate is already defined on
   random dev curves, and the low-n ideas (n-dependent GI weight, propensity feature) do not depend on any strategy.
   Only a covariate-shift correction needs strategy trajectories; keep that one as a post-P4 sensitivity applied to all
   arms. *Why:* this avoids re-running all of P3–P4 (≈ 15–25 job-h) and stops the headline predictor changing under
   the strategy comparison.
4. **Pre-register the yeast denominator N for fraction scaling (§7.3–7.4).** Under fraction scaling the HAP1 line
   count is f × N_HAP1, and f depends on whether N_yeast = 3,283 (constructed queries) or 4,458 (Q ∪ A). That choice
   alone moves the HAP1 number by 36%. Fix it now: N_yeast = 3,283 ↔ the HAP1 viable-KO query universe, a like-for-like
   pair of query universes. Show 4,458 as a sensitivity. In the subsampling experiment, N = queries in each sub-universe.
5. **Define n_80/n_90 on the interpolated curve, with CIs, and state what they mean in absolute terms (§0, §7.1).**
   The upstream curve is still rising steeply at 30% (0.061/0.083/0.102 at 10/20/30%), so n_90 probably falls in the
   1,310–1,970 grid gap. Use the existing isotonic, log-n, bootstrap machinery from n_req (add 1,640 to the grid). Also,
   every n_X in REPORT must be shown with its absolute AUPRC, lift and precision at 10% recall. "90% of a ceiling of
   ~0.15 AUPRC" must not read as "90% of the network mapped". Put next to it the fraction of all strong SLs detected by
   the screened rows at n_X (§7.2). *Why:* this is the number a consortium will quote.
6. **Make the warm-start comparator consistent.** §0 Q3 promises "warm-start gain over **Billmann-like continuation**",
   but §6.3 compares strategy vs **random** for the next 440. Add an R1-continuation arm (R1 extended from 220 to 660).
   For a consortium that already holds ~222 lines, the actionable contrast is strategy vs "keep doing what we do".
7. **Add random@490 to the m1_hgb checkpoints.** The headline rule (§3) needs the same sign with m1_hgb, and c1 is
   strategy@160 vs random@490. Neither the m1_hgb Q1 grid {30, …, 660, …} nor the Q2 checkpoints {60, 160, 330, 660}
   includes 490, so c1 cannot be checked with m1_hgb as written.

## Medium priority
8. **Power pilot: use the variance of the paired difference (§6.4, §7.1).** V_Q should be the query-bootstrap variance
   of Δ, which is paired over the same test queries. A random-only pilot cannot estimate that directly, so use the
   paired random(210) − random(160) difference as a proxy. Say so explicitly. The per-arm V_Q would overstate the
   variance and wrongly push Q2 into "bounds only".
9. **Make the scaling-law decision noise-aware (§7.3).** At 25% of the universe the dev pool has ~100 validation
   queries, so with only 3 draws the CV < 0.25 rule may be decided by noise. Compute the criterion per bootstrap
   replicate and adopt a law only if ≥ 80% of replicates favour it; otherwise report both scenarios. Add a descriptive
   check of absolute curve collapse (AUPRC_U plotted against f and against n). n_90 is relative to each universe's own
   ceiling, so it can hide level shifts.
10. **Add a batch-size scenario with a round count HAP1 can actually run (§6.3).** Batches of 50 still mean 13 rounds
    to 660, and each HAP1 round takes months. Use batch ∈ {10, 50, 110 (= 6 rounds to 660)} plus one-shot, and base
    the Q3 recommendation on the best strategy that is feasible in ≤ 6 rounds.
11. **Cut V and merge C/C-static (§6.2).** V (GP variance on the mk kernel) and F are both uncertainty-reduction
    designs on the same geometry. F has the built-in r = 1 → H ablation, so keep F. C is label-free except through the
    GI part, so run one C and report its static and adaptive parts as one arm. This leaves H, H-static, F, F-static, C,
    M, R1 on dev. That is enough to answer "does strategy help" and "does iteration help".
12. **Trim m1_hgb confirmation (§4, §6.4).** Refit only at {160, 660} for strategies and {160, 490, 660} for random,
    with 5 trajectories per arm. The "agree" rule only needs the sign of Δ(160). This saves about 6–8 job-h.

## Low priority
13. Randomised greedy (top-20 sampling): also run one temperature-0 trajectory per frozen strategy, and state that the
    recommendation is the evaluated randomised rule.
14. Round-1 external PI item 9 is only partly in v2. Add one P0 table comparing the column SMF and dev-row SL-degree
    distributions of E_AO and E_TT against all arrays (labels only, dev rows).
15. Round-1 external PI item 1 (marginal returns): report ΔAUPRC and Δ(SLs detected) per +100 lines in the Q3 table.
    These are the units a consortium will budget in.
16. Optional post-freeze items (§7.4 HAP1 checks, strict label) should run only after the REPORT draft, if compute
    remains. Say so in §9.

## No change needed
- §1 admissibility (genome excluded; SMF rules; tier 2 as sensitivity only), §2 split, §3 observation model and AB/BA
  handling, E_fixed/E_TT/E_common and the headline rule, tie-aware AP.
- §6.1 SelectorView and poisoning tests, §6.4 confirmatory family/Holm and the G/"matches" rule, §8 ledger/freeze,
  §10 one-shot final, §11 verification.
- Scope stays within the brief. Do not add more strategies or predictors; the gains will come from items 1–7.
