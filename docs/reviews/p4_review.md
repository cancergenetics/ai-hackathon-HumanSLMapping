# Review: P3 oracle ceilings + P4 strategy development (commits fb18bb0, 264ff93, 01ec9ad)

Reviewer: independent; I did not produce these results. Date 2026-10-29. Dev only: `views.dev_bundle` only, with no
test rows, test-query columns used as targets, or `final_bundle`.
Scope: `results/p3/summary.md`, `results/p4/*`, `q2dev.py`, `scripts/40_*`, `scripts/41_*`, `strategies/*`,
`SelectorView`, `harness.py`, against DESIGN v3.1 §3, §5, §6, §8, §9, §12 and D4–D9.
Review code (no production edits): `scripts/review/p4_{poison_replay,mechanism,recompute_stats,altpred,
oracle_capped,minimal_clusters}.py`, `tests/test_review_p4_replay.py` (6 pass); outputs in `results/review_p4/`.

## Verdict: **ACCEPT WITH FIXES**

- **No leakage.** Every headline statistic recomputes exactly. **The gain is legitimate** under the brief's
  observation model and survives two other predictors.
- **The fixes** cover interpretation and record-keeping: the oracle "ceiling" labels, the mechanism text, the R
  decision record, the minimal-features wording and the deterministic-run caveat. None changes the finalists
  H, C, F_r16.
- **How to read the gain:** +15 to +52% per split (mean +28–33%), with shrinkage expected on T. It is not a uniform
  +30%.

## 1. Leakage in adaptive selection: PASS

`p4_poison_replay.py` replays the selection of stored trajectory 0 to n = 160 on a poisoned copy of the dev bundle:
one poisoning per process, no predictor fit. Garbling redraws Y at random (including missingness), E and SP.

| arm (dev1, dev2) | none | (a) D_val rows | (b) D_val+T cols, all rows | (c) 1,940 unscreened pool rows | (d) D_val SMF | ctrl: 1 visible col | ctrl: all visible pool cols |
|---|---|---|---|---|---|---|---|
| H | identical | identical | identical | identical | identical | changed at pick 31 / 31 | changed at 30 / 30 |
| C | identical | identical | identical | identical | identical | changed at 10 / 10 | changed at 10 / 10 |
| F_r16 | identical | identical | identical | identical | identical | changed at 30 / 30 | changed at 30 / 30 |

- The one-column control garbles the visible candidate column with the most SLs among the first 30 screened rows.
  It changes 42–62 of the 160 picks.
- Also verified: CRN (P4 trajectory r has P1 permutation r's model seed on every split); the random arm is the
  frozen-P* P1 rows, with no predictor/harness/views/eval/features/strategy code change between 21b6070 and HEAD;
  the engine AP equals the stored AUPRC for every row loaded.

**Mechanism** (`p4_mechanism.py`, evaluation-only). n = 160, mean of dev1–3; changes are relative to R0.

| arm | SLs in screened rows | …in D_val cols | true q-degree | off-array | ΔAUPRC | on-array target rows | off-array target rows | E_VV | E_AO |
|---|---|---|---|---|---|---|---|---|---|
| H | +64% | +69% | 64 (R0 39) | 0.17 | +33% | +40% | +12% | +42% | +27% |
| C | +50% | +52% | 58 | 0.19 | +30% | +35% | +9% | +37% | +26% |
| F_r16 | +64% | +66% | 64 | 0.20 | +28% | +35% | +8% | +36% | +23% |
| M_lam0.03 | +61% | +62% | 63 | 0.02 | +29% | +44% | +5% | +34% | +25% |
| H-static | +30% | +32% | 50 | 0.15 | +13% | +15% | +9% | +15% | +13% |
| oracle_hub | **+283%** | +278% | **149** | 0.37 | +18% | +14% | **+30%** | +25% | +14% |

- **The proposed hypothesis is wrong as worded.** The hub oracle reveals about 2.5× more SLs than H, in D_val columns
  as everywhere else, and still scores lower. "More revealed information" cannot explain H > oracle.
- **What the data support:**
  - Every adaptive arm screens moderately SL-rich genes, about 1.5–1.65× random's degree. That includes C and M,
    which model no degree. They find these genes through array-side signals of candidates in visible columns (ĉ_x,
    GI column norms).
  - AUPRC(160) has Spearman 0.80–0.89 per split with SLs revealed over all non-oracle trajectories. Within the 10
    random permutations alone it is 0.36–0.58.
  - The true-degree oracle overshoots: mean degree 146–153 and 37% off-array. Its gain sits on off-array target
    rows (the D5 yeast artefact).
- **Overshoot, not the oracle's lack of cap/randomisation, is the cause.** `p4_oracle_capped.py` runs H's own
  selection machinery (top-20 randomised greedy, cluster cap, same seeds) on the TRUE degree. It still picks extreme
  hubs and reaches 0.0426 / 0.0602 / 0.0528 (dev1–3, R = 3), against H at 0.0498 / 0.0661 / 0.0571 and R0 at
  0.0416 / 0.0434 / 0.0453.
- **Channel.** B1_array_rate ranks cells by column SL rate alone. Refitted on the finalists' stored selections, it is
  **worse** than random: −11 to −15% on dev1 for all three (CIs exclude 0; R_r = 4) and −1 to −2% on dev3. SL-rich rows bias the column-rate estimate.
  The P* gain therefore comes through the target side: each on-array target's own column in screened rows, plus the
  GI fold-in kernel. That fits a gain of +35–44% on on-array target rows against +5–12% on off-array target rows.
- **Legitimate.** These are pairs that involve a screened gene, observed "regardless of side". P* sees them by design
  (§3), sees T columns the same way at the final, and every HAP1 library gene has such a column. "H beats the
  true-degree oracle" is explained. Phrase it as "adaptive arms screen SL-richer (not maximal) genes, whose rows
  reveal more SLs involving unscreened genes". Do not phrase it as "beats a ceiling".

## 2. Statistics: PASS (2 low)

- **Independent recompute** (`p4_recompute_stats.py`). It uses `eval.WeightedAP` with hand-written Welch and joint
  bootstrap; `q2dev`'s contrast code is not used.
  - H vs R0 at n = 160, all identical to summary.md: dev2 0.02276 [0.01236, 0.03316], p = 9e-6, V_Q 2.7e-5;
    dev1 0.00819 [−0.00005, 0.01643]; dev3 0.01179 [0.00477, 0.01881].

  - G pooled (H): 0.839 [0.632, 1.052] (reported 0.84 [0.63, 1.05]); with fresh independent weights 0.839
    [0.620, 1.042], per-split CIs within ±0.0003. Unit-weight AP agrees with sklearn to < 1e-6.
- **Welch SE** (`eval.py:273-294`) is §6.4 exactly. V_Q = Var_b of the trajectory-mean contrast under shared w^b.
  ν comes from the trajectory terms only.
- **Shortlist rule** (`41_p4_strategies.py:281-287`) is §8 exactly: mean over 3 splits of relative Δ(160); families
  H{H, H-static}, F{F, F-static}, C/M; ≤ 3 finalists; eligibility R ≥ 3 on all splits (M_lam0.3/M_lam3 out).
- **M λ timeline** (mtimes, job files): all trajectory-0 λ runs done by 16:18; `p4a3_M_lam0.03.json` (trajectories
  1–2 only) written 16:39; first trajectory-1 part 17:55. λ was chosen first (264ff93 already holds some r1 files).
- **No knob was tuned outside §8.** F r was chosen by the ranking itself, as §8 intends. Route calibration was set on
  off-array share, not AUPRC. Smoke runs drive no decisions.
- **L-S1.** The C/M slot is a near-tie shaped by an R = 1 screen. M_lam3 led on trajectory 0 (+32.6%) but never got
  trajectories 1–2. C wins by the rule, and M's 1.5% off-array picks (yeast artefact) make C the better choice anyway.
  Record it as rule + PI judgement, not as C > M.
- **L-S2.** The pooled G and the 3-split mean treat dev1–3 as independent, but the candidate pools overlap about 75%.
  dev2 (+52%) drives the mean, and on dev1 the H and F_r16 CIs touch 0. The finalists are the arg-max of 11
  configurations. Quote the per-split range and expect shrinkage on T.

## 3. Scenarios and diagnostics

- **Warm start: PASS.** `r1_prefix(split, r, 220)` = R1 arm[:220] = H__ws220[:220] = R0__ws220[:220] (checked on
  dev1/2 × r 0/1). The R1 arm has a 440 checkpoint, seeds are CRN, and contrasts are at 440 and 660 (§6.3).
- **Batch 50/110: PASS.** The H__b110 logs show 6 rounds of 110. n = 160 is evaluated on S[:160]: batch 1 plus the
  first 50 of batch 2 in greedy order (harness rule = §6.3).
- **Feasible pool: PASS.** It drops the bottom array-SMF decile of on-array pool genes, keeps off-array genes, and
  uses visible SMF only. The kept off-array yeast genes have unknown fitness, so the −6 to −9% loss is an optimistic
  bound for HAP1.
- **Minimal features: ISSUE (L, wording).** `<X>__min` removes only AP-MS. §6.3 says fitness + coexpression + ESM,
  but the H degree model keeps abundance, conservation and paralogs, and the cap clusters keep the AP-MS kernel.
  The clusters do not matter: I re-clustered without AP-MS (same k, dim, seed; ARI 0.38) and re-ran H__min
  (9 trajectories, n = 160). Picks diverge from pick 2–5, but AUPRC changes by −2.5 / +0.9 / +1.0% (dev1–3). Fix:
  call it "no AP-MS", or drop abundance, conservation and paralogs too. Dev only; does not block.
- **Deterministic C/F: CONFIRMED (L).** `base.calibrate_routes` maps an off-array gene's mid-rank percentile (always
  < 1) onto the remaining on-array gains, so the best off-array gene sits strictly below the on-array max (toy:
  argmax on-array even with off-array gains +3 SD). C/F recalibrate after every pick, so temperature 0 never picks
  off-array after round 0 (0.023 / 0.048 vs 0.19). H scores once per round and is unaffected (0.18). The randomised
  estimand is untouched, but C/F "det +3 / +6%" is confounded with off-array avoidance (a yeast-only effect): report
  the caveat and draw no det-vs-randomised conclusion. Optional fix: map the top percentile to 1, random ties.

## 4. Power / R: ISSUE (M, record-keeping; the substance is sound)

- **The substance is sound.** P4-measured V_Q(Δ) at |T| scale is 1.1–1.5e-5 (checked: 1.39e-5 = mean per-split
  V_Q × 526/656). That is about 25× (g/3.5)², and about 3× even the pilot's "conservative" ρ = 0 value (4.5e-6). The
  ρ_boot ≈ 0.99 sensitivity (R_s ≈ 35–39) assumed random-vs-random correlation, which fails for strategy vs random.
  At R_s = 20 / R_r = 40 the trajectory terms are < 3% of SE (0.0038 vs 0.0037 at R_s = 40), and fewer trajectories
  cannot favour the strategy.
- **A strict reader would flag it.** §5/§7.1: R is provisional from the pilot and "may only go up". The only finite R
  the pilot reported is R_s ≈ 35–39 (≈ 21 before P2), and D9's "a sensitivity, not the provisional R" was written
  after the P4 results. §12 is a second route (R_s = 39 + scenarios > 10 h → "the pilot minimum", undefined).
  **Fix:** declare D9 a protocol deviation in DESIGN before the freeze (cite §7.1 bounds-only, §12, the numbers
  above) and fix R_s = 20 / R_r = 40 for good; or run R_s = 39 / R_r = 78 if the budget allows.

## 5. Other findings

- **M-O1. The "ceilings" are mislabelled.** H, C, F_r16 and M_lam0.03 exceed the hub-first "ceiling" on every dev
  split (dev1: 0.047–0.050 vs 0.042); call it a "true-degree ranking reference" in the P3 summary, tree and
  shortlist. `oracle_greedyHO` (fold-in surrogate, 200 random candidates/round) is a heuristic reference, not an
  "attainable ceiling", and its half-B G = 0.61 is not comparable with H's all-query dev1 G = 0.60.
- **Predictor dependence: PASS** (`p4_altpred.py`; stored selections at n = 160, CRN seeds, P1 random rows of the
  same predictor). stack_lr_wr2 (§4 "unmodified predictor"), dev1/2/3: H +21/+57/+29%, C +21/+46/+27%, F_r16
  +17/+48/+27%; every CI excludes 0 except F_r16 dev1 (p = 0.04). m1_hgb (dev1, R_r = 5): H +21% [0.0012, 0.0151].
  Not a P2 artefact; the §3 m1_hgb headline condition looks likely to hold. Schedule m1_hgb on T (§4).
- **L-O3.** F-static was run at r = 8 but the finalist is F_r16: freeze its static counterpart at r = 16, or say why.
- **L-O4.** The proposal's §6.2 says "H's gain is not finding true hubs", but H's picks have 1.64× random's true
  degree. Pool-wide Spearman (0.02–0.21) is the wrong diagnostic for top-k enrichment.

## Fixes before freeze

1. Relabel the hub and held-out greedy oracles as references (M-O1); replace shortlist §6.2's mechanism with §1's.
2. Declare the R decision as a protocol deviation with its §7.1/§12 basis (§4).
3. Rename the minimal-features scenario; add the det-run caveat for C/F (§3).
4. Freeze F-static at r = 16 (L-O3); schedule m1_hgb on T for the headline rule.
5. Report per-split ranges next to means, expect shrinkage on T (L-S2); record the C/M slot as rule + judgement.
