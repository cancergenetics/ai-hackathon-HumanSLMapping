# Review: reading of the one-shot final on T (results/final, scored at 3e33655, commit 8c1818a)

Reviewer: independent, read-only. I did not produce these results. Date 2026-10-30. Inputs: PROMPT.md; DESIGN §0, §3, §6.3–6.4, §7, §13
and deviations (i)–(vii); RESEARCH_TREE.md (pre-stated readings, D8–D11, final block); results/final/*; results/p4/{summary.md,
shortlist_proposal.md}; results/p1/summary.md; results/q3_final/dev_vs_final.md; docs/reviews/p4_review.md. Pandas checks used the
result CSVs only. No jobs were run and no other file was changed.

## Verdict: the reading is CORRECT, with 4 overstatements to fix before REPORT.md

- Every number in the tree's final block matches summary.md and the CSVs (§1).
- Pre-stated reading (2) is applied faithfully: (a) is rejected for H, C and F_r16, and c1 for none. The claim is therefore "beats
  random at matched budget", and 'matches 20%' is correctly marked not allowed. Readings (1) and (3) do not trigger.
- The 4 overstatements: (1) "worth about 3× the random budget" restates c1 (490/160 = 3.06×), which was not shown; (2) "no loss
  from large batches" holds for H only, C loses 6–7%; (3) E_unobs is buried, against reading (4); (4) the secondary metrics are
  omitted, and they show the per-query gain is far smaller than the pooled gain.

## 1. Tree vs summary.md (check 1): all match

| item | tree | summary.md / CSV |
|---|---|---|
| (a) | rejected; Holm p 0.0003/0.0001/0.0014; s160 0.0555/0.0555/0.0534 vs 0.0423 | identical (+31.2/+31.0/+26.2%) |
| c1 | Δ −0.002/−0.002/−0.004 vs r490 0.0576; not rejected | −0.0020/−0.0021/−0.0042; Holm p 1.0 |
| G, G*, matches | 0.65 [0.41,0.90], 0.65 [0.43,0.87], 0.55 [0.29,0.80]; 0.85; not allowed | same; 0% non-finite reps; G* recomputes to 0.845; p ≥ 0.947 |
| headline / n_req/n | all 3 (P*, m1_hgb × all/AO/TT); 2.7 at 160, 3.4–3.9 at 330–490 | 6/6 positive signs; H 2.71/3.36/3.87, C 2.70/3.72/3.83 |
| Q1 | 0.0917; n_80 1,143 [970,1,311] (0.44); n_90 1,761 [1,571,1,899] (0.67); B1 119 | same (0.435/0.670) |
| scenarios | warm +27–35% at 440, +38–40% at 660; H-static +12–22%; F-static +1%; R1 ≈ 0 | 27.0/35.2, 39.5/39.0; 11.9–22.4% (n ≥ 160); +0.7%; −3.4…+0.5% |
| sensitivities | lockbox +32, nosamechr +33, onrows +40, offrows +2–6 | 31.9/32.7/40.5/1.9–6.2 (offrows n.s.) |
| E_unobs | 0.0521 vs 0.0457 (lift 5.2 vs 4.1); equal at 660, lift 8.8 vs 6.4 | 0.0521/0.0457, 5.18/4.13; 0.06946/0.06951, 8.80/6.38 |

## 2. Claims for REPORT.md (check 2)

Labels: [C] confirmatory (Holm, α = 0.025); [P] pre-specified rule, descriptive; [D] descriptive.

### Allowed claim sentences
- **Q1 [D].** With random selection on T, P* pooled AUPRC rose from 0.031 at 30 queries (lift 3.2) to 0.042 at 160 (4.9%; lift 4.5), 0.063
    at 660 (20%; lift 6.7) and 0.092 with the whole 2,627-gene pool screened (lift 9.8), with no plateau. 80% / 90% of the full-pool AUPRC was reached at n_80 = 1,143 [970, 1,311] and n_90 = 1,761 [1,571, 1,899]: 44% [37–50%] and
    67% [60–72%] of the pool, or 35% and 54% of all 3,283 queries. Simple baselines saturate early (B1 n_80 = 119 [112, 127]). The figures are predictor-specific: m1_hgb exceeds P* at
    n ≥ 660 (0.107 vs 0.092 at the pool max, R = 5).
- **Q2(a) [C].** At a matched budget of 160 queries (batches of 10), all three pre-specified adaptive strategies beat random
  (AUPRC 0.0423). The sign was the same on E_AO, on E_TT and with m1_hgb (headline rule met).

  | strategy | AUPRC | Δ [95% CI] | rel. | Holm p |
  |---|---|---|---|---|
  | H | 0.0555 | 0.0132 [0.0064, 0.0200] | +31% | 0.0003 |
  | C | 0.0555 | 0.0131 [0.0071, 0.0191] | +31% | 0.0001 |
  | F_r16 | 0.0534 | 0.0111 [0.0047, 0.0175] | +26% | 0.0014 |

- **Q2 c1 [C, not rejected].** That 160 strategic queries beat 490 random queries (15%) was not shown: Δ = −0.0020
  [−0.0069, 0.0028] (H), −0.0021 [−0.0065, 0.0023] (C), −0.0042 [−0.0090, 0.0007] (F_r16); Holm p = 1.0. This is not evidence of
  equivalence.
- **Brief example [P].** 160 strategic queries did not match 660 random queries: s160 − r660 = −0.0070 [−0.0116, −0.0024] (H),
  −0.0070 [−0.0114, −0.0026] (C), −0.0091 [−0.0139, −0.0042] (F_r16). The 'matches' criterion failed for all three (p ≥ 0.947).
- **Gap closure [D].** Adaptive selection closed an estimated 65% [41–90%] (H), 65% [43–87%] (C) and 55% [29–80%] (F_r16) of
  the 160→660 random gap (G* = 0.85).
- **Random-equivalent budget [D].** At n = 160, random needed 434 [317, 583] (H), 433 [329, 574] (C) and 384 [268, 516] (F_r16) queries for the same AUPRC:
    2.4–2.7× (ratio CIs ≈ 1.7–3.6×), or ~12–13% of the query set (CI ≈ 8–18%). The ratio is ≈ 1–1.6× at n = 30 and 3.4–3.9× at 330–490 (H, C). At 660 it is right-censored for 11% (H), 4% (C) and 60%
    (F_r16) of replicates.
- **Iteration [D].** Static versions recovered much less: H-static +11.9% at 160 (−14.7% vs H); F-static +0.7% n.s. (−20% vs F).
- **Billmann-like panel [D].** R1 showed no detectable difference from random at any n: rel. Δ −3.4% to +0.5%, all CIs include 0,
  n_req/n 0.89–1.03.
- **Warm start [D].** After a 220-line R1 panel, switching to adaptive gave +27% (H) and +35% (C) at 440, and +39–40% at 660, over
  continuing R1 (CIs exclude 0).
- **Batch size [D].** Batches of 50 or 110 (≤ 6 rounds) still beat random at 160 by +23% to +42% (CIs exclude 0). Against
  batch-10, H showed no loss, while C at b110 lost 6.2% [−0.0064, −0.0004] at 160 and 7.2% at 330, converging by 660.
- **Robustness [D].** The gain held with the reference predictor (+39% H at 160) and m1_hgb (+35–40%), on lockbox rows (+32%), without
    same-chromosome cells (+33%) and under the two-way bootstrap. It is larger on on-array test rows (+35–40%) than on off-array rows (+2–6%, n.s.). The finalists' off-array share
    (0.18–0.21) matches random's (0.20).
- **Where the gain is [D, must be reported].** At 160 the pooled gain (+26–31%) exceeds the per-query macro-AUPRC gain (+5–10%). AUROC is not improved (−0.009 to 0.000) and
    p@50 is ≈ 0. By SL-count tertile, the gain is detectable only in the top tertile (+22–24%); tertiles 1–2 range from −10% to +9%, all n.s. At 660, macro AUPRC is +17–19% and all tertiles gain. So at small n much of the pooled gain is better cross-query ordering (degree/calibration, cf. P2) and top-of-list precision
    on SL-rich queries.
- **E_unobs [D; next to E_fixed per reading 4].** On the arm-specific unobserved test pairs (unpaired, no CI): AUPRC at 160 / 660: H 0.052 / 0.069, C 0.057 / 0.078, F_r16 0.049 / 0.069, random 0.046 / 0.070. Adaptive arms have already observed 1.3–1.6× more SLs, so their unobserved sets are sparser (prevalence at 660: 0.0079–0.0085
    vs 0.0109) and lift is higher at equal AUPRC (H 8.8 vs 6.4). Neither AUPRC nor lift alone is like-for-like.

### Disallowed / overstated phrasings (with the fix)
- "5% strategic ≈ 15% random" / "≈ r490" / "matches 15%". c1 was not rejected, which is not equivalence. Use the n_req sentence
  (≈ 430 [320–580], ~13%).
- "Approaches / nearly matches 20%". Forbidden by reading (3). Give G with its CI and "not met".
- "Worth about 3× the random budget" (tree). This is the c1 claim (3.06×), which was not shown. Use "2.4–2.7× at 160 (CI ≈
  1.7–3.6×), descriptive".
- "Batch size doesn't matter" / "no loss from large batches on T" (tree). C at b110 is −6 to −7% (CI excludes 0).
- "+31% better SL prediction" without saying it is pooled AUPRC. Add macro AUPRC (+10%) and AUROC (no gain).
- "The gain is not the off-array artefact" as a blanket claim. Say what was shown: off-array share at base, and a larger gain on
  on-array rows. Do not cite the C/F det arms; they are not off-array-neutral.
- "Strategies also predict unobserved pairs better" (from lift). At 660, H's and F's E_unobs AUPRC equals random's.
- "H is the best strategy" / "H > C > F". There is no between-finalist test. H ≈ C at 160, and at 660 the order is F > H > C.
- "The gain grows with n". It does on T (31→44%) but not on dev (33→28%).
- "Deterministic selection is as good / better". The det arms have R = 1, so their CIs omit trajectory variance (48
  `traj_var_missing` rows), and C/F det are confounded by off-array avoidance.
- "R1 equals random" (say "no detectable difference"). Any asymptote from the saturating fit (§5).

## 3. Dev → T shrinkage (check 3)

| finalist | rel Δ(160) dev mean (splits) → T | G dev → T | n_req/n(160) dev → T | rel Δ(660) dev → T |
|---|---|---|---|---|
| H | +32.7% (19.7, 52.5, 26.0) → +31.2% | 0.84 [0.63,1.05] → 0.65 [0.41,0.90] | 3.33 → 2.71 | +28.2% → +43.5% |
| C | +29.9% (19.7, 44.6, 25.6) → +31.0% | 0.77 [0.58,0.96] → 0.65 [0.43,0.87] | 3.10 → 2.70 | +22.6% → +40.9% |
| F_r16 | +28.0% (14.7, 44.4, 24.8) → +26.2% | 0.72 [0.49,0.93] → 0.55 [0.29,0.80] | 2.91 → 2.40 | +27.1% → +47.6% |

- **No winner's curse on the effect.** The finalists' mean rel Δ(160) was 30.2% on dev and 29.5% on T, and every T value lies
  inside the dev per-split range.
- **G and n_req fell by ≈ 20% because T's random curve is steeper.** r660/r160 is 1.48 on T against 1.39 on dev. T's rel Δ with
  the dev ratio gives G ≈ 0.312/0.392 = 0.80 for H, against dev's 0.84. Budget-scale metrics depend on the evaluation set's
  random-curve slope, and REPORT.md should say so.
- **Ranking.** H > C > F is preserved nominally at 160 (31.2 ≈ 31.0 > 26.2) but is unresolved and reverses at 660. Say
  "indistinguishable at 160".

## 4. E_unobs (check 4)

- **The tree's numbers and sparser-set explanation are correct.** At 660, prevalence is 0.0079 (H) vs 0.0109 (R0), over ≈ 1.94 M
  cells per arm. The tree omits C, the only arm with an E_unobs AUPRC gain at 660 (+12%), and F (0%).
- **The brief's literal metric.** The brief asks for AUPRC on held-out queries, with pairs that involve a screened gene treated as
  observed. Read literally, that is E_unobs. E_fixed is its strategy-invariant subset, chosen for paired inference (§3).
- **Presentation.** Put E_fixed, E_TT and E_unobs side by side in the main Q2 table. State that the E_unobs advantage is +14% (H)
  at 160 and vanishes at 660 for H and F. Pair it with discovery (SLs observed: H 1.63× random at 160, 1.43× at 660) as the two
  halves of one trade-off. Follow-up (evaluation-only): a combined "test-row SLs observed + recovered at precision ≥ x" per arm.

## 5. Descriptive arms and possible bugs (checks 5, 7)

- **No descriptive arm contradicts the headline.** The only counter-sign entry is F_r16's **−10.3%, which is SL-count tertile 1**
  (the queries with the fewest SLs; AUPRC ≈ 0.015): Δ −0.0016 [−0.0047, 0.0016], p = 0.83. It is not a threat, but it belongs in
  the "where the gain is" qualifier (H and C are also n.s. in tertiles 1–2 at 160).
- **E_TT c1 matters for HAP1 translation (DESIGN §7.4).** s160 − r490 on E_TT: −0.0050 [−0.0118, 0.0017] (H), −0.0068 [−0.0128, −0.0008] (C), −0.0086 [−0.0149, −0.0023] (F). The random curve rises faster there (r160 0.047 → r490 0.068), so the budget equivalence is weaker on E_TT even though
    Δ(160) is +27–34%. Q3 should use its stratum-specific n_req, not the E_fixed 2.7×.
- **q1_satfit is degenerate.** c = 0.0103 [0.0102, 0.0105] sits on the curve_fit lower bound of 0.01
  (`scripts/score_final.py:307`), so a = 0.89 is meaningless. Dev is the same (c = 0.0100 [0.0100, 0.0100]). Report "no
  saturation detectable (fit at bound, ≈ log-linear, n ≤ 660)" and do not quote a, b or c.
- **Q1 dropped bootstrap replicates.** 9/2000 for P*: replicates whose 40 resampled trajectories include none of the 5 pool-max ones. The expected count is
    2000·(35/40)^40 ≈ 9.6, so they are immaterial. `n80_frac_censored` (0.0045) counts these dropped NaNs as censored, which is cosmetic. B1 has 18, about 2× the expected
    count; minor, and the rescorer should confirm.
- **Det arms (R = 1).** The SE omits σ_s²/R_s, so the CIs are too narrow. summary.md marks only C/F det (†); label H__det too.
- **Clean.** No NaN in the Δ, CI, p or SE columns of contrasts.csv or q2_curves.csv, and the G replicates are all finite. Units
  ran at several HEADs (0ef629d…3e33655), consistent with deviation (iv) and a passing frozen check. The Q3 stopping rule was
  frozen (51cae82) before scoring.

## 6. Q1 on T vs dev (check 6)

- **The fractions agree.** n_80 is 0.39 of the pool on dev (810/≈2,100; per split 0.40 / 0.45 / 0.31) and 0.44 [0.37, 0.50] on T,
  inside dev's split range; n_90 is 0.61 vs 0.67.
- **The absolute numbers differ because the anchors differ.** The anchor is 0.080 at ≈ 2,100 on dev and 0.0917 at 2,627 on T. At
  matched n the curves agree (660: 0.060 vs 0.0625).
- **Statement for REPORT.md.** n_80 is relative to the full-pool AUPRC of a curve that has not plateaued: "80% of what screening
  the whole pool gives", not 80% of attainable performance. It is ≈ 40–45% of the pool on both dev and T, consistent with
  fraction scaling (b ≈ 1.06, a bound). Quote T (1,143 [970, 1,311]) as primary and dev only as a consistency check.

## Fixes before REPORT.md (none changes a confirmatory result)
1. Replace "about 3×" with the n_req/ratio sentence plus its CI, and note that c1 (≥ 3.06×) was not shown.
2. Rewrite the batch line to include C's −6 to −7% at b110.
3. Put E_unobs (all three finalists plus discovery) and the E_TT c1 values in the main Q2 table, per reading (4).
4. Add the secondary-metric and tertile qualifier. 5. Drop the satfit parameters and label the det-arm CIs.
6. Present G and n_req as sensitive to the random-curve slope.
