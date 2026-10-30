# Final evaluation — DRY RUN (pseudo-final on dev1; numbers are meaningless)

> DRY RUN: dev1 stands in for T, R_s = 2, R_r = 3, reduced grids, B = 400. This file proves the pipeline runs end to end; no number here is a result.

Lock: git `c0f791a4f7`, created 2026-10-29T21:34:22+0000, plan sha `adc1b9b838f8`; units done 179/179; evaluation set `E_fixed_dev` (748157 cells, prevalence 0.00903). Predictor P* `stack_lr_wr2_lb`. Inference: B = 400 shared stratified query-bootstrap replicates; Welch SE (V_Q + σ̂²_s/R_s + σ̂²_r/R_r).

Scored with git `43af996990` (tree dirty: True); frozen check: NOT OK — config/frozen.yaml missing (scripts/make_frozen.py writes it at the freeze); code SHA-256 in summary.json `scored_with`. Headline rule: the shared sign must also be positive.

## Confirmatory family (one-sided, Holm, familywise α = 0.025)

| finalist | H | s | r | Δ | 95% CI | p | p (Holm) | reject |
|---|---|---|---|---|---|---|---|---|
| H | c1: AUPRC_H(160) > AUPRC_R0(490) | 0.0466 | 0.0517 | -0.0051 | [-0.0104, 0.0002] | 0.9706 | 1.0000 | no |
| H | a: AUPRC_H(160) > AUPRC_R0(160) | 0.0466 | 0.0396 | 0.0070 | [-0.0009, 0.0150] | 0.0401 | 0.1605 | no |
| C | c1: AUPRC_C(160) > AUPRC_R0(490) | 0.0479 | 0.0517 | -0.0038 | [-0.0079, 0.0002] | 0.9678 | 1.0000 | no |
| C | a: AUPRC_C(160) > AUPRC_R0(160) | 0.0479 | 0.0396 | 0.0083 | [0.0015, 0.0150] | 0.0082 | 0.0493 | no |
| F_r16 | c1: AUPRC_F_r16(160) > AUPRC_R0(490) | 0.0482 | 0.0517 | -0.0035 | [-0.0083, 0.0013] | 0.9243 | 1.0000 | no |
| F_r16 | a: AUPRC_F_r16(160) > AUPRC_R0(160) | 0.0482 | 0.0396 | 0.0086 | [0.0011, 0.0161] | 0.0120 | 0.0602 | no |

## Brief example (s160 vs r660), gap closure G, G*, 'matches 20% random'

| finalist | s160 − r660 [CI] | G [95% CI] (defined; non-finite reps) | G* | matches: s160 − 0.95·r660 [CI], p | claim allowed |
|---|---|---|---|---|---|
| H | -0.0083 [-0.0135, -0.0031] | 0.46 [0.05, 0.83] (yes; 0.0%) | 0.82 | -0.0056 [-0.0108, -0.0004], p = 0.982 | no |
| C | -0.0071 [-0.0106, -0.0036] | 0.54 [0.22, 0.80] (yes; 0.0%) | 0.82 | -0.0043 [-0.0079, -0.0008], p = 0.992 | no |
| F_r16 | -0.0067 [-0.0109, -0.0026] | 0.56 [0.19, 0.86] (yes; 0.0%) | 0.82 | -0.0040 [-0.0082, 0.0002], p = 0.969 | no |

'Matches' may be stated only if the matches contrast is rejected, c1 is rejected (Holm) and G* ≥ 0.5 (DESIGN §6.4).

## Headline rule (Δ(160) same sign on E_fixed, E_AO, E_TT with P* and m1_hgb)

| finalist | signs | complete | headline |
|---|---|---|---|
| H | pstar:all +, pstar:AO +, pstar:TT +, m1_hgb:all +, m1_hgb:AO +, m1_hgb:TT + | True | **headline** |
| C | pstar:all +, pstar:AO +, pstar:TT +, m1_hgb:all +, m1_hgb:AO +, m1_hgb:TT + | True | **headline** |
| F_r16 | pstar:all +, pstar:AO +, pstar:TT +, m1_hgb:all +, m1_hgb:AO +, m1_hgb:TT + | True | **headline** |

## Δ(n) vs random at matched n (P*; Bonferroni bands over each arm's checkpoints)

| arm | 30 | 160 | 490 | 660 |
|---|---|---|---|---|
| C | 0.0040 [-0.0010, 0.0091] | 0.0083 [-0.0003, 0.0169] | 0.0105 [0.0018, 0.0191] | 0.0110 [0.0032, 0.0188] |
| C__b110 | — | 0.0072 [-0.0007, 0.0152] | — | 0.0109 [0.0041, 0.0177] |
| C__b50 | — | 0.0071 [-0.0011, 0.0152] | — | 0.0112 [0.0045, 0.0180] |
| C__det † | 0.0040 [-0.0006, 0.0086] | 0.0116 [0.0017, 0.0215] | 0.0110 [0.0020, 0.0201] | 0.0115 [0.0033, 0.0197] |
| F-static_r16 | 0.0015 [-0.0090, 0.0120] | -0.0001 [-0.0061, 0.0060] | 0.0009 [-0.0037, 0.0056] | 0.0001 [-0.0047, 0.0049] |
| F_r16 | 0.0015 [-0.0090, 0.0120] | 0.0086 [-0.0009, 0.0182] | 0.0148 [0.0052, 0.0244] | 0.0140 [0.0058, 0.0221] |
| F_r16__det † | 0.0025 [-0.0022, 0.0072] | 0.0123 [0.0017, 0.0228] | 0.0151 [0.0053, 0.0248] | 0.0139 [0.0058, 0.0219] |
| H | 0.0007 [-0.0053, 0.0066] | 0.0070 [-0.0031, 0.0172] | 0.0138 [0.0045, 0.0231] | 0.0151 [0.0064, 0.0238] |
| H-static | 0.0008 [-0.0039, 0.0054] | 0.0050 [-0.0029, 0.0129] | 0.0087 [0.0015, 0.0159] | 0.0063 [0.0001, 0.0125] |
| H__b110 | — | 0.0091 [0.0007, 0.0176] | — | 0.0150 [0.0075, 0.0226] |
| H__b50 | — | 0.0091 [-0.0011, 0.0193] | — | 0.0156 [0.0076, 0.0235] |
| H__det | 0.0012 [-0.0036, 0.0060] | 0.0063 [-0.0025, 0.0150] | 0.0128 [0.0036, 0.0220] | 0.0154 [0.0065, 0.0243] |
| R1 | 0.0023 [-0.0032, 0.0079] | 0.0022 [-0.0131, 0.0174] | -0.0009 [-0.0061, 0.0043] | -0.0011 [-0.0049, 0.0027] |

† C__det: not off-array-neutral (descriptive only); F_r16__det: not off-array-neutral (descriptive only).

## Means with ranges (finalists; trajectory min–max; relative Δ range across the strata)

Per-split ranges do not exist on T (one held-out split); the dev per-split ranges are in results/p4/summary.md. Sets: all, AO, TT, onrows, offrows, tert1, tert2, tert3.

| finalist | n | s mean [traj min, max] | r mean [traj min, max] | rel Δ (E_fixed) | rel Δ range over strata |
|---|---|---|---|---|---|
| H | 160 | 0.0466 [0.0458, 0.0475] | 0.0396 [0.0377, 0.0416] | +17.8% | [-41.3%, +26.8%] |
| H | 660 | 0.0700 [0.0699, 0.0701] | 0.0550 [0.0545, 0.0554] | +27.4% | [+18.2%, +28.9%] |
| C | 160 | 0.0479 [0.0473, 0.0484] | 0.0396 [0.0377, 0.0416] | +20.9% | [-8.6%, +29.8%] |
| C | 660 | 0.0660 [0.0653, 0.0667] | 0.0550 [0.0545, 0.0554] | +20.0% | [+4.7%, +31.1%] |
| F_r16 | 160 | 0.0482 [0.0476, 0.0488] | 0.0396 [0.0377, 0.0416] | +21.8% | [-19.4%, +28.2%] |
| F_r16 | 660 | 0.0689 [0.0687, 0.0692] | 0.0550 [0.0545, 0.0554] | +25.5% | [+15.6%, +27.3%] |

## Random-equivalent budget n_req(n) (inverts the P* random curve; joint bootstrap; censored replicates coded 0/∞)

| arm | n | n_req [95% CI] | n_req / n | censored (left/right) |
|---|---|---|---|---|
| H | 30 | 33 [0, 53] | 1.10 | 0.33/0.00 |
| H | 160 | 307 [170, 500] | 1.92 | 0.00/0.00 |
| H | 490 | ∞ (right) [∞, ∞] | — | 0.00/0.97 |
| H | 660 | ∞ (right) [∞, ∞] | — | 0.00/1.00 |
| C | 30 | 54 [35, 98] | 1.82 | 0.01/0.00 |
| C | 160 | 344 [219, 493] | 2.15 | 0.00/0.00 |
| C | 490 | ∞ (right) [817, ∞] | — | 0.00/0.78 |
| C | 660 | ∞ (right) [∞, ∞] | — | 0.00/0.99 |
| F_r16 | 30 | 38 [0, 55] | 1.25 | 0.16/0.00 |
| F_r16 | 160 | 355 [208, 574] | 2.22 | 0.00/0.00 |
| F_r16 | 490 | ∞ (right) [∞, ∞] | — | 0.00/0.99 |
| F_r16 | 660 | ∞ (right) [∞, ∞] | — | 0.00/1.00 |
| R1 | 30 | 43 [33, 62] | 1.42 | 0.00/0.00 |
| R1 | 160 | 195 [137, 265] | 1.22 | 0.00/0.00 |
| R1 | 440 | 360 [295, 465] | 0.82 | 0.00/0.00 |
| R1 | 490 | 452 [344, 574] | 0.92 | 0.00/0.00 |
| R1 | 660 | 596 [456, 751] | 0.90 | 0.00/0.00 |
| H-static | 30 | 34 [0, 49] | 1.12 | 0.32/0.00 |
| H-static | 160 | 254 [163, 463] | 1.59 | 0.00/0.00 |
| H-static | 490 | 978 [735, ∞] | 2.00 | 0.00/0.46 |
| H-static | 660 | ∞ (right) [778, ∞] | — | 0.00/0.65 |
| F-static_r16 | 30 | 38 [0, 56] | 1.25 | 0.12/0.00 |
| F-static_r16 | 160 | 159 [104, 206] | 0.99 | 0.00/0.00 |
| F-static_r16 | 490 | 532 [371, 724] | 1.09 | 0.00/0.00 |
| F-static_r16 | 660 | 665 [464, 871] | 1.01 | 0.00/0.00 |
| H__det | 30 | 36 [0, 53] | 1.20 | 0.17/0.00 |
| H__det | 160 | 286 [165, 436] | 1.79 | 0.00/0.00 |
| H__det | 490 | ∞ (right) [903, ∞] | — | 0.00/0.94 |
| H__det | 660 | ∞ (right) [∞, ∞] | — | 0.00/1.00 |
| C__det | 30 | 54 [38, 81] | 1.81 | 0.00/0.00 |
| C__det | 160 | 466 [289, 687] | 2.91 | 0.00/0.00 |
| C__det | 490 | ∞ (right) [801, ∞] | — | 0.00/0.84 |
| C__det | 660 | ∞ (right) [∞, ∞] | — | 0.00/0.99 |
| F_r16__det | 30 | 44 [32, 62] | 1.46 | 0.01/0.00 |
| F_r16__det | 160 | 496 [285, 739] | 3.10 | 0.00/0.00 |
| F_r16__det | 490 | ∞ (right) [∞, ∞] | — | 0.00/0.99 |
| F_r16__det | 660 | ∞ (right) [∞, ∞] | — | 0.00/1.00 |

## Q1: n_80 / n_90 (PAVA + log-n interpolation of the random mean curve; anchor = AUPRC at the pool maximum; joint bootstrap)

| predictor | R (at anchor) | anchor n (pool max?) | AUPRC anchor | n_80 [CI] | n_90 [CI] | censored frac 80/90 | boot reps dropped |
|---|---|---|---|---|---|---|---|
| pstar | 3 (3) | 980 (NO) | 0.0604 | 359 [247, 476] | 627 [466, 714] | 0.00/0.00 | 0 |
| B1_array_rate | 3 (3) | 980 (NO) | 0.0421 | 115 [101, 129] | 223 [169, 274] | 0.00/0.00 | 0 |
| B2_product | 3 (3) | 980 (NO) | 0.0341 | 134 [107, 164] | 288 [211, 349] | 0.00/0.00 | 0 |
| stack_lr | 3 (3) | 980 (NO) | 0.0475 | 214 [150, 280] | 449 [356, 567] | 0.00/0.00 | 0 |
| m1_hgb | 2 (2) | 490 (NO) | 0.0503 | 186 [0, 243] | 301 [228, 345] | 0.22/0.00 | 0 |

## Q1 random curves (mean, SD, joint-bootstrap CI of the trajectory mean)

| predictor | n | R | AUPRC | SD | CI | lift | prec@10% recall | SLs detected (frac) |
|---|---|---|---|---|---|---|---|---|
| B1_array_rate | 30 | 3 | 0.0217 | 0.0012 | [0.0204, 0.0235] | 2.40 | 0.047 | 0.015 |
| B1_array_rate | 160 | 3 | 0.0366 | 0.0008 | [0.0344, 0.0392] | 4.06 | 0.075 | 0.078 |
| B1_array_rate | 490 | 3 | 0.0409 | 0.0002 | [0.0385, 0.0438] | 4.53 | 0.082 | 0.234 |
| B1_array_rate | 660 | 3 | 0.0412 | 0.0004 | [0.0389, 0.0440] | 4.57 | 0.084 | 0.317 |
| B1_array_rate | 980 | 3 | 0.0421 | 0.0004 | [0.0396, 0.0449] | 4.66 | 0.085 | 0.468 |
| B2_product | 30 | 3 | 0.0182 | 0.0002 | [0.0169, 0.0202] | 2.01 | 0.032 | 0.015 |
| B2_product | 160 | 3 | 0.0284 | 0.0003 | [0.0259, 0.0324] | 3.15 | 0.055 | 0.078 |
| B2_product | 490 | 3 | 0.0328 | 0.0007 | [0.0297, 0.0374] | 3.63 | 0.065 | 0.234 |
| B2_product | 660 | 3 | 0.0333 | 0.0007 | [0.0301, 0.0380] | 3.68 | 0.065 | 0.317 |
| B2_product | 980 | 3 | 0.0341 | 0.0005 | [0.0309, 0.0391] | 3.78 | 0.067 | 0.468 |
| m1_hgb | 160 | 2 | 0.0387 | 0.0025 | [0.0351, 0.0435] | 4.29 | 0.079 | 0.078 |
| m1_hgb | 490 | 2 | 0.0503 | 0.0011 | [0.0446, 0.0577] | 5.57 | 0.095 | 0.238 |
| pstar | 30 | 3 | 0.0283 | 0.0018 | [0.0257, 0.0310] | 3.14 | 0.054 | 0.015 |
| pstar | 160 | 3 | 0.0396 | 0.0019 | [0.0365, 0.0437] | 4.38 | 0.076 | 0.078 |
| pstar | 490 | 3 | 0.0517 | 0.0013 | [0.0465, 0.0591] | 5.73 | 0.099 | 0.234 |
| pstar | 660 | 3 | 0.0550 | 0.0005 | [0.0484, 0.0634] | 6.09 | 0.103 | 0.317 |
| pstar | 980 | 3 | 0.0604 | 0.0009 | [0.0525, 0.0701] | 6.69 | 0.117 | 0.468 |
| stack_lr | 30 | 3 | 0.0276 | 0.0024 | [0.0245, 0.0304] | 3.05 | 0.055 | 0.015 |
| stack_lr | 160 | 3 | 0.0361 | 0.0017 | [0.0333, 0.0394] | 4.00 | 0.068 | 0.078 |
| stack_lr | 490 | 3 | 0.0433 | 0.0007 | [0.0401, 0.0471] | 4.79 | 0.087 | 0.234 |
| stack_lr | 660 | 3 | 0.0448 | 0.0007 | [0.0413, 0.0488] | 4.96 | 0.089 | 0.317 |
| stack_lr | 980 | 3 | 0.0475 | 0.0003 | [0.0436, 0.0521] | 5.26 | 0.096 | 0.468 |

Saturating fit a − b·n^−c (P*, per permutation on the budgets every permutation has: 30|160|490|660|980; seed-cluster bootstrap): a = 0.9530 [0.9020, 1.0000], b = 0.9582, c = 0.010 [0.010, 0.010] (3 fits; descriptive).

## Strata, sensitivities, secondary predictors, descriptive arms (all descriptive)

| family | a (pred, n) | b (pred, n) | set | boot | Δ | rel | 95% CI | p | R_a/R_b |
|---|---|---|---|---|---|---|---|---|---|
| descriptive | R1 (pstar, 30) | R0 (pstar, 30) | all | query | 0.0023 | +8.3% | [-0.0011, 0.0058] | 0.064 | 2/3 |
| descriptive | R1 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0022 | +5.5% | [-0.0059, 0.0102] | 0.202 | 2/3 |
| descriptive | R1 (pstar, 490) | R0 (pstar, 490) | all | query | -0.0009 | -1.7% | [-0.0045, 0.0028] | 0.709 | 2/3 |
| descriptive | R1 (pstar, 660) | R0 (pstar, 660) | all | query | -0.0011 | -2.0% | [-0.0040, 0.0018] | 0.782 | 2/3 |
| descriptive | H-static (pstar, 30) | R0 (pstar, 30) | all | query | 0.0008 | +2.7% | [-0.0025, 0.0041] | 0.297 | 2/3 |
| descriptive | H-static (pstar, 160) | R0 (pstar, 160) | all | query | 0.0050 | +12.7% | [-0.0012, 0.0112] | 0.055 | 2/3 |
| descriptive | H-static (pstar, 490) | R0 (pstar, 490) | all | query | 0.0087 | +16.8% | [0.0031, 0.0143] | 0.001 | 2/3 |
| descriptive | H-static (pstar, 660) | R0 (pstar, 660) | all | query | 0.0063 | +11.5% | [0.0015, 0.0112] | 0.005 | 2/3 |
| descriptive | F-static_r16 (pstar, 30) | R0 (pstar, 30) | all | query | 0.0015 | +5.3% | [-0.0044, 0.0074] | 0.225 | 2/3 |
| descriptive | F-static_r16 (pstar, 160) | R0 (pstar, 160) | all | query | -0.0001 | -0.1% | [-0.0043, 0.0042] | 0.513 | 2/3 |
| descriptive | F-static_r16 (pstar, 490) | R0 (pstar, 490) | all | query | 0.0009 | +1.7% | [-0.0027, 0.0045] | 0.310 | 2/3 |
| descriptive | F-static_r16 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0001 | +0.2% | [-0.0036, 0.0038] | 0.479 | 2/3 |
| descriptive | H__det (pstar, 30) | R0 (pstar, 30) | all | query | 0.0012 | +4.2% | [-0.0021, 0.0045] | 0.199 | 1/3 |
| descriptive | H__det (pstar, 160) | R0 (pstar, 160) | all | query | 0.0063 | +15.9% | [-0.0006, 0.0131] | 0.036 | 1/3 |
| descriptive | H__det (pstar, 490) | R0 (pstar, 490) | all | query | 0.0128 | +24.8% | [0.0056, 0.0200] | 0.000 | 1/3 |
| descriptive | H__det (pstar, 660) | R0 (pstar, 660) | all | query | 0.0154 | +28.1% | [0.0085, 0.0224] | 0.000 | 1/3 |
| descriptive | C__det † (pstar, 30) | R0 (pstar, 30) | all | query | 0.0040 | +14.1% | [0.0007, 0.0073] | 0.012 | 1/3 |
| descriptive | C__det † (pstar, 160) | R0 (pstar, 160) | all | query | 0.0116 | +29.3% | [0.0038, 0.0193] | 0.002 | 1/3 |
| descriptive | C__det † (pstar, 490) | R0 (pstar, 490) | all | query | 0.0110 | +21.4% | [0.0039, 0.0182] | 0.001 | 1/3 |
| descriptive | C__det † (pstar, 660) | R0 (pstar, 660) | all | query | 0.0115 | +20.9% | [0.0051, 0.0179] | 0.000 | 1/3 |
| descriptive | F_r16__det † (pstar, 30) | R0 (pstar, 30) | all | query | 0.0025 | +9.0% | [-0.0007, 0.0058] | 0.053 | 1/3 |
| descriptive | F_r16__det † (pstar, 160) | R0 (pstar, 160) | all | query | 0.0123 | +31.0% | [0.0040, 0.0205] | 0.002 | 1/3 |
| descriptive | F_r16__det † (pstar, 490) | R0 (pstar, 490) | all | query | 0.0151 | +29.1% | [0.0074, 0.0227] | 0.000 | 1/3 |
| descriptive | F_r16__det † (pstar, 660) | R0 (pstar, 660) | all | query | 0.0139 | +25.2% | [0.0075, 0.0202] | 0.000 | 1/3 |
| descriptive | H__b50 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0091 | +22.9% | [0.0002, 0.0179] | 0.022 | 2/3 |
| descriptive | H__b50 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0156 | +28.3% | [0.0086, 0.0225] | 0.000 | 2/3 |
| descriptive | H__b110 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0091 | +23.1% | [0.0017, 0.0165] | 0.008 | 2/3 |
| descriptive | H__b110 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0150 | +27.4% | [0.0084, 0.0216] | 0.000 | 2/3 |
| descriptive | C__b50 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0071 | +17.8% | [0.0000, 0.0141] | 0.025 | 2/3 |
| descriptive | C__b50 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0112 | +20.4% | [0.0053, 0.0171] | 0.000 | 2/3 |
| descriptive | C__b110 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0072 | +18.3% | [0.0003, 0.0141] | 0.020 | 2/3 |
| descriptive | C__b110 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0109 | +19.8% | [0.0049, 0.0169] | 0.000 | 2/3 |
| full_set | H (pstar, 160) | R0 (pstar, 160) | all | query | 0.0070 | +17.8% | [-0.0009, 0.0150] | 0.040 | 2/3 |
| full_set | H (pstar, 660) | R0 (pstar, 660) | all | query | 0.0151 | +27.4% | [0.0082, 0.0219] | 0.000 | 2/3 |
| strata | H (pstar, 160) | R0 (pstar, 160) | AO | query | 0.0057 | +15.1% | [-0.0010, 0.0125] | 0.047 | 2/3 |
| strata | H (pstar, 660) | R0 (pstar, 660) | AO | query | 0.0151 | +28.9% | [0.0080, 0.0222] | 0.000 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | AO | query | -0.0057 | -11.5% | [-0.0109, -0.0004] | 0.983 | 2/3 |
| strata | H (pstar, 160) | R0 (pstar, 160) | TT | query | 0.0090 | +20.9% | [-0.0019, 0.0199] | 0.053 | 2/3 |
| strata | H (pstar, 660) | R0 (pstar, 660) | TT | query | 0.0152 | +25.2% | [0.0071, 0.0234] | 0.000 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | TT | query | -0.0047 | -8.3% | [-0.0118, 0.0025] | 0.902 | 2/3 |
| strata | H (pstar, 160) | R0 (pstar, 160) | onrows | query | 0.0105 | +26.8% | [0.0015, 0.0196] | 0.011 | 2/3 |
| strata | H (pstar, 660) | R0 (pstar, 660) | onrows | query | 0.0159 | +27.7% | [0.0077, 0.0240] | 0.000 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | onrows | query | -0.0038 | -7.0% | [-0.0094, 0.0019] | 0.906 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | offrows | query | -0.0003 | -0.6% | [-0.0089, 0.0083] | 0.524 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | offrows | query | 0.0130 | +26.8% | [0.0014, 0.0247] | 0.014 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | offrows | query | -0.0059 | -12.4% | [-0.0136, 0.0018] | 0.935 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | tert1 | query | -0.0077 | -41.3% | [-0.0117, -0.0037] | 1.000 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | tert1 | query | 0.0050 | +21.6% | [-0.0076, 0.0176] | 0.219 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | tert1 | query | -0.0100 | -47.6% | [-0.0143, -0.0057] | 1.000 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | tert2 | query | -0.0048 | -12.2% | [-0.0148, 0.0052] | 0.828 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | tert2 | query | 0.0090 | +18.2% | [-0.0049, 0.0229] | 0.103 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | tert2 | query | -0.0135 | -28.1% | [-0.0224, -0.0045] | 0.998 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | tert3 | query | 0.0123 | +20.5% | [0.0026, 0.0220] | 0.007 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | tert3 | query | 0.0163 | +20.2% | [0.0079, 0.0247] | 0.000 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | tert3 | query | -0.0044 | -5.7% | [-0.0102, 0.0015] | 0.928 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | lockbox | query | -0.0137 | -66.5% | [-0.0670, 0.0397] | 0.692 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | lockbox | query | -0.0115 | -49.5% | [-0.0494, 0.0263] | 0.725 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | lockbox | query | -0.0155 | -69.2% | [-0.0445, 0.0136] | 0.853 | 2/3 |
| sensitivity | H (pstar, 160) | R0 (pstar, 160) | nosamechr | query | 0.0064 | +16.4% | [-0.0010, 0.0138] | 0.045 | 2/3 |
| sensitivity | H (pstar, 660) | R0 (pstar, 660) | nosamechr | query | 0.0150 | +28.1% | [0.0080, 0.0220] | 0.000 | 2/3 |
| sensitivity_c1 | H (pstar, 160) | R0 (pstar, 490) | nosamechr | query | -0.0052 | -10.3% | [-0.0101, -0.0002] | 0.980 | 2/3 |
| sensitivity_two_way | H (pstar, 160) | R0 (pstar, 490) | all | two_way | -0.0051 | -9.8% | [-0.0134, 0.0032] | 0.885 | 2/3 |
| sensitivity_two_way | H (pstar, 160) | R0 (pstar, 160) | all | two_way | 0.0070 | +17.8% | [-0.0050, 0.0191] | 0.126 | 2/3 |
| m1_hgb | H (m1_hgb, 160) | R0 (m1_hgb, 160) | all | query | 0.0065 | +16.7% | [-0.0022, 0.0152] | 0.070 | 1/2 |
| m1_hgb | H (m1_hgb, 160) | R0 (m1_hgb, 160) | AO | query | 0.0066 | +18.4% | [-0.0009, 0.0142] | 0.041 | 1/2 |
| m1_hgb | H (m1_hgb, 160) | R0 (m1_hgb, 160) | TT | query | 0.0062 | +14.3% | [-0.0056, 0.0180] | 0.146 | 1/2 |
| m1_hgb_c1 | H (m1_hgb, 160) | R0 (m1_hgb, 490) | all | query | -0.0051 | -10.2% | [-0.0107, 0.0005] | 0.963 | 1/2 |
| full_set | C (pstar, 160) | R0 (pstar, 160) | all | query | 0.0083 | +20.9% | [0.0015, 0.0150] | 0.008 | 2/3 |
| full_set | C (pstar, 660) | R0 (pstar, 660) | all | query | 0.0110 | +20.0% | [0.0049, 0.0171] | 0.000 | 2/3 |
| strata | C (pstar, 160) | R0 (pstar, 160) | AO | query | 0.0068 | +18.1% | [0.0013, 0.0124] | 0.009 | 2/3 |
| strata | C (pstar, 660) | R0 (pstar, 660) | AO | query | 0.0113 | +21.7% | [0.0049, 0.0178] | 0.000 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | AO | query | -0.0046 | -9.3% | [-0.0088, -0.0003] | 0.982 | 2/3 |
| strata | C (pstar, 160) | R0 (pstar, 160) | TT | query | 0.0109 | +25.5% | [0.0011, 0.0208] | 0.015 | 2/3 |
| strata | C (pstar, 660) | R0 (pstar, 660) | TT | query | 0.0107 | +17.7% | [0.0038, 0.0176] | 0.001 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | TT | query | -0.0027 | -4.8% | [-0.0088, 0.0033] | 0.816 | 2/3 |
| strata | C (pstar, 160) | R0 (pstar, 160) | onrows | query | 0.0117 | +29.8% | [0.0034, 0.0200] | 0.003 | 2/3 |
| strata | C (pstar, 660) | R0 (pstar, 660) | onrows | query | 0.0126 | +22.0% | [0.0049, 0.0202] | 0.001 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | onrows | query | -0.0026 | -4.8% | [-0.0073, 0.0021] | 0.860 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | offrows | query | -0.0023 | -5.5% | [-0.0055, 0.0009] | 0.922 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | offrows | query | 0.0045 | +9.3% | [-0.0019, 0.0110] | 0.084 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | offrows | query | -0.0079 | -16.7% | [-0.0124, -0.0035] | 1.000 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | tert1 | query | -0.0016 | -8.6% | [-0.0106, 0.0074] | 0.642 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | tert1 | query | 0.0072 | +31.1% | [-0.0052, 0.0196] | 0.127 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | tert1 | query | -0.0039 | -18.4% | [-0.0125, 0.0048] | 0.822 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | tert2 | query | -0.0018 | -4.5% | [-0.0090, 0.0055] | 0.687 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | tert2 | query | 0.0023 | +4.7% | [-0.0072, 0.0118] | 0.315 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | tert2 | query | -0.0104 | -21.8% | [-0.0164, -0.0045] | 1.000 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | tert3 | query | 0.0112 | +18.8% | [0.0029, 0.0196] | 0.004 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | tert3 | query | 0.0119 | +14.8% | [0.0044, 0.0194] | 0.001 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | tert3 | query | -0.0054 | -7.1% | [-0.0097, -0.0011] | 0.993 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | lockbox | query | -0.0126 | -61.6% | [-0.0648, 0.0395] | 0.683 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | lockbox | query | -0.0116 | -49.7% | [-0.0431, 0.0200] | 0.766 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | lockbox | query | -0.0145 | -64.7% | [-0.0415, 0.0125] | 0.854 | 2/3 |
| sensitivity | C (pstar, 160) | R0 (pstar, 160) | nosamechr | query | 0.0079 | +20.3% | [0.0017, 0.0142] | 0.007 | 2/3 |
| sensitivity | C (pstar, 660) | R0 (pstar, 660) | nosamechr | query | 0.0109 | +20.3% | [0.0046, 0.0171] | 0.000 | 2/3 |
| sensitivity_c1 | C (pstar, 160) | R0 (pstar, 490) | nosamechr | query | -0.0037 | -7.3% | [-0.0074, 0.0000] | 0.974 | 2/3 |
| sensitivity_two_way | C (pstar, 160) | R0 (pstar, 490) | all | two_way | -0.0038 | -7.4% | [-0.0107, 0.0030] | 0.865 | 2/3 |
| sensitivity_two_way | C (pstar, 160) | R0 (pstar, 160) | all | two_way | 0.0083 | +20.9% | [-0.0023, 0.0189] | 0.063 | 2/3 |
| m1_hgb | C (m1_hgb, 160) | R0 (m1_hgb, 160) | all | query | 0.0044 | +11.4% | [-0.0035, 0.0123] | 0.128 | 1/2 |
| m1_hgb | C (m1_hgb, 160) | R0 (m1_hgb, 160) | AO | query | 0.0038 | +10.6% | [-0.0030, 0.0107] | 0.131 | 1/2 |
| m1_hgb | C (m1_hgb, 160) | R0 (m1_hgb, 160) | TT | query | 0.0055 | +12.7% | [-0.0051, 0.0162] | 0.147 | 1/2 |
| m1_hgb_c1 | C (m1_hgb, 160) | R0 (m1_hgb, 490) | all | query | -0.0072 | -14.3% | [-0.0124, -0.0020] | 0.996 | 1/2 |
| full_set | F_r16 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0086 | +21.8% | [0.0011, 0.0161] | 0.012 | 2/3 |
| full_set | F_r16 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0140 | +25.5% | [0.0076, 0.0204] | 0.000 | 2/3 |
| strata | F_r16 (pstar, 160) | R0 (pstar, 160) | AO | query | 0.0075 | +19.8% | [0.0013, 0.0137] | 0.009 | 2/3 |
| strata | F_r16 (pstar, 660) | R0 (pstar, 660) | AO | query | 0.0142 | +27.3% | [0.0077, 0.0208] | 0.000 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | AO | query | -0.0039 | -8.0% | [-0.0087, 0.0009] | 0.946 | 2/3 |
| strata | F_r16 (pstar, 160) | R0 (pstar, 160) | TT | query | 0.0106 | +24.7% | [-0.0003, 0.0215] | 0.028 | 2/3 |
| strata | F_r16 (pstar, 660) | R0 (pstar, 660) | TT | query | 0.0136 | +22.4% | [0.0058, 0.0213] | 0.000 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | TT | query | -0.0031 | -5.5% | [-0.0100, 0.0038] | 0.813 | 2/3 |
| strata | F_r16 (pstar, 160) | R0 (pstar, 160) | onrows | query | 0.0111 | +28.2% | [0.0015, 0.0207] | 0.012 | 2/3 |
| strata | F_r16 (pstar, 660) | R0 (pstar, 660) | onrows | query | 0.0150 | +26.2% | [0.0072, 0.0227] | 0.000 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | onrows | query | -0.0032 | -6.0% | [-0.0097, 0.0033] | 0.840 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | offrows | query | 0.0008 | +2.0% | [-0.0088, 0.0105] | 0.414 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | offrows | query | 0.0104 | +21.3% | [0.0005, 0.0203] | 0.020 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | offrows | query | -0.0048 | -10.1% | [-0.0151, 0.0055] | 0.880 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | tert1 | query | -0.0036 | -19.4% | [-0.0083, 0.0010] | 0.938 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | tert1 | query | 0.0062 | +26.9% | [-0.0067, 0.0191] | 0.172 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | tert1 | query | -0.0059 | -28.1% | [-0.0104, -0.0013] | 0.994 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | tert2 | query | -0.0029 | -7.3% | [-0.0109, 0.0052] | 0.758 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | tert2 | query | 0.0077 | +15.6% | [-0.0045, 0.0199] | 0.108 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | tert2 | query | -0.0115 | -24.1% | [-0.0180, -0.0050] | 1.000 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | tert3 | query | 0.0123 | +20.5% | [0.0031, 0.0214] | 0.004 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | tert3 | query | 0.0159 | +19.8% | [0.0076, 0.0242] | 0.000 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | tert3 | query | -0.0044 | -5.7% | [-0.0096, 0.0009] | 0.950 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | lockbox | query | -0.0093 | -45.1% | [-0.0553, 0.0367] | 0.653 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | lockbox | query | -0.0121 | -51.9% | [-0.0515, 0.0273] | 0.727 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | lockbox | query | -0.0111 | -49.6% | [-0.0389, 0.0168] | 0.783 | 2/3 |
| sensitivity | F_r16 (pstar, 160) | R0 (pstar, 160) | nosamechr | query | 0.0081 | +20.9% | [0.0011, 0.0151] | 0.012 | 2/3 |
| sensitivity | F_r16 (pstar, 660) | R0 (pstar, 660) | nosamechr | query | 0.0139 | +26.0% | [0.0074, 0.0205] | 0.000 | 2/3 |
| sensitivity_c1 | F_r16 (pstar, 160) | R0 (pstar, 490) | nosamechr | query | -0.0034 | -6.8% | [-0.0080, 0.0011] | 0.932 | 2/3 |
| sensitivity_two_way | F_r16 (pstar, 160) | R0 (pstar, 490) | all | two_way | -0.0035 | -6.8% | [-0.0116, 0.0046] | 0.802 | 2/3 |
| sensitivity_two_way | F_r16 (pstar, 160) | R0 (pstar, 160) | all | two_way | 0.0086 | +21.8% | [-0.0036, 0.0208] | 0.082 | 2/3 |
| m1_hgb | F_r16 (m1_hgb, 160) | R0 (m1_hgb, 160) | all | query | 0.0069 | +17.9% | [-0.0019, 0.0157] | 0.059 | 1/2 |
| m1_hgb | F_r16 (m1_hgb, 160) | R0 (m1_hgb, 160) | AO | query | 0.0068 | +18.9% | [-0.0005, 0.0142] | 0.034 | 1/2 |
| m1_hgb | F_r16 (m1_hgb, 160) | R0 (m1_hgb, 160) | TT | query | 0.0074 | +17.1% | [-0.0045, 0.0194] | 0.108 | 1/2 |
| m1_hgb_c1 | F_r16 (m1_hgb, 160) | R0 (m1_hgb, 490) | all | query | -0.0047 | -9.3% | [-0.0097, 0.0003] | 0.966 | 1/2 |
| reference_predictor | H (ref, 160) | R0 (ref, 160) | all | query | 0.0087 | +24.0% | [0.0001, 0.0173] | 0.023 | 2/2 |
| reference_predictor | H (ref, 660) | R0 (ref, 660) | all | query | 0.0145 | +27.0% | [0.0077, 0.0213] | 0.000 | 2/2 |
| reference_predictor | C (ref, 160) | R0 (ref, 160) | all | query | 0.0093 | +25.6% | [0.0016, 0.0169] | 0.010 | 2/2 |
| reference_predictor | C (ref, 660) | R0 (ref, 660) | all | query | 0.0107 | +20.0% | [0.0047, 0.0168] | 0.000 | 2/2 |
| reference_predictor | F_r16 (ref, 160) | R0 (ref, 160) | all | query | 0.0095 | +26.4% | [0.0013, 0.0178] | 0.012 | 2/2 |
| reference_predictor | F_r16 (ref, 660) | R0 (ref, 660) | all | query | 0.0135 | +25.2% | [0.0072, 0.0199] | 0.000 | 2/2 |
| descriptive:same_arm_adaptive_randomised | H-static (pstar, 30) | H (pstar, 30) | all | query | 0.0001 | +0.4% | [-0.0059, 0.0062] | 0.463 | 2/2 |
| descriptive:same_arm_adaptive_randomised | H-static (pstar, 160) | H (pstar, 160) | all | query | -0.0020 | -4.4% | [-0.0096, 0.0055] | 0.702 | 2/2 |
| descriptive:same_arm_adaptive_randomised | H-static (pstar, 490) | H (pstar, 490) | all | query | -0.0051 | -7.8% | [-0.0090, -0.0012] | 0.995 | 2/2 |
| descriptive:same_arm_adaptive_randomised | H-static (pstar, 660) | H (pstar, 660) | all | query | -0.0088 | -12.5% | [-0.0135, -0.0040] | 1.000 | 2/2 |
| descriptive:same_arm_adaptive_randomised | F-static_r16 (pstar, 30) | F_r16 (pstar, 30) | all | query | 0.0000 | +0.0% | [-0.0078, 0.0078] | 0.500 | 2/2 |
| descriptive:same_arm_adaptive_randomised | F-static_r16 (pstar, 160) | F_r16 (pstar, 160) | all | query | -0.0087 | -18.0% | [-0.0164, -0.0011] | 0.987 | 2/2 |
| descriptive:same_arm_adaptive_randomised | F-static_r16 (pstar, 490) | F_r16 (pstar, 490) | all | query | -0.0139 | -20.9% | [-0.0207, -0.0071] | 1.000 | 2/2 |
| descriptive:same_arm_adaptive_randomised | F-static_r16 (pstar, 660) | F_r16 (pstar, 660) | all | query | -0.0139 | -20.1% | [-0.0205, -0.0073] | 1.000 | 2/2 |
| descriptive:same_arm_adaptive_randomised | H__det (pstar, 30) | H (pstar, 30) | all | query | 0.0005 | +1.9% | [-0.0051, 0.0062] | 0.332 | 1/2 |
| descriptive:same_arm_adaptive_randomised | H__det (pstar, 160) | H (pstar, 160) | all | query | -0.0008 | -1.6% | [-0.0050, 0.0035] | 0.641 | 1/2 |
| descriptive:same_arm_adaptive_randomised | H__det (pstar, 490) | H (pstar, 490) | all | query | -0.0010 | -1.5% | [-0.0018, -0.0002] | 0.990 | 1/2 |
| descriptive:same_arm_adaptive_randomised | H__det (pstar, 660) | H (pstar, 660) | all | query | 0.0003 | +0.5% | [-0.0004, 0.0010] | 0.173 | 1/2 |
| descriptive:same_arm_adaptive_randomised | C__det † (pstar, 30) | C (pstar, 30) | all | query | -0.0000 | -0.0% | [-0.0031, 0.0031] | 0.504 | 1/2 |
| descriptive:same_arm_adaptive_randomised | C__det † (pstar, 160) | C (pstar, 160) | all | query | 0.0033 | +6.9% | [-0.0000, 0.0066] | 0.026 | 1/2 |
| descriptive:same_arm_adaptive_randomised | C__det † (pstar, 490) | C (pstar, 490) | all | query | 0.0006 | +0.9% | [-0.0029, 0.0041] | 0.355 | 1/2 |
| descriptive:same_arm_adaptive_randomised | C__det † (pstar, 660) | C (pstar, 660) | all | query | 0.0005 | +0.7% | [-0.0026, 0.0035] | 0.377 | 1/2 |
| descriptive:same_arm_adaptive_randomised | F_r16__det † (pstar, 30) | F_r16 (pstar, 30) | all | query | 0.0010 | +3.5% | [-0.0049, 0.0070] | 0.286 | 1/2 |
| descriptive:same_arm_adaptive_randomised | F_r16__det † (pstar, 160) | F_r16 (pstar, 160) | all | query | 0.0036 | +7.5% | [0.0008, 0.0064] | 0.007 | 1/2 |
| descriptive:same_arm_adaptive_randomised | F_r16__det † (pstar, 490) | F_r16 (pstar, 490) | all | query | 0.0003 | +0.4% | [-0.0019, 0.0024] | 0.404 | 1/2 |
| descriptive:same_arm_adaptive_randomised | F_r16__det † (pstar, 660) | F_r16 (pstar, 660) | all | query | -0.0001 | -0.2% | [-0.0022, 0.0019] | 0.553 | 1/2 |
| descriptive:R1_continuation | H__ws220 (pstar, 440) | R1 (pstar, 440) | all | query | 0.0104 | +21.4% | [0.0036, 0.0171] | 0.001 | 2/2 |
| descriptive:R1_continuation | H__ws220 (pstar, 660) | R1 (pstar, 660) | all | query | 0.0137 | +25.4% | [0.0069, 0.0204] | 0.000 | 2/2 |
| descriptive:random_continuation | H__ws220 (pstar, 330) | R0__ws220 (pstar, 330) | all | query | 0.0041 | +9.0% | [-0.0011, 0.0093] | 0.057 | 2/2 |
| descriptive:random_continuation | H__ws220 (pstar, 440) | R0__ws220 (pstar, 440) | all | query | 0.0101 | +20.7% | [0.0030, 0.0171] | 0.003 | 2/2 |
| descriptive:random_continuation | H__ws220 (pstar, 660) | R0__ws220 (pstar, 660) | all | query | 0.0125 | +22.7% | [0.0063, 0.0187] | 0.000 | 2/2 |
| descriptive:R1_continuation | C__ws220 (pstar, 440) | R1 (pstar, 440) | all | query | 0.0082 | +16.9% | [0.0026, 0.0138] | 0.002 | 2/2 |
| descriptive:R1_continuation | C__ws220 (pstar, 660) | R1 (pstar, 660) | all | query | 0.0096 | +17.8% | [0.0031, 0.0161] | 0.003 | 2/2 |
| descriptive:random_continuation | C__ws220 (pstar, 330) | R0__ws220 (pstar, 330) | all | query | 0.0068 | +14.9% | [0.0010, 0.0126] | 0.011 | 2/2 |
| descriptive:random_continuation | C__ws220 (pstar, 440) | R0__ws220 (pstar, 440) | all | query | 0.0079 | +16.2% | [0.0017, 0.0141] | 0.007 | 2/2 |
| descriptive:random_continuation | C__ws220 (pstar, 660) | R0__ws220 (pstar, 660) | all | query | 0.0084 | +15.3% | [0.0025, 0.0143] | 0.004 | 2/2 |
| descriptive:same_arm_batch10 | H__b50 (pstar, 160) | H (pstar, 160) | all | query | 0.0020 | +4.4% | [-0.0060, 0.0101] | 0.233 | 2/2 |
| descriptive:same_arm_batch10 | H__b50 (pstar, 660) | H (pstar, 660) | all | query | 0.0005 | +0.7% | [-0.0006, 0.0016] | 0.202 | 2/2 |
| descriptive:random | H__b50 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0091 | +22.9% | [0.0002, 0.0179] | 0.022 | 2/3 |
| descriptive:random | H__b50 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0156 | +28.3% | [0.0086, 0.0225] | 0.000 | 2/3 |
| descriptive:same_arm_batch10 | H__b110 (pstar, 160) | H (pstar, 160) | all | query | 0.0021 | +4.5% | [-0.0028, 0.0070] | 0.198 | 2/2 |
| descriptive:same_arm_batch10 | H__b110 (pstar, 660) | H (pstar, 660) | all | query | -0.0000 | -0.1% | [-0.0018, 0.0017] | 0.523 | 2/2 |
| descriptive:random | H__b110 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0091 | +23.1% | [0.0017, 0.0165] | 0.008 | 2/3 |
| descriptive:random | H__b110 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0150 | +27.4% | [0.0084, 0.0216] | 0.000 | 2/3 |
| descriptive:same_arm_batch10 | C__b50 (pstar, 160) | C (pstar, 160) | all | query | -0.0012 | -2.6% | [-0.0079, 0.0054] | 0.724 | 2/2 |
| descriptive:same_arm_batch10 | C__b50 (pstar, 660) | C (pstar, 660) | all | query | 0.0002 | +0.3% | [-0.0031, 0.0035] | 0.413 | 2/2 |
| descriptive:random | C__b50 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0071 | +17.8% | [0.0000, 0.0141] | 0.025 | 2/3 |
| descriptive:random | C__b50 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0112 | +20.4% | [0.0053, 0.0171] | 0.000 | 2/3 |
| descriptive:same_arm_batch10 | C__b110 (pstar, 160) | C (pstar, 160) | all | query | -0.0011 | -2.2% | [-0.0052, 0.0031] | 0.723 | 2/2 |
| descriptive:same_arm_batch10 | C__b110 (pstar, 660) | C (pstar, 660) | all | query | -0.0001 | -0.2% | [-0.0028, 0.0026] | 0.540 | 2/2 |
| descriptive:random | C__b110 (pstar, 160) | R0 (pstar, 160) | all | query | 0.0072 | +18.3% | [0.0003, 0.0141] | 0.020 | 2/3 |
| descriptive:random | C__b110 (pstar, 660) | R0 (pstar, 660) | all | query | 0.0109 | +19.8% | [0.0049, 0.0169] | 0.000 | 2/3 |

Paired P* − reference (stack_lr_wr2) on the same trajectories:

| arm | n | Δ | rel | 95% CI | R |
|---|---|---|---|---|---|
| H | 160 | 0.0017 | +3.8% | [0.0007, 0.0027] | 2 |
| H | 660 | 0.0019 | +2.9% | [0.0013, 0.0026] | 2 |
| C | 160 | 0.0024 | +5.2% | [0.0017, 0.0031] | 2 |
| C | 660 | 0.0016 | +2.5% | [0.0010, 0.0023] | 2 |
| F_r16 | 160 | 0.0024 | +5.3% | [0.0016, 0.0033] | 2 |
| F_r16 | 660 | 0.0018 | +2.7% | [0.0012, 0.0025] | 2 |
| R0 | 160 | 0.0023 | +6.5% | [-0.0019, 0.0066] | 2 |
| R0 | 660 | 0.0011 | +2.1% | [0.0006, 0.0016] | 2 |

## E_unobs per arm (descriptive, P*; not paired)

| arm | n | R | AUPRC | prevalence | lift |
|---|---|---|---|---|---|
| C | 30 | 2 | 0.0348 | 0.01085 | 3.21 |
| C | 160 | 2 | 0.0469 | 0.01014 | 4.62 |
| C | 490 | 2 | 0.0525 | 0.00875 | 6.00 |
| C | 660 | 2 | 0.0566 | 0.00845 | 6.70 |
| C__b110 | 160 | 2 | 0.0478 | 0.01043 | 4.59 |
| C__b110 | 660 | 2 | 0.0559 | 0.00850 | 6.59 |
| C__b50 | 160 | 2 | 0.0465 | 0.01015 | 4.58 |
| C__b50 | 660 | 2 | 0.0562 | 0.00845 | 6.65 |
| C__det | 30 | 1 | 0.0339 | 0.01079 | 3.14 |
| C__det | 160 | 1 | 0.0482 | 0.01003 | 4.81 |
| C__det | 490 | 1 | 0.0536 | 0.00859 | 6.24 |
| C__det | 660 | 1 | 0.0560 | 0.00827 | 6.77 |
| C__ws220 | 330 | 2 | 0.0503 | 0.00997 | 5.05 |
| C__ws220 | 440 | 2 | 0.0507 | 0.00934 | 5.43 |
| C__ws220 | 660 | 2 | 0.0553 | 0.00863 | 6.41 |
| F-static_r16 | 30 | 2 | 0.0321 | 0.01096 | 2.92 |
| F-static_r16 | 160 | 2 | 0.0415 | 0.01074 | 3.87 |
| F-static_r16 | 490 | 2 | 0.0567 | 0.01026 | 5.52 |
| F-static_r16 | 660 | 2 | 0.0574 | 0.01013 | 5.67 |
| F_r16 | 30 | 2 | 0.0321 | 0.01096 | 2.92 |
| F_r16 | 160 | 2 | 0.0415 | 0.00989 | 4.20 |
| F_r16 | 490 | 2 | 0.0500 | 0.00856 | 5.85 |
| F_r16 | 660 | 2 | 0.0527 | 0.00821 | 6.42 |
| F_r16__det | 30 | 1 | 0.0324 | 0.01095 | 2.96 |
| F_r16__det | 160 | 1 | 0.0426 | 0.00970 | 4.39 |
| F_r16__det | 490 | 1 | 0.0508 | 0.00831 | 6.11 |
| F_r16__det | 660 | 1 | 0.0536 | 0.00803 | 6.67 |
| H | 30 | 2 | 0.0308 | 0.01080 | 2.85 |
| H | 160 | 2 | 0.0395 | 0.00987 | 4.01 |
| H | 490 | 2 | 0.0496 | 0.00827 | 6.00 |
| H | 660 | 2 | 0.0534 | 0.00790 | 6.76 |
| H-static | 30 | 2 | 0.0305 | 0.01079 | 2.82 |
| H-static | 160 | 2 | 0.0448 | 0.01014 | 4.42 |
| H-static | 490 | 2 | 0.0561 | 0.00931 | 6.02 |
| H-static | 660 | 2 | 0.0553 | 0.00897 | 6.17 |
| H__b110 | 160 | 2 | 0.0450 | 0.01000 | 4.50 |
| H__b110 | 660 | 2 | 0.0542 | 0.00798 | 6.80 |
| H__b50 | 160 | 2 | 0.0437 | 0.00979 | 4.46 |
| H__b50 | 660 | 2 | 0.0535 | 0.00795 | 6.72 |
| H__det | 30 | 1 | 0.0310 | 0.01078 | 2.88 |
| H__det | 160 | 1 | 0.0416 | 0.00989 | 4.20 |
| H__det | 490 | 1 | 0.0494 | 0.00828 | 5.97 |
| H__det | 660 | 1 | 0.0537 | 0.00789 | 6.81 |
| H__ws220 | 330 | 2 | 0.0458 | 0.00989 | 4.63 |
| H__ws220 | 440 | 2 | 0.0478 | 0.00918 | 5.20 |
| H__ws220 | 660 | 2 | 0.0521 | 0.00829 | 6.28 |
| R0 | 30 | 3 | 0.0304 | 0.01099 | 2.76 |
| R0 | 160 | 3 | 0.0414 | 0.01092 | 3.79 |
| R0 | 490 | 3 | 0.0560 | 0.01068 | 5.25 |
| R0 | 660 | 3 | 0.0606 | 0.01059 | 5.72 |
| R0__ws220 | 330 | 2 | 0.0496 | 0.01089 | 4.56 |
| R0__ws220 | 440 | 2 | 0.0530 | 0.01086 | 4.88 |
| R0__ws220 | 660 | 2 | 0.0609 | 0.01076 | 5.66 |
| R1 | 30 | 2 | 0.0324 | 0.01100 | 2.94 |
| R1 | 160 | 2 | 0.0438 | 0.01097 | 4.00 |
| R1 | 440 | 2 | 0.0517 | 0.01089 | 4.75 |
| R1 | 490 | 2 | 0.0551 | 0.01088 | 5.06 |
| R1 | 660 | 2 | 0.0588 | 0.01083 | 5.43 |

## Discovery, off-array share and DESIGN §3 secondary metrics (P* rows; means over trajectories)

| arm | n | SLs detected | frac of pool SLs | off-array share | AUPRC | macro AUPRC | AUROC | p@10 | p@50 |
|---|---|---|---|---|---|---|---|---|---|
| C | 30 | 1392 | 0.017 | 0.133 | 0.0323 | 0.0604 | 0.718 | 0.055 | 0.051 |
| C | 160 | 8948 | 0.107 | 0.144 | 0.0479 | 0.0839 | 0.749 | 0.081 | 0.057 |
| C | 490 | 25678 | 0.307 | 0.184 | 0.0622 | 0.1060 | 0.763 | 0.106 | 0.065 |
| C | 660 | 32380 | 0.387 | 0.190 | 0.0660 | 0.1127 | 0.766 | 0.107 | 0.067 |
| C__b110 | 160 | 8668 | 0.104 | 0.162 | 0.0468 | 0.0850 | 0.748 | 0.083 | 0.059 |
| C__b110 | 660 | 33730 | 0.403 | 0.191 | 0.0659 | 0.1098 | 0.767 | 0.105 | 0.067 |
| C__b50 | 160 | 9476 | 0.113 | 0.172 | 0.0466 | 0.0831 | 0.754 | 0.083 | 0.058 |
| C__b50 | 660 | 33054 | 0.395 | 0.191 | 0.0662 | 0.1101 | 0.768 | 0.108 | 0.067 |
| C__det | 30 | 1624 | 0.019 | 0.133 | 0.0323 | 0.0670 | 0.719 | 0.060 | 0.051 |
| C__det | 160 | 9033 | 0.108 | 0.025 | 0.0512 | 0.0861 | 0.749 | 0.085 | 0.059 |
| C__det | 490 | 24294 | 0.290 | 0.008 | 0.0628 | 0.1033 | 0.761 | 0.103 | 0.065 |
| C__det | 660 | 30231 | 0.361 | 0.006 | 0.0664 | 0.1099 | 0.767 | 0.107 | 0.066 |
| C__ws220 | 330 | 15614 | 0.187 | 0.183 | 0.0527 | 0.0914 | 0.767 | 0.094 | 0.062 |
| C__ws220 | 440 | 22091 | 0.264 | 0.176 | 0.0566 | 0.0946 | 0.769 | 0.096 | 0.064 |
| C__ws220 | 660 | 32134 | 0.384 | 0.187 | 0.0634 | 0.1068 | 0.770 | 0.105 | 0.067 |
| F-static_r16 | 30 | 1256 | 0.015 | 0.250 | 0.0298 | 0.0594 | 0.715 | 0.055 | 0.053 |
| F-static_r16 | 160 | 7308 | 0.087 | 0.219 | 0.0395 | 0.0801 | 0.750 | 0.083 | 0.057 |
| F-static_r16 | 490 | 20779 | 0.248 | 0.209 | 0.0526 | 0.0916 | 0.768 | 0.096 | 0.062 |
| F-static_r16 | 660 | 28248 | 0.338 | 0.198 | 0.0551 | 0.0960 | 0.772 | 0.098 | 0.065 |
| F_r16 | 30 | 1256 | 0.015 | 0.250 | 0.0298 | 0.0594 | 0.715 | 0.055 | 0.053 |
| F_r16 | 160 | 10812 | 0.129 | 0.162 | 0.0482 | 0.0812 | 0.748 | 0.080 | 0.055 |
| F_r16 | 490 | 29264 | 0.350 | 0.191 | 0.0665 | 0.1030 | 0.761 | 0.106 | 0.063 |
| F_r16 | 660 | 36193 | 0.432 | 0.199 | 0.0689 | 0.1088 | 0.766 | 0.109 | 0.065 |
| F_r16__det | 30 | 1131 | 0.014 | 0.200 | 0.0309 | 0.0609 | 0.719 | 0.060 | 0.054 |
| F_r16__det | 160 | 10567 | 0.126 | 0.037 | 0.0518 | 0.0821 | 0.748 | 0.080 | 0.056 |
| F_r16__det | 490 | 26883 | 0.321 | 0.012 | 0.0668 | 0.1036 | 0.763 | 0.102 | 0.065 |
| F_r16__det | 660 | 32697 | 0.391 | 0.009 | 0.0688 | 0.1089 | 0.768 | 0.106 | 0.066 |
| H | 30 | 1259 | 0.015 | 0.183 | 0.0290 | 0.0589 | 0.702 | 0.049 | 0.049 |
| H | 160 | 11244 | 0.134 | 0.194 | 0.0466 | 0.0757 | 0.743 | 0.079 | 0.054 |
| H | 490 | 29016 | 0.347 | 0.187 | 0.0655 | 0.1046 | 0.758 | 0.106 | 0.063 |
| H | 660 | 36006 | 0.430 | 0.186 | 0.0700 | 0.1118 | 0.764 | 0.110 | 0.066 |
| H-static | 30 | 1056 | 0.013 | 0.100 | 0.0291 | 0.0607 | 0.701 | 0.052 | 0.050 |
| H-static | 160 | 8234 | 0.098 | 0.134 | 0.0446 | 0.0807 | 0.753 | 0.081 | 0.055 |
| H-static | 490 | 23282 | 0.278 | 0.149 | 0.0604 | 0.0967 | 0.769 | 0.099 | 0.063 |
| H-static | 660 | 32242 | 0.385 | 0.152 | 0.0613 | 0.1002 | 0.772 | 0.101 | 0.065 |
| H__b110 | 160 | 9892 | 0.118 | 0.125 | 0.0487 | 0.0862 | 0.751 | 0.085 | 0.059 |
| H__b110 | 660 | 35168 | 0.420 | 0.167 | 0.0700 | 0.1133 | 0.768 | 0.110 | 0.066 |
| H__b50 | 160 | 10568 | 0.126 | 0.138 | 0.0487 | 0.0883 | 0.748 | 0.085 | 0.056 |
| H__b50 | 660 | 36326 | 0.434 | 0.182 | 0.0705 | 0.1125 | 0.764 | 0.112 | 0.066 |
| H__det | 30 | 1079 | 0.013 | 0.100 | 0.0295 | 0.0610 | 0.710 | 0.055 | 0.049 |
| H__det | 160 | 10788 | 0.129 | 0.181 | 0.0459 | 0.0770 | 0.745 | 0.078 | 0.053 |
| H__det | 490 | 28670 | 0.343 | 0.194 | 0.0645 | 0.1004 | 0.757 | 0.103 | 0.062 |
| H__det | 660 | 36333 | 0.434 | 0.195 | 0.0704 | 0.1113 | 0.763 | 0.113 | 0.065 |
| H__ws220 | 330 | 16748 | 0.200 | 0.198 | 0.0500 | 0.0893 | 0.764 | 0.091 | 0.060 |
| H__ws220 | 440 | 23983 | 0.287 | 0.194 | 0.0587 | 0.0942 | 0.766 | 0.097 | 0.062 |
| H__ws220 | 660 | 35054 | 0.419 | 0.190 | 0.0675 | 0.1065 | 0.769 | 0.108 | 0.065 |
| R0 | 30 | 1296 | 0.015 | 0.167 | 0.0283 | 0.0588 | 0.711 | 0.049 | 0.051 |
| R0 | 160 | 6546 | 0.078 | 0.175 | 0.0396 | 0.0804 | 0.752 | 0.080 | 0.059 |
| R0 | 490 | 19571 | 0.234 | 0.178 | 0.0517 | 0.0959 | 0.769 | 0.098 | 0.064 |
| R0 | 660 | 26557 | 0.317 | 0.184 | 0.0550 | 0.0979 | 0.771 | 0.102 | 0.065 |
| R0 | 980 | 39162 | 0.468 | 0.185 | 0.0604 | 0.1043 | 0.776 | 0.107 | 0.066 |
| R0__ws220 | 330 | 12802 | 0.153 | 0.198 | 0.0459 | 0.0896 | 0.766 | 0.090 | 0.061 |
| R0__ws220 | 440 | 16558 | 0.198 | 0.191 | 0.0487 | 0.0939 | 0.769 | 0.094 | 0.063 |
| R0__ws220 | 660 | 25662 | 0.307 | 0.192 | 0.0550 | 0.0985 | 0.775 | 0.101 | 0.065 |
| R1 | 30 | 912 | 0.011 | 0.150 | 0.0307 | 0.0603 | 0.717 | 0.053 | 0.053 |
| R1 | 160 | 6086 | 0.073 | 0.216 | 0.0417 | 0.0826 | 0.756 | 0.082 | 0.059 |
| R1 | 440 | 16916 | 0.202 | 0.190 | 0.0484 | 0.0921 | 0.768 | 0.094 | 0.062 |
| R1 | 490 | 19152 | 0.229 | 0.183 | 0.0508 | 0.0953 | 0.769 | 0.097 | 0.063 |
| R1 | 660 | 26366 | 0.315 | 0.192 | 0.0538 | 0.0978 | 0.773 | 0.099 | 0.064 |

## Runtime

See `runtime.csv` (per arm × predictor × n: fit+predict, unit wall, peak RSS). Bootstrap vectors computed this run: 0.

Files: `summary.json`, `confirmatory.csv`, `contrasts.csv`, `q2_curves.csv`, `n_req.csv`, `q1_curves.csv`, `q1_nx.csv`, `q1_satfit.csv`, `e_unobs.csv`, `discovery.csv`, `runtime.csv`, `runs.csv`, `selections.json.gz`, `lock_record.json`, `run_ledger.jsonl`, `fig_final_q1.png`, `fig_final_q2_delta.png`, `fig_final_frontier.png`

