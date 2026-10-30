# P1 — dev random curves (R0), DEV DATA ONLY

DESIGN v3.1 §7.1/§8. Splits dev1–3 (D_val = inner fold k, candidates = D \ fold k); evaluation set E_fixed_dev = D_val × (A \ D_cand \ T) minus identity exclusions (no test-gene columns). Pooled AUPRC (tie-grouped AP, frozen `hsl.eval`). Trajectories: seeded permutations of the candidate pool (`dev_traj_base_seed`, harness R0 streams), nested prefixes; model seed = trajectory evaluation stream (at the pool maximum all permutations coincide, so the SD there is model-seed variation only). Rows: `results/p1/curves.csv`; per-(split, predictor, n) table: `results/p1/curve_table.csv`.

Figures: `fig_q1_dev_curves.png` (AUPRC vs n, mean ± SD band), `fig_q1_dev_strata.png` (strata).

## AUPRC at key budgets (mean ± SD over permutations; E_fixed_dev)

**dev1** (prevalence 0.0090, pool max 2101)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 | 1310 | 2101 |
|---|---|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2_lb | 10 | 0.0288 ± 0.0036 | 0.0385 ± 0.0014 | 0.0416 ± 0.0021 | 0.0435 ± 0.0024 | 0.0467 ± 0.0027 | 0.0509 ± 0.0028 | 0.0553 ± 0.0022 | 0.0656 ± 0.0014 | 0.0725 ± 0.0000 |
| stack_lr_wr2 (pre-P2 P*) | 10 | 0.0274 ± 0.0032 | 0.0353 ± 0.0012 | 0.0389 ± 0.0021 | 0.0414 ± 0.0024 | 0.0452 ± 0.0025 | 0.0498 ± 0.0027 | 0.0544 ± 0.0022 | 0.0649 ± 0.0014 | — |
| stack_lr | 10 | 0.0279 ± 0.0032 | 0.0345 ± 0.0008 | 0.0376 ± 0.0007 | 0.0391 ± 0.0008 | 0.0412 ± 0.0005 | 0.0433 ± 0.0006 | 0.0450 ± 0.0007 | 0.0492 ± 0.0006 | 0.0524 ± 0.0001 |
| B1 array rate | 10 | 0.0236 ± 0.0011 | 0.0338 ± 0.0008 | 0.0372 ± 0.0007 | 0.0386 ± 0.0006 | 0.0400 ± 0.0004 | 0.0411 ± 0.0003 | 0.0417 ± 0.0003 | 0.0425 ± 0.0001 | 0.0428 ± 0.0000 |
| B2 product | 10 | 0.0193 ± 0.0009 | 0.0263 ± 0.0009 | 0.0290 ± 0.0008 | 0.0302 ± 0.0009 | 0.0318 ± 0.0006 | 0.0329 ± 0.0005 | 0.0336 ± 0.0005 | 0.0347 ± 0.0002 | 0.0352 ± 0.0000 |
| m1_hgb (no genome) | 5 | 0.0299 ± 0.0039 | — | 0.0395 ± 0.0020 | — | — | 0.0521 ± 0.0028 | 0.0561 ± 0.0021 | 0.0694 ± 0.0009 | 0.0783 ± 0.0004 |

**dev2** (prevalence 0.0098, pool max 2100)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 | 1310 | 2100 |
|---|---|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2_lb | 10 | 0.0309 ± 0.0022 | 0.0394 ± 0.0012 | 0.0434 ± 0.0012 | 0.0466 ± 0.0025 | 0.0522 ± 0.0040 | 0.0586 ± 0.0040 | 0.0647 ± 0.0037 | 0.0810 ± 0.0018 | 0.0913 ± 0.0001 |
| stack_lr_wr2 (pre-P2 P*) | 10 | 0.0298 ± 0.0020 | 0.0368 ± 0.0010 | 0.0409 ± 0.0011 | 0.0447 ± 0.0025 | 0.0505 ± 0.0039 | 0.0575 ± 0.0039 | 0.0638 ± 0.0035 | 0.0805 ± 0.0018 | — |
| stack_lr | 10 | 0.0300 ± 0.0021 | 0.0365 ± 0.0011 | 0.0389 ± 0.0011 | 0.0409 ± 0.0012 | 0.0428 ± 0.0011 | 0.0453 ± 0.0012 | 0.0472 ± 0.0013 | 0.0528 ± 0.0009 | 0.0569 ± 0.0001 |
| B1 array rate | 10 | 0.0228 ± 0.0017 | 0.0335 ± 0.0012 | 0.0369 ± 0.0008 | 0.0384 ± 0.0005 | 0.0400 ± 0.0006 | 0.0413 ± 0.0004 | 0.0420 ± 0.0004 | 0.0431 ± 0.0001 | 0.0436 ± 0.0000 |
| B2 product | 10 | 0.0209 ± 0.0014 | 0.0287 ± 0.0011 | 0.0314 ± 0.0009 | 0.0333 ± 0.0010 | 0.0348 ± 0.0011 | 0.0361 ± 0.0009 | 0.0370 ± 0.0006 | 0.0383 ± 0.0002 | 0.0388 ± 0.0000 |

**dev3** (prevalence 0.0096, pool max 2103)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 | 1310 | 2103 |
|---|---|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2_lb | 10 | 0.0328 ± 0.0024 | 0.0405 ± 0.0022 | 0.0453 ± 0.0020 | 0.0489 ± 0.0027 | 0.0534 ± 0.0025 | 0.0577 ± 0.0023 | 0.0612 ± 0.0025 | 0.0700 ± 0.0011 | 0.0763 ± 0.0001 |
| stack_lr_wr2 (pre-P2 P*) | 10 | 0.0310 ± 0.0019 | 0.0380 ± 0.0024 | 0.0435 ± 0.0023 | 0.0475 ± 0.0029 | 0.0526 ± 0.0026 | 0.0573 ± 0.0022 | 0.0610 ± 0.0025 | 0.0701 ± 0.0011 | — |
| stack_lr | 10 | 0.0313 ± 0.0017 | 0.0369 ± 0.0020 | 0.0406 ± 0.0016 | 0.0428 ± 0.0018 | 0.0457 ± 0.0016 | 0.0482 ± 0.0013 | 0.0502 ± 0.0011 | 0.0550 ± 0.0007 | 0.0585 ± 0.0002 |
| B1 array rate | 10 | 0.0250 ± 0.0021 | 0.0342 ± 0.0010 | 0.0373 ± 0.0008 | 0.0387 ± 0.0007 | 0.0402 ± 0.0005 | 0.0412 ± 0.0004 | 0.0417 ± 0.0003 | 0.0428 ± 0.0001 | 0.0431 ± 0.0000 |
| B2 product | 10 | 0.0247 ± 0.0030 | 0.0338 ± 0.0028 | 0.0372 ± 0.0022 | 0.0394 ± 0.0023 | 0.0412 ± 0.0013 | 0.0429 ± 0.0005 | 0.0438 ± 0.0006 | 0.0448 ± 0.0004 | 0.0451 ± 0.0000 |

**Mean over splits** (mean of the per-split means; ± SD of the 3 split means; pool max = each split's)

| predictor | 30 | 100 | 160 | 220 | 330 | 490 | 660 | 1310 | max |
|---|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2_lb | 0.0308 ± 0.0020 | 0.0395 ± 0.0010 | 0.0434 ± 0.0018 | 0.0463 ± 0.0027 | 0.0508 ± 0.0036 | 0.0558 ± 0.0042 | 0.0604 ± 0.0047 | 0.0722 ± 0.0079 | 0.0800 ± 0.0099 |
| stack_lr_wr2 (pre-P2 P*) | 0.0294 ± 0.0018 | 0.0367 ± 0.0013 | 0.0411 ± 0.0023 | 0.0445 ± 0.0031 | 0.0494 ± 0.0038 | 0.0549 ± 0.0044 | 0.0597 ± 0.0048 | 0.0718 ± 0.0079 | 0.0759 ± 0.0091 (n = 1640; no pool-max run) |
| stack_lr | 0.0297 ± 0.0017 | 0.0360 ± 0.0013 | 0.0390 ± 0.0015 | 0.0409 ± 0.0018 | 0.0432 ± 0.0023 | 0.0456 ± 0.0024 | 0.0475 ± 0.0026 | 0.0523 ± 0.0029 | 0.0559 ± 0.0032 |
| B1 array rate | 0.0238 ± 0.0011 | 0.0338 ± 0.0004 | 0.0371 ± 0.0003 | 0.0386 ± 0.0001 | 0.0401 ± 0.0001 | 0.0412 ± 0.0001 | 0.0418 ± 0.0001 | 0.0428 ± 0.0003 | 0.0432 ± 0.0004 |
| B2 product | 0.0216 ± 0.0028 | 0.0296 ± 0.0038 | 0.0325 ± 0.0042 | 0.0343 ± 0.0047 | 0.0359 ± 0.0048 | 0.0373 ± 0.0051 | 0.0381 ± 0.0052 | 0.0392 ± 0.0051 | 0.0397 ± 0.0050 |

**P* vs baselines on the same trajectories** (mean over splits of the per-split mean relative gain; range over splits):

| n | P* vs pre-P2 P* | P* vs stack_lr | P* vs B1 | P* vs B2_product |
|---|---|---|---|---|
| 30 | +4.9% [+3.8%, +6.1%] | +3.7% [+3.0%, +5.0%] | +29.7% [+21.9%, +35.4%] | +43.7% [+33.8%, +49.1%] |
| 100 | +7.5% [+6.6%, +8.9%] | +9.8% [+8.1%, +11.5%] | +16.7% [+13.8%, +18.5%] | +34.7% [+20.7%, +46.1%] |
| 160 | +5.7% [+4.2%, +6.9%] | +11.3% [+10.6%, +11.7%] | +16.9% [+11.9%, +21.2%] | +34.5% [+21.9%, +43.3%] |
| 220 | +4.1% [+2.8%, +5.2%] | +13.2% [+11.2%, +14.1%] | +20.1% [+12.6%, +26.3%] | +36.1% [+24.2%, +44.1%] |
| 330 | +2.8% [+1.5%, +3.4%] | +17.4% [+13.4%, +21.9%] | +26.7% [+16.7%, +32.7%] | +42.3% [+29.5%, +50.3%] |
| 490 | +1.6% [+0.7%, +2.2%] | +22.3% [+17.6%, +29.3%] | +35.3% [+23.9%, +41.8%] | +50.6% [+34.5%, +62.3%] |
| 660 | +1.1% [+0.4%, +1.7%] | +27.3% [+21.9%, +37.0%] | +44.5% [+32.6%, +54.1%] | +59.8% [+39.8%, +74.9%] |
| 1310 | +0.5% [-0.2%, +1.1%] | +38.0% [+27.2%, +53.4%] | +68.5% [+54.4%, +87.8%] | +85.7% [+56.3%, +111.7%] |

## Strata (P*, mean over permutations, split-averaged)

| n | all | E_AO | E_VV | on-array val rows | off-array val rows | lift | macro AUPRC | AUROC | prec@10% recall | SLs detected (frac of dev rows' SLs) |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 0.0308 | 0.0290 | 0.0345 | 0.0311 | 0.0347 | 3.26 | 0.0637 | 0.709 | 0.058 | 0.011 |
| 100 | 0.0395 | 0.0376 | 0.0432 | 0.0398 | 0.0415 | 4.17 | 0.0786 | 0.746 | 0.074 | 0.037 |
| 160 | 0.0434 | 0.0411 | 0.0480 | 0.0440 | 0.0436 | 4.59 | 0.0834 | 0.755 | 0.082 | 0.060 |
| 220 | 0.0463 | 0.0436 | 0.0517 | 0.0473 | 0.0451 | 4.89 | 0.0870 | 0.760 | 0.088 | 0.083 |
| 330 | 0.0508 | 0.0473 | 0.0574 | 0.0525 | 0.0476 | 5.36 | 0.0912 | 0.766 | 0.098 | 0.125 |
| 490 | 0.0558 | 0.0516 | 0.0639 | 0.0583 | 0.0498 | 5.89 | 0.0968 | 0.772 | 0.107 | 0.185 |
| 660 | 0.0604 | 0.0556 | 0.0694 | 0.0637 | 0.0520 | 6.38 | 0.1011 | 0.775 | 0.118 | 0.250 |
| 1310 | 0.0722 | 0.0662 | 0.0833 | 0.0768 | 0.0595 | 7.61 | 0.1119 | 0.782 | 0.148 | 0.498 |

## n_80 / n_90 (budget at which the PAVA, log-n interpolated mean random curve reaches 80/90% of its value at the dev pool maximum; frozen `eval.n_req`; P*: joint bootstrap 95% CI, B = 2,000)

| split | predictor | R | AUPRC at pool max (PAVA) | n_80 [95% CI] | n_90 [95% CI] |
|---|---|---|---|---|---|
| dev1 | P* stack_lr_wr2_lb | 10 | 0.0725 [0.0613, 0.0854] | 836 [624, 945] | 1277 [1124, 1403] |
| dev2 | P* stack_lr_wr2_lb | 10 | 0.0913 [0.0781, 0.1048] | 941 [802, 1054] | 1374 [1237, 1471] |
| dev3 | P* stack_lr_wr2_lb | 10 | 0.0763 [0.0662, 0.0881] | 651 [487, 832] | 1185 [974, 1373] |
| dev1 | stack_lr | 10 | 0.0524 | 384 | 957 |
| dev2 | stack_lr | 10 | 0.0569 | 506 | 1092 |
| dev3 | stack_lr | 10 | 0.0585 | 385 | 914 |
| dev1 | B1 array rate | 10 | 0.0428 | 106 | 216 |
| dev2 | B1 array rate | 10 | 0.0436 | 118 | 270 |
| dev3 | B1 array rate | 10 | 0.0431 | 104 | 229 |
| dev1 | B2 product | 10 | 0.0352 | 137 | 324 |
| dev2 | B2 product | 10 | 0.0388 | 148 | 342 |
| dev3 | B2 product | 10 | 0.0451 | 134 | 290 |
| dev1 | m1_hgb (no genome) | 5 | 0.0783 | 926 | 1386 |

P* mean over splits: n_80 = 810, n_90 = 1279 (pool ≈ 2,100; as fractions of the pool: 0.385 / 0.609).

## Saturating fit a − b·n^(−c) per permutation (P*, descriptive; permutation-bootstrap 95% CI)

| split | fits | a | b | c |
|---|---|---|---|---|
| dev1 | 10 | 0.9302 [0.9032, 0.9539] | 0.9346 [0.9065, 0.9595] | 0.0100 [0.0100, 0.0100] |
| dev1 | 10 | 0.9620 [0.9375, 0.9824] | 0.9663 [0.9409, 0.9880] | 0.0100 [0.0100, 0.0101] |
| dev2 | 10 | 1.0000 [1.0000, 1.0000] | 1.0164 [1.0149, 1.0179] | 0.0127 [0.0123, 0.0131] |
| dev2 | 10 | 1.0000 [1.0000, 1.0000] | 1.0189 [1.0174, 1.0205] | 0.0136 [0.0133, 0.0139] |
| dev3 | 10 | 0.9917 [0.9750, 1.0000] | 0.9979 [0.9794, 1.0087] | 0.0107 [0.0103, 0.0111] |
| dev3 | 10 | 0.9958 [0.9923, 0.9992] | 1.0000 [0.9947, 1.0050] | 0.0105 [0.0102, 0.0108] |

## Runtime per fit (fit+predict s, mean over splits/permutations; peak RSS sampled during the fit, max)

| n | P* stack_lr_wr2_lb | stack_lr_wr2 (pre-P2 P*) | stack_lr | B1 array rate | B2 product | m1_hgb (no genome) |
|---|---|---|---|---|---|---|
| 10 | 49.4 s (1.8 GB, 2 thr) | 32.2 s (1.2 GB, 2 thr) | 0.4 s (1.2 GB, 2 thr) | 0.0 s (1.2 GB, 2 thr) | 0.0 s (1.2 GB, 2 thr) | — |
| 30 | 7.2 s (2.0 GB, 2 thr) | 32.0 s (1.6 GB, 2 thr) | 0.6 s (1.4 GB, 2 thr) | 0.0 s (1.4 GB, 2 thr) | 0.0 s (1.4 GB, 2 thr) | 236.4 s (1.9 GB, 2 thr) |
| 100 | 46.3 s (1.8 GB, 2 thr) | 34.0 s (1.8 GB, 2 thr) | 1.3 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | — |
| 160 | 7.0 s (2.0 GB, 2 thr) | 37.2 s (1.9 GB, 2 thr) | 2.1 s (1.3 GB, 2 thr) | 0.0 s (1.1 GB, 2 thr) | 0.0 s (1.1 GB, 2 thr) | 394.9 s (1.9 GB, 2 thr) |
| 330 | 9.4 s (2.5 GB, 2 thr) | 47.5 s (2.0 GB, 2 thr) | 4.2 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | — |
| 660 | 15.6 s (3.2 GB, 2 thr) | 76.0 s (2.7 GB, 2 thr) | 9.1 s (1.7 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 424.5 s (2.7 GB, 2 thr) |
| 980 | 75.9 s (3.7 GB, 3 thr) | 99.7 s (3.3 GB, 2 thr) | 12.8 s (2.1 GB, 2 thr) | 0.0 s (1.9 GB, 2 thr) | 0.1 s (1.9 GB, 2 thr) | — |
| 1310 | 123.2 s (4.4 GB, 3 thr) | 168.1 s (3.7 GB, 2 thr) | 20.8 s (2.9 GB, 2 thr) | 0.0 s (2.4 GB, 2 thr) | 0.1 s (2.4 GB, 2 thr) | 575.3 s (4.5 GB, 2 thr) |
| 1640 | 183.9 s (4.8 GB, 3 thr) | 255.5 s (4.0 GB, 2 thr) | 27.9 s (2.8 GB, 2 thr) | 0.0 s (2.5 GB, 2 thr) | 0.1 s (2.5 GB, 2 thr) | — |
| 1970 | 233.2 s (5.5 GB, 3 thr) | — | 22.0 s (3.5 GB, 3 thr) | 0.0 s (3.0 GB, 3 thr) | 0.1 s (3.0 GB, 3 thr) | — |
| 2100 | 256.1 s (5.8 GB, 3 thr) | — | 26.0 s (3.6 GB, 3 thr) | 0.0 s (3.2 GB, 3 thr) | 0.1 s (3.2 GB, 3 thr) | — |
| 2101 | 275.5 s (5.9 GB, 3 thr) | — | 24.5 s (3.7 GB, 3 thr) | 0.0 s (3.2 GB, 3 thr) | 0.1 s (3.2 GB, 3 thr) | 564.3 s (5.3 GB, 3 thr) |
| 2103 | 292.5 s (5.8 GB, 3 thr) | — | 24.3 s (3.7 GB, 3 thr) | 0.0 s (3.1 GB, 3 thr) | 0.1 s (3.1 GB, 3 thr) | — |


Timings are wall-clock under contention (≤ 4 P1 jobs × 2 threads, n ≤ 660; plus other agents' jobs on the same 8-core M2; load average 38–120 at times), so they are upper bounds; single-job timings (3 threads): docs/recon.md §6 (stack_lr_wr2 12–25 s at n ≤ 650).
Total fit+predict hours by predictor: P* stack_lr_wr2_lb 10.82 h, stack_lr_wr2 (pre-P2 P*) 9.35 h, stack_lr 1.20 h, B1 array rate 0.00 h, B2 product 0.01 h, m1_hgb (no genome) 3.61 h.

## Upstream-protocol sanity (dev split 1, n = 660 ≈ 20% of 3,283; R = 10 permutations; mean ± SD)

Upstream protocol = every measured cell of the held-out rows incl. cells whose reciprocal is known (array screened); here minus test-gene columns (dev rule) and self pairs. Upstream reference = YEAST_SLs REPORT §3 (S1 dev, 20% screened, ~2,600 held-out rows, seeds 0–4).

| set / stratum | P* stack_lr_wr2_lb ours | upstream | stack_lr_wr2 (pre-P2 P*) ours | upstream | stack_lr ours | upstream | B1 array rate ours | upstream |
|---|---|---|---|---|---|---|---|---|
| E_fixed_dev | 0.0553 ± 0.0022 | — | 0.0544 ± 0.0022 | — | 0.0450 ± 0.0007 | — | 0.0417 ± 0.0003 | — |
| all measured D_val cells (upstream) | 0.0663 ± 0.0028 | — | 0.0657 ± 0.0027 | 0.0754 (lift 6.6) | 0.0543 ± 0.0017 | 0.0590 (lift 5.1) | 0.0417 ± 0.0001 | 0.0430 (lift 3.7) |
| reciprocal available | 0.1003 ± 0.0101 | — | 0.0999 ± 0.0098 | 0.126 | 0.0959 ± 0.0101 | 0.117 | 0.0392 ± 0.0026 | 0.038 |
| on-array, no recip. | 0.0623 ± 0.0029 | — | 0.0613 ± 0.0028 | 0.073 | 0.0458 ± 0.0008 | 0.051 | 0.0405 ± 0.0006 | 0.043 |
| off-array | 0.0578 ± 0.0033 | — | 0.0579 ± 0.0034 | 0.058 | 0.0477 ± 0.0004 | 0.047 | 0.0479 ± 0.0006 | 0.048 |
| lift (upstream set) | 6.0267 ± 0.2566 | — | 5.9698 ± 0.2470 | — | 4.9330 ± 0.1584 | — | 3.7853 ± 0.0122 | — |
| prevalence (upstream set) | 0.0110 ± 0.0000 | — | 0.0110 ± 0.0000 | 0.0115 | 0.0110 ± 0.0000 | — | 0.0110 ± 0.0000 | — |

Ratio P*/stack_lr (paired by permutation): upstream-protocol set 1.210 ± 0.035, E_fixed_dev 1.209 ± 0.037; upstream dev 1.28 (REPORT §4.3).

Gap diagnosis (`upstream_protocol_diag.md`): our harness reproduces upstream's own dev split 0 (P* 0.0774 vs 0.0775); dev1's target rows are the hardest fold (upstream's own predictions on dev1 rows: 0.0673); dev2/dev3 upstream-protocol P* = 0.0796 / 0.0757 at n = 660.

Cells (perm 0): array_screened 252,198, recip_available 193,301, x_off_array 291,754, x_on_array_no_recip 1,068,137; set 1,553,192 cells.

## End-to-end controls (cited, not re-run)

`tests/test_harness.py::test_end_to_end_controls_dev1_n100` (slow; P0): shuffled screened labels (ε, P, labels permuted jointly over the screened block, array SMF permuted) → AUPRC within 0.7–1.3× prevalence; oracle back door → AUPRC > 0.99, both on dev1 E_fixed_dev through `harness.run_predictor`. P1 code paths use the same `run_predictor`; every P1 engine AP (unit weights) equals the stored harness AUPRC to 1e-9.

## Power pilot

See `power_pilot.md`: g = 0.0017, V_Q(Δ) = 4.46e-06, σ_r(160) = 0.0021 → V_Q(Δ) alone exceeds (g/3.5)^2: Q2 claims phrased as bounds (DESIGN §7.1).

