# P1 — INTERIM: dev random curves n ≤ 660 (R0), DEV DATA ONLY

DESIGN v3.1 §7.1/§8. Splits dev1–3 (D_val = inner fold k, candidates = D \ fold k); evaluation set E_fixed_dev = D_val × (A \ D_cand \ T) minus identity exclusions (no test-gene columns). Pooled AUPRC (tie-grouped AP, frozen `hsl.eval`). Trajectories: seeded permutations of the candidate pool (`dev_traj_base_seed`, harness R0 streams), nested prefixes; model seed = trajectory evaluation stream (at the pool maximum all permutations coincide, so the SD there is model-seed variation only). Rows: `results/p1/curves.csv`; per-(split, predictor, n) table: `results/p1/curve_table.csv`.

Figures: `fig_q1_dev_curves_interim.png` (AUPRC vs n, mean ± SD band), `fig_q1_dev_strata_interim.png` (strata).

**Interim (PI decision D6):** n ≥ 820, the pool maximum, m1_hgb at 1310/max and universe subsampling are on hold until P2 is done, so n_80/n_90 (which need the pool-max anchor) are not reported here.

## AUPRC at key budgets (mean ± SD over permutations; E_fixed_dev)

**dev1** (prevalence 0.0090, pool max 2101)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 |
|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2 | 10 | 0.0274 ± 0.0032 | 0.0353 ± 0.0012 | 0.0389 ± 0.0021 | 0.0414 ± 0.0024 | 0.0452 ± 0.0025 | 0.0498 ± 0.0027 | 0.0544 ± 0.0022 |
| stack_lr | 10 | 0.0279 ± 0.0032 | 0.0345 ± 0.0008 | 0.0376 ± 0.0007 | 0.0391 ± 0.0008 | 0.0412 ± 0.0005 | 0.0433 ± 0.0006 | 0.0450 ± 0.0007 |
| B1 array rate | 10 | 0.0236 ± 0.0011 | 0.0338 ± 0.0008 | 0.0372 ± 0.0007 | 0.0386 ± 0.0006 | 0.0400 ± 0.0004 | 0.0411 ± 0.0003 | 0.0417 ± 0.0003 |
| B2 product | 10 | 0.0193 ± 0.0009 | 0.0263 ± 0.0009 | 0.0290 ± 0.0008 | 0.0302 ± 0.0009 | 0.0318 ± 0.0006 | 0.0329 ± 0.0005 | 0.0336 ± 0.0005 |
| m1_hgb (no genome) | 5 | 0.0299 ± 0.0039 | — | 0.0395 ± 0.0020 | — | — | 0.0521 ± 0.0028 | 0.0561 ± 0.0021 |

**dev2** (prevalence 0.0098, pool max 2100)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 |
|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2 | 10 | 0.0298 ± 0.0020 | 0.0368 ± 0.0010 | 0.0409 ± 0.0011 | 0.0447 ± 0.0025 | 0.0505 ± 0.0039 | 0.0575 ± 0.0039 | 0.0638 ± 0.0035 |
| stack_lr | 10 | 0.0300 ± 0.0021 | 0.0365 ± 0.0011 | 0.0389 ± 0.0011 | 0.0409 ± 0.0012 | 0.0428 ± 0.0011 | 0.0453 ± 0.0012 | 0.0472 ± 0.0013 |
| B1 array rate | 10 | 0.0228 ± 0.0017 | 0.0335 ± 0.0012 | 0.0369 ± 0.0008 | 0.0384 ± 0.0005 | 0.0400 ± 0.0006 | 0.0413 ± 0.0004 | 0.0420 ± 0.0004 |
| B2 product | 10 | 0.0209 ± 0.0014 | 0.0287 ± 0.0011 | 0.0314 ± 0.0009 | 0.0333 ± 0.0010 | 0.0348 ± 0.0011 | 0.0361 ± 0.0009 | 0.0370 ± 0.0006 |

**dev3** (prevalence 0.0096, pool max 2103)

| predictor | R | 30 | 100 | 160 | 220 | 330 | 490 | 660 |
|---|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2 | 10 | 0.0310 ± 0.0019 | 0.0380 ± 0.0024 | 0.0435 ± 0.0023 | 0.0475 ± 0.0029 | 0.0526 ± 0.0026 | 0.0573 ± 0.0022 | 0.0610 ± 0.0025 |
| stack_lr | 10 | 0.0313 ± 0.0017 | 0.0369 ± 0.0020 | 0.0406 ± 0.0016 | 0.0428 ± 0.0018 | 0.0457 ± 0.0016 | 0.0482 ± 0.0013 | 0.0502 ± 0.0011 |
| B1 array rate | 10 | 0.0250 ± 0.0021 | 0.0342 ± 0.0010 | 0.0373 ± 0.0008 | 0.0387 ± 0.0007 | 0.0402 ± 0.0005 | 0.0412 ± 0.0004 | 0.0417 ± 0.0003 |
| B2 product | 10 | 0.0247 ± 0.0030 | 0.0338 ± 0.0028 | 0.0372 ± 0.0022 | 0.0394 ± 0.0023 | 0.0412 ± 0.0013 | 0.0429 ± 0.0005 | 0.0438 ± 0.0006 |

**Mean over splits** (mean of the per-split means; ± SD of the 3 split means; pool max = each split's)

| predictor | 30 | 100 | 160 | 220 | 330 | 490 | 660 |
|---|---|---|---|---|---|---|---|
| P* stack_lr_wr2 | 0.0294 ± 0.0018 | 0.0367 ± 0.0013 | 0.0411 ± 0.0023 | 0.0445 ± 0.0031 | 0.0494 ± 0.0038 | 0.0549 ± 0.0044 | 0.0597 ± 0.0048 |
| stack_lr | 0.0297 ± 0.0017 | 0.0360 ± 0.0013 | 0.0390 ± 0.0015 | 0.0409 ± 0.0018 | 0.0432 ± 0.0023 | 0.0456 ± 0.0024 | 0.0475 ± 0.0026 |
| B1 array rate | 0.0238 ± 0.0011 | 0.0338 ± 0.0004 | 0.0371 ± 0.0003 | 0.0386 ± 0.0001 | 0.0401 ± 0.0001 | 0.0412 ± 0.0001 | 0.0418 ± 0.0001 |
| B2 product | 0.0216 ± 0.0028 | 0.0296 ± 0.0038 | 0.0325 ± 0.0042 | 0.0343 ± 0.0047 | 0.0359 ± 0.0048 | 0.0373 ± 0.0051 | 0.0381 ± 0.0052 |

**P* vs baselines on the same trajectories** (mean over splits of the per-split mean relative gain; range over splits):

| n | P* vs stack_lr | P* vs B1 | P* vs B2_product |
|---|---|---|---|
| 30 | -1.2% [-1.7%, -0.7%] | +23.8% [+16.3%, +30.7%] | +37.2% [+26.3%, +42.9%] |
| 100 | +2.1% [+0.9%, +2.9%] | +8.6% [+4.5%, +11.1%] | +25.2% [+13.1%, +34.2%] |
| 160 | +5.3% [+3.5%, +7.0%] | +10.7% [+4.7%, +16.3%] | +27.2% [+17.0%, +34.1%] |
| 220 | +8.7% [+5.8%, +11.1%] | +15.4% [+7.1%, +22.9%] | +30.7% [+20.8%, +37.1%] |
| 330 | +14.2% [+9.8%, +17.9%] | +23.4% [+13.0%, +30.7%] | +38.4% [+27.6%, +45.3%] |
| 490 | +20.3% [+15.1%, +26.8%] | +33.2% [+21.2%, +39.3%] | +48.1% [+33.6%, +59.2%] |
| 660 | +25.8% [+20.9%, +35.1%] | +42.9% [+30.4%, +52.1%] | +57.9% [+39.3%, +72.6%] |

## Strata (P*, mean over permutations, split-averaged)

| n | all | E_AO | E_VV | on-array val rows | off-array val rows | lift | macro AUPRC | AUROC | prec@10% recall | SLs detected (frac of dev rows' SLs) |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 0.0294 | 0.0276 | 0.0335 | 0.0294 | 0.0345 | 3.10 | 0.0590 | 0.696 | 0.057 | 0.011 |
| 100 | 0.0367 | 0.0346 | 0.0414 | 0.0366 | 0.0392 | 3.88 | 0.0708 | 0.731 | 0.069 | 0.037 |
| 160 | 0.0411 | 0.0385 | 0.0465 | 0.0413 | 0.0419 | 4.34 | 0.0774 | 0.741 | 0.077 | 0.060 |
| 220 | 0.0445 | 0.0416 | 0.0504 | 0.0453 | 0.0437 | 4.70 | 0.0821 | 0.748 | 0.084 | 0.083 |
| 330 | 0.0494 | 0.0459 | 0.0562 | 0.0510 | 0.0467 | 5.22 | 0.0886 | 0.756 | 0.095 | 0.125 |
| 490 | 0.0549 | 0.0506 | 0.0631 | 0.0572 | 0.0495 | 5.79 | 0.0955 | 0.763 | 0.106 | 0.185 |
| 660 | 0.0597 | 0.0549 | 0.0688 | 0.0628 | 0.0520 | 6.30 | 0.0997 | 0.769 | 0.117 | 0.250 |


## Runtime per fit (fit+predict s, mean over splits/permutations; peak RSS sampled during the fit, max)

| n | P* stack_lr_wr2 | stack_lr | B1 array rate | B2 product | m1_hgb (no genome) |
|---|---|---|---|---|---|
| 10 | 32.2 s (1.2 GB, 2 thr) | 0.4 s (1.2 GB, 2 thr) | 0.0 s (1.2 GB, 2 thr) | 0.0 s (1.2 GB, 2 thr) | — |
| 30 | 32.0 s (1.6 GB, 2 thr) | 0.6 s (1.4 GB, 2 thr) | 0.0 s (1.4 GB, 2 thr) | 0.0 s (1.4 GB, 2 thr) | 236.4 s (1.9 GB, 2 thr) |
| 100 | 34.0 s (1.8 GB, 2 thr) | 1.3 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | — |
| 160 | 37.2 s (1.9 GB, 2 thr) | 2.1 s (1.3 GB, 2 thr) | 0.0 s (1.1 GB, 2 thr) | 0.0 s (1.1 GB, 2 thr) | 394.9 s (1.9 GB, 2 thr) |
| 330 | 47.5 s (2.0 GB, 2 thr) | 4.2 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | — |
| 660 | 76.0 s (2.7 GB, 2 thr) | 9.1 s (1.7 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 0.0 s (1.6 GB, 2 thr) | 424.5 s (2.7 GB, 2 thr) |


Timings are wall-clock under contention (≤ 4 P1 jobs × 2 threads, n ≤ 660; plus other agents' jobs on the same 8-core M2; load average 38–120 at times), so they are upper bounds; single-job timings (3 threads): docs/recon.md §6 (stack_lr_wr2 12–25 s at n ≤ 650).
Total fit+predict hours by predictor: P* stack_lr_wr2 4.32 h, stack_lr 0.31 h, B1 array rate 0.00 h, B2 product 0.00 h, m1_hgb (no genome) 2.03 h.

## Upstream-protocol sanity (dev split 1, n = 660 ≈ 20% of 3,283; R = 10 permutations; mean ± SD)

Upstream protocol = every measured cell of the held-out rows incl. cells whose reciprocal is known (array screened); here minus test-gene columns (dev rule) and self pairs. Upstream reference = YEAST_SLs REPORT §3 (S1 dev, 20% screened, ~2,600 held-out rows, seeds 0–4).

| set / stratum | P* stack_lr_wr2 ours | upstream | stack_lr ours | upstream | B1 array rate ours | upstream |
|---|---|---|---|---|---|---|
| E_fixed_dev | 0.0544 ± 0.0022 | — | 0.0450 ± 0.0007 | — | 0.0417 ± 0.0003 | — |
| all measured D_val cells (upstream) | 0.0657 ± 0.0027 | 0.0754 (lift 6.6) | 0.0543 ± 0.0017 | 0.0590 (lift 5.1) | 0.0417 ± 0.0001 | 0.0430 (lift 3.7) |
| reciprocal available | 0.0999 ± 0.0098 | 0.126 | 0.0959 ± 0.0101 | 0.117 | 0.0392 ± 0.0026 | 0.038 |
| on-array, no recip. | 0.0613 ± 0.0028 | 0.073 | 0.0458 ± 0.0008 | 0.051 | 0.0405 ± 0.0006 | 0.043 |
| off-array | 0.0579 ± 0.0034 | 0.058 | 0.0477 ± 0.0004 | 0.047 | 0.0479 ± 0.0006 | 0.048 |
| lift (upstream set) | 5.9698 ± 0.2470 | — | 4.9330 ± 0.1584 | — | 3.7853 ± 0.0122 | — |
| prevalence (upstream set) | 0.0110 ± 0.0000 | 0.0115 | 0.0110 ± 0.0000 | — | 0.0110 ± 0.0000 | — |

Ratio P*/stack_lr (paired by permutation): upstream-protocol set 1.210 ± 0.035, E_fixed_dev 1.209 ± 0.037; upstream dev 1.28 (REPORT §4.3).

Gap diagnosis (`upstream_protocol_diag.md`): our harness reproduces upstream's own dev split 0 (P* 0.0774 vs 0.0775); dev1's target rows are the hardest fold (upstream's own predictions on dev1 rows: 0.0673); dev2/dev3 upstream-protocol P* = 0.0796 / 0.0757 at n = 660.

Cells (perm 0): array_screened 252,198, recip_available 193,301, x_off_array 291,754, x_on_array_no_recip 1,068,137; set 1,553,192 cells.

## End-to-end controls (cited, not re-run)

`tests/test_harness.py::test_end_to_end_controls_dev1_n100` (slow; P0): shuffled screened labels (ε, P, labels permuted jointly over the screened block, array SMF permuted) → AUPRC within 0.7–1.3× prevalence; oracle back door → AUPRC > 0.99, both on dev1 E_fixed_dev through `harness.run_predictor`. P1 code paths use the same `run_predictor`; every P1 engine AP (unit weights) equals the stored harness AUPRC to 1e-9.

## Power pilot

See `power_pilot.md`: g = 0.0021, V_Q(Δ) = 3.94e-06, σ_r(160) = 0.0021 → V_Q(Δ) alone exceeds (g/3.5)^2: Q2 claims phrased as bounds (DESIGN §7.1).

