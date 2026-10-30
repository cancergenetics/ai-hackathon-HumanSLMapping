# Cost projection for the REAL final (plan of mode `final`)

Inputs: per-unit wall time (fit + predict + E_fixed scoring + save) = median of measured runs per (predictor, n): dry-run rows on dev1 (this machine, 1 job, 2 threads, shared with a running batch) and P1 random-curve rows (up to 4 concurrent jobs, loaded machine), interpolated in log n; above the largest measured n: power-law fit (time) / linear fit (RSS) over n >= max/3 (flagged). Scaled to |T| = 656 targets per predictor (pstar ×1.00, ref ×1.00, m1_hgb ×1.08, stack_lr ×1.00, B1_array_rate ×1.00, B2_product ×1.00; from same-day calibration pairs, re-split 1050 vs dev1 526 targets, fit+predict ratios pstar 0.94, m1_hgb 1.32; predictors without a pair use P*'s ratio) for the larger prediction grid. Selection per trajectory = dry-run median by arm role (full trajectories to 660): batch 6 s, det 20 s, ws 8 s, finalist 15 s, static 9 s, ws_random 0 s, random 0 s, R1 13 s. Per-job overhead 15 s. Job memory = the launch.py `rss_gb` hint (finalists.yaml compute: fitted on measured peak RSS + calibrated per-target term + 0.5 GB), the same number the gate uses.

Plan: 4019 units in 393 jobs (+ labels job).

Dry-run / P1 speed ratios per predictor: B1_array_rate 0.72, B2_product 0.74, m1_hgb 0.70, pstar 0.77, ref 0.49, stack_lr 0.71.

## Units and job-hours by arm role × predictor (lower / upper bound)

| role | predictor | units | job-h (1-job speed) | job-h (loaded) |
|---|---|---|---|---|
| R1 | pstar | 180 | 2.3 | 3.0 |
| batch | pstar | 240 | 3.2 | 4.2 |
| det | pstar | 24 | 0.3 | 0.4 |
| finalist | m1_hgb | 45 | 3.1 | 4.4 |
| finalist | pstar | 480 | 6.0 | 7.8 |
| finalist | ref | 240 | 1.5 | 3.1 |
| random | B1_array_rate | 535 | 0.1 | 0.1 |
| random | B2_product | 535 | 0.1 | 0.1 |
| random | m1_hgb | 30 | 2.5 | 3.6 |
| random | pstar | 535 | 7.8 | 10.2 |
| random | ref | 80 | 0.5 | 1.0 |
| random | stack_lr | 535 | 0.6 | 0.9 |
| static | pstar | 320 | 4.0 | 5.2 |
| ws | pstar | 160 | 2.5 | 3.2 |
| ws_random | pstar | 80 | 1.2 | 1.6 |

Totals in units: 35.8 (1-job speed) to 48.7 (loaded) job-h, + 1.6 h job overhead + selection 0.9 h.

## Wall-clock (event simulation of the launch.py queue in plan order)

`gate` = D10 memory gate: ≤ 6 jobs, ≤ 2 jobs at n ≥ 1970, start only if 18 GB (memory available to us: PI check 2026-10-29, free + inactive ≈ 15 of 24 GB, other non-project use ~0.7 GB; 12 / 22 GB shown as sensitivities) minus the running jobs' hints ≥ hint + 2 GB. `old` = previous policy for comparison (4 jobs, n ≥ 1970 alone, 20 GB). CPU contention at 6 jobs is not modelled beyond the loaded-speed bound (P1 rows ran up to 4 jobs on a shared machine).

| scenario | jobs | job-h | wall h, gate 18 GB (mean jobs running) | gate 12 GB | gate 22 GB | old, 4 jobs | > 10 h at 18 GB? |
|---|---|---|---|---|---|---|---|
| 1-job speed (lower bound): D10 plan as configured | 393 | 38.3 | 10.2 (3.7) | 20.2 | 8.3 | 11.0 | YES |
| 1-job speed (lower bound): without the reference predictor | 353 | 36.0 | 9.8 (3.7) | 19.4 | 7.9 | 10.4 | no |
| 1-job speed (lower bound): §12 trims: m1_hgb Q1 {30, 1310} and batch 50 dropped | 348 | 35.5 | 9.5 (3.7) | 18.8 | 7.8 | 10.3 | no |
| 1-job speed (lower bound): both of the above | 308 | 33.2 | 9.1 (3.7) | 18.0 | 7.4 | 9.7 | no |
| P1 loaded speed (upper bound): D10 plan as configured | 393 | 51.2 | 13.6 (3.8) | 26.7 | 11.0 | 14.7 | YES |
| P1 loaded speed (upper bound): without the reference predictor | 353 | 46.8 | 12.8 (3.7) | 25.3 | 10.3 | 13.6 | YES |
| P1 loaded speed (upper bound): §12 trims: m1_hgb Q1 {30, 1310} and batch 50 dropped | 348 | 47.6 | 12.7 (3.7) | 24.9 | 10.3 | 13.8 | YES |
| P1 loaded speed (upper bound): both of the above | 308 | 43.2 | 11.8 (3.7) | 23.5 | 9.6 | 12.7 | YES |

## Calibration runs (run_final.py calibrate; DEV ONLY; one job at a time, 2 threads)

| predictor | split | targets | n | peak RSS GB | fit+predict s | hint for T at this n |
|---|---|---|---|---|---|---|
| pstar | resplit | 1050 | 660 | 2.58 | 36 | 3.53 |
| pstar | dev1 | 526 | 660 | 2.56 | 37 | 3.53 |
| pstar | resplit | 1050 | 1577 | 4.45 | 117 | 5.34 |
| pstar | dev1 | 526 | 1577 | 4.24 | 128 | 5.34 |
| ref | resplit | 1050 | 660 | 2.22 | 37 | 3.11 |
| ref | resplit | 1050 | 1577 | 3.85 | 124 | 4.59 |
| stack_lr | resplit | 1050 | 660 | 1.78 | 7 | 2.53 |
| stack_lr | resplit | 1050 | 1577 | 2.68 | 12 | 3.55 |
| m1_hgb | resplit | 1050 | 660 | 2.94 | 291 | 3.54 |
| m1_hgb | dev1 | 526 | 660 | 3.04 | 221 | 3.54 |
| m1_hgb | dev1 | 526 | 2101 | 5.21 | 703 | 6.25 |

## Measured per-unit times used (s, before × T scaling) and peak RSS (GB)

| source | predictor | n | median s | × T | max RSS | k |
|---|---|---|---|---|---|---|
| dryrun | m1_hgb | 160 | 180.7 | 194.8 | 2.33 | 5 |
| dryrun | m1_hgb | 490 | 226.3 | 244.0 | 2.58 | 2 |
| dryrun | m1_hgb | 660 | 222.4 | 239.9 | 3.04 | 1 |
| dryrun | m1_hgb | 2101 | 703.9 | 759.0 | 5.21 | 1 |
| dryrun | pstar | 30 | 22.2 | 22.2 | 1.84 | 18 |
| dryrun | pstar | 160 | 20.6 | 20.6 | 2.05 | 26 |
| dryrun | pstar | 330 | 26.1 | 26.1 | 2.03 | 6 |
| dryrun | pstar | 440 | 27.1 | 27.1 | 2.66 | 8 |
| dryrun | pstar | 490 | 34.9 | 34.9 | 2.85 | 18 |
| dryrun | pstar | 660 | 37.0 | 37.0 | 3.37 | 33 |
| dryrun | pstar | 980 | 83.2 | 83.2 | 3.49 | 3 |
| dryrun | pstar | 1577 | 128.9 | 128.9 | 4.24 | 1 |
| dryrun | ref | 160 | 18.9 | 18.9 | 1.50 | 8 |
| dryrun | ref | 660 | 34.8 | 34.8 | 2.50 | 8 |
| dryrun | stack_lr | 30 | 1.1 | 1.1 | 1.33 | 3 |
| dryrun | stack_lr | 160 | 2.0 | 2.0 | 1.32 | 3 |
| dryrun | stack_lr | 490 | 4.9 | 4.9 | 1.71 | 3 |
| dryrun | stack_lr | 660 | 7.1 | 7.1 | 1.98 | 3 |
| dryrun | stack_lr | 980 | 11.3 | 11.3 | 2.27 | 3 |
| p1 | m1_hgb | 30 | 242.9 | 261.9 | 1.91 | 5 |
| p1 | m1_hgb | 160 | 334.6 | 360.8 | 1.91 | 5 |
| p1 | m1_hgb | 490 | 314.8 | 339.4 | 2.43 | 5 |
| p1 | m1_hgb | 660 | 323.5 | 348.8 | 2.68 | 5 |
| p1 | m1_hgb | 1310 | 577.2 | 622.4 | 4.46 | 5 |
| p1 | m1_hgb | 2101 | 594.1 | 640.6 | 5.27 | 5 |
| p1 | pstar | 10 | 52.4 | 52.4 | 1.76 | 30 |
| p1 | pstar | 60 | 45.9 | 45.9 | 1.78 | 30 |
| p1 | pstar | 100 | 48.1 | 48.1 | 1.85 | 30 |
| p1 | pstar | 120 | 50.2 | 50.2 | 1.82 | 30 |
| p1 | pstar | 210 | 56.1 | 56.1 | 1.90 | 30 |
| p1 | pstar | 220 | 54.8 | 54.8 | 1.99 | 30 |
| p1 | pstar | 410 | 71.1 | 71.1 | 2.30 | 30 |
| p1 | pstar | 490 | 79.4 | 79.4 | 2.33 | 30 |
| p1 | pstar | 820 | 65.0 | 65.0 | 3.22 | 30 |
| p1 | pstar | 980 | 75.6 | 75.6 | 3.75 | 30 |
| p1 | pstar | 1310 | 122.0 | 122.0 | 4.42 | 30 |
| p1 | pstar | 1640 | 182.0 | 182.0 | 4.82 | 30 |
| p1 | pstar | 1970 | 241.8 | 241.8 | 5.50 | 30 |
| p1 | pstar | 2100 | 244.8 | 244.8 | 5.80 | 5 |
| p1 | pstar | 2101 | 246.7 | 246.7 | 5.90 | 5 |
| p1 | pstar | 2103 | 306.9 | 306.9 | 5.79 | 5 |
| p1 | ref | 10 | 32.8 | 32.8 | 1.18 | 30 |
| p1 | ref | 30 | 30.3 | 30.3 | 1.63 | 30 |
| p1 | ref | 60 | 29.3 | 29.3 | 1.68 | 30 |
| p1 | ref | 100 | 33.7 | 33.7 | 1.83 | 30 |
| p1 | ref | 120 | 32.8 | 32.8 | 1.68 | 30 |
| p1 | ref | 160 | 32.1 | 32.1 | 1.92 | 30 |
| p1 | ref | 210 | 38.5 | 38.5 | 1.83 | 30 |
| p1 | ref | 220 | 40.3 | 40.3 | 1.68 | 30 |
| p1 | ref | 330 | 51.9 | 51.9 | 2.00 | 30 |
| p1 | ref | 410 | 55.9 | 55.9 | 2.22 | 30 |
| p1 | ref | 490 | 61.1 | 61.1 | 2.20 | 30 |
| p1 | ref | 660 | 88.4 | 88.4 | 2.74 | 30 |
| p1 | ref | 820 | 71.1 | 71.1 | 2.71 | 30 |
| p1 | ref | 980 | 88.8 | 88.8 | 3.28 | 30 |
| p1 | ref | 1310 | 146.3 | 146.3 | 3.66 | 30 |
| p1 | ref | 1640 | 225.4 | 225.4 | 4.00 | 30 |
| p1 | stack_lr | 10 | 1.2 | 1.2 | 1.19 | 30 |
| p1 | stack_lr | 30 | 1.6 | 1.6 | 1.36 | 30 |
| p1 | stack_lr | 60 | 1.9 | 1.9 | 1.60 | 30 |
| p1 | stack_lr | 100 | 2.4 | 2.4 | 1.61 | 30 |
| p1 | stack_lr | 120 | 2.5 | 2.5 | 1.59 | 30 |
| p1 | stack_lr | 160 | 2.8 | 2.8 | 1.28 | 30 |
| p1 | stack_lr | 210 | 3.9 | 3.9 | 1.36 | 30 |
| p1 | stack_lr | 220 | 3.6 | 3.6 | 1.36 | 30 |
| p1 | stack_lr | 330 | 5.8 | 5.8 | 1.60 | 30 |
| p1 | stack_lr | 410 | 7.1 | 7.1 | 1.57 | 30 |
| p1 | stack_lr | 490 | 8.2 | 8.2 | 1.71 | 30 |
| p1 | stack_lr | 660 | 11.2 | 11.2 | 1.75 | 30 |
| p1 | stack_lr | 820 | 10.4 | 10.4 | 1.93 | 30 |
| p1 | stack_lr | 980 | 12.1 | 12.1 | 2.14 | 30 |
| p1 | stack_lr | 1310 | 19.2 | 19.2 | 2.88 | 30 |
| p1 | stack_lr | 1640 | 25.8 | 25.8 | 2.83 | 30 |
| p1 | stack_lr | 1970 | 21.3 | 21.3 | 3.45 | 30 |
| p1 | stack_lr | 2100 | 28.2 | 28.2 | 3.62 | 5 |
| p1 | stack_lr | 2101 | 25.9 | 25.9 | 3.70 | 5 |
| p1 | stack_lr | 2103 | 21.6 | 21.6 | 3.68 | 5 |

Largest job memory hint (pool max; RSS extrapolated beyond n = 2103): 7.4 GB.
Score step (separate, 1 process, no bundle): ≈ 5219 bootstrap vectors + 200 two-way vectors ≈ 2.4 h (measured dry-run engine speed scaled to |E_fixed|).

Extrapolated (no measurement at this n): B1_array_rate@2627, B2_product@2627, m1_hgb@2627, pstar@2627, stack_lr@2627

