# Cost projection for the REAL final (plan of mode `final`)

Inputs: per-unit wall time (fit + predict + E_fixed scoring + save) = median of measured runs per (predictor, n): dry-run rows on dev1 (this machine, 1 job, 2 threads, shared with a running batch) and P1 random-curve rows (up to 4 concurrent jobs, loaded machine), interpolated in log n; above the largest measured n: power-law fit (time) / linear fit (RSS) over n >= max/3 (flagged). Scaled × 1.25 (= |T| 656 / |dev1 targets| 526) for the larger prediction grid. Selection per trajectory = dry-run median by arm role (full trajectories to 660): batch 6 s, det 20 s, ws 8 s, finalist 15 s, static 9 s, ws_random 0 s, random 0 s, R1 13 s. Per-job overhead 15 s. Peak RSS per job = max over its units (+1 GB process; +10% if extrapolated).

Plan: 4619 units in 453 jobs (+ labels job).

Dry-run / P1 speed ratios per predictor: B1_array_rate 0.72, B2_product 0.74, m1_hgb 0.63, pstar 0.77, ref 0.49, stack_lr 0.71.

## Units and job-hours by arm role × predictor (lower / upper bound)

| role | predictor | units | job-h (1-job speed) | job-h (loaded) |
|---|---|---|---|---|
| R1 | pstar | 180 | 2.9 | 3.7 |
| batch | pstar | 240 | 4.0 | 5.3 |
| det | pstar | 24 | 0.4 | 0.5 |
| finalist | m1_hgb | 45 | 3.2 | 5.1 |
| finalist | pstar | 480 | 7.5 | 9.7 |
| finalist | ref | 240 | 1.9 | 3.8 |
| random | B1_array_rate | 685 | 0.1 | 0.2 |
| random | B2_product | 685 | 0.1 | 0.2 |
| random | m1_hgb | 30 | 2.9 | 4.6 |
| random | pstar | 685 | 15.3 | 19.8 |
| random | ref | 80 | 0.6 | 1.3 |
| random | stack_lr | 685 | 1.4 | 2.0 |
| static | pstar | 320 | 5.0 | 6.5 |
| ws | pstar | 160 | 3.1 | 4.0 |
| ws_random | pstar | 80 | 1.5 | 2.0 |

Totals in units: 50.0 (1-job speed) to 68.7 (loaded) job-h, + 1.9 h job overhead + selection 0.9 h.

## Wall-clock (event simulation of the launch.py queue; 20 GB RAM cap)

`alone` = protocol rule (n ≥ 1970 runs alone); `ram` = n ≥ 1970 jobs share the machine within the RAM cap (as P1 batch6 did, D-log 2026-10-29).

| scenario | jobs | job-h | wall h, 3 jobs (alone) | 3 jobs (ram) | 4 jobs (alone) | 4 jobs (ram) | > 10 h? |
|---|---|---|---|---|---|---|---|
| 1-job speed (lower bound): as configured | 453 | 52.7 | 20.6 | 18.0 | 17.1 | 14.5 | YES |
| 1-job speed (lower bound): random n > 660 at R = 10 | 393 | 46.3 | 17.0 | 15.8 | 13.5 | 12.3 | YES |
| 1-job speed (lower bound): random n > 660 at R = 5 | 383 | 45.2 | 16.4 | 15.4 | 12.9 | 11.9 | YES |
| 1-job speed (lower bound): random n > 660 at R = 5 and no reference predictor | 343 | 42.4 | 15.4 | 14.5 | 12.2 | 11.2 | YES |
| P1 loaded speed (upper bound): as configured | 453 | 71.4 | 28.0 | 24.4 | 23.2 | 19.7 | YES |
| P1 loaded speed (upper bound): random n > 660 at R = 10 | 393 | 63.0 | 23.3 | 21.6 | 18.5 | 16.9 | YES |
| P1 loaded speed (upper bound): random n > 660 at R = 5 | 383 | 61.6 | 22.4 | 21.1 | 17.7 | 16.4 | YES |
| P1 loaded speed (upper bound): random n > 660 at R = 5 and no reference predictor | 343 | 56.2 | 20.6 | 19.3 | 16.4 | 15.0 | YES |

## Measured per-unit times used (s, before × T scaling) and peak RSS (GB)

| source | predictor | n | median s | × T | max RSS | k |
|---|---|---|---|---|---|---|
| dryrun | m1_hgb | 160 | 180.7 | 225.4 | 2.33 | 5 |
| dryrun | m1_hgb | 490 | 226.3 | 282.2 | 2.58 | 2 |
| dryrun | pstar | 30 | 22.2 | 27.7 | 1.84 | 18 |
| dryrun | pstar | 160 | 20.6 | 25.7 | 2.05 | 26 |
| dryrun | pstar | 330 | 26.1 | 32.5 | 2.03 | 6 |
| dryrun | pstar | 440 | 27.1 | 33.8 | 2.66 | 8 |
| dryrun | pstar | 490 | 34.9 | 43.5 | 2.85 | 18 |
| dryrun | pstar | 660 | 36.8 | 45.9 | 3.37 | 32 |
| dryrun | pstar | 980 | 83.2 | 103.8 | 3.49 | 3 |
| dryrun | ref | 160 | 18.9 | 23.5 | 1.50 | 8 |
| dryrun | ref | 660 | 34.8 | 43.4 | 2.50 | 8 |
| dryrun | stack_lr | 30 | 1.1 | 1.4 | 1.33 | 3 |
| dryrun | stack_lr | 160 | 2.0 | 2.6 | 1.32 | 3 |
| dryrun | stack_lr | 490 | 4.9 | 6.1 | 1.71 | 3 |
| dryrun | stack_lr | 660 | 7.1 | 8.9 | 1.98 | 3 |
| dryrun | stack_lr | 980 | 11.3 | 14.1 | 2.27 | 3 |
| p1 | m1_hgb | 30 | 242.9 | 303.0 | 1.91 | 5 |
| p1 | m1_hgb | 160 | 334.6 | 417.3 | 1.91 | 5 |
| p1 | m1_hgb | 490 | 314.8 | 392.6 | 2.43 | 5 |
| p1 | m1_hgb | 660 | 323.5 | 403.4 | 2.68 | 5 |
| p1 | m1_hgb | 1310 | 577.2 | 719.9 | 4.46 | 5 |
| p1 | pstar | 10 | 52.4 | 65.4 | 1.76 | 30 |
| p1 | pstar | 60 | 45.9 | 57.2 | 1.78 | 30 |
| p1 | pstar | 100 | 48.1 | 60.0 | 1.85 | 30 |
| p1 | pstar | 120 | 50.2 | 62.6 | 1.82 | 30 |
| p1 | pstar | 210 | 56.1 | 69.9 | 1.90 | 30 |
| p1 | pstar | 220 | 54.8 | 68.4 | 1.99 | 30 |
| p1 | pstar | 410 | 71.1 | 88.7 | 2.30 | 30 |
| p1 | pstar | 490 | 79.4 | 99.1 | 2.33 | 30 |
| p1 | pstar | 820 | 65.0 | 81.0 | 3.22 | 30 |
| p1 | pstar | 980 | 75.6 | 94.3 | 3.75 | 30 |
| p1 | pstar | 1310 | 122.0 | 152.2 | 4.42 | 30 |
| p1 | pstar | 1640 | 182.0 | 227.0 | 4.82 | 30 |
| p1 | pstar | 1970 | 241.9 | 301.6 | 5.50 | 29 |
| p1 | pstar | 2100 | 244.8 | 305.3 | 5.80 | 5 |
| p1 | pstar | 2101 | 246.7 | 307.7 | 5.90 | 5 |
| p1 | pstar | 2103 | 306.9 | 382.7 | 5.79 | 5 |
| p1 | ref | 10 | 32.8 | 40.9 | 1.18 | 30 |
| p1 | ref | 30 | 30.3 | 37.8 | 1.63 | 30 |
| p1 | ref | 60 | 29.3 | 36.6 | 1.68 | 30 |
| p1 | ref | 100 | 33.7 | 42.0 | 1.83 | 30 |
| p1 | ref | 120 | 32.8 | 40.9 | 1.68 | 30 |
| p1 | ref | 160 | 32.1 | 40.0 | 1.92 | 30 |
| p1 | ref | 210 | 38.5 | 48.0 | 1.83 | 30 |
| p1 | ref | 220 | 40.3 | 50.2 | 1.68 | 30 |
| p1 | ref | 330 | 51.9 | 64.8 | 2.00 | 30 |
| p1 | ref | 410 | 55.9 | 69.7 | 2.22 | 30 |
| p1 | ref | 490 | 61.1 | 76.1 | 2.20 | 30 |
| p1 | ref | 660 | 88.4 | 110.2 | 2.74 | 30 |
| p1 | ref | 820 | 71.1 | 88.7 | 2.71 | 30 |
| p1 | ref | 980 | 88.8 | 110.8 | 3.28 | 30 |
| p1 | ref | 1310 | 146.3 | 182.5 | 3.66 | 30 |
| p1 | ref | 1640 | 225.4 | 281.1 | 4.00 | 30 |
| p1 | stack_lr | 10 | 1.2 | 1.5 | 1.19 | 30 |
| p1 | stack_lr | 30 | 1.6 | 2.0 | 1.36 | 30 |
| p1 | stack_lr | 60 | 1.9 | 2.3 | 1.60 | 30 |
| p1 | stack_lr | 100 | 2.4 | 3.0 | 1.61 | 30 |
| p1 | stack_lr | 120 | 2.5 | 3.2 | 1.59 | 30 |
| p1 | stack_lr | 160 | 2.8 | 3.5 | 1.28 | 30 |
| p1 | stack_lr | 210 | 3.9 | 4.8 | 1.36 | 30 |
| p1 | stack_lr | 220 | 3.6 | 4.5 | 1.36 | 30 |
| p1 | stack_lr | 330 | 5.8 | 7.2 | 1.60 | 30 |
| p1 | stack_lr | 410 | 7.1 | 8.8 | 1.57 | 30 |
| p1 | stack_lr | 490 | 8.2 | 10.3 | 1.71 | 30 |
| p1 | stack_lr | 660 | 11.2 | 14.0 | 1.75 | 30 |
| p1 | stack_lr | 820 | 10.4 | 13.0 | 1.93 | 30 |
| p1 | stack_lr | 980 | 12.1 | 15.0 | 2.14 | 30 |
| p1 | stack_lr | 1310 | 19.2 | 23.9 | 2.88 | 30 |
| p1 | stack_lr | 1640 | 25.8 | 32.2 | 2.83 | 30 |
| p1 | stack_lr | 1970 | 21.5 | 26.8 | 3.45 | 28 |
| p1 | stack_lr | 2100 | 28.2 | 35.2 | 3.62 | 5 |
| p1 | stack_lr | 2101 | 25.9 | 32.3 | 3.70 | 5 |
| p1 | stack_lr | 2103 | 21.6 | 27.0 | 3.68 | 5 |

Largest projected job RSS (pool max, extrapolated): 9.6 GB.
Score step (separate, 1 process, no bundle): ≈ 5819 bootstrap vectors + 200 two-way vectors ≈ 2.6 h (measured dry-run engine speed scaled to |E_fixed|).

Extrapolated (no measurement at this n): B1_array_rate@2627, B2_product@2627, m1_hgb@2627, pstar@2627, stack_lr@2627

