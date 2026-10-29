# P1 power pilot (dev split 1, R0 random, predictor P* stack_lr_wr2_lb)

R = 10 permutations; B = 2000 stratified query-bootstrap replicates (frozen ReplicateEngine, seed `bootstrap_seed` child 'p1_power_pilot_dev1'); |D_val| = 526, |T| = 656. DESIGN §7.1.

| n | R | mean AUPRC | σ_r(n) | V_Q,single(n) (traj-mean) | SD_Q single traj |
|---|---|---|---|---|---|
| 10 | 10 | 0.0253 | 0.0016 | 5.18e-07 (SD 0.0007) | 0.0010 |
| 30 | 10 | 0.0288 | 0.0036 | 8.22e-07 (SD 0.0009) | 0.0011 |
| 60 | 10 | 0.0343 | 0.0008 | 1.28e-06 (SD 0.0011) | 0.0014 |
| 100 | 10 | 0.0385 | 0.0014 | 1.97e-06 (SD 0.0014) | 0.0016 |
| 120 | 10 | 0.0394 | 0.0015 | 2.07e-06 (SD 0.0014) | 0.0017 |
| 160 | 10 | 0.0416 | 0.0021 | 2.78e-06 (SD 0.0017) | 0.0021 |
| 210 | 10 | 0.0433 | 0.0024 | 3.32e-06 (SD 0.0018) | 0.0022 |
| 220 | 10 | 0.0435 | 0.0024 | 3.45e-06 (SD 0.0019) | 0.0023 |
| 330 | 10 | 0.0467 | 0.0027 | 5.54e-06 (SD 0.0024) | 0.0028 |
| 410 | 10 | 0.0492 | 0.0030 | 8.06e-06 (SD 0.0028) | 0.0032 |
| 490 | 10 | 0.0509 | 0.0028 | 9.73e-06 (SD 0.0031) | 0.0034 |
| 660 | 10 | 0.0553 | 0.0022 | 1.59e-05 (SD 0.0040) | 0.0043 |

- g = r210 − r160 = 0.0433 − 0.0416 = **0.0017** (paired SE 0.0004, 95% CI 0.0009 to 0.0025)
- V_Q(Δ) := 2·V_Q,single(160)·|D_val|/|T| = 2·2.78e-06·526/656 = **4.46e-06** (SD 0.0021)
- σ_r(160) = 0.0021; target SE(Δ̂(160)) ≤ g/3.5 = 0.0005 (SE² 2.3e-07)
- **V_Q(Δ) alone exceeds (g/3.5)^2: Q2 claims phrased as bounds (DESIGN §7.1)**
- SE(Δ̂(160)) by R_s (R_r = 2R_s): R_s=3: 0.0026, R_s=5: 0.0024, R_s=10: 0.0023, R_s=20: 0.0022, R_s=40: 0.0022

Sensitivity (NOT the frozen §7.1 rule): strategy and random arms are scored on the same queries with shared bootstrap weights, so the real V_Q(Δ) is reduced by their replicate correlation. Proxy = two random budgets on the same queries: ρ_boot(r160, r210) = 0.993; V_Q of the paired contrast (scaled to |T|) = 5.32e-08 (SD 0.0002); 2·V_Q,single·(1−ρ)·|D_val|/|T| = 3.1e-08 → R_s = 39 (paired proxy) / 35 (ρ formula), R_r = 2R_s. A strategy whose selections differ from random decorrelates more than two nested random prefixes do, so this is optimistic; P4 re-estimates V_Q(Δ) from real strategy-vs-random replicates.

Caveats: σ_s := σ_r (strategy trajectory variance unknown until P4); V_Q(Δ) assumes ρ = 0 (conservative); g is estimated from 10 nested-prefix trajectories (its CI above) — a smaller true g raises R steeply (R ∝ σ²/((g/3.5)² − V_Q)). R is re-estimated from P4 dev strategy runs and fixed at the freeze (may only go up).
