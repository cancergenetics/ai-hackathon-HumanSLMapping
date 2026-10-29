# P1 power pilot (dev split 1, R0 random, predictor stack_lr_wr2 (pre-P2 P*))

R = 10 permutations; B = 2000 stratified query-bootstrap replicates (frozen ReplicateEngine, seed `bootstrap_seed` child 'p1_power_pilot_dev1'); |D_val| = 526, |T| = 656. DESIGN §7.1.

| n | R | mean AUPRC | σ_r(n) | V_Q,single(n) (traj-mean) | SD_Q single traj |
|---|---|---|---|---|---|
| 10 | 10 | 0.0252 | 0.0018 | 5.27e-07 (SD 0.0007) | 0.0010 |
| 30 | 10 | 0.0274 | 0.0032 | 7.57e-07 (SD 0.0009) | 0.0011 |
| 60 | 10 | 0.0318 | 0.0007 | 1.07e-06 (SD 0.0010) | 0.0012 |
| 100 | 10 | 0.0353 | 0.0012 | 1.62e-06 (SD 0.0013) | 0.0015 |
| 120 | 10 | 0.0364 | 0.0013 | 1.72e-06 (SD 0.0013) | 0.0015 |
| 160 | 10 | 0.0389 | 0.0021 | 2.46e-06 (SD 0.0016) | 0.0019 |
| 210 | 10 | 0.0411 | 0.0024 | 3.05e-06 (SD 0.0017) | 0.0021 |
| 220 | 10 | 0.0414 | 0.0024 | 3.21e-06 (SD 0.0018) | 0.0021 |
| 330 | 10 | 0.0452 | 0.0025 | 5.55e-06 (SD 0.0024) | 0.0028 |
| 410 | 10 | 0.0479 | 0.0029 | 8.06e-06 (SD 0.0028) | 0.0032 |
| 490 | 10 | 0.0498 | 0.0027 | 9.67e-06 (SD 0.0031) | 0.0034 |
| 660 | 10 | 0.0544 | 0.0022 | 1.58e-05 (SD 0.0040) | 0.0042 |

- g = r210 − r160 = 0.0411 − 0.0389 = **0.0021** (paired SE 0.0004, 95% CI 0.0012 to 0.0030)
- V_Q(Δ) := 2·V_Q,single(160)·|D_val|/|T| = 2·2.46e-06·526/656 = **3.94e-06** (SD 0.0020)
- σ_r(160) = 0.0021; target SE(Δ̂(160)) ≤ g/3.5 = 0.0006 (SE² 3.73e-07)
- **V_Q(Δ) alone exceeds (g/3.5)^2: Q2 claims phrased as bounds (DESIGN §7.1)**
- SE(Δ̂(160)) by R_s (R_r = 2R_s): R_s=3: 0.0025, R_s=5: 0.0023, R_s=10: 0.0021, R_s=20: 0.0021, R_s=40: 0.0020

Sensitivity (NOT the frozen §7.1 rule): strategy and random arms are scored on the same queries with shared bootstrap weights, so the real V_Q(Δ) is reduced by their replicate correlation. Proxy = two random budgets on the same queries: ρ_boot(r160, r210) = 0.994; V_Q of the paired contrast (scaled to |T|) = 5.2e-08 (SD 0.0002); 2·V_Q,single·(1−ρ)·|D_val|/|T| = 2.34e-08 → R_s = 21 (paired proxy) / 19 (ρ formula), R_r = 2R_s. A strategy whose selections differ from random decorrelates more than two nested random prefixes do, so this is optimistic; P4 re-estimates V_Q(Δ) from real strategy-vs-random replicates.

Caveats: σ_s := σ_r (strategy trajectory variance unknown until P4); V_Q(Δ) assumes ρ = 0 (conservative); g is estimated from 10 nested-prefix trajectories (its CI above) — a smaller true g raises R steeply (R ∝ σ²/((g/3.5)² − V_Q)). R is re-estimated from P4 dev strategy runs and fixed at the freeze (may only go up).
