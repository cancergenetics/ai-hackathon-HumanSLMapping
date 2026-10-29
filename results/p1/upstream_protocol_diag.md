# Upstream-protocol gap diagnosis (dev only)

Question: P* (stack_lr_wr2) scores 0.0657 on the upstream-protocol set of dev1 at n = 660 vs upstream's dev 0.0754 at 20%. Set = every measured cell of the held-out rows incl. reciprocal-known cells, minus test-gene (T) columns (dev rule) and self pairs. Strata as upstream's evaluator (`recip_available`: target on array, array's ORF screened, reciprocal cell observed).

## dev1, n = 660, our R0 random sets (same sets for all predictors)

| predictor | R | AUPRC | prev | lift | recip: AUPRC / prev / lift | on-array no-recip | off-array | upstream REPORT §3: AUPRC (lift); recip / on-no-recip / off |
|---|---|---|---|---|---|---|---|---|
| stack_lr_wr2 | 10 | 0.0657 ± 0.0027 | 0.0110 | 5.97 | 0.0999 / 0.0116 / 8.58 | 0.0613 / 0.0098 / 6.25 | 0.0579 / 0.0150 / 3.87 | 0.0754 (lift 6.6); 0.126 / 0.073 / 0.058 |
| stack_lr | 10 | 0.0543 ± 0.0017 | 0.0110 | 4.93 | 0.0959 / 0.0116 / 8.24 | 0.0458 / 0.0098 / 4.67 | 0.0477 / 0.0150 / 3.18 | 0.0590 (lift 5.1); 0.117 / 0.051 / 0.047 |
| B1_array_rate | 10 | 0.0417 ± 0.0001 | 0.0110 | 3.79 | 0.0392 / 0.0116 / 3.37 | 0.0405 / 0.0098 / 4.13 | 0.0479 / 0.0150 / 3.20 | 0.0430 (lift 3.7); 0.038 / 0.043 / 0.048 |

## dev2, n = 660, our R0 random sets (same sets for all predictors)

| predictor | R | AUPRC | prev | lift | recip: AUPRC / prev / lift | on-array no-recip | off-array | upstream REPORT §3: AUPRC (lift); recip / on-no-recip / off |
|---|---|---|---|---|---|---|---|---|
| stack_lr_wr2 | 3 | 0.0796 ± 0.0038 | 0.0115 | 6.92 | 0.1308 / 0.0111 / 11.88 | 0.0754 / 0.0100 / 7.57 | 0.0677 / 0.0171 / 3.96 | 0.0754 (lift 6.6); 0.126 / 0.073 / 0.058 |
| stack_lr | 3 | 0.0609 ± 0.0005 | 0.0115 | 5.30 | 0.1228 / 0.0111 / 11.15 | 0.0509 / 0.0100 / 5.11 | 0.0516 / 0.0171 / 3.02 | 0.0590 (lift 5.1); 0.117 / 0.051 / 0.047 |
| B1_array_rate | 3 | 0.0448 ± 0.0001 | 0.0115 | 3.90 | 0.0423 / 0.0111 / 3.81 | 0.0435 / 0.0100 / 4.37 | 0.0523 / 0.0171 / 3.06 | 0.0430 (lift 3.7); 0.038 / 0.043 / 0.048 |

## dev3, n = 660, our R0 random sets (same sets for all predictors)

| predictor | R | AUPRC | prev | lift | recip: AUPRC / prev / lift | on-array no-recip | off-array | upstream REPORT §3: AUPRC (lift); recip / on-no-recip / off |
|---|---|---|---|---|---|---|---|---|
| stack_lr_wr2 | 3 | 0.0757 ± 0.0053 | 0.0116 | 6.55 | 0.1318 / 0.0129 / 10.20 | 0.0712 / 0.0107 / 6.68 | 0.0533 / 0.0138 / 3.88 | 0.0754 (lift 6.6); 0.126 / 0.073 / 0.058 |
| stack_lr | 3 | 0.0628 ± 0.0037 | 0.0116 | 5.44 | 0.1241 / 0.0129 / 9.60 | 0.0535 / 0.0107 / 5.02 | 0.0447 / 0.0138 / 3.25 | 0.0590 (lift 5.1); 0.117 / 0.051 / 0.047 |
| B1_array_rate | 3 | 0.0434 ± 0.0002 | 0.0116 | 3.76 | 0.0435 / 0.0129 / 3.36 | 0.0429 / 0.0107 / 4.02 | 0.0452 / 0.0138 / 3.29 | 0.0430 (lift 3.7); 0.038 / 0.043 / 0.048 |

## Upstream's own dev splits re-run here (make_split('s1','dev',seed); model_seed 17; T genes dropped)

Screened 653 of 657 (T genes have no rows in the dev bundle), targets 1951 of 1970; T columns excluded.

| predictor | R | AUPRC | prev | lift | recip: AUPRC / prev / lift | on-array no-recip | off-array | upstream REPORT §3: AUPRC (lift); recip / on-no-recip / off |
|---|---|---|---|---|---|---|---|---|
| stack_lr_wr2 | 1 | 0.0774 ± 0.0000 | 0.0114 | 6.81 | 0.1312 / 0.0115 / 11.41 | 0.0715 / 0.0105 / 6.84 | 0.0612 / 0.0146 / 4.19 | 0.0754 (lift 6.6); 0.126 / 0.073 / 0.058 |
| stack_lr | 1 | 0.0615 ± 0.0000 | 0.0114 | 5.40 | 0.1236 / 0.0115 / 10.75 | 0.0510 / 0.0105 / 4.88 | 0.0483 / 0.0146 / 3.30 | 0.0590 (lift 5.1); 0.117 / 0.051 / 0.047 |
| B1_array_rate | 1 | 0.0433 ± 0.0000 | 0.0114 | 3.81 | 0.0406 / 0.0115 / 3.53 | 0.0424 / 0.0105 / 4.05 | 0.0487 / 0.0146 / 3.33 | 0.0430 (lift 3.7); 0.038 / 0.043 / 0.048 |

Upstream m4b5 stack_lr_wr2 split 0 (all columns incl. lockbox): 0.0775; here 0.0774.

Same upstream-split predictions, restricted to the target rows that belong to each of our dev folds (AUPRC / prevalence / lift):

| predictor | dev1 rows | dev2 rows | dev3 rows | folds 4–5 rows |
|---|---|---|---|---|
| stack_lr_wr2 | 0.0673 / 0.0109 / 6.19 | 0.0864 / 0.0116 / 7.47 | 0.0750 / 0.0116 / 6.47 | 0.0796 / 0.0114 / 6.97 |
| stack_lr | 0.0563 / 0.0109 / 5.18 | 0.0656 / 0.0116 / 5.66 | 0.0614 / 0.0116 / 5.30 | 0.0622 / 0.0114 / 5.44 |
| B1_array_rate | 0.0409 / 0.0109 / 3.77 | 0.0458 / 0.0116 / 3.95 | 0.0431 / 0.0116 / 3.72 | 0.0435 / 0.0114 / 3.81 |

## Ratios to upstream (ours / upstream REPORT §3 AUPRC)

- dev1: stack_lr_wr2 0.871, stack_lr 0.920, B1_array_rate 0.969
- dev2: stack_lr_wr2 1.056, stack_lr 1.033, B1_array_rate 1.043
- dev3: stack_lr_wr2 1.004, stack_lr 1.064, B1_array_rate 1.010
- upstream dev splits re-run: stack_lr_wr2 1.027, stack_lr 1.042, B1_array_rate 1.007

## Conclusion

- **Pipeline/stores/params are correct:** upstream's own dev split 0 re-run through our harness path (`hsl.upstream.construct_view` + `run_method`, registered defaults = upstream recipe `params: {}`, model seed 17) gives P* 0.0774 vs upstream's recorded 0.0775 (lift 6.8), stack_lr 0.0615, B1 0.0433. All six external pair stores used by `mk`/`ext_off` (apms_cosine, apms_edge, coex_z, dom_ipr_jaccard, esm_cos, para_identity) cover dev1 D_val and the screened genes like the dev pool as a whole (coverage 0.70–1.00 per store, every gene in the 4,426-gene store).
- **Cause = dev1's target rows (query composition), not the predictor:** with upstream's own split-0 predictions, the dev1 rows alone score 0.0673 (P*), 0.0563 (stack_lr), 0.0409 (B1) — the same ~0.87/0.92/0.97 shortfall we see — while dev2 rows score 0.0864 and dev3 rows 0.0750. Our dev2/dev3 upstream-protocol runs at n = 660 give P* 0.0796 / 0.0757 (ratio to upstream 1.06 / 1.00). dev1 has lower prevalence (0.0109) and its queries are harder for every GI-based model; the gap grows with model strength because the GI-driven gain (recip / on-array strata) is what varies between query sets. Not a composition shift of the cell set (B1 lift 3.79 ≈ upstream 3.7) and not a P*-only defect.
- Implication: absolute dev AUPRCs differ by split by up to ~±10% (query sampling); within-split paired contrasts (same queries) are unaffected. The power pilot (split 1) is on the hardest fold.

