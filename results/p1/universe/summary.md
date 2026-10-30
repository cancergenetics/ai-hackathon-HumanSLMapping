# Universe subsampling — scaling of n_80 / n_90 with the candidate-pool size (DESIGN §7.3, dev only)

Predictor stack_lr_wr2_lb (frozen P*); R0 random permutations; level 100% = dev splits 1–3 evaluated on the same pool-fraction grid (30, 60, …, 1640, max); sub-levels: grid scaled to the sub-pool (`hsl.universe.scaled_grid`). n_X relative to the curve's own (sub-)pool maximum (frozen PAVA/n_req).

| level | draw | N (pool) | |D_val| | E_fixed_dev positives | R | AUPRC at max (PAVA) | n_80 | n_90 | n_80/N | n_90/N |
|---|---|---|---|---|---|---|---|---|---|---|
| 100% | 1 | 2101 | 526 | 6755 | 10 | 0.0725 | 791 | 1277 | 0.377 | 0.608 |
| 100% | 2 | 2100 | 527 | 7323 | 10 | 0.0913 | 945 | 1374 | 0.450 | 0.654 |
| 100% | 3 | 2103 | 524 | 7106 | 10 | 0.0763 | 651 | 1185 | 0.310 | 0.563 |
| 50% | 0 | 1043 | 261 | 1646 | 10 | 0.0575 | 378 | 635 | 0.363 | 0.609 |
| 50% | 1 | 1045 | 260 | 1809 | 10 | 0.0588 | 334 | 582 | 0.320 | 0.557 |
| 50% | 2 | 1035 | 259 | 1942 | 10 | 0.0779 | 429 | 698 | 0.415 | 0.674 |
| 25% | 0 | 534 | 134 | 424 | 10 | 0.0398 | 132 | 228 | 0.247 | 0.427 |
| 25% | 1 | 520 | 131 | 462 | 10 | 0.0662 | 203 | 367 | 0.390 | 0.706 |
| 25% | 2 | 528 | 132 | 504 | 10 | 0.0683 | 221 | 299 | 0.419 | 0.566 |

Slope of log n_X on log N (OLS over draws; 95% CI = draw bootstrap stratified by level, B = 2,000):

- **n80: b = 1.06** [0.85, 1.28] (9 draws). b ≈ 0 → count scaling, b ≈ 1 → fraction scaling.
- **n90: b = 1.06** [0.89, 1.24] (9 draws). b ≈ 0 → count scaling, b ≈ 1 → fraction scaling.

## What was and was not restricted to the sub-universe

- Restricted (recomputed inside U): query rows and array columns (same genes removed from both sides; T columns dropped), inner folds (rebuilt, val = fold 1), candidate pool, E_fixed_dev, and every statistic the predictor derives from the view — array rates (incl. LOO), query propensity, kNN similarities, `tied` low-rank, 2-hop scores, GI-profile kernel, external-kernel normalisation over the training genes (mean diagonal, PSD projection), λ scaling, LR feature standardisation; array SMF is per gene.
- Not restricted: the upstream pair-feature store used by `mk`/`ext_off` (coexpression Fisher z, AP-MS edge score and neighbourhood cosine, ESM-2 standardised cosine, paralog identity, InterPro Jaccard). The upstream builders compute these over each SOURCE's own gene set (all genes of each expression study; all AP-MS network genes, deliberately 'so degrees are not truncated'; all ESM-embedded proteins), not over the S1 universe, so rebuilding the store on U with the upstream code would return identical values. Restricting them would require re-deriving the sources. The sub-universes therefore keep full-genome external context: **b is reported as a bound** (direction not determined; external features are per pair, so the effect is expected to be small).
- Caveat (DESIGN §7.3): subsampling mimics fewer genes per module, not more modules. At 25% E_fixed_dev has ~420–500 positives (noise check above).

