# Universe subsampling — scaling of n_80 / n_90 with the candidate-pool size (DESIGN §7.3, dev only)

Predictor stack_lr_wr2 (pre-P2 P*, reference); R0 random permutations; level 100% = dev splits 1–3 evaluated on the same pool-fraction grid (30, 60, …, 1640, max); sub-levels: grid scaled to the sub-pool (`hsl.universe.scaled_grid`). n_X relative to the curve's own (sub-)pool maximum (frozen PAVA/n_req).

| level | draw | N (pool) | |D_val| | E_fixed_dev positives | R | AUPRC at max (PAVA) | n_80 | n_90 | n_80/N | n_90/N |
|---|---|---|---|---|---|---|---|---|---|---|
| 50% | 0 | 1043 | 261 | 1646 | 10 | 0.0574 | 400 | 642 | 0.383 | 0.616 |
| 50% | 1 | 1045 | 260 | 1809 | 10 | 0.0586 | 348 | 608 | 0.333 | 0.581 |
| 50% | 2 | 1035 | 259 | 1942 | 10 | 0.0774 | 449 | 709 | 0.434 | 0.685 |

Slope of log n_X on log N (OLS over draws; 95% CI = draw bootstrap stratified by level, B = 2,000):

- Slope not estimable: fewer than two universe levels with a (sub-)pool-maximum anchor for this predictor (the pre-P2 reference has no pool-maximum run at 100%, PI decision for batch6).

## What was and was not restricted to the sub-universe

- Restricted (recomputed inside U): query rows and array columns (same genes removed from both sides; T columns dropped), inner folds (rebuilt, val = fold 1), candidate pool, E_fixed_dev, and every statistic the predictor derives from the view — array rates (incl. LOO), query propensity, kNN similarities, `tied` low-rank, 2-hop scores, GI-profile kernel, external-kernel normalisation over the training genes (mean diagonal, PSD projection), λ scaling, LR feature standardisation; array SMF is per gene.
- Not restricted: the upstream pair-feature store used by `mk`/`ext_off` (coexpression Fisher z, AP-MS edge score and neighbourhood cosine, ESM-2 standardised cosine, paralog identity, InterPro Jaccard). The upstream builders compute these over each SOURCE's own gene set (all genes of each expression study; all AP-MS network genes, deliberately 'so degrees are not truncated'; all ESM-embedded proteins), not over the S1 universe, so rebuilding the store on U with the upstream code would return identical values. Restricting them would require re-deriving the sources. The sub-universes therefore keep full-genome external context: **b is reported as a bound** (direction not determined; external features are per pair, so the effect is expected to be small).
- Caveat (DESIGN §7.3): subsampling mimics fewer genes per module, not more modules. At 25% E_fixed_dev has ~420–500 positives (noise check above).

