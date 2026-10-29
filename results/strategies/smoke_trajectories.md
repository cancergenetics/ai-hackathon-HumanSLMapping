# Strategy smoke trajectories (dev split 1) — SMOKE, NOT FOR DECISIONS

One trajectory per arm, stack_lr, E_fixed_dev pooled AUPRC. Proves the pipeline end to end (SelectorView -> selector -> harness -> predictor -> metrics); single trajectories carry no comparative information (no random arm here by design; comparisons are P4, after the P* freeze).

| arm | AUPRC@30 | AUPRC@60 | AUPRC@100 | select s/round (mean / max) | off-array picked |
|---|---|---|---|---|---|
| R1 | 0.0286 | 0.0302 | 0.0328 | 1.46 / 14.56 | 21 |
| H | 0.0291 | 0.0320 | 0.0356 | 0.02 / 0.04 | 15 |
| H-static | 0.0297 | 0.0324 | 0.0355 | 0.01 / 0.03 | 12 |
| F | 0.0266 | 0.0333 | 0.0371 | 2.22 / 14.55 | 17 |
| F-static | 0.0266 | 0.0305 | 0.0334 | 0.02 / 0.03 | 18 |
| C | 0.0265 | 0.0355 | 0.0357 | 0.16 / 0.17 | 14 |
| M | 0.0279 | 0.0286 | 0.0339 | 0.01 / 0.02 | 1 |
