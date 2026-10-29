# Strategy smoke trajectories (dev split 1) — SMOKE, NOT FOR DECISIONS

One trajectory per arm, stack_lr, E_fixed_dev pooled AUPRC. Proves the pipeline end to end (SelectorView -> selector -> harness -> predictor -> metrics); single trajectories carry no comparative information (no random arm here by design; comparisons are P4, after the P* freeze).

| arm | AUPRC@30 | AUPRC@60 | AUPRC@100 | select s/round (mean / max) | off-array picked |
|---|---|---|---|---|---|
| R1 | 0.0282 | 0.0306 | 0.0358 | 1.65 / 16.51 | 19 |
| H | 0.0306 | 0.0327 | 0.0373 | 0.01 / 0.04 | 12 |
| H-static | 0.0302 | 0.0316 | 0.0357 | 0.01 / 0.03 | 16 |
| F | 0.0266 | 0.0319 | 0.0350 | 2.08 / 13.19 | 5 |
| F-static | 0.0266 | 0.0305 | 0.0334 | 0.02 / 0.04 | 18 |
| C | 0.0280 | 0.0338 | 0.0349 | 0.18 / 0.23 | 5 |
| M | 0.0291 | 0.0319 | 0.0349 | 0.01 / 0.02 | 1 |
