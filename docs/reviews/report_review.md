# Review: REPORT.md (draft at 6ec6c2b)

Independent, read-only reviewer; wrote none of REPORT.md. 2026-10-30. The harness blocked the reviewer's own file write, so the PI saved this file. It is condensed from the reviewer's returned text; the full text is in the session transcript.

## Verdict: READY WITH FIXES

### Numbers
About 110 numbers in §0 and §2–§6 were checked against committed sources. 9 were mismatches or misleading roundings; none changed a conclusion.

| # | REPORT said | Source says | Fix applied |
|---|---|---|---|
| 1 | dev n_80 per split 0.38–0.45 | 0.31–0.45 | yes |
| 2 | E_unobs @660 0.069 / 0.070 | identical: 0.0695 | yes |
| 3 | 3.2–3.8× at 330–490 | 3.2–4.3× (3.2–3.8× for H, C) | yes |
| 4 | runtimes were eval_total | P\* fit_predict 28–31 / 63 / 548 s; stack_lr 1.8 / 34 s | yes |
| 5 | static rule "(partial at 660)" | 3 fits at 660 | deleted |
| 6 | "Pearson 0.97–0.99" | pooled over 3 H + 3 random trajectories | reworded |
| 7 | threshold 0.12–0.18 | holds for n ≥ 100 only; 0.075 at n = 30 | yes |
| 8 | b CI 1.28 | the translation's re-bootstrap gives 1.30 | yes |
| 9 | "p@50 ≈ 0" | the difference is ≈ 0 | yes |

Minor: "2–2.7%" → "1.7–2.7%".

### Claims
Labels and pre-stated readings were honoured. Still too strong:
- "H lost nothing at b110";
- "match the precision";
- "R1 gives the same as random";
- "most of the benefit comes from screening the right kind of gene".

All four were softened.

### Static-rule provenance
The coefficients were fixed at bde1622, before scoring. The code that evaluated the rule on T was written after scoring (f623c5e). Two rules were run on T: k1 recovered 77–97% and k3 recovered 91–108%. The report now states both, plus the dev values (72–81%).

### Q3
- **Unit error:** "rounds of up to about 110 lines" should be about 300 HAP1 lines (110 yeast queries). Fixed.
- **Missing caveats:** the CI omitting trajectory variance and the ±27% denominator lever are now carried into §0 and §6.
- **b wording:** b is "an approximation of undetermined direction", not "a bound". Fixed.
- **Additions:**
  - on-array top-10 values;
  - the "what 80% buys" sentence;
  - line feasibility for low-fitness lines.

### Figures
- **Q3 figures:** none were embedded. `fig_hap1_curves.png` and `fig_choice_profile.png` are now added as Figures 5 and 6.
- **Numbering:** Figures 3 and 4 were renumbered into order of appearance.

### Clarity
- Added a §0 setup preamble defining the held-out set, AUPRC, lift, pool, H/C/F/R1 and n_80.
- Optional cuts were not applied: the report is about 30% longer than necessary (§7 deviations duplicate DESIGN, and §2 Freeze duplicates §9).
