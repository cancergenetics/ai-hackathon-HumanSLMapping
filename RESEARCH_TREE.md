# RESEARCH TREE — Designing human HAP1 SL-mapping screens from the yeast GI network

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-29, M0 setup)

Project started. PROMPT.md saved verbatim. Upstream method (YEAST_SLs, commit 9b5dce9, read-only) is a kernel-ridge
fold-in + stacked combiner that predicts stringent negative GIs (ε < −0.12, P < 0.05) for unscreened yeast queries
from the Costanzo 2016 NxN matrix (lockbox AUPRC 0.083 at 20% screened, prevalence 0.011; 10/20/30% → 0.061/0.083/0.102).
It already reads Costanzo S1/CellMap bundles and tier-1 external features (coexpression, AP-MS, ESM-2, paralogs, etc.).

Next: literature review (LITERATURE.md) and technical reconnaissance of the upstream code (how to call it with an
arbitrary screened set, runtime vs. number screened, data characterisation), in parallel; then the fixed test split,
then DESIGN.md and up to three review rounds.

## Decision log
- 2026-10-29 D0: Do not copy upstream data (11 GB; 38 GB free disk). Read upstream processed data in place; any code
  that could write into YEAST_SLs is vendored/redirected. A stamp file (`.yeast_readonly_stamp`) plus upstream
  `git status` are checked at each milestone to prove YEAST_SLs was not modified.

## Tree
- [~] M0 Setup: repo, PROMPT.md, tree, LESSONS, DATA.md skeleton
- [ ] M1 Literature + reconnaissance
  - [ ] LITERATURE.md (yeast SGA/GI maps, HAP1 GI screens incl. Billmann et al. 2026, SL prediction, active learning / experimental design for matrix completion, human-analogue feature datasets)
  - [ ] Upstream method recon: API, runtime scaling, data characterisation
- [ ] M2 Fixed test split (20% query genes, stratified by #strong SLs) — before any development
- [ ] M3 DESIGN.md + review rounds (≤3)
- [ ] M4 Q1: performance vs. fraction screened (random)
- [ ] M5 Q2: strategic iterative selection (batches of 10) vs. random at matched budgets
- [ ] M6 Q3: translation to HAP1 (query numbers, choice of query lines, role of prediction)
- [ ] M7 REPORT.md, figures, final verification
