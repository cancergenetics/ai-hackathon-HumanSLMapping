This is an exercise in assessing the bioinformatics capabilities of agentic AI. This particular exercise was run by the Ryan lab at University College Dublin as part of a larger hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects. All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.

Models and machine: the project was run in Claude Code. The PI agent and most sub-agents used Claude Opus 5.5 (`claude-opus-5-5`); one retrospective reviewer used Claude Sonnet 5.5. The sub-agents covered literature review, design reviewers, data and evaluation infrastructure, predictor development, strategy evaluation, the final-evaluation pipeline, an independent rescorer, the Q3 analyses, report writing, and code and results reviewers. All computation ran locally on an Apple MacBook Air (M2) with an 8-core CPU (4 performance + 4 efficiency cores) and 24 GB unified memory. The integrated 10-core GPU was not used; all model fitting was CPU-only.

Currently none of the results have been verified by a human. We plan to verify the agent’s methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It is a simulation study of experimental design. It includes code for the SL predictor, the query-selection strategies and the one-shot held-out evaluation, and an AI-written report (REPORT.md) on the design, results and recommendations. The report describes the final approach and the path followed to reach it, including the pre-registered design, independent reviews and declared deviations. It benchmarks the query-selection strategies against random and diversity-based query panels, and the predictor against simpler baselines.

---

# Human SL mapping: designing HAP1 synthetic-lethal screens from the yeast network

We used the essentially complete yeast genetic-interaction network (Costanzo et al., *Science* 2016) to ask how human query-by-genome CRISPR screens in HAP1 cells (Billmann, Costanzo et al., *Cell* 2026) should be designed to find strong negative / synthetic-lethal interactions (ε < −0.12, P < 0.05). The questions were:

- **Q1.** How does SL prediction for unscreened genes scale with the fraction of queries screened at random?
- **Q2.** Can a smaller, iteratively chosen set of queries (batches of 10) match a larger random set at matched budgets?
- **Q3.** What does this imply for the number and choice of HAP1 query lines, and what can prediction contribute?

## Start here

| File | What it is |
|---|---|
| [BRIEFING.md](BRIEFING.md) | Research briefing for a broad scientific audience: the biology and the agentic process, with figures |
| [REPORT.md](REPORT.md) | Full report: approach, results on the held-out test set, the HAP1 recommendation, limitations, follow-ups |
| [docs/explainer/pairs_that_kill.html](docs/explainer/pairs_that_kill.html) | Short illustrated summary for systems-biology readers |
| [RESEARCH_TREE.md](RESEARCH_TREE.md) | Project log: status, decision log, research tree |
| [DESIGN.md](DESIGN.md) | Pre-registered design (v3.1) and declared deviations |
| [PROMPT.md](PROMPT.md) | The original brief |
| [LITERATURE.md](LITERATURE.md), [DATA.md](DATA.md), [LESSONS.md](LESSONS.md) | Background, data provenance, process lessons |
| `docs/reviews/` | Every independent review (design, code, predictor, strategies, final pipeline, results, report) |

## Headline results (held-out test queries, one-shot evaluation)

- **Random sampling.** AUPRC for strong negatives on held-out queries rises with no plateau: 0.042 at 5% of queries screened, 0.063 at 20%, 0.092 with the whole pool screened.
- **Adaptive selection.** Adaptive batch selection (hub-first, coverage, Fisher design) beats random at 160 queries by +26–31% (Holm p ≤ 0.0014). That is equivalent to about 2.4–2.7× as many random queries. It does not match 20% random.
- **HAP1 projection.** About 1,000 adaptively chosen query lines in about 4 rounds of 300 reach 80% of full-pool predictive performance. Random selection would need about 3,300. Both figures carry scaling and denominator uncertainty (see REPORT §5).

## Data

No datasets are included in this repository. The yeast network and external feature stores are read in place from a separate local checkout of the upstream YEAST_SLs project, at the commit recorded in `config/yeast_sls_commit.txt`. `DATA.md` lists the source, version, checksum and licence of each dataset. Paths are configured in `config/paths.yaml`, and `/data/`, prediction vectors and bootstrap caches are git-ignored.

## Layout

- `src/hsl/` — library: evaluation engine, views and leakage guard, harness, strategies, predictors, final-evaluation runner.
- `scripts/` — numbered pipeline: split → P1 random curves → P2 predictor → P3 oracles → P4 strategies → final (`run_final.py`, `score_final.py`).
- `analysis/q3/` — HAP1 translation, choice of queries, prediction contribution, stopping rule.
- `analysis/rescore/` — independent rescorer. It imports nothing from `src/`.
- `analysis/report/` — report figures.
- `config/` — split, seeds, budgets, protocol constants, frozen predictor, finalists and frozen hashes.
- `results/` — committed summaries, CSVs and figures for each phase; `results/final/` holds the scored test-set outputs.
- `tests/` — unit, leakage-poisoning and review tests.
- `third_party/ysl/` — vendored upstream method code (read-only).

## Reproducing

```bash
uv sync                                   # environment (pyproject.toml, uv.lock)
# point config/paths.yaml at a local YEAST_SLs checkout (commit in config/yeast_sls_commit.txt)
.venv/bin/python -m pytest -q             # tests
.venv/bin/python scripts/run_final.py check --mode final
.venv/bin/python scripts/score_final.py --mode final                      # re-score saved predictions
.venv/bin/python analysis/rescore/rescore.py --run-dir results/final --out analysis/rescore/out/final
bash analysis/q3/run_all.sh final
.venv/bin/python analysis/report/make_figures.py
```

The final evaluation is locked (`results/final/lock_record.json`). Re-running it would break the pre-registered one-shot design. REPORT.md §9 gives the full command sequence and the freeze/lock commits.
