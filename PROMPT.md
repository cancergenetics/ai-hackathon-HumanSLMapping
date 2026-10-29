/goal You are the PI of a computational biology project and work autonomously. You lead a team of agents (Opus 5.5 or Sonnet 5.5, as you see fit). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

In budding yeast, large-scale genetic interaction maps were assembled by screening query mutants against an array of deletion mutants (and temperature-sensitive alleles of essential genes). Analogous efforts are now under way in human HAP1 cells, in which isogenic query knockout lines are screened against genome-wide CRISPR libraries (Billmann et al., Cell 2026). Your task is to use the existing yeast genetic interaction network to determine how best to design these experiments in human cells.

Focus: identification of strong synthetic lethal (SL) interactions. Positive interactions and quantitative interaction scores are out of scope. Define "strong SL" as [ε < -0.12, p < 0.05]

Resources:
- A method for predicting SL interactions from incomplete yeast screens. For example, given screens for 20% of query genes, it predicts interactions for the remaining 80%. It is located at /Users/colmclaude/claude_experiments/YEAST_SLs. This directory is read-only: never modify, add, or delete files in it. 
- The full yeast genetic interaction matrix, for evaluation only. No part of your approach (query selection, hyperparameter tuning, or stopping rules) may use information from unscreened queries. You can use other datasets as input to your model (e.g. gene expression, protein-protein interaction networks and features) but you should avoid datasets that might introduce circularity (no genetic interactions from any source, no text mining or resources that derive from text mining such as STRING, no gene ontology as the ontology annotations may have arisen through interaction experiments). You should restrict your feature sources to ones where there is an obvious analogue for humans (e.g. yeast has comprehensive knockouts profiled with proteomics but this is not available for humans so should be avoided).

Write all new code, files, and outputs to /Users/colmclaude/claude_experiments/humanSLmapping

Questions to address:
1. How does SL prediction performance scale with the fraction of queries screened (random selection)?
2. Can a smaller, strategically chosen query set match the performance of a larger random set? For example, can 5% chosen genes match 20% random? Assume screening is iterative: you may screen batches of 10 queries and use the results to choose the next batch. Compare against random selection at matched budgets.
3. What do these results imply for the number and choice of query lines needed to map the HAP1 SL network, and what can computational prediction realistically contribute?

Before any development, split query genes into a fixed held-out test set (20%, stratified by number of SLs) and a development pool. Test queries are never screened, selected, or used for tuning; all strategies are evaluated on the same test set. Treat pairs involving a screened gene as observed regardless of which side (query or array) it was screened on, and state how symmetric pairs are handled. Evaluate performance using AUPRC on held-out queries.

Deliverables: a written summary with figures, recommended human screening strategy and estimated query numbers, reproducible code.

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and "what worked well?". 

Auto-research: Your job is to build the best possible method. Report runtime alongside accuracy. Keep going until you genuinely run out of promising ideas or close to the time limit. Then write REPORT.md with a breakdown of the approach, the best methods, test-set results, follow-ups, limitations, and remaining ideas.
