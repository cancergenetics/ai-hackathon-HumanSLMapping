# LESSONS — process rules for running this project efficiently

Each rule: one line, with its evidence. Pruned at each milestone.

Pruned 2026-10-30 from the fresh-agent retrospective (docs/reviews/retro_m_final.md). Earlier lessons 1–9 are merged into
the list below.
1. YEAST_SLs is read-only. Read upstream data in place. Prefix git with GIT_OPTIONAL_LOCKS=0. Run scripts/check_readonly.sh at each milestone.
2. Repo hygiene with several agents active:
   - anchored .gitignore patterns;
   - stage explicit paths only;
   - edit RESEARCH_TREE.md with anchored inserts and verify each edit with grep (chained replaces silently no-op'd once).
3. Every policy change states whether it applies retroactively.
4. Use one global compute budget across all agents. Size it from measured CPU idle and RSS, not thread caps (jobs use about 1 core each).
   - Set the BLAS thread env vars in every job.
   - Grant each agent a total budget that includes draining jobs.
   - Enforce the budget in the launcher's memory/concurrency gate.
   - Background waits track exactly one PID.
5. Dry runs must exercise every final-only code path. The first final `start` failed on an untested deep-preflight path.
6. Evaluate every pre-registered rule and trigger on pilot output before freezing it. The ρ = 0 power rule and the §12 R action were both undefined when needed (D9, D10).
7. Label dev-selected winners and oracles for what they are: best-of-k, in-sample, reference. Expect shrinkage. Never quote them as ceilings or headlines.
8. Decisions live in DESIGN (deviations) and in the status summary, not only in the decision log. Rewrite the summary at every milestone.
9. Subagents cannot write report-type files (the harness blocks it). Brief writer and reviewer agents to return text; the PI saves it. Also, write the final-evaluation code for any post-hoc analysis before scoring, not only its coefficients (the static-rule provenance caveat) (2026-10-30).
