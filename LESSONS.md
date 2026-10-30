# LESSONS — process rules for running this project efficiently

Each rule: one line, with its evidence. Pruned at each milestone.

1. Read large upstream data in place behind a write guard (provenance + no copies); re-check disk when constraints change. (Evidence: 38 GB free at start; user later freed ~500 GB, 2026-10-29.)
2. Anchor gitignore data patterns (`/data/`, not `data/`): the bare pattern silently dropped the vendored `third_party/ysl/data/` package from the first commit (2026-10-29).
3. Run `git` in YEAST_SLs only with `GIT_OPTIONAL_LOCKS=0`: a plain `git status` there touched `.git` (opportunistic index lock), which the stamp check flags (2026-10-29).
4. While sub-agents are writing in the repo, the PI stages explicit paths only (never `git add -A`). (Evidence: 2026-10-29 tree commit could sweep in P0 agent's in-progress files.)
5. When changing a storage/output policy mid-run, say explicitly whether it applies retroactively and who may delete what. (Evidence: 2026-10-29 P1 agent bulk-deleted its own regenerable baseline .npz files to apply the new disk policy; flagged by the safety monitor, cost an investigation.)
6. Set BLAS thread env vars (OMP/OPENBLAS/VECLIB/MKL) in every job and count all agents' jobs against one core budget; parallel agents each 'within limits' oversubscribed 8 cores to load 120 and slowed everything ~1.5–3×. (2026-10-29) Grant each agent a total thread budget covering all its concurrent jobs (incl. ones still draining), not a per-batch 'k × t': a per-batch grant stacked on a draining tranche gave load 40 (2026-10-29). But size budgets by measured CPU/RAM, not thread caps: jobs are mostly single-threaded (~1 core each), so 'k × 2 threads' left the machine 44% idle; check `top` idle % and RSS after each reallocation (2026-10-29).
7. Edit RESEARCH_TREE.md with anchored inserts and verify with grep after each edit; chained str.replace edits silently no-op'd and lost two tree entries (2026-10-29).
8. Background waits must track exactly one PID (have the wrapper `exec` the launcher): a chained launcher returned two PIDs, `kill -0` waits broke, and the job budget was silently exceeded for ~20 min (2026-10-29).
9. Dry runs must exercise every final-only code path (deep preflight, label-gate interactions); the first final `start` failed on an untested final-only checksum path (2026-10-30).
