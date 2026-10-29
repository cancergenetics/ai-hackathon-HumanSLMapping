# LESSONS — process rules for running this project efficiently

Each rule: one line, with its evidence. Pruned at each milestone.

1. Read large upstream data in place behind a write guard (provenance + no copies); re-check disk when constraints change. (Evidence: 38 GB free at start; user later freed ~500 GB, 2026-10-29.)
2. Anchor gitignore data patterns (`/data/`, not `data/`): the bare pattern silently dropped the vendored `third_party/ysl/data/` package from the first commit (2026-10-29).
3. Run `git` in YEAST_SLs only with `GIT_OPTIONAL_LOCKS=0`: a plain `git status` there touched `.git` (opportunistic index lock), which the stamp check flags (2026-10-29).
4. While sub-agents are writing in the repo, the PI stages explicit paths only (never `git add -A`). (Evidence: 2026-10-29 tree commit could sweep in P0 agent's in-progress files.)
5. When changing a storage/output policy mid-run, say explicitly whether it applies retroactively and who may delete what. (Evidence: 2026-10-29 P1 agent bulk-deleted its own regenerable baseline .npz files to apply the new disk policy; flagged by the safety monitor, cost an investigation.)
6. Set BLAS thread env vars (OMP/OPENBLAS/VECLIB/MKL) in every job and count all agents' jobs against one core budget; parallel agents each 'within limits' oversubscribed 8 cores to load 120 and slowed everything ~1.5–3×. (2026-10-29)
