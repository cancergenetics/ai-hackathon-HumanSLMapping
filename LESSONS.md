# LESSONS — process rules for running this project efficiently

Each rule: one line, with its evidence. Pruned at each milestone.

1. Check disk before copying data; read large upstream data in place behind a write guard. (Evidence: 38 GB free vs 11 GB upstream processed data, 2026-10-29.)
2. Anchor gitignore data patterns (`/data/`, not `data/`): the bare pattern silently dropped the vendored `third_party/ysl/data/` package from the first commit (2026-10-29).
3. Run `git` in YEAST_SLs only with `GIT_OPTIONAL_LOCKS=0`: a plain `git status` there touched `.git` (opportunistic index lock), which the stamp check flags (2026-10-29).
4. While sub-agents are writing in the repo, the PI stages explicit paths only (never `git add -A`). (Evidence: 2026-10-29 tree commit could sweep in P0 agent's in-progress files.)
