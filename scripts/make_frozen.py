#!/usr/bin/env python
"""Write config/frozen.yaml at the freeze (DESIGN v3.1 §8): SHA-256 of every tracked file matched by the freeze list.

Groups (hsl.final.FROZEN_SPEC):
  design_s8   the DESIGN §8 list: hsl/eval.py, views.py, guard.py, harness.py, protocol.yaml, split files, seeds,
              budgets, config/predictor.yaml, src/hsl/predictors/*.py, third_party/ysl/PROVENANCE.md,
              config/feature_provenance.yaml, DATA.md, yeast_sls_commit.txt, the data checksum table
  final_plan  config/finalists.yaml (must be committed first)
  extended    code that produces the final numbers but §8 does not name (selectors, features, upstream adapter,
              hsl.final, the final scripts, launch.py, the vendored predictor code); PI decision to keep
Only TRACKED files are hashed (git ls-files pathspec globs), and only from a clean tree, so the hashes are those of
the committed blobs. scripts/run_final.py refuses to start unless every hash still matches.

    .venv/bin/python scripts/make_frozen.py [--out config/frozen.yaml] [--force] [--allow-dirty]   (not run before
    the PI's freeze decision; writes nothing else)
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import argparse
import sys
import time
from pathlib import Path

import yaml

from hsl import final as F


def make(root: Path, out: Path, force: bool = False, allow_dirty: bool = False) -> dict:
    root = Path(root)
    if out.exists() and not force:
        raise SystemExit(f"{out} exists (use --force to regenerate before committing)")
    c = F.check_finalists_committed(root)
    if not c.ok:
        raise SystemExit(f"refusing: {c.detail}")
    clean = F.check_tree_clean(root)
    if not clean.ok and not allow_dirty:
        raise SystemExit(f"refusing: tree not clean ({clean.detail})")
    files = F.frozen_hashes(root, F.FROZEN_SPEC)
    rec = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "generator": "scripts/make_frozen.py",
           "git_head_at_generation": F._git(root, "rev-parse", "HEAD").stdout.strip(),
           "design": "v3.1 §8", "spec": F.FROZEN_SPEC, "files": files,
           "n_files": {k: len(v) for k, v in files.items()}}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("# FROZEN at the freeze (DESIGN §8). Generated; never hand-edited. Checked by scripts/run_final.py.\n"
                   + yaml.safe_dump(rec, sort_keys=False, width=120))
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", type=Path, default=F.ROOT)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--allow-dirty", action="store_true", help="tests only")
    a = ap.parse_args(argv)
    out = a.out or (a.root / F.FROZEN)
    rec = make(a.root, out, a.force, a.allow_dirty)
    print(f"{out}: " + ", ".join(f"{k} {v} files" for k, v in rec["n_files"].items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
