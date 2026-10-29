#!/usr/bin/env python
"""Evaluation-set manifest (DESIGN v3 §3, before any scoring): sizes / prevalence / strata of E_fixed_dev on every dev
split and the re-split check (dev rows), gene-identity groups and removals, inner-fold sizes, and LABEL-ONLY counts of
the final E_fixed (hsl.views.efixed_label_summary: aggregate counts inside the gate; no scores, no per-gene output).

    .venv/bin/python scripts/04_eset_manifest.py   ->  results/esets/manifest.json
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import json

from hsl import paths, views


def main() -> dict:
    man = views.eset_manifest(include_final=True)
    man["git"] = views.git_state()
    out = paths.output_root() / "results" / "esets"
    out.mkdir(parents=True, exist_ok=True)
    (out / "manifest.json").write_text(json.dumps(man, indent=1, default=str) + "\n")
    fin = man["final"]["E_fixed_label_counts"]
    print(json.dumps({"final_E_fixed": {k: fin[k] for k in ("n_cells", "n_pos", "prevalence", "AO", "TT",
                                                             "n_rows_zero_pos")},
                      "identity_excluded_arrays_final": fin["identity_excluded_arrays"],
                      "dev": {k: v["E_fixed_dev"] for k, v in man["splits"].items()}}, indent=1, default=str))
    return man


if __name__ == "__main__":
    main()
