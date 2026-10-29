"""Project predictors built on the vendored upstream components (third_party/ysl is never edited).

Importing this package registers the project methods in the upstream registry (`ysl.models.registry`), so the frozen
harness (`hsl.harness.run_predictor` -> `hsl.upstream.run_method` -> `registry.get`) can run them by name.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard; must be first)

from hsl.predictors import lowbudget  # noqa: F401  (registers stack_lr_wr2x)
