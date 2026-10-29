"""[humanSLmapping PATCH — not part of upstream YEAST_SLs] Path roots for the vendored `ysl` package.

Upstream derives ROOT from the package location (`Path(__file__).parents[N]` = the YEAST_SLs checkout). In the
vendored copy that would point at humanSLmapping, so every upstream ROOT constant is replaced by one of:

  UPSTREAM_ROOT  read-only root holding upstream `data/` and `config/` (the YEAST_SLs checkout). Resolved from the
                 env var YSL_UPSTREAM_ROOT, else from humanSLmapping `config/paths.yaml` via `hsl.paths`.
  OUTPUT_ROOT    writable root for anything the upstream runner would write (results, ledger, predictions);
                 env var YSL_OUTPUT_ROOT, else `hsl.paths`. Never inside UPSTREAM_ROOT.

Resolving through `hsl.paths` also installs the YEAST_SLs write guard (`hsl.guard`), so any use of the vendored
package is guarded.
"""
from __future__ import annotations

import os
from pathlib import Path


def _resolve() -> tuple[Path, Path]:
    up, out = os.environ.get("YSL_UPSTREAM_ROOT"), os.environ.get("YSL_OUTPUT_ROOT")
    try:  # always install the guard when hsl is importable
        from hsl import paths as _hp  # noqa: WPS433
        up = up or str(_hp.yeast_sls_root())
        out = out or str(_hp.output_root())
    except ImportError:
        if not (up and out):
            raise RuntimeError("set YSL_UPSTREAM_ROOT and YSL_OUTPUT_ROOT, or make the `hsl` package importable")
    up_p, out_p = Path(up).resolve(), Path(out).resolve()
    if out_p == up_p or up_p in out_p.parents:
        raise RuntimeError(f"OUTPUT_ROOT {out_p} must not lie inside the read-only UPSTREAM_ROOT {up_p}")
    return up_p, out_p


UPSTREAM_ROOT, OUTPUT_ROOT = _resolve()
