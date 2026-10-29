"""Project paths from config/paths.yaml (with env-var overrides). Importing this module installs the YEAST_SLs
write guard (hsl.guard) for the resolved upstream root."""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "paths.yaml"


@lru_cache(maxsize=1)
def config() -> dict:
    with open(CONFIG_FILE, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def _resolve(value: str) -> Path:
    p = Path(os.path.expanduser(value))
    return (p if p.is_absolute() else PROJECT_ROOT / p).resolve()


def yeast_sls_root() -> Path:
    v = os.environ.get("HSL_YEAST_SLS_ROOT") or os.environ.get("YSL_UPSTREAM_ROOT") or config()["yeast_sls_root"]
    return _resolve(v)


def output_root() -> Path:
    v = os.environ.get("HSL_OUTPUT_ROOT") or os.environ.get("YSL_OUTPUT_ROOT") or config().get("output_root", ".")
    return _resolve(v)


def costanzo_dataset() -> str:
    return str(config().get("costanzo_dataset", "s1"))


def lockbox_file() -> Path:
    """Upstream S1 lockbox list (read-only)."""
    return yeast_sls_root() / "config" / f"lockbox_{costanzo_dataset()}.txt"


def configure_ysl_env() -> None:
    """Export the roots for the vendored ysl package (third_party/ysl/_paths.py)."""
    os.environ.setdefault("YSL_UPSTREAM_ROOT", str(yeast_sls_root()))
    os.environ.setdefault("YSL_OUTPUT_ROOT", str(output_root()))


# side effects on import: guard first, then env for ysl
from hsl import guard as _guard  # noqa: E402

_guard.install([yeast_sls_root()])
configure_ysl_env()
