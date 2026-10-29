"""Read-only access to the project configuration files (config/*.yaml). Single source for seeds, budgets and
protocol constants; code never hardcodes these values."""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import yaml

from hsl.paths import PROJECT_ROOT

CONFIG_DIR = PROJECT_ROOT / "config"


@lru_cache(maxsize=None)
def _load(name: str) -> dict:
    with open(CONFIG_DIR / f"{name}.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def seeds() -> dict:
    return _load("seeds")


def budgets() -> dict:
    return _load("budgets")


def protocol() -> dict:
    return _load("protocol")


def seed(name: str) -> int:
    s = seeds()
    if name not in s:
        raise KeyError(f"seed {name!r} is not listed in config/seeds.yaml (append it before first use)")
    return int(s[name])


def trajectory_seeds(kind: str, R: int) -> list[np.random.SeedSequence]:
    """SeedSequence(base).spawn(R); base = dev_traj_base_seed or final_traj_base_seed."""
    if kind not in ("dev", "final"):
        raise ValueError(kind)
    return np.random.SeedSequence(seed(f"{kind}_traj_base_seed")).spawn(int(R))
