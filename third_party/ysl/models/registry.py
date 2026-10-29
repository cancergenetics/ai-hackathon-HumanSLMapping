"""Method registry. A method is `fit_predict(view: ScreenView, ctx: MethodContext) -> float32 [n_targets x n_arrays]`.

    from ysl.models.registry import register

    @register("my_method", defaults={"k": 50}, description="...")
    def my_method(view, ctx):
        ...

Every registered method is run by the real-pipeline leakage tests (tests/frozen/test_leakage_real.py), so a new
method cannot skip them. The runner resolves methods by name from this registry.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

METHODS: dict = {}


@dataclass(frozen=True)
class MethodSpec:
    name: str
    fn: Callable
    defaults: dict = field(default_factory=dict)
    description: str = ""
    logit_scale: bool = False   # True if scores are log-odds (oracle logit-sum diagnostic is meaningful)
    test_only: bool = False     # deliberately broken mutants used by tests; never run by the runner


def register(name: str, defaults: dict | None = None, description: str = "", logit_scale: bool = False,
             test_only: bool = False):
    def deco(fn):
        if name in METHODS and METHODS[name].fn is not fn:
            raise ValueError(f"method {name!r} registered twice")
        METHODS[name] = MethodSpec(name, fn, dict(defaults or {}), description, logit_scale, test_only)
        return fn
    return deco


def get(name: str) -> MethodSpec:
    load_all()
    if name not in METHODS:
        raise KeyError(f"unknown method {name!r}; registered: {sorted(METHODS)}")
    return METHODS[name]


def load_all() -> dict:
    """Import every module that registers methods."""
    import ysl.models.baselines  # noqa: F401
    import ysl.models.external  # noqa: F401  (B5 external-only; needs data/processed/features/pairs/)
    import ysl.models.propensity  # noqa: F401  (M4 query-propensity sub-model; stack_lr_prop)
    import ysl.models.completion  # noqa: F401  (M4 within-row completion; stack_lr_wr)
    import ysl.models.foldin2  # noqa: F401  (M4b fold-in v2: multi-kernel / external-only / neural; stack_lr_wr2)
    import ysl.models.combiner  # noqa: F401  (M4b integrated combiners: stack_lr_all, m1_hgb)
    import ysl.models.tier2  # noqa: F401  (M4c TIER 2 co-reported: stack_lr_wr2_t2, Hoepfner co-fitness kernel)
    import ysl.models.m4d  # noqa: F401  (M4d pre-stated blends: blend_lr_hgb, m1_hgb_t2, blend_lr_hgb_t2)
    return METHODS


def production_methods() -> list[str]:
    load_all()
    return sorted(n for n, s in METHODS.items() if not s.test_only)
