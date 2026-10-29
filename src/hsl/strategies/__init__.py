"""Query-selection strategies (DESIGN v3.1 §6.1–§6.3). Every selector follows the harness protocol
(`needs_view`, `static`, `reset(candidates, rng)`, `select(view: SelectorView, batch_size, rng) -> ordered genes`) and
receives ONLY a SelectorView. Knobs: config/strategies.yaml (dev defaults); protocol constants: config/protocol.yaml.

    from hsl import strategies
    sel = strategies.make("F", r=8)          # registry name + knob overrides
    harness.Trajectory(sel, "dev1", seed, ...).run()

Registry: R1 (Billmann-like, static), H, H-static, F, F-static, C, C-static (diagnostic only, not an arm), M.
Evaluation-only oracles live in hsl.oracles and must never be imported from this package (tested).
"""
from __future__ import annotations

from hsl.strategies.base import Selector, knobs  # noqa: F401
from hsl.strategies.coverage import Coverage
from hsl.strategies.fisher import Fisher
from hsl.strategies.hub import Hub
from hsl.strategies.maxvar import MaxVar
from hsl.strategies.r1 import Billmann

REGISTRY = {
    "R1": (Billmann, {}),
    "H": (Hub, {"static": False}),
    "H-static": (Hub, {"static": True}),
    "F": (Fisher, {"static": False}),
    "F-static": (Fisher, {"static": True}),
    "C": (Coverage, {"static": False}),
    "C-static": (Coverage, {"static": True}),
    "M": (MaxVar, {}),
}
ARMS = ("R1", "H", "H-static", "F", "F-static", "C", "M")   # DESIGN §9 P4 strategy set
STATIC = ("R1", "H-static", "F-static", "C-static")


def make(name: str, **overrides) -> Selector:
    cls, fixed = REGISTRY[name]
    sel = cls(**{**fixed, **overrides})
    sel.name = name
    return sel


def run_trajectory(traj) -> dict:
    """`harness.Trajectory(...).run()` plus the selector's per-round diagnostics (phase, picks, their first-pick gains
    and ranks, predicted logit q̂ of the picks, off-array flags, cluster cap, GCV λ's, fills at random, ...) written to
    `<trajectory dir>/<traj_id>.selector.json` next to the trajectory log (harness.py is frozen; review L-1). Feeds the
    Q3 'choice of queries' profiling."""
    import json

    from hsl import harness
    log = traj.run()
    sel = traj.selector
    diag = [{("round" if k == "_r" else k): v for k, v in d.items()} for d in getattr(sel, "diagnostics", [])]
    if traj.write_log:
        d = harness.TRAJ_DIR / traj.split
        d.mkdir(parents=True, exist_ok=True)
        p = d / f"{log['traj_id']}.selector.json"
        p.write_text(json.dumps({"traj_id": log["traj_id"], "selector": getattr(sel, "name", type(sel).__name__),
                                 "params": getattr(sel, "params", {}), "rounds": diag}, indent=1,
                                default=harness._json_default))
        log["selector_diagnostics_file"] = str(p)
    log["selector_diagnostics"] = diag
    return log
