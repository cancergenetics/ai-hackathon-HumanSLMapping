"""M4d (last dev round before the final): pre-stated, untuned combinations of the M4b combiners and the M4c tier-2
co-fitness kernel. Recipes fixed before any result (PI brief M4d); nothing here is tuned.

  blend_lr_hgb     0.5 * log-odds(m1_hgb) + 0.5 * log-odds(stack_lr_all). Both components are registered
                   logit_scale=True (decision_function of a logistic model / a binomial HGB), so the log-odds are
                   averaged directly (no z-scoring). One fit of each component per split; each component runs with
                   its own registered defaults (m1_hgb: K-fold additive stats, matched test views; stack_lr_all:
                   exact-LOO additive stats, full test view) and the same MethodContext seeds, so each component is
                   bit-identical to its standalone run.
  m1_hgb_t2        TIER 2 (co-reported): m1_hgb with the Hoepfner co-fitness kernel added to the within-row features
                   mk (w = tier2.COFIT_W_SELECTED = 1.0) and ext_off (2 w), exactly as tier2.stack_lr_wr2_t2 does:
                   tier2.t2_params (mk__kern / ext_off__kern) + the scoped kernel dispatcher tier2.tier2_kernels.
                   Everything else is the frozen m1_hgb recipe. cofit_w = 0 gives m1_hgb bit for bit (tested).
  blend_lr_hgb_t2  TIER 2: 0.5/0.5 log-odds blend of m1_hgb_t2 and stack_lr_all with the same co-fitness kernel
                   (stack_lr_all_t2, not registered on its own).
Component overrides (tests only): params "hgb__<key>" / "lr__<key>" are passed to the HGB / LR component.
No changes to combiner.py or tier2.py (reused by import).
"""
from __future__ import annotations

import dataclasses

import numpy as np

from ysl.models import combiner as CB
from ysl.models import tier2 as T2
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

W_HGB = 0.5
T2_DEFAULTS = {"cofit_w": T2.COFIT_W_SELECTED, "cofit_ext_ratio": T2.COFIT_EXT_RATIO, "cofit_store": T2.DEFAULT_STORE}


def _sub_ctx(ctx: MethodContext, defaults: dict, prefix: str | None = None, extra: dict | None = None) -> MethodContext:
    """Component context: the component's registered defaults, + blend-level shared settings, + prefixed overrides."""
    p = dict(defaults) | (extra or {})
    if prefix:
        p |= {k[len(prefix):]: v for k, v in ctx.params.items() if k.startswith(prefix)}
    return dataclasses.replace(ctx, params=p, info={})


def _with_t2(fn, v: ScreenView, ctx: MethodContext) -> np.ndarray:
    """Run a combiner with the tier-2 co-fitness kernel added to mk / ext_off (as tier2.stack_t2_predict)."""
    prm = T2.t2_params(dict(ctx.params))
    ctx2 = dataclasses.replace(ctx, params=prm, info={})
    with T2.tier2_kernels(prm.get("cofit_store", T2.DEFAULT_STORE)):
        S = fn(v, ctx2)
    ctx.info.update(ctx2.info)
    return S


def m1_hgb_t2(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return _with_t2(CB.m1_hgb, v, ctx)


def stack_lr_all_t2(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    return _with_t2(CB.stack_lr_all, v, ctx)


HGB_T2_DEFAULTS = dict(CB.HGB_DEFAULTS) | T2_DEFAULTS
LR_DEFAULTS = CB.COMBINER_DEFAULTS | {"test_view": "full"}       # = registered stack_lr_all defaults
LR_T2_DEFAULTS = LR_DEFAULTS | T2_DEFAULTS


def blend(v: ScreenView, ctx: MethodContext, t2: bool) -> np.ndarray:
    w = float(ctx.get("w_hgb", W_HGB))
    hgb_fn, lr_fn = (m1_hgb_t2, stack_lr_all_t2) if t2 else (CB.m1_hgb, CB.stack_lr_all)
    shared = {k: ctx.params[k] for k in T2_DEFAULTS if t2 and k in ctx.params}   # blend-level cofit settings
    ch = _sub_ctx(ctx, HGB_T2_DEFAULTS if t2 else CB.HGB_DEFAULTS, "hgb__", shared)
    cl = _sub_ctx(ctx, LR_T2_DEFAULTS if t2 else LR_DEFAULTS, "lr__", shared)
    Sh = hgb_fn(v, ch).astype(np.float64)
    Sl = lr_fn(v, cl).astype(np.float64)
    ctx.info["components"] = {"m1_hgb": ch.info, "stack_lr_all": cl.info}
    ctx.info["blend"] = {"w_hgb": w, "scale": "log-odds (both components logit_scale)",
                         "corr_hgb_lr": float(np.corrcoef(Sh.ravel()[::97], Sl.ravel()[::97])[0, 1])}
    return (w * Sh + (1.0 - w) * Sl).astype(np.float32)


register("blend_lr_hgb", defaults={"w_hgb": W_HGB}, logit_scale=True,
         description="0.5/0.5 log-odds blend of m1_hgb and stack_lr_all (pre-stated, untuned)")(
    lambda v, ctx: blend(v, ctx, t2=False))
register("m1_hgb_t2", defaults=HGB_T2_DEFAULTS, logit_scale=True,
         description="TIER 2 (co-reported): m1_hgb with the Hoepfner co-fitness kernel in mk / ext_off (tier2)")(
    m1_hgb_t2)
register("blend_lr_hgb_t2", defaults={"w_hgb": W_HGB} | T2_DEFAULTS, logit_scale=True,
         description="TIER 2: 0.5/0.5 log-odds blend of m1_hgb_t2 and stack_lr_all + co-fitness kernel")(
    lambda v, ctx: blend(v, ctx, t2=True))
