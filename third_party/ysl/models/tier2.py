"""Tier-2 co-reported variant (DESIGN.md §7): `stack_lr_wr2_t2` = stack_lr_wr2 with a chemogenomic co-fitness
kernel (Hoepfner 2014 HIP/HOP; `ysl.features.cofit`) added to
  * the multi-kernel fold-in `mk` (weight w = `cofit_w`), and
  * the external-only fold-in used for off-array targets `ext_off` (weight `cofit_ext_ratio` * w; 2x mirrors the
    ext/mk weight ratio of the tier-1 kernels: AP-MS 1.0 vs 0.5, coexpression 0.5 vs 0.3).
Everything else (features, cross-fitting, LR) is foldin2.stack_wr2_predict unchanged: the co-fitness kernel enters
only through foldin2's own `kern` parameterisation (params `mk__kern` / `ext_off__kern`).

Mechanism: foldin2.mk_predict looks up external kernels by name via the module-level `foldin2.ext_kernel`. The
context manager `tier2_kernels(store)` swaps in a dispatcher for the duration of one prediction call: names in
TIER2_KERNELS are served by `cofit_kernel` (same construction as foldin2.ext_kernel: transformed similarity, unit
diagonal for covered training genes, PSD projection, mean-diagonal normalisation), every other name goes to the
original function. It is always restored (try/finally; tested). foldin2.py itself is not modified.

The kernel is label-free: it reads only the co-fitness store (chemogenomic scores, no Costanzo data) and the views'
ORF names. With cofit_w = 0 the method is bit-identical to stack_lr_wr2 (zero-weight kernels are dropped by mk_predict).

Kernel weight: pre-stated grid COFIT_GRID, to be chosen on S1 dev seed 0 only, then frozen in COFIT_W_SELECTED.
Until it is chosen the method refuses to run with default parameters.
"""
from __future__ import annotations

import contextlib
import dataclasses
from pathlib import Path

import numpy as np

from ysl.features import cofit as CF
from ysl.models import baselines as B
from ysl.models import foldin2 as F2
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
# name -> (transform, scale) as in foldin2._transform: 'pos' = max(r, 0) / scale clipped to 1 (anti-correlated
# co-fitness is not used as similarity, as for coexpression)
TIER2_KERNELS = {"cofit": ("pos", 1.0)}
COFIT_GRID = (0.1, 0.3, 0.5, 1.0)        # pre-stated (m4c brief), mk weight w; ext_off weight = 2 w
# Selected on S1 dev seed 0 only (2026-10-27; results/m4c_sel_w*): pooled AUPRC 0.07855 / 0.08051 / 0.08187 / 0.08407
# for w = 0.1 / 0.3 / 0.5 / 1.0 (stack_lr_wr2 seed 0: 0.07748). Best = grid edge 1.0; the grid was NOT extended.
COFIT_W_SELECTED: float | None = 1.0
COFIT_EXT_RATIO = 2.0
DEFAULT_STORE = "data/processed/features/hoepfner"


def _store_path(store) -> str:
    p = Path(store)
    return str(p if p.is_absolute() else ROOT / p)


def cofit_kernel(v: ScreenView, name: str, store=DEFAULT_STORE) -> tuple[np.ndarray, np.ndarray]:
    """(K_TT PSD with unit diagonal for covered training genes, K_xT) / mean diagonal, as foldin2.ext_kernel, from the
    co-fitness store. Label-free."""
    tf, scale = TIER2_KERNELS[name]
    st = _store_path(store)
    S_tt = CF.block(v.train_orfs_current, v.train_orfs_current, st)
    cov = ~np.isnan(S_tt).all(1)
    K = F2._transform(S_tt, tf, scale)
    np.fill_diagonal(K, cov.astype(np.float64))
    K = F2._psd(K)
    Kx = F2._transform(CF.block(v.target_orfs_current, v.train_orfs_current, st), tf, scale)
    s = float(np.trace(K)) / max(len(K), 1)
    s = s if s > 0 else 1.0
    return K / s, Kx / s


@contextlib.contextmanager
def tier2_kernels(store=DEFAULT_STORE):
    """Make foldin2's kernel lookup serve TIER2_KERNELS names from the co-fitness store (restored on exit)."""
    orig = F2.ext_kernel

    def ext_kernel(v, name):
        if name in TIER2_KERNELS:
            return cofit_kernel(v, name, store)
        return orig(v, name)
    F2.ext_kernel = ext_kernel
    try:
        yield
    finally:
        F2.ext_kernel = orig


def t2_params(params: dict) -> dict:
    """stack_lr_wr2 params with the co-fitness kernel added to the mk and ext_off kernel sets."""
    w = params.get("cofit_w", COFIT_W_SELECTED)
    if w is None:
        raise ValueError("stack_lr_wr2_t2: cofit_w not selected yet (pre-stated grid COFIT_GRID; choose on dev seed 0)")
    w = float(w)
    r = float(params.get("cofit_ext_ratio", COFIT_EXT_RATIO))
    mk_kern = dict(F2.MK_RECIPE["kern"]) | {"cofit": w}
    ext_kern = dict(F2.EXT_RECIPE["kern"]) | {"cofit": r * w}
    return dict(params) | {"mk__kern": mk_kern, "ext_off__kern": ext_kern}


def stack_t2_predict(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    prm = t2_params(dict(ctx.params))
    ctx2 = dataclasses.replace(ctx, params=prm, info={})
    with tier2_kernels(prm.get("cofit_store", DEFAULT_STORE)):
        S = F2.stack_wr2_predict(v, ctx2)
    ctx.info.update(ctx2.info)
    return S


register("stack_lr_wr2_t2",
         defaults=B._LR_DEFAULTS | B.KNN_DEFAULTS | {"wr": list(F2.WR2_STACK_NAMES), "cofit_w": COFIT_W_SELECTED,
                                                     "cofit_ext_ratio": COFIT_EXT_RATIO, "cofit_store": DEFAULT_STORE},
         logit_scale=True,
         description="TIER 2 (co-reported, never the headline): stack_lr_wr2 with a chemogenomic co-fitness kernel "
                     "(Hoepfner 2014 HIP/HOP) added to the multi-kernel and external-only fold-ins")(stack_t2_predict)
