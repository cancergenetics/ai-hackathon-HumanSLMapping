"""Thin adapter to run any registered upstream (`ysl`) method on an ARBITRARY screened set and target set
(docs/recon.md). No evaluation happens here: callers get the score matrix and do their own scoring.

    import hsl
    from hsl import upstream as U
    b = U.load_s1()                                   # unblinded S1 bundle (all 3,283 query rows), read in place
    view = U.make_view(b, screened_orfs, target_orfs) # ScreenView: only the screened rows are physically present
    S, info = U.run_method(view, "stack_lr_wr2", n_threads=3)   # float32 [n_targets x n_arrays], rows = view.target_orfs

Why not `ysl.predict.costanzo_view`: it loads the BLINDED bundle and refuses upstream-lockbox ORFs as screened genes
or targets; our test split deliberately draws from that lockbox. `make_view` instead builds a final-kind
`ysl.eval.splits.Split` (which may contain lockbox ORFs) and calls the same `ScreenView.from_split` the upstream
runner uses, on the unblinded bundle. Leakage guarantees are unchanged: the view holds only the screened rows; the
method runs inside upstream's `guarded_call` (Costanzo files cannot be opened while it runs).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import time

import numpy as np

from hsl import paths


def load_s1(dataset: str | None = None):
    from ysl.data.dense import load_bundle
    return load_bundle(dataset or paths.costanzo_dataset(), unblind=True)


def make_view(bundle, screened, targets, split_seed: int = 0):
    """ScreenView with training rows = `screened` and targets = `targets` (gene-level ORF names, disjoint)."""
    from ysl.eval.splits import Split, read_lockbox
    from ysl.models.view import ScreenView
    scr = np.unique(np.asarray(list(screened), dtype=str))
    tg = np.unique(np.asarray(list(targets), dtype=str))
    if len(scr) != len(list(screened)) or len(tg) != len(list(targets)):
        raise ValueError("duplicate ORFs in screened/targets")
    if np.isin(scr, tg).any():
        raise ValueError("an ORF is both screened and a target")
    unknown = np.setdiff1d(np.concatenate([scr, tg]), bundle.all_query_orfs)
    if len(unknown):
        raise KeyError(f"not query ORFs of {bundle.dataset}: {unknown[:5].tolist()}")
    split = Split("final", int(split_seed), scr, tg, read_lockbox(bundle.dataset), bundle.dataset,
                  len(bundle.all_query_orfs))
    return ScreenView.from_split(bundle, split)


def run_method(view, method: str, model_seed: int = 17, split_seed: int = 0, n_threads: int = 3,
               n_folds: int = 10, chunk_targets: int = 250, params: dict | None = None):
    """Run registered `method` exactly as upstream's runner/predict do (registered defaults <- params; MethodContext;
    guarded_call). Unlike ysl.predict.predict, seeds are not restricted to upstream's seeds.yaml (our seeds live in
    config/seeds.yaml). Returns (scores float32 [n_targets x n_arrays], info dict)."""
    from ysl.models import registry
    from ysl.models.view import MethodContext
    from ysl.runner import guarded_call, set_thread_limits
    set_thread_limits(int(n_threads))
    spec = registry.get(method)
    if spec.test_only:
        raise ValueError(f"{method} is a test-only mutant")
    prm = dict(spec.defaults) | dict(params or {})
    ctx = MethodContext(model_seed=int(model_seed), split_seed=int(split_seed), params=prm, n_folds=int(n_folds),
                        n_threads=int(n_threads), chunk_targets=int(chunk_targets))
    t0 = time.perf_counter()
    S = np.asarray(guarded_call(spec.fn, view, ctx))
    dt = time.perf_counter() - t0
    if S.shape != (view.n_targets, view.n_arrays):
        raise AssertionError(f"{method} returned {S.shape}, expected {(view.n_targets, view.n_arrays)}")
    info = {"method": method, "params": prm, "model_seed": int(model_seed), "split_seed": int(split_seed),
            "n_folds": int(n_folds), "n_threads": int(n_threads), "fit_predict_s": dt,
            "n_screened": int(view.n), "n_targets": int(view.n_targets), "n_arrays": int(view.n_arrays),
            "finite": bool(np.isfinite(S).all()), "logit_scale": bool(spec.logit_scale), "method_info": ctx.info}
    return S, info


def construct_view(bundle, screened, targets):
    """ScreenView via the upstream private constructor `ScreenView._construct(_BundleSource(bundle), ...)` — the same
    constructor from_split/without_queries use — WITHOUT the upstream lockbox refusal. Used for dev views on
    hsl.views.dev_bundle() (Bundle.blind(T): unblinded=False, yet 23 upstream-lockbox genes sit in D). Training rows
    must exist in `bundle` (gene_rows raises otherwise); targets enter by name only."""
    from ysl.models.view import ScreenView, _BundleSource
    scr = np.unique(np.asarray(list(screened), dtype=str))
    tg = np.unique(np.asarray(list(targets), dtype=str))
    if len(scr) != len(list(screened)) or len(tg) != len(list(targets)):
        raise ValueError("duplicate ORFs in screened/targets")
    if np.isin(scr, tg).any():
        raise ValueError("an ORF is both screened and a target")
    unknown = np.setdiff1d(np.concatenate([scr, tg]), bundle.all_query_orfs)
    if len(unknown):
        raise KeyError(f"not query ORFs of {bundle.dataset}: {unknown[:5].tolist()}")
    return ScreenView._construct(_BundleSource(bundle), scr, tg)
