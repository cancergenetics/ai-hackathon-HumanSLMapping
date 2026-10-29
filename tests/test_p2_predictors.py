"""P2 predictors (src/hsl/predictors): stack_lr_wr2x with every switch off is stack_lr_wr2 bit for bit; the feature
cache never changes results; the frozen P* recipe (config/predictor.yaml) passes prediction poisoning — garbling the
D_val rows or the unscreened pool rows leaves predictions identical, garbling screened rows changes them."""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import numpy as np
import pytest
import yaml

import hsl.predictors  # noqa: F401  (registers stack_lr_wr2x)
from hsl import paths, upstream as U, views
from hsl.predictors import lowbudget as LB


def _poison_all(b, orfs, seed):
    from test_review_p0 import _poison_all as pa
    return pa(b, orfs, seed)


def _frozen() -> dict:
    return yaml.safe_load((paths.PROJECT_ROOT / "config" / "predictor.yaml").read_text())["primary"]


def _S(n, seed=42):
    sp = views.dev_split(1)
    return sp, np.sort(np.random.default_rng(seed).choice(sp.pool, n, replace=False))


def _score(b, sp, S, method, params=None, seed=7):
    return U.run_method(U.construct_view(b, S, sp.targets), method, model_seed=seed, n_threads=2,
                        params=params or {})[0]


@pytest.mark.slow
def test_wr2x_all_off_is_stack_lr_wr2(dev_b):
    sp, S = _S(30)
    assert np.array_equal(_score(dev_b, sp, S, "stack_lr_wr2x"), _score(dev_b, sp, S, "stack_lr_wr2"))


@pytest.mark.slow
def test_feature_cache_is_transparent(dev_b):
    sp, S = _S(30, seed=3)
    prm = {"extra": ["dense", "prop"], "gi_n0": 100.0}
    ref = _score(dev_b, sp, S, "stack_lr_wr2x", prm)
    with LB.feature_cache():
        _score(dev_b, sp, S, "stack_lr_wr2x")                       # fills the cache with the parent's features
        assert np.array_equal(_score(dev_b, sp, S, "stack_lr_wr2x", prm), ref)


def test_dense_loo_matches_recomputation(dev_b):
    """Exact LOO dense statistics of training row t == dense_test of the view without t (array side) and t's
    array-side column over the other rows (query side)."""
    sp, S = _S(25, seed=9)
    v = U.construct_view(dev_b, S, sp.targets)
    d = LB.dense_loo(v, 1.0)
    for t in (0, 7, 24):
        fv = v.without_queries(v.train_orfs[[t]])
        dt = LB.dense_test(fv, 1.0)
        np.testing.assert_allclose(d["arr_mild"][t], dt["arr_mild"], rtol=1e-10)
        np.testing.assert_allclose(d["arr_negeps"][t], dt["arr_negeps"], rtol=1e-10)
        np.testing.assert_allclose(d["q_mild"][t], dt["q_mild"][0], rtol=1e-10)
        np.testing.assert_allclose(d["q_negeps"][t], dt["q_negeps"][0], rtol=1e-10)


@pytest.mark.slow
def test_frozen_recipe_prediction_poisoning(dev_b):
    fr = _frozen()
    sp, S = _S(60)
    base = _score(dev_b, sp, S, fr["method"], fr["params"])
    assert np.array_equal(_score(_poison_all(dev_b, sp.targets, 1), sp, S, fr["method"], fr["params"]), base)
    unscreened = np.setdiff1d(sp.pool, S)
    assert np.array_equal(_score(_poison_all(dev_b, unscreened, 2), sp, S, fr["method"], fr["params"]), base)
    assert not np.array_equal(_score(_poison_all(dev_b, S[:5], 3), sp, S, fr["method"], fr["params"]), base)


def test_frozen_yaml_matches_registered_method():
    from ysl.models.registry import get
    fr = _frozen()
    assert fr["method"] == LB.PSTAR_NAME and fr["params"] == LB.PSTAR_PARAMS
    d = get(fr["method"]).defaults
    assert all(d[k] == v for k, v in fr["params"].items())
    assert {k: v for k, v in d.items() if k not in fr["params"]} == \
        {k: v for k, v in get("stack_lr_wr2x").defaults.items() if k not in fr["params"]}
