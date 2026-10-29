"""Independent review of the P2 frozen P* (stack_lr_wr2_lb; docs/reviews/p2_predictor_review.md). DEV DATA ONLY.

Covers: store provenance of the new aprior inputs, own-label exclusion of the fold-view array prior, poisoning at small
n (targets / unscreened rows / one screened cell), and degenerate inputs (n = 10 with zero strong SLs; off-array
targets)."""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import numpy as np
import pytest
import yaml

import hsl.predictors  # noqa: F401
from hsl import guard, harness as H, paths, upstream as U, views
from hsl.predictors import lowbudget as LB


def _frozen() -> dict:
    return yaml.safe_load((paths.PROJECT_ROOT / "config" / "predictor.yaml").read_text())["primary"]


def _poison_all(b, orfs, seed):
    from test_review_p0 import _poison_all as pa
    return pa(b, orfs, seed)


def _score(b, sp, S, seed=7):
    fr = _frozen()
    return U.run_method(U.construct_view(b, S, sp.targets), fr["method"], model_seed=seed, n_threads=2,
                        params=fr["params"])[0]


def _S(n, seed):
    sp = views.dev_split(1)
    return sp, np.sort(np.random.default_rng(seed).choice(sp.pool, n, replace=False))


# ----------------------------------------------------------------------------- provenance
def _reads_of_pstar():
    fr = _frozen()
    _, S = _S(30, 0)
    LB._ADESIGN.clear()                                   # force the gene-scalar stores to be (re)opened
    from ysl.features import genes as G
    for fn in (G.gene_table_current, G._esm_current, G._coexp):
        fn.cache_clear()
    with guard.record_reads() as ps:
        H.run_predictor(H.RunSpec("dev1", tuple(S)), fr["method"], fr["params"], seeds={"model_seed": 1},
                        save=False, ledger=False)
    return ps


@pytest.mark.slow
def test_pstar_opens_only_admissible_stores():
    """Every store P* opens is an admissible (label-free, GI/text/GO-free) store already cleared for m1_hgb or
    stack_lr_wr2."""
    from test_provenance import _uncovered
    ps = _reads_of_pstar()
    bad = [p for p in _uncovered(ps, "m1_hgb") if p in _uncovered(ps, "stack_lr_wr2")]
    assert not bad, bad


@pytest.mark.slow
# REVIEW P2 finding (fixed by PI): used_by now lists stack_lr_wr2_lb.
def test_pstar_listed_as_consumer_in_provenance():
    from test_provenance import _uncovered
    assert not _uncovered(_reads_of_pstar(), LB.PSTAR_NAME)


# ----------------------------------------------------------------------------- aprior own-label exclusion
def test_fold_view_array_prior_ignores_fold_rows(dev_b):
    """The aprior column of training rows in fold k comes from the fold view (rows outside k): garbling the fold
    rows' data must not change it; garbling one non-fold screened row must."""
    sp, S = _S(40, 11)
    fold = S[:4]
    v = U.construct_view(dev_b, S, sp.targets)
    ap = LB.array_prior(v.without_queries(fold))
    vb = U.construct_view(_poison_all(dev_b, fold, 5), S, sp.targets)
    assert np.array_equal(LB.array_prior(vb.without_queries(fold)), ap)
    assert not np.array_equal(LB.array_prior(vb), LB.array_prior(v))        # full view sees the fold rows
    vc = U.construct_view(_poison_all(dev_b, S[10:11], 6), S, sp.targets)
    assert not np.array_equal(LB.array_prior(vc.without_queries(fold)), ap)


def test_one_visible_cell_changes_array_side_features(dev_b):
    """Positive control at cell level: flipping one screened cell (label + epsilon) changes that array's
    dense statistic and the GLM prior."""
    sp, S = _S(40, 12)
    v = U.construct_view(dev_b, S, sp.targets)
    r = dev_b.gene_rows(S[:1])[0]
    j = int(np.flatnonzero(~np.isnan(dev_b.E[r]) & (dev_b.Y[r] == 0))[0])
    M = {k: np.array(getattr(dev_b, k), copy=True) for k in ("E", "SP", "Y", "DMF", "DMF_sd")}
    M["E"][r, j], M["Y"][r, j] = -0.9, 1
    v2 = U.construct_view(dev_b.with_matrices(**M), S, sp.targets)
    d1, d2 = LB.dense_test(v, 1.0), LB.dense_test(v2, 1.0)
    assert d1["arr_mild"][j] != d2["arr_mild"][j]
    dj = abs(d1["arr_mild"][j] - d2["arr_mild"][j])
    assert dj > 100 * np.abs(np.delete(d1["arr_mild"] - d2["arr_mild"], j)).max()   # others: global prior only
    assert not np.array_equal(LB.array_prior(v), LB.array_prior(v2))


# ----------------------------------------------------------------------------- poisoning at n = 30
@pytest.mark.slow
def test_pstar_poisoning_n30(dev_b):
    sp, S = _S(30, 13)
    base = _score(dev_b, sp, S)
    assert np.array_equal(_score(_poison_all(dev_b, sp.targets, 1), sp, S), base)                  # D_val rows
    assert np.array_equal(_score(_poison_all(dev_b, np.setdiff1d(sp.pool, S), 2), sp, S), base)    # unscreened
    assert not np.array_equal(_score(_poison_all(dev_b, S[:1], 3), sp, S), base)                   # one screened row


# ----------------------------------------------------------------------------- degenerate inputs
def test_array_prior_finite_on_zero_sl_view(dev_b):
    sp = views.dev_split(1)
    zero = np.sort(sp.pool[views.query_sl_counts(dev_b, sp.pool) == 0])
    v = U.construct_view(dev_b, zero[:10], sp.targets)
    assert np.nansum(v.Y) == 0 and np.isfinite(LB.array_prior(v)).all()


@pytest.mark.slow
def test_pstar_n10_zero_sl_fails_like_parent(dev_b):
    """All-zero-SL screened set: the PARENT already fails (upstream baselines.shrunk_rate divides by the global
    rate 0); P* must not add a new failure mode -> same exception type."""
    sp = views.dev_split(1)
    zero = np.sort(sp.pool[views.query_sl_counts(dev_b, sp.pool) == 0])
    S = np.sort(np.random.default_rng(14).choice(zero, 10, replace=False))
    v = U.construct_view(dev_b, S, sp.targets)
    with pytest.raises(ZeroDivisionError):
        U.run_method(v, "stack_lr_wr2", model_seed=7, n_threads=2)
    with pytest.raises(ZeroDivisionError):
        _score(dev_b, sp, S)


@pytest.mark.slow
def test_pstar_n10_random_offarray_targets_finite(dev_b):
    sp, S = _S(10, 15)
    v = U.construct_view(dev_b, S, sp.targets)
    s = _score(dev_b, sp, S)
    off = v.target_col < 0
    assert off.any() and np.isfinite(s).all()
