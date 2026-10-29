"""SelectorView leakage tests (DESIGN v3 §6.1): per-channel poisoning on a PSEUDO-test fixture drawn from D (no
unblinded data). Each channel alone must leave selections bit-identical; each test has a broken-twin view that leaks
that channel and must be caught. Positive control (adaptive selector): poisoning a visible column of a screened row
changes the selection. Static selectors: poisoning any screened row leaves selections identical."""
from __future__ import annotations

import hsl  # noqa: F401

import numpy as np
import pytest

from hsl import views

K = 10


# ----------------------------------------------------------------------------- toy selectors
def adaptive_toy(view, k=K):
    """Data-using toy: candidate score = Σ_s 1[SL(s, c)]·(1 + rate_s) + 1e-3 Σ_s ε(s, c)·v_s, rate_s = row SL rate over
    the visible columns, v_s = Σ_j ε(s, j)·SMF_j, plus (SMF_c − mean visible SMF)² (uses every visible column and
    every visible SMF value)."""
    Y = view.Y.astype(float)
    Y[Y < 0] = np.nan
    rate = np.nan_to_num(np.nanmean(Y, axis=1))
    E = np.nan_to_num(view.E.astype(np.float64))
    v = E @ np.nan_to_num(view.array_smf.astype(np.float64))
    col = view.candidate_col
    on = col >= 0
    sc = np.full(len(col), -np.inf)
    Yc = np.nan_to_num(Y[:, col[on]])
    m = float(np.nanmean(view.array_smf))                  # uses every visible SMF value
    smf_c = view.array_smf[col[on]].astype(np.float64)
    sc[on] = (Yc * (1 + rate[:, None])).sum(0) + 1e-3 * (E[:, col[on]] * v[:, None]).sum(0) \
        + np.nan_to_num((smf_c - m) ** 2)
    order = np.lexsort((np.arange(len(sc)), -sc))
    return [str(view.candidates[i]) for i in order[:k]]


def static_toy(view, k=K):
    """Static: external features + array SMF of candidates only."""
    f = view.gene_features(view.candidates)
    smf = np.where(view.candidate_col >= 0, view.array_smf[np.maximum(view.candidate_col, 0)], np.nan)
    sc = np.nan_to_num(f["log10_ppm_mueller2020"].to_numpy()) - np.nan_to_num(smf, nan=1.0)
    order = np.lexsort((np.arange(len(sc)), -sc))
    return [str(view.candidates[i]) for i in order[:k]]


def fake_static_toy(view, k=K):
    """Broken twin of a static selector: secretly uses screened rows."""
    sc = np.nan_to_num(view.E.astype(float)).sum(0)
    col = view.candidate_col
    s = np.where(col >= 0, sc[np.maximum(col, 0)], -np.inf)
    order = np.lexsort((np.arange(len(s)), -s))
    return [str(view.candidates[i]) for i in order[:k]]


# ----------------------------------------------------------------------------- broken-twin views
class LeakRows(views.SelectorView):            # (a)/(d): also loads excluded / unscreened rows
    __slots__ = ()
    EXTRA = None

    @classmethod
    def build(cls, b, screened, split):
        v = super().build(b, screened, split)
        r = b.gene_rows(cls.EXTRA)
        for k in ("E", "SP", "Y"):
            extra = np.array(getattr(b, k)[r])
            object.__setattr__(v, k, np.vstack([getattr(v, k), extra]))
        return v


class LeakCols(views.SelectorView):            # (b): no column masking
    __slots__ = ()

    @staticmethod
    def _masked_cols(b, split):
        return np.zeros(b.n_arrays, dtype=bool)


class LeakSMF(views.SelectorView):             # (c): SMF of excluded genes not masked
    __slots__ = ()

    @staticmethod
    def _smf(b, masked):
        return np.array(b.array_smf, dtype=np.float32)


# ----------------------------------------------------------------------------- poisoning helpers
def poison_rows(b, orfs, seed=0):
    rng = np.random.default_rng(seed)
    r = b.gene_rows(orfs)
    E, SP, Y = (np.array(getattr(b, k), copy=True) for k in ("E", "SP", "Y"))
    obs = Y[r] >= 0                                  # keep the measurement pattern; garble values and labels
    E[r] = np.where(obs, rng.normal(-1, 2, E[r].shape), np.nan).astype(np.float32)
    SP[r] = np.where(obs, rng.normal(0, 20, SP[r].shape), np.nan).astype(np.float32)
    Y[r] = np.where(obs, rng.integers(0, 2, Y[r].shape), -1).astype(np.int8)
    return b.with_matrices(E=E, SP=SP, Y=Y)


def poison_cols(b, rows_orfs, col_genes, seed=0):
    rng = np.random.default_rng(seed)
    r = b.gene_rows(rows_orfs)
    c = np.flatnonzero(np.isin(b.array_orfs, col_genes))
    E, SP, Y = (np.array(getattr(b, k), copy=True) for k in ("E", "SP", "Y"))
    E[np.ix_(r, c)] = rng.normal(-3, 3, (len(r), len(c))).astype(np.float32)
    SP[np.ix_(r, c)] = rng.normal(0, 30, (len(r), len(c))).astype(np.float32)
    Y[np.ix_(r, c)] = rng.integers(0, 2, (len(r), len(c))).astype(np.int8)
    return b.with_matrices(E=E, SP=SP, Y=Y)


def poison_smf(b, genes, seed=0):
    nb = b.with_matrices()
    smf = np.array(b.array_smf, copy=True)
    m = np.isin(b.array_orfs, genes)
    smf[m] = np.random.default_rng(seed).normal(0, 1e4, m.sum()).astype(np.float32)
    nb.array_smf = smf
    return nb


def select(cls, b, fx, sel=adaptive_toy):
    return sel(cls.build(b, fx["S"], fx["split"]))


# ----------------------------------------------------------------------------- the tests
def test_view_exposes_only_whitelist(dev_b, pseudo):
    v = views.SelectorView.build(dev_b, pseudo["S"], pseudo["split"])
    for forbidden in ("DMF", "DMF_sd", "strain_ids", "query_smf", "_b", "bundle", "strain_temperature"):
        assert not hasattr(v, forbidden)
    with pytest.raises(AttributeError):
        v.new_attr = 1
    assert not v.E.flags.writeable and not v.Y.flags.writeable
    assert v.Y.dtype == np.int8
    masked = np.setdiff1d(np.arange(v.n_arrays), v.visible_cols)
    assert (v.Y[:, masked] == -1).all() and np.isnan(v.E[:, masked]).all() and np.isnan(v.array_smf[masked]).all()
    assert views.cols_matching(dev_b, pseudo["split"].excluded)[masked].all()
    with pytest.raises(KeyError):
        v.gene_features(pseudo["PT"][:2])
    with pytest.raises(KeyError):
        v.similarity("coexp", pseudo["PV"][:2], v.candidates[:2])
    assert not np.isin(v.candidates, pseudo["split"].excluded).any()
    sv = v.screen_view()
    assert sv.n == len(pseudo["S"]) and np.isnan(sv.DMF).all()


def test_a_test_rows(dev_b, pseudo):
    base = select(views.SelectorView, dev_b, pseudo)
    pb = poison_rows(dev_b, np.union1d(pseudo["PT"], pseudo["PV"]))
    assert select(views.SelectorView, pb, pseudo) == base
    LeakRows.EXTRA = pseudo["PT"]
    assert select(LeakRows, pb, pseudo) != select(LeakRows, dev_b, pseudo)        # broken twin caught


def test_b_test_columns_of_screened_rows(dev_b, pseudo):
    base = select(views.SelectorView, dev_b, pseudo)
    pb = poison_cols(dev_b, pseudo["S"], np.union1d(pseudo["PT"], pseudo["PV"]))
    assert select(views.SelectorView, pb, pseudo) == base
    assert select(LeakCols, pb, pseudo) != select(LeakCols, dev_b, pseudo)


def test_c_test_smf(dev_b, pseudo):
    base = select(views.SelectorView, dev_b, pseudo)
    pb = poison_smf(dev_b, np.union1d(pseudo["PT"], pseudo["PV"]))
    assert select(views.SelectorView, pb, pseudo) == base
    assert select(LeakSMF, pb, pseudo) != select(LeakSMF, dev_b, pseudo)


def test_d_unscreened_pool_rows(dev_b, pseudo):
    base = select(views.SelectorView, dev_b, pseudo)
    unscreened = np.setdiff1d(pseudo["pool"], pseudo["S"])
    pb = poison_rows(dev_b, unscreened)
    assert select(views.SelectorView, pb, pseudo) == base
    LeakRows.EXTRA = unscreened[:200]
    assert select(LeakRows, pb, pseudo) != select(LeakRows, dev_b, pseudo)


def test_positive_control_visible_column(dev_b, pseudo):
    v = views.SelectorView.build(dev_b, pseudo["S"], pseudo["split"])
    base = adaptive_toy(v)
    on = v.candidates[(v.candidate_col >= 0)]
    target = [g for g in on if g not in base][-1]
    pb = poison_cols(dev_b, pseudo["S"], [target])
    r = pb.gene_rows(pseudo["S"])
    c = int(np.flatnonzero(pb.array_orfs == target)[0])
    Y = np.array(pb.Y, copy=True)
    Y[r, c] = 1
    pb = pb.with_matrices(Y=Y)
    new = select(views.SelectorView, pb, pseudo)
    assert new != base and target in new


def test_static_selector_ignores_screened_rows(dev_b, pseudo):
    base = select(views.SelectorView, dev_b, pseudo, static_toy)
    pb = poison_rows(dev_b, pseudo["S"])
    assert select(views.SelectorView, pb, pseudo, static_toy) == base
    # broken twin: a "static" selector that reads rows is caught by the same test
    assert select(views.SelectorView, pb, pseudo, fake_static_toy) != \
        select(views.SelectorView, dev_b, pseudo, fake_static_toy)
