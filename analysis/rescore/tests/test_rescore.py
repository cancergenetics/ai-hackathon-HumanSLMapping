"""Toy tests with known answers for the independent rescorer (DESIGN §11). Run:
    .venv/bin/python -m pytest -p no:cacheprovider analysis/rescore/tests -q
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats
from sklearn.metrics import average_precision_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import rescore as RS  # noqa: E402


# ------------------------------------------------------------------ independence
def test_no_project_imports():
    code = ("import sys; sys.path.insert(0, %r); import rescore; "
            "bad = [m for m in sys.modules if m.split('.')[0] in ('hsl', 'ysl')]; "
            "assert not bad, bad; print('ok')") % str(HERE.parent)
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr
    src = (HERE.parent / "rescore.py").read_text()
    imports = [ln.strip() for ln in src.splitlines() if ln.strip().startswith(("import ", "from "))]
    roots = {ln.split()[1].split(".")[0] for ln in imports}
    assert roots <= {"__future__", "os", "argparse", "glob", "json", "re", "resource", "sys", "time", "zlib", "pathlib",
                     "numpy", "pandas", "scipy", "sklearn"}, roots
    assert "sys.path" not in src


# ------------------------------------------------------------------ AP
def test_ap_literals():
    assert RS.ap_point([1, 0, 1, 0], [4, 3, 2, 1]) == pytest.approx((1 + 2 / 3) / 2)
    assert RS.ap_point([0, 1], [0.5, 0.5]) == pytest.approx(0.5)          # one tie group -> prevalence
    y = np.r_[np.ones(3), np.zeros(17)]
    assert RS.ap_point(y, np.zeros(20)) == pytest.approx(0.15)           # all tied -> prevalence


@pytest.mark.parametrize("seed", range(20))
def test_ap_equals_sklearn_with_ties(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(50, 3000))
    y = (rng.random(n) < rng.uniform(0.01, 0.4)).astype(int)
    y[0] = 1
    s = np.round(rng.normal(size=n) + y * rng.uniform(0, 2), int(rng.integers(0, 3)))   # heavy ties
    assert abs(RS.ap_point(y, s) - average_precision_score(y, s)) < 1e-12
    e = RS.ReplicateAP(y, s, rng.integers(0, 7, n), 7)
    assert abs(e.ap(np.ones(7))[0] - average_precision_score(y, s)) < 1e-12


def test_ap_float32_scores_like_saved_predictions():
    rng = np.random.default_rng(5)
    y = (rng.random(20000) < 0.02).astype(int)
    s = (rng.normal(size=20000) + y).astype(np.float32)
    s[::7] = s[0]
    assert abs(RS.ap_point(y, s) - average_precision_score(y, s)) < 1e-12


def test_ap_cell_order_invariant():
    rng = np.random.default_rng(1)
    y = (rng.random(500) < 0.1).astype(int)
    s = np.round(rng.normal(size=500), 1)
    p = rng.permutation(500)
    assert RS.ap_point(y, s) == pytest.approx(RS.ap_point(y[p], s[p]), abs=1e-15)


@pytest.mark.parametrize("seed", range(10))
def test_weighted_ap_equals_duplicated_cells(seed):
    """Integer query weights == every cell of query q duplicated w_q times (sklearn on the expanded data)."""
    rng = np.random.default_rng(100 + seed)
    nq = 9
    n = 800
    rows = rng.integers(0, nq, n)
    y = (rng.random(n) < 0.08).astype(int)
    y[:nq] = 1
    rows[:nq] = np.arange(nq)
    s = np.round(rng.normal(size=n) + y, 1)
    W = rng.integers(0, 4, size=(5, nq))
    W[:, 0] = 1
    e = RS.ReplicateAP(y, s, rows, nq)
    got = e.ap(W, chunk=2)
    for b in range(5):
        rep = W[b][rows]
        exp = average_precision_score(np.repeat(y, rep), np.repeat(s, rep))
        assert abs(got[b] - exp) < 1e-12
    # sklearn's sample_weight path agrees too
    assert abs(got[0] - average_precision_score(y, s, sample_weight=W[0][rows])) < 1e-12


def test_stratified_weights():
    strata = np.repeat(np.arange(11), [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
    W = RS.stratified_query_weights(strata, 300, np.random.default_rng(0))
    for h in range(11):
        assert np.all(W[:, strata == h].sum(axis=1) == (strata == h).sum())     # stratum sizes preserved
    assert W.min() >= 0 and np.allclose(W.mean(), 1.0, atol=0.05)


# ------------------------------------------------------------------ Holm
def test_holm_known():
    p = np.array([0.01, 0.04, 0.03, 0.005, 0.5, 0.0001])
    adj, rej = RS.holm(p, 0.025)
    # sorted: 0.0001*6, 0.005*5, 0.01*4, 0.03*3, 0.04*2, 0.5*1 -> cummax
    exp_sorted = np.maximum.accumulate([0.0006, 0.025, 0.04, 0.09, 0.09, 0.5])
    order = np.argsort(p)
    assert np.allclose(adj[order], exp_sorted)
    assert list(rej) == [False, False, False, True, False, True]


def test_holm_caps_and_ties():
    adj, rej = RS.holm([0.3, 0.3, 0.3], 0.025)
    assert np.allclose(adj, 0.9)
    adj, _ = RS.holm([0.6, 0.7], 0.025)
    assert np.allclose(adj, [1.0, 1.0])


# ------------------------------------------------------------------ Welch
def test_welch_toy_hand_computed():
    pa = np.array([0.050, 0.052, 0.054])           # var 4e-6
    pb = np.array([0.040, 0.041, 0.042, 0.043])    # var 1.6667e-6
    ca, cb = 0.002, 0.001
    t = np.arange(400) % 2 * 2 - 1.0                # +-1, exact variance known
    ra = 0.052 + ca * t[:, None] + np.zeros((1, 3))
    rb = 0.0415 + cb * (-t)[:, None] + np.zeros((1, 4))
    w = RS.welch(ra, pa, rb, pb)
    VQ = np.var((ca + cb) * t, ddof=1)
    ta, tb = 4e-6 / 3, (np.var(pb, ddof=1)) / 4
    se2 = VQ + ta + tb
    nu = se2 ** 2 / (ta ** 2 / 2 + tb ** 2 / 3)
    d = pa.mean() - pb.mean()
    assert w["delta"] == pytest.approx(d)
    assert w["V_Q"] == pytest.approx(VQ)
    assert w["se"] == pytest.approx(np.sqrt(se2))
    assert w["df"] == pytest.approx(nu)
    assert w["p"] == pytest.approx(1 - stats.t.cdf(d / np.sqrt(se2), nu))
    q = stats.t.ppf(0.975, nu)
    assert (w["lo"], w["hi"]) == pytest.approx((d - q * np.sqrt(se2), d + q * np.sqrt(se2)))


def test_welch_scaled_matches_contrast():
    pa = np.array([0.05, 0.06])
    pb = np.array([0.05, 0.055, 0.06])
    z = np.zeros((10, 1))
    w = RS.welch(z + pa, pa, z + pb, pb, scale_b=0.95)
    assert w["V_Q"] == pytest.approx(0.0, abs=1e-20)
    ta, tb = np.var(pa, ddof=1) / 2, 0.95 ** 2 * np.var(pb, ddof=1) / 3
    assert w["delta"] == pytest.approx(pa.mean() - 0.95 * pb.mean())
    assert w["se"] == pytest.approx(np.sqrt(ta + tb))
    assert w["df"] == pytest.approx((ta + tb) ** 2 / (ta ** 2 / 1 + tb ** 2 / 2))


def test_welch_large_R_reduces_to_normal():
    rng = np.random.default_rng(3)
    pa, pb = rng.normal(1.0, 0.1, 5000), rng.normal(0.9, 0.1, 5000)
    w = RS.welch(np.zeros((50, 1)), pa, np.zeros((50, 1)), pb)
    z = (pa.mean() - pb.mean()) / np.sqrt(pa.var(ddof=1) / 5000 + pb.var(ddof=1) / 5000)
    assert w["p"] == pytest.approx(stats.norm.sf(z), rel=1e-3, abs=1e-300)


# ------------------------------------------------------------------ PAVA and inversion (n_80 / n_req)
def test_pava():
    assert np.allclose(RS.pava([1, 3, 2, 4]), [1, 2.5, 2.5, 4])
    assert np.allclose(RS.pava([3, 2, 1]), [2, 2, 2])
    assert np.allclose(RS.pava([1, 2, 3]), [1, 2, 3])
    v = np.random.default_rng(0).normal(size=30)
    f = RS.pava(v)
    assert np.all(np.diff(f) >= -1e-15) and f.sum() == pytest.approx(v.sum())


def test_invert_known_curve():
    grid = [10, 100, 1000]
    f = [0.1, 0.2, 0.3]
    v, c = RS.invert_curve(0.15, grid, f)
    assert c == "" and v == pytest.approx(np.sqrt(10 * 100))          # midpoint in log n
    assert RS.invert_curve(0.2, grid, f) == (pytest.approx(100.0), "")
    assert RS.invert_curve(0.05, grid, f) == (0.0, "left")
    assert RS.invert_curve(0.1, grid, f) == (0.0, "left")               # reached at the first grid point
    assert RS.invert_curve(0.31, grid, f)[1] == "right"
    assert np.isinf(RS.invert_curve(0.31, grid, f)[0])


def test_invert_flat_segment_left_endpoint():
    grid = [10, 20, 40, 80]
    fit = RS.pava([0.1, 0.3, 0.2, 0.4])         # -> 0.1, 0.25, 0.25, 0.4
    v, c = RS.invert_curve(0.25, grid, fit)
    assert v == pytest.approx(20.0) and c == ""


def test_pct_ci_with_censoring():
    lo, hi = RS.pct_ci(np.r_[np.full(90, 100.0), np.full(10, np.inf)])
    assert lo == 100.0 and np.isinf(hi)
