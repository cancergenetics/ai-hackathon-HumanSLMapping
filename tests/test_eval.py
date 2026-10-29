"""Known-answer tests of the frozen evaluation module (hsl.eval). Toy data only."""
from __future__ import annotations

import hsl  # noqa: F401

import numpy as np
import pytest
from scipy import stats
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score as sk_ap

from hsl import eval as E

TOL = 1e-12


# ----------------------------------------------------------------------------- broken twins (must be caught)
def ap_stable_sort_mutant(y, s, w=None):
    """WRONG: per-cell precision after a stable descending sort (ties broken by input order)."""
    y = np.asarray(y, bool)
    o = np.argsort(-np.asarray(s, float), kind="stable")
    yy = y[o]
    tp = np.cumsum(yy)
    prec = tp / np.arange(1, len(yy) + 1)
    return float(prec[yy].sum() / yy.sum())


def ap_desc_tiebreak_mutant(y, s, w=None):
    """WRONG: ties broken with positives first (optimistic)."""
    y = np.asarray(y, bool)
    o = np.lexsort((~y, -np.asarray(s, float)))
    yy = y[o]
    tp = np.cumsum(yy)
    return float((tp / np.arange(1, len(yy) + 1))[yy].sum() / yy.sum())


def ap_battery(fn) -> list[str]:
    """The AP property checks; returns the names of failed checks."""
    fails = []
    lit = [([1, 0, 1, 0], [4, 3, 2, 1], 5 / 6), ([1, 0], [1, 1], 0.5), ([1, 0, 1], [3, 2, 2], 5 / 6),
           ([0, 1, 1, 0], [2, 2, 2, 2], 0.5), ([0, 0, 1], [3, 2, 1], 1 / 3)]
    for y, s, want in lit:
        if abs(fn(np.array(y), np.array(s, float)) - want) > TOL:
            fails.append(f"literal {y} {s}")
    rng = np.random.default_rng(1)
    y = rng.random(2000) < 0.05
    y[0] = True
    if abs(fn(y, np.zeros(2000)) - y.mean()) > TOL:
        fails.append("all tied -> prevalence")
    s = np.round(rng.normal(size=2000) + y, 1)
    p = rng.permutation(2000)
    if fn(y, s) != fn(y[p], s[p]):
        fails.append("cell-order permutation")
    if abs(fn(y, s) - sk_ap(y, s)) > TOL:
        fails.append("sklearn agreement")
    return fails


def test_ap_battery_passes_on_eval_module():
    assert ap_battery(E.average_precision) == []


@pytest.mark.parametrize("mutant", [ap_stable_sort_mutant, ap_desc_tiebreak_mutant])
def test_ap_battery_catches_broken_twins(mutant):
    assert ap_battery(mutant), "a tie-breaking AP implementation passed the battery"


def test_sklearn_agreement_1e5_cells_heavy_ties():
    rng = np.random.default_rng(2)
    n = 100_000
    y = rng.random(n) < 0.0114
    s = np.round(rng.normal(size=n) + 0.8 * y, 2)          # many ties
    assert abs(E.average_precision(y, s) - sk_ap(y, s)) < TOL
    s2 = rng.normal(size=n).astype(np.float32).astype(float) + y
    assert abs(E.average_precision(y, s2) - sk_ap(y, s2)) < TOL


def test_weighted_ap_engine_exact():
    rng = np.random.default_rng(3)
    n, Q = 30_000, 300
    y = rng.random(n) < 0.02
    s = np.round(rng.normal(size=n) + y, 2)
    q = rng.integers(0, Q, n)
    eng = E.WeightedAP(y, s, q, Q)
    assert abs(eng.ap(np.ones(Q))[0] - sk_ap(y, s)) < TOL                 # unweighted == sklearn
    for seed in range(3):
        w = np.random.default_rng(seed).integers(0, 4, Q)
        idx = np.repeat(np.arange(n), w[q])                               # duplicated cells
        want = sk_ap(y[idx], s[idx])
        assert abs(eng.ap(w)[0] - want) < TOL
        assert abs(eng.ap_gather(w) - want) < TOL
        assert abs(E.average_precision(y, s, w[q]) - want) < TOL
        assert abs(E.SortedAP(y, s).ap(w[q]) - want) < TOL
    W = E.bootstrap_weights(rng.integers(0, 5, Q), 50, 7)
    batch = eng.ap(W)
    assert np.allclose(batch, [eng.ap_gather(w) for w in W], atol=TOL, rtol=0)


def test_bootstrap_weights_stratified():
    strata = np.repeat([0, 1, 2], [5, 10, 20])
    W = E.bootstrap_weights(strata, 200, 0)
    for h, n in zip(range(3), (5, 10, 20)):
        assert (W[:, strata == h].sum(1) == n).all()
    assert np.array_equal(W, E.bootstrap_weights(strata, 200, 0))


def test_macro_precision_at_k_auroc():
    y = np.array([1, 0, 0, 0, 1, 0, 0, 0])
    s = np.array([0.9, 0.8, 0.1, 0.0, 0.5, 0.5, 0.2, 0.1])
    q = np.array([0, 0, 0, 1, 1, 1, 2, 2])
    m = E.macro_ap(y, s, q)
    assert m["n_zero_pos"] == 1 and m["n_used"] == 2
    assert abs(m["macro"] - (1.0 + 0.5) / 2) < TOL                         # q1: tie of pos and neg -> 0.5
    p1 = E.precision_at_k(y, s, q, 1)
    assert p1[0] == 1.0 and abs(p1[1] - 0.5) < TOL and p1[2] == 0.0          # tie at cut-off counted in expectation
    assert abs(E.auroc(y, s) - __import__("sklearn.metrics").metrics.roc_auc_score(y, s)) < TOL


def _toy_population(rng, Q, m, base=0.04):
    """Queries with heterogeneous SL rates; scores informative."""
    rate = np.clip(rng.gamma(2.0, base / 2.0, size=Q), 0.001, 0.5)
    y = rng.random((Q, m)) < rate[:, None]
    s = np.round(rng.normal(size=(Q, m)) + 1.2 * y + 3 * rate[:, None], 2)
    return y.ravel(), s.ravel(), np.repeat(np.arange(Q), m), rate


def test_bootstrap_ci_coverage_toy():
    """Wald CI from the stratified query bootstrap covers the population AP ≈ 95% of the time."""
    rng = np.random.default_rng(11)
    yb, sb, _, _ = _toy_population(rng, 40_000, 60)
    truth = E.average_precision(yb, sb)
    hits, n_draw = 0, 200
    for d in range(n_draw):
        r = np.random.default_rng(1000 + d)
        y, s, q, rate = _toy_population(r, 200, 60)
        strata = (rate > np.median(rate)).astype(int)
        eng = E.WeightedAP(y, s, q, 200)
        est = eng.ap(np.ones(200))[0]
        reps = eng.ap(E.bootstrap_weights(strata, 300, r))
        se = reps.std(ddof=1)
        hits += abs(est - truth) <= 1.96 * se
    cov = hits / n_draw
    assert 0.89 <= cov <= 0.99, cov


def test_welch_reduces_to_welch_t_test_without_query_variance():
    rng = np.random.default_rng(4)
    a, b = rng.normal(0.30, 0.01, 12), rng.normal(0.29, 0.02, 20)
    boot_a = np.tile(a, (50, 1))          # identical replicates -> V_Q = 0
    boot_b = np.tile(b, (50, 1))
    r = E.welch(boot_a, a, boot_b, b)
    t = stats.ttest_ind(a, b, equal_var=False, alternative="greater")
    assert abs(r["V_Q"]) < 1e-20
    assert abs(r["t"] - t.statistic) < 1e-9 and abs(r["df"] - t.df) < 1e-6 and abs(r["p_one_sided"] - t.pvalue) < 1e-9


def test_welch_adds_query_variance():
    rng = np.random.default_rng(5)
    a, b = rng.normal(0.3, 0.01, 10), rng.normal(0.29, 0.01, 20)
    noise = rng.normal(0, 0.02, (400, 1))
    r = E.welch(a + noise, a, b + noise * 0.5, b)
    assert abs(r["V_Q"] - np.var(noise[:, 0] * 0.5, ddof=1)) < 1e-12
    assert abs(r["se"] ** 2 - (r["V_Q"] + r["a"] + r["b"])) < 1e-15


def test_paired_reduces_to_one_sample_t():
    rng = np.random.default_rng(6)
    p = rng.normal(0.3, 0.02, 10)
    c = p + rng.normal(0.005, 0.003, 10)
    r = E.paired(np.tile(c, (30, 1)), c, np.tile(p, (30, 1)), p)
    t = stats.ttest_1samp(c - p, 0, alternative="greater")
    assert abs(r["p_one_sided"] - t.pvalue) < 1e-9 and abs(r["df"] - 9) < 1e-9


def test_holm_known_values():
    r = E.holm([0.01, 0.04, 0.03, 0.005], alpha=0.05)
    assert np.allclose(r["p_adj"], [0.03, 0.06, 0.06, 0.02])
    assert r["reject"].tolist() == [True, False, False, True]


def test_pava_matches_sklearn():
    rng = np.random.default_rng(7)
    for _ in range(20):
        y = rng.normal(size=15).cumsum() * 0.1 + rng.normal(size=15)
        w = rng.random(15) + 0.1
        want = IsotonicRegression().fit(np.arange(15), y, sample_weight=w).predict(np.arange(15))
        assert np.allclose(E.pava(y, w), want)


def test_n_req_known_curve_and_censoring():
    grid = [10, 30, 100, 300, 1000]
    curve = np.log10(grid)                                    # linear in log n -> exact inversion
    assert abs(E.n_req(np.log10(50), grid, curve).n - 50) < 1e-9
    assert E.n_req(np.log10(10), grid, curve) == E.NReq(10.0, "")
    assert E.n_req(0.5, grid, curve) == E.NReq(0.0, "left")
    assert E.n_req(3.5, grid, curve) == E.NReq(float("inf"), "right")
    # non-monotone curve: PAVA pools [2, 1.5] -> 1.75 flat on (30, 100); reaching it gives the LEFT endpoint
    g2, c2 = [10, 30, 100, 300], [1.0, 2.0, 1.5, 3.0]
    assert E.n_req(1.75, g2, c2) == E.NReq(30.0, "")
    r = E.n_req(1.375, g2, c2)                               # between f(10)=1 and f(30)=1.75: halfway in log n
    assert abs(r.n - np.sqrt(10 * 30)) < 1e-9 and r.censored == ""


def test_gap_closure_and_g_star():
    assert E.gap_closure(0.25, 0.2, 0.3) == pytest.approx(0.5)
    assert E.g_star(0.2, 0.3, 0.05) == pytest.approx(1 - 0.05 * 0.3 / 0.1)


def test_joint_resample_and_bands():
    A = np.random.default_rng(8).normal(size=(100, 6, 4))
    M1, M2 = E.joint_resample(A, 3), E.joint_resample(A, 3)
    assert M1.shape == (100, 4) and np.array_equal(M1, M2)
    b = E.bonferroni_bands([0.1, 0.2], [0.01, 0.01])
    assert (b["crit"] > 1.96).all()


def test_replicate_engine_shared_weights():
    rng = np.random.default_rng(9)
    Q = 50
    W = E.bootstrap_weights(np.zeros(Q, int), 40, 1)
    eng = E.ReplicateEngine(W, (10, 20), {"rand": 2, "strat": 1})
    for arm, R in (("rand", 2), ("strat", 1)):
        for r in range(R):
            for n in (10, 20):
                y = rng.random(Q * 20) < 0.05
                eng.add(arm, r, n, y, rng.normal(size=Q * 20) + y, np.repeat(np.arange(Q), 20))
    boot, est = eng.at("rand", 20)
    assert boot.shape == (40, 2) and np.isfinite(boot).all() and est.shape == (2,)
