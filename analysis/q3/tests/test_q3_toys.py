"""Toy tests with known answers for the Q3 analyses (no project data, no model fits).

    .venv/bin/python -m pytest analysis/q3/tests -q
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from analysis.q3 import choice, frontier, hap1, lbo, predcontrib as PC, stopping


# ============================================================================ hap1
def test_slope_recovers_planted_exponent():
    rng = np.random.default_rng(0)
    rows = []
    for lvl, N in ((25, 500), (50, 1000), (100, 2000)):
        for d in range(3):
            rows.append({"level": lvl, "pool_n": N * (1 + 0.01 * d), "n80": 0.4 * (N * (1 + 0.01 * d)) ** 0.7})
    b0, reps = hap1.slope_boot(pd.DataFrame(rows), "n80", 200, rng)
    assert b0 == pytest.approx(0.7, abs=1e-9)
    assert np.nanmax(np.abs(reps - 0.7)) < 1e-6           # noiseless -> every replicate exact


def test_translate_bounds():
    assert hap1.translate(800, 8000, 3200, 0.0) == pytest.approx(800)        # count scaling
    assert hap1.translate(800, 8000, 3200, 1.0) == pytest.approx(2000)       # fraction scaling


def test_nyeast_boot_matches_ci():
    rng = np.random.default_rng(1)
    nx = pd.DataFrame({"n80": [800.0], "n80_lo": [600.0], "n80_hi": [800.0 ** 2 / 600.0]})
    pt, reps = hap1.nyeast_boot(nx, "n80", 20000, rng)
    lo, hi = np.percentile(reps, [2.5, 97.5])
    assert pt == 800 and lo == pytest.approx(600, rel=0.03) and hi == pytest.approx(800 ** 2 / 600, rel=0.03)


def test_arm_budget_for_target_inverts_curve():
    cur = pd.DataFrame({"n": [30, 60, 120, 240], "auprc": [0.1, 0.2, 0.3, 0.4]})
    n, c = hap1.arm_budget_for_target(cur, 0.25)
    assert c == "" and n == pytest.approx(np.exp((np.log(60) + np.log(120)) / 2))
    assert hap1.arm_budget_for_target(cur, 0.5)[1] == "right"


# ============================================================================ frontier
def test_network_detected_vs_touched():
    # rows a, b (queries); columns a, b, c, t (t = test gene, excluded)
    Y = np.array([[-1, 1, 0, 1],      # a: SL with b (called in row a), SL with t (excluded)
                  [0, -1, 1, 0]])     # b: SL with c; reciprocal b->a not called
    net = frontier.Network(Y, ["a", "b"], ["a", "b", "c", "t"], ["t"])
    assert net.n_cells == 2 and net.n_pairs == 2           # {a,b}, {b,c}
    c = net.counts(["a"])
    assert c["cells_detected"] == 1 and c["pairs_detected"] == 1 and c["pairs_touched"] == 1
    c = net.counts(["c"])                                  # c never screened as a row: touches {b,c}, detects 0
    assert c["pairs_detected"] == 0 and c["pairs_touched"] == 1
    c = net.counts(["a", "b"])
    assert c["pairs_detected"] == 2 and c["pairs_touched"] == 2


# ============================================================================ lbo
class ToyBundle:
    """Minimal bundle: 12 query genes g0..g11 (all also arrays) + 3 array-only + 2 test genes as arrays."""

    def __init__(self, seed=0):
        rng = np.random.default_rng(seed)
        q = np.array([f"g{i:02d}" for i in range(12)])
        self.query_orfs = q
        self.all_query_orfs = np.r_[q, ["t0", "t1"]]
        self.all_query_orfs_current = self.all_query_orfs.copy()
        self.array_orfs = np.r_[q, ["a0", "a1", "a2", "t0", "t1"]]
        self.array_orfs_current = self.array_orfs.copy()
        Y = (rng.random((12, len(self.array_orfs))) < 0.3).astype(np.int8)
        Y[rng.random(Y.shape) < 0.1] = -1
        self.Y = Y

    def gene_rows(self, orfs):
        idx = {g: i for i, g in enumerate(self.query_orfs)}
        return np.array([idx[o] for o in orfs])


def test_fold_of_batches():
    f = lbo.fold_of_batches(30, 10, 5)
    assert list(np.unique(f)) == [0, 1, 2] and (f[:10] == 0).all() and (f[20:] == 2).all()
    f = lbo.fold_of_batches(70, 10, 5)
    assert f[60] == 1                                       # batch 6 -> fold 1 (interleaved)


def test_lbo_cells_and_controls():
    b = ToyBundle()
    sel = list(b.query_orfs)
    truth = {(g, a): b.Y[i, j] for i, g in enumerate(b.query_orfs) for j, a in enumerate(b.array_orfs)}
    seen = []

    def oracle(train, test, nf):            # positive control: returns the true labels of the held-out rows
        seen.append((set(train), set(test)))
        tg = np.sort(test)
        return np.array([[float(truth[(g, a)] == 1) for a in b.array_orfs] for g in tg])

    r = lbo.lbo_at(b, sel, 8, 0, batch=2, k_folds=4, fit_fn=oracle, excluded_cols=["t0", "t1"])
    assert r["cv_auprc"] == pytest.approx(1.0)
    cells = r["_cells"]
    S = set(sel[:8])
    cols = b.array_orfs[cells["col"]]
    assert not set(cols) & S                                # no screened column (no reciprocal available)
    assert not set(cols) & {"t0", "t1"}                     # never a test-gene column
    assert set(cells["row_gene"]) <= S                      # only screened rows
    assert all(not (tr & te) for tr, te in seen)            # held-out fold never in training
    assert set().union(*[te for _, te in seen]) == S

    rng = np.random.default_rng(3)
    r0 = lbo.lbo_at(b, sel, 12, 0, batch=2, k_folds=3,
                    fit_fn=lambda tr, te, nf: rng.random((len(te), len(b.array_orfs))), excluded_cols=["t0", "t1"])
    assert abs(r0["cv_auprc"] - r0["prevalence"]) < 0.25    # negative control: ~ prevalence

    with pytest.raises(PermissionError):
        lbo.lbo_at(b, ["t0"] + sel, 4, 0, fit_fn=oracle, excluded_cols=["t0", "t1"])


# ============================================================================ predcontrib
def test_pr_helpers():
    y = np.array([1, 0, 1, 0, 0, 1])
    s = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.1])
    assert PC.precision_at_recall(y, s, 0.5) == pytest.approx(2 / 3)
    assert PC.recall_at_precision(y, s, 0.6) == pytest.approx(2 / 3)
    tau = PC.threshold_for_recall(y, s, 0.6)
    assert tau == pytest.approx(0.7) and PC.apply_threshold(y, s, tau)[:2] == pytest.approx((2 / 3, 2 / 3))


def test_reciprocal_reference():
    # targets x, z both on-array; cells (x,z), (z,x), (x,a) ; a is array-only
    L = {"targets": np.array(["x", "z"]), "array_orfs": np.array(["x", "z", "a"]),
         "rows": np.array([0, 1, 0, 1]), "cols": np.array([1, 0, 2, 2]), "y": np.array([1, 0, 1, 0], np.int8),
         "OA": np.array([True, True, False, False])}
    i, j = PC.reciprocal_index(L)
    assert sorted(zip(i, j)) == [(0, 1), (1, 0)]
    ref = PC.reference(L)
    assert ref["n_cells"] == 2 and ref["ref_precision"] == 0 and ref["ref_recall"] == 0
    L["y"] = np.array([1, 1, 0, 0], np.int8)
    assert PC.reference(L)["ref_precision"] == 1.0


def test_topk():
    rows = np.array([0, 0, 0, 1, 1, 1])
    y = np.array([1, 0, 0, 0, 0, 1])
    s = np.array([0.9, 0.1, 0.2, 0.9, 0.8, 0.1])
    t = PC.topk(y, s, rows, ks=(1, 3))
    assert t[1]["true_sl_per_gene"] == 0.5 and t[3]["true_sl_per_gene"] == 1.0
    assert t[1]["random_expect"] == pytest.approx(1 / 3)


# ============================================================================ stopping
def test_plateau_rule_and_knee():
    ns = np.array([30, 60, 100, 160, 220, 330, 490, 660])
    cv = np.array([0.02, 0.03, 0.04, 0.05, 0.0505, 0.051, 0.0512, 0.0513])
    n, trig = stopping.plateau_stop(ns, cv, 0.05, 1)
    assert trig and n == 220
    n, trig = stopping.plateau_stop(ns, cv, 0.05, 2)
    assert n == 330
    n, trig = stopping.plateau_stop(ns, np.linspace(0.02, 0.2, 8), 0.01, 1)
    assert not trig and n == 660
    kn, kg = stopping.knee(ns, cv, 0.95)
    assert kg == 160 and 100 < kn <= 160


# ============================================================================ choice
def test_static_pick_and_rule_route_neutral():
    sc = pd.Series([5.0, 4, 3, 2, 1], index=list("abcde"))
    assert choice.static_pick(sc, 3, -1, 2) == ["a", "b", "c"]
    picks = choice.static_pick(sc, 5, 0, 2)
    assert sorted(picks) == list("abcde")
    # null features: off-array genes (no SMF) must be picked at their base rate under the route-stratified rule
    rng = np.random.default_rng(0)
    n = 2000
    F = pd.DataFrame({k: rng.standard_normal(n) for k in choice.FEATURES}, index=[f"g{i}" for i in range(n)])
    off = rng.random(n) < 0.2
    F.loc[off, "smf"] = np.nan
    s = pd.Series(choice.rule_scores(F, F.index.to_numpy(), {"smf": -1.0, "conservation": 0.5}), index=F.index)
    top = s.sort_values(ascending=False).index[:400]
    assert abs(np.isin(top, F.index[off]).mean() - off.mean()) < 0.05


def test_forward_select_finds_planted_feature():
    rng = np.random.default_rng(2)
    n = 3000
    Z = pd.DataFrame({k: rng.standard_normal(n) for k in ("a", "b", "c")})
    p = 1 / (1 + np.exp(-(-2 + 1.5 * Z["b"].to_numpy())))
    y = rng.binomial(9, p) / 9
    path = choice.forward_select(Z, y, np.full(n, 9.0), k_max=2)
    assert path[0]["features"] == ["b"] and path[0]["coef"]["b"] == pytest.approx(1.5, abs=0.15)


def test_projection_rule():
    ns = np.array([30, 60, 100, 160, 220, 330, 490, 660])
    x = 0.01 * np.log(ns)                    # unbounded log-linear growth: x(660) = 0.0649
    n, trig = stopping.projection_stop(ns, x, 0.93, w=3, n_max=660)
    # stop when log(n) >= 0.93 log(660) -> n >= 660^0.93 = 418 -> first checkpoint 490
    assert trig and n == 490
    flat = np.r_[0.01, 0.02, 0.03, np.full(5, 0.031)]
    n, trig = stopping.projection_stop(ns, flat, 0.95, w=3)
    assert trig and n == 220
    assert stopping.stop_n(ns, x, "plateau", 0.5, 1)[0] == 60
