"""Independent P0 code-review checks (docs/reviews/p0_code_review.md). Written by a reviewer who did not write the
code under test: brute-force reference implementations, adversarial toy cases, and dev-only leakage checks. Known
issues are pinned as strict xfails (they flip to failures once fixed, so the xfail marker must then be removed).
Dev data only: no test-query rows or labels are read (dev_bundle has none; final E_fixed checked by column names)."""
from __future__ import annotations

import hsl  # noqa: F401

import warnings

import numpy as np
import pytest
from scipy import stats
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score as sk_ap

from hsl import config, features, harness as H, views
from hsl import eval as E

TOL = 1e-12


# ============================================================================ reference implementations
def ref_ap(y, s, w=None) -> float:
    """AP from first principles: loop over DISTINCT thresholds (descending), ΔRecall × Precision."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    P = float((w * y).sum())
    if P <= 0:
        return float("nan")
    out, prev = 0.0, 0.0
    for t in np.unique(s)[::-1]:
        m = s >= t
        tp, fp = float((w * (y & m)).sum()), float((w * (~y & m)).sum())
        out += (tp / P - prev) * (tp / (tp + fp) if tp + fp > 0 else 0.0)
        prev = tp / P
    return out


def ref_holm(p, alpha):
    p = np.asarray(p, float)
    rej = np.zeros(len(p), bool)
    for i, j in enumerate(np.argsort(p, kind="mergesort")):
        if p[j] <= alpha / (len(p) - i):
            rej[j] = True
        else:
            break
    return rej


def ref_n_req(t, n, c):
    """Isotonic fit (sklearn), piecewise linear in log n, first point of a fine grid reaching t."""
    f = IsotonicRegression().fit(np.arange(len(c)), c).predict(np.arange(len(c)))
    if t < f[0] - 1e-15:
        return 0.0, "left"
    if t > f[-1] + 1e-15:
        return float("inf"), "right"
    x = np.log(np.asarray(n, float))
    fine = np.linspace(x[0], x[-1], 200_001)
    k = int(np.argmax(np.interp(fine, x, f) >= t - 1e-13))
    return float(np.exp(fine[k])), ""


# ============================================================================ 1. AP engine
def test_pooled_and_weighted_ap_vs_reference_and_sklearn():
    rng = np.random.default_rng(101)
    for it in range(300):
        n = int(rng.integers(1, 60))
        y = rng.random(n) < rng.random()
        s = rng.integers(0, int(rng.integers(1, 6)), n).astype(float)        # heavy ties
        w = rng.random(n) * 3
        if not y.any():
            assert np.isnan(E.average_precision(y, s))
            continue
        a = E.average_precision(y, s)
        assert abs(a - ref_ap(y, s)) < TOL
        assert abs(E.average_precision(y, s, w) - ref_ap(y, s, w)) < TOL
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert abs(a - sk_ap(y, s)) < TOL
            assert abs(E.average_precision(y, s, w) - sk_ap(y, s, sample_weight=w)) < TOL


def test_ap_adversarial_cases():
    rng = np.random.default_rng(102)
    y = rng.random(1000) < 0.03
    y[0] = True
    assert abs(E.average_precision(y, np.ones(1000)) - y.mean()) < 1e-14       # all tied -> prevalence
    assert E.average_precision([1, 1, 1], [1, 2, 3]) == 1.0
    assert E.average_precision([1], [0.3]) == 1.0
    assert np.isnan(E.average_precision([0], [0.3]))
    assert np.isnan(E.average_precision(np.zeros(0), np.zeros(0)))
    assert E.average_precision([1, 0], [0.0, -0.0]) == 0.5                     # signed zeros are one tie group
    with pytest.raises(ValueError):
        E.average_precision([1, 0], [np.nan, 1.0])                             # NaN scores refused, not dropped
    with pytest.raises(ValueError):
        E.average_precision([1, -1], [1.0, 2.0])                               # missing labels refused, not negatives


def test_weighted_ap_engine_equals_duplicated_cells_adversarial():
    rng = np.random.default_rng(103)
    for it in range(200):
        Q, n = int(rng.integers(1, 8)), int(rng.integers(1, 80))
        q, y = rng.integers(0, Q, n), rng.random(n) < 0.3
        s = rng.integers(0, 4, n).astype(float)
        eng = E.WeightedAP(y, s, q, Q)
        W = rng.integers(0, 4, (5, Q))                                         # includes zero-weight queries
        got = eng.ap(W)
        for b in range(5):
            idx = np.repeat(np.arange(n), W[b][q])
            want = ref_ap(y[idx], s[idx]) if y[idx].any() else np.nan
            assert (np.isnan(want) and np.isnan(got[b])) or abs(got[b] - want) < TOL
        w = rng.random(Q)                                                      # fractional weights
        if (w[q] * y).sum() > 0:
            assert abs(eng.ap(w)[0] - ref_ap(y, s, w[q])) < TOL
    assert abs(E.WeightedAP([1, 0, 0, 1], [1, 1, 1, 1], [0] * 4, 1).ap([3])[0] - 0.5) < 1e-15
    assert np.isnan(E.WeightedAP([0, 0], [1, 2], [0, 1], 2).ap([1, 1])[0])
    assert np.isnan(E.WeightedAP([1, 0], [1, 2], [0, 1], 2).ap([0, 1])[0])


def test_bootstrap_stratified_and_shared_across_arms():
    rng = np.random.default_rng(104)
    strata = np.repeat(np.arange(11), rng.integers(1, 30, 11))
    W = E.bootstrap_weights(strata, 4000, 5)
    for h in range(11):
        m = strata == h
        assert (W[:, m].sum(1) == m.sum()).all()                                # n_h draws within each stratum
        assert np.allclose(W[:, m].var(0), 1 - 1 / m.sum(), atol=0.15)        # multinomial variance
    Q = 20
    Wb = E.bootstrap_weights(np.zeros(Q, int), 30, 2)
    eng = E.ReplicateEngine(Wb, (10, 20), {"a": 1, "b": 2}, chunk=7)
    y = rng.random(400) < 0.1
    s, q = rng.normal(size=400) + y, np.repeat(np.arange(Q), 20)
    eng.add("a", 0, 10, y, s, q)
    eng.add("b", 1, 20, y, s, q)
    assert np.array_equal(eng.A["a"][:, 0, 0], eng.A["b"][:, 1, 1])           # same w^b for every arm and budget
    assert np.allclose(eng.A["a"][:, 0, 0], [E.average_precision(y, s, w[q]) for w in Wb], atol=TOL, rtol=0)
    assert abs(eng.est["a"][0, 0] - E.average_precision(y, s)) < TOL
    assert np.isnan(eng.A["b"][:, 0, :]).all()                                 # missing (r, n) stay NaN


# ============================================================================ 1b. inference
def test_welch_and_paired_by_hand():
    rng = np.random.default_rng(105)
    a, b = rng.normal(0.3, 0.01, 7), rng.normal(0.29, 0.02, 14)
    qn = rng.normal(0, 0.01, (300, 1))
    ba, bb = a[None] + qn, b[None] + 0.7 * qn + rng.normal(0, 0.003, (300, 1))
    r = E.welch(ba, a, bb, b, margin=0.001)
    VQ = np.var(ba.mean(1) - bb.mean(1), ddof=1)
    sa, sb = a.var(ddof=1) / 7, b.var(ddof=1) / 14
    se = np.sqrt(VQ + sa + sb)
    nu = se ** 4 / (sa ** 2 / 6 + sb ** 2 / 13)
    assert abs(r["se"] - se) < 1e-14 and abs(r["df"] - nu) < 1e-8
    assert abs(r["p_one_sided"] - (1 - stats.t.cdf((a.mean() - b.mean() - 0.001) / se, nu))) < 1e-12
    c = a + rng.normal(0.002, 0.001, 7)
    rp = E.paired(ba + 0.002, c, ba, a)
    tr = (c - a).var(ddof=1) / 7
    assert abs(rp["se"] ** 2 - (tr + np.var((ba + 0.002).mean(1) - ba.mean(1), ddof=1))) < 1e-15
    assert abs(rp["df"] - rp["se"] ** 4 / (tr ** 2 / 6)) < 1e-6


def test_holm_vs_step_down_reference():
    rng = np.random.default_rng(106)
    for it in range(500):
        p = rng.random(int(rng.integers(1, 7))) * 0.05
        if it % 3 == 0:
            p = np.round(p, 3)                                                  # ties
        assert np.array_equal(E.holm(p, 0.025)["reject"], ref_holm(p, 0.025))


def test_n_req_vs_brute_force_random_non_monotone_curves():
    rng = np.random.default_rng(107)
    for it in range(400):
        m = int(rng.integers(2, 10))
        n = np.sort(rng.choice(np.arange(10, 3000), m, replace=False)).astype(float)
        c = np.round(rng.normal(size=m).cumsum() * 0.05 + rng.normal(size=m) * 0.05, 2)
        f = E.pava(c)
        t = float(rng.uniform(f.min() - 0.05, f.max() + 0.05) if it % 2 else rng.choice(f))
        u = np.unique(f)
        if (np.diff(u) < 1e-12).any() and np.min(np.abs(u - t)) < 1e-12:
            continue      # adjacent PAVA blocks equal up to 1 ulp: "left endpoint" is float-ambiguous (not a defect)
        got, (wn, wc) = E.n_req(t, n, c), ref_n_req(t, n, c)
        assert got.censored == wc, (t, n, c, got)
        assert (np.isinf(got.n) and np.isinf(wn)) or (got.n == wn == 0) or abs(np.log(got.n / wn)) < 2e-3
    assert E.n_req(0.5, [10, 30, 100], [0.5, 0.5, 0.5]) == E.NReq(10.0, "")     # flat: left endpoint
    assert E.n_req(0.6, [10, 30, 100], [0.4, 0.6, 0.6]) == E.NReq(30.0, "")
    with pytest.raises(ValueError):
        E.n_req(0.5, [30, 10, 100], [0.1, 0.2, 0.3])                            # unsorted budgets refused


def test_precision_at_k_equals_expectation_over_tie_breaks():
    rng = np.random.default_rng(108)
    for it in range(60):
        n = int(rng.integers(1, 12))
        y, s, k = rng.random(n) < 0.4, rng.integers(0, 3, n).astype(float), int(rng.integers(1, 6))
        got = E.precision_at_k(y, s, np.zeros(n, int), k)[0]
        sims = [y[np.lexsort((rng.random(n), -s))[:min(k, n)]].mean() for _ in range(3000)]
        assert abs(got - np.mean(sims)) < 0.03


def test_bonferroni_crit():
    b = E.bonferroni_bands([0.1] * 8, [0.01] * 8, df=[20] * 8, level=0.95)
    assert abs(b["crit"][0] - stats.t.ppf(1 - 0.05 / 16, 20)) < 1e-12


def test_percentile_ci_keeps_censored_bounds():
    x = np.r_[np.full(950, 100.0), np.full(50, np.inf)]                         # 5% right-censored n_req replicates
    lo, hi = E.percentile_ci(x)
    assert lo == 100.0 and hi == np.inf


# ============================================================================ 2. leakage (dev data only)
def test_dev_bundle_every_row_array_blinded(dev_b):
    T = views.test_queries()
    for k in ("E", "SP", "Y", "DMF", "DMF_sd", "strain_ids", "strain_orf", "strain_canonical", "strain_is_supp",
              "strain_temperature"):
        assert getattr(dev_b, k).shape[0] == 2627, k
    assert not np.isin(dev_b.strain_orf, T).any() and set(dev_b.blinded_lockbox) == set(T)
    assert not dev_b.unblinded


@pytest.mark.parametrize("name", ["dev1", "dev2", "dev3", "resplit"])
def test_e_fixed_dev_complete_and_clean(dev_b, name):
    sp = views.get_split(name)
    T = views.test_queries()
    es = views.e_fixed(dev_b, sp)
    cur = dev_b.array_orfs_current[es.cols]
    assert not np.isin(dev_b.array_orfs[es.cols], T).any() and not np.isin(cur, views.current_of(dev_b, T)).any()
    assert not np.isin(cur, views.current_of(dev_b, sp.cand_full)).any()
    assert set(es.targets) == set(sp.targets) and not np.isin(sp.targets, sp.pool).any()
    rows = dev_b.gene_rows(es.targets[es.rows])
    assert np.array_equal(dev_b.Y[rows, es.cols], es.y) and set(np.unique(es.y)) <= {0, 1}
    ok = ~views.cols_matching(dev_b, np.union1d(sp.cand_full, T))              # independent recount
    r = dev_b.gene_rows(sp.targets)
    selfm = views.current_of(dev_b, sp.targets)[:, None] == dev_b.array_orfs_current[None, :]
    assert int(((dev_b.Y[r] >= 0) & ok[None, :] & ~selfm).sum()) == es.n


def test_final_e_fixed_columns_hold_no_dev_pool_array(dev_b):
    D = views.dev_pool()
    ok, _ = views._fixed_col_ok(dev_b, views.final_split())
    assert not np.isin(dev_b.array_orfs[ok], D).any()
    assert not np.isin(dev_b.array_orfs_current[ok], views.current_of(dev_b, D)).any()


def test_selector_view_masks_dev_split1(dev_b):
    sp = views.dev_split(1)
    T = views.test_queries()
    S = np.sort(np.random.default_rng(42).choice(sp.pool, 60, replace=False))
    v = views.SelectorView.build(dev_b, S, sp)
    masked = views.cols_matching(dev_b, np.union1d(sp.targets, T))
    assert np.array_equal(np.setdiff1d(np.arange(dev_b.n_arrays), v.visible_cols), np.flatnonzero(masked))
    assert np.isnan(v.E[:, masked]).all() and np.isnan(v.SP[:, masked]).all() and (v.Y[:, masked] == -1).all()
    r = dev_b.gene_rows(S)
    assert np.array_equal(v.E[:, ~masked], dev_b.E[r][:, ~masked], equal_nan=True)
    assert np.isnan(v.array_smf[masked]).all()
    assert set(v.candidates) == set(np.setdiff1d(sp.pool, S))
    assert not (set(v._visible_genes) & (set(sp.targets) | set(T)))
    public = {a for a in dir(v) if not a.startswith("_")}
    assert public == {"E", "SP", "Y", "array_orfs", "array_orfs_current", "array_smf", "build", "candidate_col",
                      "candidates", "clusters", "gene_features", "n_arrays", "n_screened", "screen_view",
                      "screened_col", "screened_orfs", "similarity", "visible_cols"}
    for g in (sp.targets[:3], T[:3]):
        for f in (v.gene_features, v.clusters):
            with pytest.raises(KeyError):
                f(g)
    sv = v.screen_view()
    assert np.isnan(sv.DMF).all() and np.isnan(sv.E[:, masked]).all() and np.isnan(sv.array_smf[masked]).all()
    assert not any(c.startswith("has_") for c in features.feature_columns())


def _poison_all(b, orfs, seed):
    """Garble every row matrix (incl. DMF/DMF_sd) AND the measurement pattern of `orfs`."""
    rng = np.random.default_rng(seed)
    r = b.gene_rows(orfs)
    M = {k: np.array(getattr(b, k), copy=True) for k in ("E", "SP", "Y", "DMF", "DMF_sd")}
    sh = M["E"][r].shape
    obs = rng.random(sh) < 0.9
    M["E"][r] = np.where(obs, rng.normal(-1, 2, sh), np.nan).astype(np.float32)
    M["SP"][r] = np.where(obs, rng.normal(0, 20, sh), np.nan).astype(np.float32)
    M["Y"][r] = np.where(obs, rng.integers(0, 2, sh), -1).astype(np.int8)
    M["DMF"][r] = np.where(obs, rng.random(sh), np.nan).astype(M["DMF"].dtype)
    M["DMF_sd"][r] = np.where(obs, rng.random(sh), np.nan).astype(M["DMF_sd"].dtype)
    return b.with_matrices(**M)


def test_prediction_poisoning_stack_lr_n60_dev1(dev_b):
    from hsl import upstream as U
    sp = views.dev_split(1)
    S = np.sort(np.random.default_rng(42).choice(sp.pool, 60, replace=False))

    def score(b):
        return U.run_method(U.construct_view(b, S, sp.targets), "stack_lr", model_seed=7, n_threads=3)[0]
    base = score(dev_b)
    assert np.array_equal(score(_poison_all(dev_b, sp.targets, 1)), base)               # D_val rows
    assert np.array_equal(score(_poison_all(dev_b, np.setdiff1d(sp.pool, S), 2)), base)  # unscreened pool rows
    assert not np.array_equal(score(_poison_all(dev_b, S[:5], 3)), base)                 # positive control
    r = H.run_predictor(H.RunSpec("dev1", tuple(S)), "stack_lr", seeds={"model_seed": 7}, bundle=dev_b,
                        save=False, ledger=False)
    es = r.esets["E_fixed"]
    assert np.array_equal(r.scores["E_fixed"], base[es.rows, es.cols].astype(np.float32))   # row/col alignment


@pytest.mark.xfail(strict=True, reason="ACCEPTED DEVIATION D4 (label-free external features only): k-means cluster universe includes D_val genes in dev "
                                       "(DESIGN §6.2: 'dev: also excluding D_val')")
def test_cluster_universe_excludes_dval_in_dev():
    sp = views.dev_split(1)
    assert not np.isin(sp.targets, features.cluster_universe()).any()


# ============================================================================ 3. harness / trajectories
class _RG:
    needs_view, static, name = True, False, "review_rg"

    def reset(self, c, rng):
        pass

    def select(self, view, k, rng):
        Y = view.Y.astype(float)
        Y[Y < 0] = 0
        col = view.candidate_col
        g = np.where(col >= 0, Y[:, np.maximum(col, 0)].sum(0), 0.0) if view.n_screened else np.zeros(len(col))
        return [str(view.candidates[i]) for i in H.randomised_greedy(g, k, rng)]


def _traj(sel, max_n, cps=(), preds=(), r=0, bs=None):
    return H.Trajectory(sel, "dev1", config.trajectory_seeds("dev", r + 1)[r], r, checkpoints=cps, max_n=max_n,
                        predictors=preds, batch_size=bs, save=False, ledger=False, write_log=False).run()


def test_checkpoints_do_not_change_randomised_greedy(monkeypatch, tmp_path):
    monkeypatch.setattr(H, "LEDGER", tmp_path / "L.jsonl")
    a = _traj(_RG(), 50)["selected"]
    b = _traj(_RG(), 50, cps=(13, 30, 47), preds=(("stack_lr", {}),))
    assert a == b["selected"] and set(b["checkpoints"]) == {"13", "30", "47"}


def test_random_selector_nested_across_batch_sizes_and_lengths():
    r10 = _traj(H.RandomSelector(), 100, bs=10)["selected"]
    r50 = _traj(H.RandomSelector(), 100, bs=50)["selected"]
    r300 = _traj(H.RandomSelector(), 300, bs=10)["selected"]
    assert r10 == r50 == r300[:100]


def test_randomised_greedy_uniform_top_k():
    g = np.arange(30, dtype=float)[::-1]
    cnt = np.bincount([H.randomised_greedy(g, 1, np.random.default_rng(s))[0] for s in range(20000)], minlength=30)
    assert cnt[20:].sum() == 0
    assert stats.chisquare(cnt[:20]).pvalue > 1e-3


def test_trajectory_does_not_overshoot_max_n():
    assert len(_traj(H.RandomSelector(), 655, bs=50)["selected"]) == 655


def test_pool_max_checkpoint_not_silently_dropped():
    log = _traj(H.RandomSelector(), 2627, cps=(2627,), bs=100)
    assert len(log["checkpoints"]) == 1
