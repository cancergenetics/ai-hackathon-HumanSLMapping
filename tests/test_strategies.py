"""Strategy selectors (DESIGN v3.1 §6.1–§6.2, §11): toy positive/negative controls, the F(r=1) ≡ H ablation,
per-channel poisoning on the P0 pseudo-test fixture for every selector (each check with a broken twin that must be
caught), determinism, uniform top-k sampling, batch cap, valid picks, and the strategies/oracles separation."""
from __future__ import annotations

import hsl  # noqa: F401

import ast
import subprocess
import sys

import numpy as np
import pytest
from scipy.stats import chisquare

from hsl import harness, paths, strategies, views
from hsl.strategies import kernels
from strategy_toys import ToyWorld, run_toy
from test_selector_view import LeakCols, LeakSMF, poison_cols, poison_rows, poison_smf

ALL = tuple(strategies.REGISTRY)
ADAPTIVE = tuple(n for n in ALL if not strategies.make(n).static)
STATIC = tuple(n for n in ALL if strategies.make(n).static)


def lexorder(g):
    g = np.asarray(g, dtype=np.float64)
    return np.lexsort((np.arange(len(g)), -g))


# ============================================================================ toy positive / negative controls
def test_planted_hub_H_beats_random():
    diff = []
    for s in range(6):
        w = ToyWorld(100 + s, "hub")
        h = run_toy(w, strategies.make("H"), s, 60)
        r = run_toy(w, harness.RandomSelector(), s, 60)
        diff.append(np.mean([w.a[g] for g in h]) - np.mean([w.a[g] for g in r]))
    diff = np.array(diff)
    assert diff.mean() > 3 * diff.std(ddof=1) / np.sqrt(len(diff)) and diff.mean() > 0.3, diff


def test_planted_cluster_C_covers_all_clusters_first():
    def n_cover(w, S):
        seen = set()
        for i, g in enumerate(S):
            seen.add(int(w.cl[w.gi[g]]))
            if len(seen) == 10:
                return i + 1
        return 10 ** 6
    nc, nr = [], []
    for s in range(5):
        w = ToyWorld(200 + s, "cluster")
        nc.append(n_cover(w, run_toy(w, strategies.make("C"), s, 60)))
        nr.append(n_cover(w, run_toy(w, harness.RandomSelector(), s, 60)))
    assert all(n == 10 for n in nc), nc                 # one gene per cluster in the first 10 picks
    assert np.mean(nr) > 15, nr


@pytest.mark.parametrize("n", [0, 40, 80])
def test_F_rank1_orders_exactly_like_H(n):
    w = ToyWorld(300 + n, "hub")
    S = list(np.random.default_rng(n).choice(w.queries, n, replace=False)) if n else []
    v = w.view(S)
    h, f = strategies.make("H"), strategies.make("F", r=1)
    for sel in (h, f):
        sel.reset(v.candidates, np.random.default_rng(0))
    assert np.array_equal(lexorder(h.gains(v)), lexorder(f.gains(v)))
    for det in (True, False):
        for seed in range(3):
            a = strategies.make("H", deterministic=det)
            b = strategies.make("F", r=1, deterministic=det)
            ra, rb = np.random.default_rng(seed), np.random.default_rng(seed)
            a.reset(v.candidates, ra)
            b.reset(v.candidates, rb)
            assert a.select(v, 10, ra) == b.select(v, 10, rb)


def test_null_toy_no_strategy_beats_random():
    """Negative control: with no structure, no selector's picks hold more true SLs than random beyond noise."""
    seeds = range(8)
    worlds = [ToyWorld(400 + s, "null") for s in seeds]
    rnd = np.array([w.true_degree(run_toy(w, harness.RandomSelector(), s, 40)).sum() for s, w in zip(seeds, worlds)])
    for name in strategies.ARMS:
        st = np.array([w.true_degree(run_toy(w, strategies.make(name), s, 40)).sum() for s, w in zip(seeds, worlds)])
        d = st - rnd
        se = np.sqrt(st.var(ddof=1) / len(st) + rnd.var(ddof=1) / len(rnd))
        assert abs(d.mean()) < 4 * se, (name, d.mean(), se)


# ============================================================================ determinism / sampling / cap / validity
@pytest.mark.parametrize("name", ALL)
def test_determinism_and_seed_dependence(name):
    w = ToyWorld(500, "hub")
    a = run_toy(w, strategies.make(name), 7, 50)
    b = run_toy(w, strategies.make(name), 7, 50)
    c = run_toy(w, strategies.make(name), 8, 50)
    assert a == b
    assert a != c


def test_uniform_top_k_sampling():
    w = ToyWorld(600, "hub")
    v = w.view([])
    g = strategies.make("H").gains(v)
    top = set(np.asarray(v.candidates)[lexorder(g)[:20]])
    picks = []
    for s in range(1200):
        sel = strategies.make("H")
        rng = np.random.default_rng(s)
        sel.reset(v.candidates, rng)
        picks += sel.select(v, 1, rng)
    assert set(picks) == top
    counts = np.array([picks.count(t) for t in sorted(top)])
    assert chisquare(counts).pvalue > 1e-3


@pytest.mark.parametrize("name", ["H", "F", "H-static", "F-static"])
@pytest.mark.parametrize("n", [0, 40])
def test_batch_cap(name, n):
    w = ToyWorld(700 + n, "hub", cluster_of_top_hubs=True)
    S = list(np.random.default_rng(n).choice(w.queries, n, replace=False)) if n else []
    v = w.view(S)
    cap = 2
    for det in (True, False):
        sel = strategies.make(name, deterministic=det)
        rng = np.random.default_rng(1)
        sel.reset(v.candidates, rng)
        picks = sel.select(v, 10, rng)
        cl = v.clusters(picks)
        assert np.bincount(cl).max() <= cap, (name, det, np.bincount(cl))
    # without the cap the planted hub cluster would dominate H's first batch (the test has power)
    if name == "H":
        g = sel.gains(v)
        assert (v.clusters(np.asarray(v.candidates)[lexorder(g)[:10]]) == 0).sum() > cap


# ============================================================================ real fixture: poisoning per channel
class LeakRowsAsScreened(views.SelectorView):
    """Broken twin for channels (a)/(d): extra (excluded or unscreened) genes' rows are exposed as if screened."""
    __slots__ = ()
    EXTRA = None

    @classmethod
    def build(cls, b, screened, split):
        v = super().build(b, screened, split)
        extra = np.sort(np.asarray(cls.EXTRA, dtype=str))
        r = b.gene_rows(extra)
        masked = np.setdiff1d(np.arange(b.n_arrays), v.visible_cols)
        S = np.concatenate([v.screened_orfs, extra])
        o = np.argsort(S)
        for k, fill in (("E", np.nan), ("SP", np.nan), ("Y", -1)):
            x = np.array(getattr(b, k)[r])
            x[:, masked] = fill
            object.__setattr__(v, k, np.vstack([getattr(v, k), x])[o])
        cmap = {g: j for j, g in enumerate(b.array_orfs)}
        v.screened_orfs = S[o]
        v.screened_col = np.array([cmap.get(g, -1) for g in S[o]], dtype=np.int64)
        keep = ~np.isin(v.candidates, extra)
        v.candidates, v.candidate_col = v.candidates[keep], v.candidate_col[keep]
        v._visible_genes = v._visible_genes | frozenset(extra)
        return v


class _Trap(Exception):
    pass


def _trap(name):
    def get(self):
        raise _Trap(f"selector read {name}")
    return property(get)


class TrapRows(views.SelectorView):          # any access to screened-row data raises
    __slots__ = ()
    E, SP, Y = _trap("E"), _trap("SP"), _trap("Y")


class TrapSMF(views.SelectorView):           # any access to array SMF raises
    __slots__ = ()
    array_smf = _trap("array_smf")


def build(view_cls, b, fx):
    if view_cls in (TrapRows, TrapSMF):
        v = views.SelectorView.build(b, fx["S"], fx["split"])
        v.__class__ = view_cls
        return v
    return view_cls.build(b, fx["S"], fx["split"])


def run_sel(name, view_cls, b, fx, seed=0, sel=None):
    sel = sel or strategies.make(name)
    rng = np.random.default_rng(seed)
    sel.reset(np.setdiff1d(fx["pool"], fx["S"]), rng)
    v = build(view_cls, b, fx)
    picks = sel.select(v, 10, rng)
    assert len(picks) == 10 and len(set(picks)) == 10
    assert np.isin(picks, fx["pool"]).all() and not np.isin(picks, fx["S"]).any()
    assert not np.isin(picks, fx["split"].excluded).any()
    return picks


@pytest.fixture(scope="module")
def poisoned(dev_b, pseudo):
    PTV = np.union1d(pseudo["PT"], pseudo["PV"])
    unscreened = np.setdiff1d(pseudo["pool"], pseudo["S"])
    return {"a": poison_rows(dev_b, PTV), "b": poison_cols(dev_b, pseudo["S"], PTV), "c": poison_smf(dev_b, PTV),
            "d": poison_rows(dev_b, unscreened), "S": poison_rows(dev_b, pseudo["S"]),
            "d_extra": unscreened[:100], "a_extra": pseudo["PT"]}


@pytest.fixture(scope="module")
def pseudo160(pseudo):
    """Same pseudo-test split, 160 screened genes: the degree model is in its informative regime (at n = 40 GCV often
    shrinks every coefficient to ~0, so an SMF leak cannot move H/F there — no power for the twin)."""
    from hsl import config
    rng = np.random.default_rng([config.seed("test_fixture_seed"), 160])
    return dict(pseudo, S=np.sort(rng.choice(pseudo["pool"], 160, replace=False)))


@pytest.fixture(scope="module")
def base_picks(dev_b, pseudo):
    return {n: run_sel(n, views.SelectorView, dev_b, pseudo) for n in ALL}


@pytest.mark.parametrize("name", ALL)
@pytest.mark.parametrize("ch", ["a", "b", "c", "d"])
def test_channel_poisoning_bit_identical(name, ch, dev_b, pseudo, poisoned, base_picks):
    assert run_sel(name, views.SelectorView, poisoned[ch], pseudo) == base_picks[name]


@pytest.mark.parametrize("name", ALL)
def test_broken_twins_caught(name, dev_b, pseudo, pseudo160, poisoned):
    """For each channel the selector consumes, a view that leaks that channel makes the poisoning visible; for each
    channel it does not consume, a trap view proves non-consumption (and the trap fires on a consuming selector)."""
    reads = strategies.make(name).reads
    if "rows" in reads:
        for ch in ("a", "d"):
            LeakRowsAsScreened.EXTRA = poisoned[f"{ch}_extra"]
            assert run_sel(name, LeakRowsAsScreened, poisoned[ch], pseudo) != \
                run_sel(name, LeakRowsAsScreened, dev_b, pseudo), (name, ch)
        assert run_sel(name, LeakCols, poisoned["b"], pseudo) != run_sel(name, LeakCols, dev_b, pseudo), name
    else:
        assert run_sel(name, TrapRows, dev_b, pseudo) == run_sel(name, views.SelectorView, dev_b, pseudo)
    if "smf" in reads:
        assert run_sel(name, views.SelectorView, poisoned["c"], pseudo160) == \
            run_sel(name, views.SelectorView, dev_b, pseudo160), name
        assert run_sel(name, LeakSMF, poisoned["c"], pseudo160) != run_sel(name, LeakSMF, dev_b, pseudo160), name
    else:
        assert run_sel(name, TrapSMF, dev_b, pseudo) == run_sel(name, views.SelectorView, dev_b, pseudo)


def test_traps_fire_on_consumers(dev_b, pseudo):
    """Broken twin of the non-consumption checks: a selector that does read the channel is caught by the trap."""
    for name, trap in (("H", TrapRows), ("C", TrapRows), ("M", TrapRows), ("F", TrapRows), ("H", TrapSMF),
                       ("R1", TrapSMF), ("H-static", TrapSMF)):
        with pytest.raises(_Trap):
            run_sel(name, trap, dev_b, pseudo)


@pytest.mark.parametrize("name", ADAPTIVE)
def test_adaptive_positive_control(name, dev_b, pseudo, poisoned, base_picks):
    """Poisoning the visible columns of the screened rows changes the selection of every adaptive selector."""
    assert run_sel(name, views.SelectorView, poisoned["S"], pseudo) != base_picks[name]


@pytest.mark.parametrize("name", STATIC)
def test_static_ignores_screened_rows(name, dev_b, pseudo, poisoned, base_picks):
    assert run_sel(name, views.SelectorView, poisoned["S"], pseudo) == base_picks[name]
    assert run_sel(name, TrapRows, dev_b, pseudo) == base_picks[name]


def test_static_broken_twin_caught(dev_b, pseudo, poisoned):
    """A selector labelled static that reads screened rows (adaptive H relabelled) fails the static test."""
    fake = strategies.make("H")
    fake.static = True
    assert run_sel("H", views.SelectorView, poisoned["S"], pseudo, sel=fake) != \
        run_sel("H", views.SelectorView, dev_b, pseudo)


def test_cached_kernel_parts_never_read_rows(dev_b, pseudo):
    """Everything cached across rounds (external kernel, kernel PCA, R1 clusters) is built without screened rows."""
    kernels.clear_cache()
    v = build(TrapRows, dev_b, pseudo)
    g = np.sort(v.candidates[:300])
    from hsl.strategies.base import FeatureCache
    fc = FeatureCache()
    K = kernels.external_kernel(v, g, fc.get(v, g))
    assert np.allclose(K, K.T) and np.linalg.eigvalsh(K).min() > -1e-8
    kernels.external_kpca(v, g, fc.get(v, g), 7)
    kernels.clear_cache()


def test_gi_inputs_match_upstream_mk(dev_b, pseudo):
    """The GI part equals upstream foldin_inputs('sym') on the upstream ScreenView built from the same view."""
    from ysl.models import completion as CM
    v = views.SelectorView.build(dev_b, pseudo["S"], pseudo["split"])
    Zs, Zc = kernels.gi_inputs(v, 1.0)
    sv = v.screen_view()
    Z, Zx, _, _ = CM.foldin_inputs(sv, "sym", "negeps", "row")
    c = np.sqrt(1.0 / (np.einsum("ij,ij->", Z, Z) / len(Z)))
    assert np.allclose(Zs, Z * c)
    assert np.allclose(Zc[v.candidate_col >= 0], Zx * c)


# ============================================================================ separation from the oracles
def test_strategies_never_import_oracles():
    pkg = paths.PROJECT_ROOT / "src" / "hsl" / "strategies"
    for f in pkg.glob("*.py"):
        tree = ast.parse(f.read_text())
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""] + [a.name for a in node.names]
            assert not any("oracle" in n for n in names), (f.name, names)
    code = "import hsl.strategies, sys; assert 'hsl.oracles' not in sys.modules"
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_registry_and_knobs():
    assert set(strategies.ARMS) <= set(strategies.REGISTRY)
    for n in ALL:
        s = strategies.make(n)
        assert s.needs_view and s.name == n and isinstance(s.params, dict)
    assert strategies.make("F").r == strategies.knobs()["F"]["r"]


# ============================================================================ numerical correctness of the greedy updates
def test_F_gains_and_sherman_morrison_match_bruteforce_logdet():
    w = ToyWorld(800, "hub")
    S = list(np.random.default_rng(1).choice(w.queries, 50, replace=False))
    v = w.view(S)
    f = strategies.make("F", r=4)
    f.reset(v.candidates, np.random.default_rng(0))
    st = f._fitted_state(v)
    M, Wc, uC = st["M"], st["W_c"], st["uC"]

    def obj(extra):
        A = M.copy()
        for j in extra:
            A = A + Wc[j][:, None, None] * np.outer(uC[j], uC[j])[None]
        return np.linalg.slogdet(A)[1].sum()
    g0 = st["gains"]()
    base_obj = obj([])
    js = [0, 5, 17, int(np.argmax(g0))]
    assert np.allclose([obj([j]) - base_obj for j in js], g0[js], rtol=1e-8, atol=1e-8)
    i = int(np.argmax(g0))
    g1 = st["update"](g0, i)
    js = [0, 5, 17, 33]
    assert np.allclose([obj([i, j]) - obj([i]) for j in js], g1[js], rtol=1e-7, atol=1e-8)


def test_M_batch_update_matches_direct_posterior():
    w = ToyWorld(801, "hub")
    S = list(np.random.default_rng(2).choice(w.queries, 30, replace=False))
    v = w.view(S)
    m = strategies.make("M")
    m.reset(w.queries, np.random.default_rng(0))
    gains, update = m._state(v)
    Zs, Zc = kernels.gi_inputs(v)
    lam = m.lam * strategies.knobs()["kernels"]["w_gi"]
    g = gains()
    i = int(np.nanargmax(np.where(np.isfinite(g), g, -1)))
    g1 = update(g, i)
    j2 = int(np.nanargmax(np.where(np.isfinite(g1) & (np.arange(len(g1)) != i), g1, -1)))
    g2 = update(g1, j2)
    Za = np.vstack([Zs, Zc[i], Zc[j2]])
    A = np.linalg.inv(Za @ Za.T + lam * np.eye(len(Za)))
    k = Zc @ Za.T
    direct = np.einsum("ij,ij->i", Zc, Zc) - np.einsum("ij,jk,ik->i", k, A, k)
    on = v.candidate_col >= 0
    assert np.allclose(g2[on], np.maximum(direct[on], 0), atol=1e-9)
    assert np.isneginf(g2[~on]).all()


def test_M_round0_is_R0_first_batch():
    w = ToyWorld(802, "hub")
    for s in range(3):
        assert run_toy(w, strategies.make("M"), s, 10) == run_toy(w, harness.RandomSelector(), s, 10)


def test_C_gains_match_bruteforce_facility_location():
    w = ToyWorld(803, "cluster")
    S = list(np.random.default_rng(3).choice(w.queries, 20, replace=False))
    v = w.view(S)
    c = strategies.make("C", offarray_calibration="none")       # raw objective; calibration checked below
    c.reset(v.candidates, np.random.default_rng(0))
    Kcc, Kcs = c.kernel_blocks(v)
    gains, update = c._state(v)
    cal = strategies.make("C")
    cal.reset(v.candidates, np.random.default_rng(0))
    on = v.candidate_col >= 0
    gc = cal.gains(v)
    assert np.array_equal(gc[on], gains()[on]) and not np.array_equal(gc[~on], gains()[~on])

    def F(extra):
        cols = [Kcs] + [Kcc[:, [j]] for j in extra]
        return np.maximum(np.hstack(cols).max(1), 0).sum()
    g = gains()
    assert np.allclose([F([j]) - F([]) for j in range(0, 200, 13)], g[0:200:13])
    i = int(np.argmax(g))
    g1 = update(g, i)
    assert np.allclose([F([i, j]) - F([i]) for j in range(0, 200, 13)], g1[0:200:13])


_PROV_SCRIPT = r"""
import hsl, json, sys, numpy as np
from hsl import config, guard, strategies, views
b = views.dev_bundle()
D = views.dev_pool(); rng = np.random.default_rng(config.seed("test_fixture_seed")); perm = D[rng.permutation(len(D))]
pool = np.sort(perm[300:]); excl = np.union1d(perm[:300], views.test_queries())
sp = views.Split("fixture", "dev", np.sort(perm[150:300]), pool, excl, pool)
out = {}
for n in (0, 40):
    S = np.sort(pool[:n])
    v = views.SelectorView.build(b, S, sp)
    with guard.record_reads() as ps:
        for name in strategies.REGISTRY:
            sel = strategies.make(name); r = np.random.default_rng(0); sel.reset(v.candidates, r); sel.select(v, 10, r)
    out[n] = sorted(ps)
print(json.dumps(out))
"""


def test_strategies_only_open_selector_stores():
    """Fresh process (no warm upstream caches): every file the selectors open is a store whitelisted for `selector`."""
    import json
    import os
    from test_provenance import _uncovered
    env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", VECLIB_MAXIMUM_THREADS="2")
    r = subprocess.run([sys.executable, "-c", _PROV_SCRIPT], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr[-2000:]
    reads = json.loads(r.stdout.strip().splitlines()[-1])
    assert reads["0"], "nothing recorded: the check would be vacuous"
    for n, ps in reads.items():
        assert not _uncovered(ps, "selector"), (n, _uncovered(ps, "selector"))


# ============================================================================ off-array routing is not a score term
def test_quantile_map_properties():
    from hsl.strategies import degree
    rng = np.random.default_rng(0)
    main_on, route_on, route_off = rng.normal(0, 1, 300), rng.normal(5, 3, 300), rng.normal(5, 3, 50)
    out = degree.quantile_map(main_on, route_on, route_off)
    assert out.min() >= main_on.min() and out.max() <= main_on.max()
    assert (np.diff(out[np.argsort(route_off)]) >= 0).all()          # order by the route score preserved


def test_offarray_share_of_H_picks_matches_base_rate():
    """No structure tied to off-array status: H's picks hold off-array genes at the pool's base rate (routing to the
    external-only route must not become an off-array score term, DESIGN §1)."""
    shares = []
    for s in range(8):
        w = ToyWorld(900 + s, "null")
        S = run_toy(w, strategies.make("H"), s, 60)
        shares.append(np.mean([g not in w.col for g in S]))
    base_rate = 40 / 240
    se = np.sqrt(base_rate * (1 - base_rate) / (60 * 8))
    assert abs(np.mean(shares) - base_rate) < 4 * se, (np.mean(shares), base_rate)


def test_F_calibration_touches_only_offarray_gains():
    w = ToyWorld(804, "hub")
    v = w.view(list(np.random.default_rng(4).choice(w.queries, 50, replace=False)))
    on = v.candidate_col >= 0
    raw, cal = strategies.make("F", offarray_calibration="none"), strategies.make("F")
    for s_ in (raw, cal):
        s_.reset(v.candidates, np.random.default_rng(0))
    gr, gc = raw.gains(v), cal.gains(v)
    assert np.allclose(gr[on], gc[on]) and not np.allclose(gr[~on], gc[~on])


@pytest.mark.parametrize("name", ["F", "C", "R1", "H", "H-static", "F-static", "C-static"])
def test_null_world_offarray_share_all(name):
    """Acceptance (PI, DESIGN §1/D5): exchangeable genes -> off-array share of picks within ±5 pp of the base rate."""
    sh = []
    for s in range(10):
        w = ToyWorld(7300 + s, "null")
        sh += [g not in w.col for g in run_toy(w, strategies.make(name), s, 80)]
    assert abs(np.mean(sh) - 40 / 240) < 0.05, (name, np.mean(sh))


def test_run_trajectory_writes_selector_diagnostics(tmp_path, monkeypatch):
    """The wrapper writes per-round diagnostics (picks, gains, q̂) next to the trajectory log (harness untouched)."""
    import json
    calls = {}

    class FakeTraj:
        split, write_log = "devX", True
        selector = strategies.make("H")

        def run(self):
            w = ToyWorld(805, "hub")
            run_toy(w, self.selector, 0, 40)
            calls["ran"] = True
            return {"traj_id": "devX_H_r0_test"}
    monkeypatch.setattr(harness, "TRAJ_DIR", tmp_path)
    log = strategies.run_trajectory(FakeTraj())
    d = json.loads((tmp_path / "devX" / "devX_H_r0_test.selector.json").read_text())
    assert calls["ran"] and len(d["rounds"]) == 4 and log["selector_diagnostics_file"]
    r3 = d["rounds"][3]
    assert r3["phase"] == "model" and len(r3["picks_qhat_logit"]) == 10 and len(r3["picks_gain0"]) == 10
    assert all(isinstance(x, bool) for x in r3["picks_offarray"])
