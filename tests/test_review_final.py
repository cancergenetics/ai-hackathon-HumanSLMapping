"""Independent review of the final pipeline (docs/reviews/final_pipeline_review.md). No T access: every test uses the
dry-run outputs on dev1 (results/final_dryrun/run, gitignored working files) or label-free names/columns.

1. Recompute one Welch contrast, one G and one n_80 of the dry run from the SAVED predictions with hsl.eval primitives
   plus a from-scratch numpy/sklearn cross-check (AP, Welch SE/ν).
2. Document the censoring defect of q2dev.ci (drops +inf right-censored replicates; used by n_req_joint and G).
3. Label-free structural checks of the final split: pool ∩ T = ∅, E_fixed column mask excludes every D column.
"""
from __future__ import annotations

import hsl  # noqa: F401

import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from hsl import config, eval as E, q2dev

ROOT = Path(__file__).resolve().parents[1]
DRY = ROOT / "results" / "final_dryrun"
RUN = DRY / "run"
need_run = pytest.mark.skipif(not (RUN / "labels.npz").exists(), reason="dry-run working files absent")


def _rows():
    rows = {}
    for f in sorted((RUN / "parts").glob("*.jsonl")):
        for ln in f.read_text().splitlines():
            if ln.strip():
                x = json.loads(ln)
                rows[x["uid"]] = x
    return list(rows.values())


def _sel(rows, cid, pred, n):
    return sorted([x for x in rows if x["cid"] == cid and x["pred"] == pred and int(x["n"]) == n],
                  key=lambda x: int(x["traj"]))


@pytest.fixture(scope="module")
def ctx():
    z = np.load(RUN / "labels.npz")
    cfg = yaml.safe_load((ROOT / "config" / "finalists.yaml").read_text())
    B = int(cfg["modes"]["dryrun"]["overrides"]["B_query_bootstrap"])
    base = config.seed(cfg["seeds"]["dryrun"]["bootstrap"])
    ss = lambda name: np.random.SeedSequence(base, spawn_key=(zlib.crc32(name.encode()),))  # noqa: E731
    W = E.bootstrap_weights(z["query_strata"], B, np.random.default_rng(ss("query_weights")))
    return {"y": z["y"].astype(bool), "rows": z["rows"], "nt": len(z["targets"]), "W": W, "ss": ss,
            "strata": z["query_strata"], "rws": _rows()}


def _vec(ctx, pred_file):
    s = np.load(ROOT / pred_file)["E_fixed"].astype(np.float32)
    e = E.WeightedAP(ctx["y"], s, ctx["rows"], ctx["nt"])
    return float(e.ap(np.ones(ctx["nt"]))[0]), e.ap(ctx["W"]), s


def _arm(ctx, cid, n, pred="pstar"):
    xs = _sel(ctx["rws"], cid, pred, n)
    out = [_vec(ctx, x["pred_file"]) for x in xs]
    return np.array([o[0] for o in out]), np.stack([o[1] for o in out], axis=1), out


@need_run
def test_strata_are_the_11_split_strata(ctx):
    assert sorted(np.unique(ctx["strata"]).tolist()) == list(range(11))


@need_run
def test_unit_ap_matches_sklearn_and_inline(ctx):
    from sklearn.metrics import average_precision_score
    x = _sel(ctx["rws"], "H", "pstar", 160)[0]
    est, _, s = _vec(ctx, x["pred_file"])
    assert abs(est - average_precision_score(ctx["y"], s.astype(np.float64))) < 1e-9
    assert abs(est - float(x["auprc_inline"])) < 1e-9


@need_run
def test_recompute_welch_family_H_a(ctx):
    es, As, _ = _arm(ctx, "H", 160)
    er, Ar, _ = _arm(ctx, "R0", 160)
    w = E.welch(As, es, Ar, er)
    # from-scratch Welch (DESIGN §6.4)
    VQ = np.var(As.mean(1) - Ar.mean(1), ddof=1)
    a, b = np.var(es, ddof=1) / len(es), np.var(er, ddof=1) / len(er)
    se = np.sqrt(VQ + a + b)
    nu = se ** 4 / (a ** 2 / (len(es) - 1) + b ** 2 / (len(er) - 1))
    from scipy import stats
    p = stats.t.sf((es.mean() - er.mean()) / se, nu)
    assert np.isclose(w["se"], se) and np.isclose(w["df"], nu) and np.isclose(w["p_one_sided"], p)
    conf = pd.read_csv(DRY / "confirmatory.csv").set_index(["finalist", "hypothesis"]).loc[("H", "a")]
    assert np.isclose(conf.delta, w["delta"], atol=1e-12)
    assert np.isclose(conf.p, w["p_one_sided"], rtol=1e-9)
    assert np.isclose(conf.se, w["se"], rtol=1e-9) and np.isclose(conf.df, w["df"], rtol=1e-9)


@need_run
def test_recompute_G_for_H(ctx):
    es, As, _ = _arm(ctx, "H", 160)
    e1, A1, _ = _arm(ctx, "R0", 160)
    e2, A2, _ = _arm(ctx, "R0", 660)
    G = (es.mean() - e1.mean()) / (e2.mean() - e1.mean())
    bg = {r["finalist"]: r for r in json.loads((DRY / "summary.json").read_text())["brief_G"]}["H"]
    assert np.isclose(G, bg["G"], rtol=1e-10)
    rng = np.random.default_rng(ctx["ss"]("G_H"))
    Ms = E.joint_resample(As[:, :, None], rng)[:, 0]
    Mr = E.joint_resample(np.stack([A1, A2], axis=2), rng)
    Gb = (Ms - Mr[:, 0]) / (Mr[:, 1] - Mr[:, 0])
    lo, hi = E.percentile_ci(Gb)
    assert np.isclose(lo, bg["G_lo"]) and np.isclose(hi, bg["G_hi"])
    gs = 1 - 0.05 * e2.mean() / (e2.mean() - e1.mean())
    assert np.isclose(gs, bg["G_star"])


@need_run
def test_recompute_n80_pstar(ctx):
    grid = [30, 160, 490, 660, 980]
    curve = np.array([_arm(ctx, "R0", n)[0].mean() for n in grid])
    # from-scratch: curve is checked monotone here so PAVA is the identity; log-n interpolation
    f = E.pava(curve)
    t = 0.8 * f[-1]
    k = int(np.argmax(f >= t))
    n80 = float(np.exp(np.log(grid[k - 1]) + (t - f[k - 1]) / (f[k] - f[k - 1]) * (np.log(grid[k]) - np.log(grid[k - 1]))))
    nx = pd.read_csv(DRY / "q1_nx.csv").set_index("predictor").loc["pstar"]
    assert np.isclose(n80, nx.n80, rtol=1e-9)
    assert np.isclose(E.n_req(t, grid, curve).n, n80)


def test_q2dev_ci_drops_right_censored_replicates():
    """REVIEW FINDING: q2dev.ci keeps only finite values, so +inf (right-censored) n_req replicates are dropped and the
    upper CI bound is finite even when most replicates are censored (dry run: H@490 97% right-censored, CI [624, 974]).
    E.percentile_ci keeps ±inf by design (P0 review M1). This test documents the defect (xfail until fixed)."""
    v = np.r_[np.full(90, np.inf), np.linspace(500, 900, 10)]
    lo_e, hi_e = E.percentile_ci(v)
    lo_q, hi_q = q2dev.ci(v)
    assert hi_e == np.inf
    if np.isfinite(hi_q):
        pytest.xfail(f"q2dev.ci drops +inf: CI ({lo_q:.0f}, {hi_q:.0f}) vs percentile_ci ({lo_e}, {hi_e})")


def test_final_split_pool_excludes_T_and_efixed_columns_exclude_D(dev_b):
    """Label-free: names/columns only (dev bundle has no T rows)."""
    from hsl import views
    sp = views.final_split()
    T, D = views.test_queries(), views.dev_pool()
    assert not np.isin(sp.pool, T).any() and np.isin(sp.pool, D).all()
    assert np.array_equal(np.sort(sp.excluded), np.sort(T))
    ok, ident = views._fixed_col_ok(dev_b, sp)
    assert not np.isin(dev_b.array_orfs[ok], D).any()
    # every on-array T gene is an E_fixed column unless identity-excluded
    on_T = np.isin(dev_b.array_orfs, T)
    assert ok[on_T].sum() >= on_T.sum() - len(ident)
