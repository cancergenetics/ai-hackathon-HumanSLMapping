"""P3/P4 helpers (hsl.q2dev) on toys with known answers: config-id parsing, feasible-pool restriction (picks stay in
the allowed set; unrestricted selectors unchanged), HAP1-minimal knob map, joint-bootstrap G and the §6.4 contrast
suite on synthetic replicate arrays, n_req on a known curve. No dev data needed."""
from __future__ import annotations

import hsl  # noqa: F401

import numpy as np
import pytest

from hsl import q2dev, strategies
from strategy_toys import ToyWorld, run_toy


def test_parse_config_ids():
    s = q2dev.parse("F_r16__b50__det")
    assert (s.arm, s.knobs, s.batch_size, s.deterministic) == ("F", {"r": 16}, 50, True)
    s = q2dev.parse("M_lam0.03__ws220__feas__min")
    assert (s.arm, s.knobs, s.warm, s.feasible, s.minimal) == ("M", {"lam": 0.03}, 220, True, True)
    assert q2dev.parse("oracle_greedy").arm == "oracle_greedy"
    assert q2dev.parse("F-static_r8").arm == "F-static"
    assert [q2dev.family(a) for a in ("H", "H-static", "F", "F-static", "C", "M", "R1")] == \
        ["H", "H", "F", "F", "C/M", "C/M", "R1"]
    with pytest.raises(ValueError):
        q2dev.parse("H_x3")
    with pytest.raises(ValueError):
        q2dev.parse("H__b")


@pytest.mark.parametrize("name", ["H", "H-static", "F", "C", "R1"])
def test_restrict_pool_keeps_picks_inside_allowed(name):
    w = ToyWorld(3, "hub")
    allowed = np.sort(w.queries[::2])
    sel = q2dev.RestrictPool(strategies.make(name), allowed)
    S = run_toy(w, sel, 7, 60)
    assert set(S) <= set(allowed) and len(set(S)) == 60


def test_restrict_pool_with_everything_allowed_is_identity():
    w = ToyWorld(4, "hub")
    a = run_toy(w, strategies.make("H"), 11, 50)
    b = run_toy(w, q2dev.RestrictPool(strategies.make("H"), w.queries), 11, 50)
    assert a == b


def test_minimal_knobs_drop_apms_only():
    from hsl.strategies import base
    k = base.knobs()
    m = q2dev.minimal_knobs(k)
    assert m["degree_model"]["features"]["apms_degree"] is None
    assert set(m["kernels"]["external"]) == set(k["kernels"]["external"]) - {"apms", "apmse"}
    assert k["degree_model"]["features"]["apms_degree"] is not None          # original untouched


def _synthetic(rng, means, R, B=400, sd_traj=0.0005, sd_q=0.001):
    """Replicate arrays [B x R x J] with a shared query effect per replicate (common to arms) + trajectory noise."""
    J = len(means)
    q = rng.normal(0, sd_q, size=(B, 1, 1))
    tr = rng.normal(0, sd_traj, size=(1, R, J))
    A = np.asarray(means)[None, None, :] + q + tr
    est = np.asarray(means)[None, :] + tr[0]
    return A, est, q


def test_gap_closure_joint_known_answer():
    rng = np.random.default_rng(0)
    # random 0.040 @160, 0.060 @660; strategy 0.050 @160 -> G = 0.5
    Ar, er, q = _synthetic(rng, [0.040, 0.060], 10)
    tr_s = rng.normal(0, 0.0005, size=(1, 5, 1))
    As = 0.050 + q[:, :, :1] + tr_s
    es = 0.050 + tr_s[0]                                                    # unit-weight estimates (no query shift)
    g = q2dev.gap_closure_joint([(As, Ar, es, er)], np.random.default_rng(1))
    assert abs(g["G"] - 0.5) < 0.06 and g["defined"]
    assert g["G_ci"][0] < g["G"] < g["G_ci"][1]
    # pooling two identical splits leaves G unchanged
    g2 = q2dev.gap_closure_joint([(As, Ar, es, er)] * 2, np.random.default_rng(1))
    assert abs(g2["G"] - g["G"]) < 1e-12


def test_contrast_suite_signs_and_matches():
    rng = np.random.default_rng(2)
    grid = [160, 490, 660]
    Ar, er, q = _synthetic(rng, [0.040, 0.052, 0.060], 10)
    # strong strategy: 0.062 at 160 (beats r490, r660, and 0.95*r660)
    As = np.full((len(q), 4, 3), np.nan)
    As[:, :, 0] = 0.062 + q[:, :, 0] + rng.normal(0, 0.0003, size=(1, 4))
    es = np.full((4, 3), np.nan)
    es[:, 0] = 0.062 + rng.normal(0, 0.0003, size=4)
    cs = q2dev.contrast_suite(As, es, Ar, er, grid, "toy")
    assert set(cs["delta"]) == {160}                                        # strategy has no 490/660 checkpoints
    assert cs["delta"][160]["delta"] > 0.02 and cs["delta"][160]["p_one_sided"] < 1e-3
    assert cs["c1"]["delta"] > 0 and cs["c1"]["p_one_sided"] < 0.025
    assert cs["matches"]["p_one_sided"] < 0.025
    assert abs(cs["matches"]["delta"] - (np.mean(es[:, 0]) - 0.95 * np.mean(er[:, 2]))) < 1e-12
    assert cs["G"]["G"] > 1 and abs(cs["G_star"] - (1 - 0.05 * er[:, 2].mean() / (er[:, 2].mean() - er[:, 0].mean()))) < 1e-12
    # null strategy = random at 160 -> Δ CI covers 0
    es0 = es.copy()
    es0[:, 0] = 0.040 + rng.normal(0, 0.0005, size=4)
    As0 = As.copy()
    As0[:, :, 0] = es0[None, :, 0] + q[:, :, 0]
    cs0 = q2dev.contrast_suite(As0, es0, Ar, er, grid, "toy0")
    lo, hi = cs0["delta"][160]["ci95"]
    assert lo < 0 < hi


def test_n_req_joint_known_curve():
    rng = np.random.default_rng(3)
    grid = [100, 200, 400, 800]
    Ar, er, q = _synthetic(rng, [0.02, 0.03, 0.04, 0.05], 10, sd_traj=1e-5, sd_q=1e-5)
    As = np.full((len(q), 3, 1), 0.035) + q[:, :, :1]
    es = np.full((3, 1), 0.035)
    out = q2dev.n_req_joint(As, es, Ar, er, grid, 0, np.random.default_rng(0))
    assert abs(out["n_req"] - np.exp((np.log(200) + np.log(400)) / 2)) < 2.0   # halfway in log n
    assert out["ci"][0] <= out["n_req"] <= out["ci"][1]


def test_minimal_patch_is_scoped():
    from hsl.strategies import base
    full = base.knobs()
    sha0 = base.knobs_sha()
    q2dev._patch_minimal()
    assert base.knobs()["degree_model"]["features"]["apms_degree"] is None
    assert base.knobs_sha() != sha0                                         # trajectory ids differ (knob hash)
    q2dev._unpatch_minimal()
    assert base.knobs() is full and base.knobs_sha() == sha0
