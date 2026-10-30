"""Cross-checks of the Q3 code against stored project results (dev data; needs data/q3_cache, built by
`python -m analysis.q3.cache`). No model fits.

    .venv/bin/python -m pytest analysis/q3/tests -q
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from analysis.q3 import common as C

pytestmark = pytest.mark.skipif(not (C.ROOT / "data" / "q3_cache" / "dev_cache.npz").exists(),
                                reason="dev label cache not built")


def test_stratum_nx_reproduces_p1():
    from analysis.q3 import hap1
    nx = C.nx_table("dev")
    x = hap1.nx_from_rows(C.rows("dev"), "auprc", nx)
    assert np.allclose(x.n80, nx.n80, rtol=1e-9) and np.allclose(x.n90, nx.n90, rtol=1e-9)


def test_slope_point_reproduces_p1():
    from analysis.q3 import hap1
    s = json.loads((C.ROOT / "results" / "p1" / "universe" / "slope.json").read_text())
    b0, _ = hap1.slope_boot(C.universe_draws(), "n80", 10, np.random.default_rng(0))
    assert b0 == pytest.approx(s["n80"]["b"], abs=1e-9)


def test_frontier_cells_match_stored_counts():
    """Recomputed cell-level counts exclude T columns; P1/P4 stored counts include them -> ours <= stored, close."""
    from analysis.q3 import frontier
    net = frontier.build_network()
    rows = C.rows("dev")
    sels = C.selections("dev", arms=["R0", "H"])
    for (arm, sp, r), s in list(sels.items())[:6]:
        x = rows[(rows.arm == arm) & (rows.split == sp) & (rows.traj == r) & (rows.n == 160)]
        if x.empty:
            continue
        c = net.counts(s["selected"][:160])["cells_detected"]
        assert c <= x.sl_detected.iloc[0] and c >= 0.75 * x.sl_detected.iloc[0]


def test_eval_labels_dev_matches_dryrun_labels():
    """The scorer's dry-run label file (dev1 stands in for T) equals our dev1 cache cell for cell."""
    p = C.ROOT / "results" / "final_dryrun" / "run" / "labels.npz"
    if not p.exists():
        pytest.skip("no dry-run labels")
    a, b = C.eval_labels("dev", "dev1"), C.eval_labels("dryrun", "dev1")
    for k in ("rows", "cols", "y", "OA", "AO"):
        assert np.array_equal(a[k], b[k]), k


def test_random_profile_is_null():
    p = C.ROOT / "results" / "q3_dev" / "choice_profile.csv"
    if not p.exists():
        pytest.skip("profile not run")
    d = pd.read_csv(p)
    r = d[(d.arm == "R0") & (d.feature != "offarray_share")]
    assert (r.smd.abs() < 0.1).all()
