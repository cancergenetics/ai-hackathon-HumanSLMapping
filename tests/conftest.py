"""Shared fixtures. Fast suite: `.venv/bin/python -m pytest -m "not slow"`; full: `.venv/bin/python -m pytest`.
No fixture reads test-query rows: the dev bundle has none, and leakage fixtures use a PSEUDO-test set drawn from D."""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import numpy as np
import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: end-to-end / long-running tests (excluded by -m 'not slow')")
    config.addinivalue_line("markers", "final_only: needs HSL_FINAL=1 (unblinded data); skipped otherwise")


@pytest.fixture(scope="session")
def dev_b():
    from hsl import views
    return views.dev_bundle()


@pytest.fixture(scope="session")
def pseudo():
    """Pseudo-test split inside D (seeded by test_fixture_seed): PT = pseudo-test genes (rows present in the dev
    bundle, masked like T), PV = pseudo-validation targets, pool = the rest of D. Excluded = PT ∪ PV ∪ T."""
    from hsl import config, views
    D = views.dev_pool()
    rng = np.random.default_rng(config.seed("test_fixture_seed"))
    perm = D[rng.permutation(len(D))]
    PT, PV = np.sort(perm[:150]), np.sort(perm[150:300])
    pool = np.sort(perm[300:])
    excl = np.union1d(np.union1d(PT, PV), views.test_queries())
    sp = views.Split("fixture", "dev", PV, pool, excl, pool)
    S = np.sort(pool[rng.permutation(len(pool))[:40]])
    return {"PT": PT, "PV": PV, "pool": pool, "split": sp, "S": S}
