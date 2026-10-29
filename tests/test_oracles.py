"""P3 oracle diagnostics (hsl.oracles) on toys: hub-first ordering, greedy-oracle round rule, dev-only refusal.
Not run on dev splits here (P3)."""
from __future__ import annotations

import hsl  # noqa: F401

import numpy as np
import pytest

from hsl import oracles

GENES = np.array([f"g{i:03d}" for i in range(60)])


def test_oracle_hub_first_orders_by_true_degree_with_seeded_ties():
    deg = {g: (i * 7) % 13 for i, g in enumerate(GENES)}          # many ties
    o = oracles.OracleHubFirst(counts_fn=lambda gs: np.array([deg[g] for g in gs]))
    o.reset(GENES[::-1], np.random.default_rng(1))
    picks = o.select(None, 10, None) + o.select(None, 10, None)
    d = [deg[g] for g in picks]
    assert d == sorted(d, reverse=True) and d[0] == 12 and len(set(picks)) == 20
    o2 = oracles.OracleHubFirst(counts_fn=lambda gs: np.array([deg[g] for g in gs]))
    o2.reset(GENES, np.random.default_rng(1))
    assert o2.select(None, 20, None) == picks                      # input order irrelevant; ties seeded
    full = {}
    for seed in (1, 2):
        o3 = oracles.OracleHubFirst(counts_fn=lambda gs: np.array([deg[g] for g in gs]))
        o3.reset(GENES, np.random.default_rng(seed))
        full[seed] = o3.select(None, 60, None)
    assert full[1] != full[2] and full[1][:20] == picks              # only the tie order depends on the seed
    assert [deg[g] for g in full[1]] == [deg[g] for g in full[2]]


def _coverage_score(weights, groups):
    def f(S):
        return float(sum(weights[c] for c in {groups[g] for g in S}))
    return f


def test_greedy_oracle_takes_top_gains_of_the_sampled_candidates():
    rng = np.random.default_rng(0)
    groups = {g: int(i % 15) for i, g in enumerate(GENES)}
    w = rng.random(15)
    o = oracles.GreedyOracle(n_candidates=200, score_fn=_coverage_score(w, groups))
    o.reset(GENES, np.random.default_rng(5))
    first = o.select(None, 10, np.random.default_rng(5))
    # all 60 candidates scored (pool < 200): picks = the 10 best singleton gains, in gain order (ties: index)
    gains = np.array([w[groups[g]] for g in GENES])
    expect = [str(GENES[i]) for i in np.lexsort((np.arange(60), -gains))[:10]]
    assert first == expect
    assert o.log[0]["base"] == 0.0 and o.log[0]["gains"] == sorted(o.log[0]["gains"], reverse=True)
    second = o.select(None, 10, np.random.default_rng(6))
    assert not set(first) & set(second)
    covered = {groups[g] for g in first}
    g2 = [w[groups[g]] if groups[g] not in covered else 0.0 for g in second]
    assert g2 == sorted(g2, reverse=True)                          # gains are marginal to the current S


def test_greedy_oracle_samples_only_n_candidates():
    seen = []

    def f(S):
        if S:
            seen.append(S[-1])
        return float(len(S))
    o = oracles.GreedyOracle(n_candidates=25, score_fn=f)
    o.reset(GENES, np.random.default_rng(0))
    picks = o.select(None, 10, np.random.default_rng(3))
    scored = set(seen)
    assert len(scored) == 25 and set(picks) <= scored


def test_oracles_refuse_final_split():
    with pytest.raises(PermissionError):
        oracles.true_query_degree("final", GENES[:2])
    o = oracles.GreedyOracle(split="final")
    with pytest.raises(PermissionError):
        o._auprc(["YAL001C"])
