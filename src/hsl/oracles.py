"""P3 ceilings (DESIGN v3.1 §9 P3): EVALUATION-ONLY oracle selectors. They read pool / D_val labels, so they are
NOT strategies: they live outside hsl.strategies, which must never import this module (tests/test_strategies.py
checks it), and they refuse the final split.

  OracleHubFirst(split)   rank the candidate pool by TRUE query-row strong-SL count (dev rows; columns of the split's
                          excluded genes not counted), ties broken by the trajectory's selection stream; static.
  GreedyOracle(split)     each round: `n_candidates` (200) random remaining candidates; each scored by the realised
                          E_fixed_dev AUPRC of the closed-form fold-in `F2_mk_foldin` fitted on S ∪ {c}; the top
                          batch (10) by gain over AUPRC(S) is taken in gain order (ties: candidate index).
  degree_diagnostic(view, scores)  Spearman of a selector's q̂ with the true pool degree (dev diagnostic of H).
Both oracles follow the harness selector protocol with needs_view = False (they track S themselves). Scoring
functions are injectable for toy tests.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import numpy as np

from hsl import views


def _dev_split(split):
    sp = views.get_split(split) if isinstance(split, str) else split
    if sp.kind != "dev":
        raise PermissionError("oracles are dev-only evaluation diagnostics (never the final split)")
    return sp


def true_query_degree(split, genes, bundle=None) -> np.ndarray:
    """Strong-SL count of each gene's query row over measured cells, excluding the split's masked columns."""
    sp = _dev_split(split)
    b = views.dev_bundle() if bundle is None else bundle
    keep = ~views.cols_matching(b, sp.excluded)
    Y = b.Y[b.gene_rows(np.asarray(genes, dtype=str))]
    return ((Y == 1) & keep[None, :]).sum(1).astype(np.int64)


class OracleHubFirst:
    needs_view = False
    static = True
    name = "oracle_hub_first"

    def __init__(self, split=None, counts_fn=None, bundle=None):
        self.split, self.bundle = split, bundle
        self.counts_fn = counts_fn or (lambda genes: true_query_degree(self.split, genes, self.bundle))
        self.params = {"oracle": "hub_first"}

    def reset(self, candidates, rng) -> None:
        c = np.sort(np.asarray(candidates, dtype=str))
        deg = np.asarray(self.counts_fn(c), dtype=np.float64)
        tie = rng.permutation(len(c))
        self._order = list(c[np.lexsort((tie, -deg))])
        self._pos = 0

    def select(self, view, k: int, rng) -> list:
        out = self._order[self._pos:self._pos + int(k)]
        self._pos += len(out)
        return out


class GreedyOracle:
    needs_view = False
    static = False
    name = "oracle_greedy"

    def __init__(self, split=None, n_candidates: int = 200, method: str = "F2_mk_foldin", score_fn=None,
                 bundle=None, n_threads: int = 2):
        self.split, self.bundle, self.method, self.n_threads = split, bundle, method, int(n_threads)
        self.n_candidates = int(n_candidates)
        self.score_fn = score_fn or self._auprc
        self.params = {"oracle": "greedy", "n_candidates": self.n_candidates, "method": method}
        self.log: list = []

    def _auprc(self, S) -> float:
        """Realised E_fixed_dev AUPRC of `method` fitted on screened set S (0 for S = ∅)."""
        from hsl import eval as E
        from hsl import harness as H
        from hsl import upstream as U
        if not len(S):
            return 0.0
        sp = _dev_split(self.split)
        b = views.dev_bundle() if self.bundle is None else self.bundle
        v = U.construct_view(b, S, sp.targets)
        Smat, _ = U.run_method(v, self.method, n_threads=self.n_threads, n_folds=min(10, max(len(S), 2)))
        es = H.get_eset(b, sp, "E_fixed")
        return float(E.average_precision(es.y, es.scores(Smat)))

    def reset(self, candidates, rng) -> None:
        self._remaining = np.sort(np.asarray(candidates, dtype=str))
        self._S: list = []
        self._base = None
        self.log = []

    def set_s0(self, s0) -> None:
        """Warm start: genes already screened before the first round (the harness passes only candidates)."""
        self._S = [str(g) for g in s0]

    def select(self, view, k: int, rng) -> list:
        rem = np.setdiff1d(self._remaining, self._S)
        m = min(self.n_candidates, len(rem))
        pick = np.sort(rng.choice(len(rem), size=m, replace=False))
        cands = rem[pick]
        base = self.score_fn(list(self._S)) if self._base is None else self._base
        gains = np.array([self.score_fn(list(self._S) + [str(c)]) - base for c in cands])
        order = np.lexsort((np.arange(m), -gains))[: int(k)]
        out = [str(cands[i]) for i in order]
        self._S.extend(out)
        self._base = None                                   # AUPRC(S) of the new S is computed next round
        self.log.append({"n_before": len(self._S) - len(out), "base": base, "gains": gains[order].tolist()})
        return out


def degree_diagnostic(view, scores, split, bundle=None) -> float:
    """Spearman(selector q̂ over view.candidates, true query degree) — DEV diagnostic (reads pool labels)."""
    from scipy.stats import spearmanr
    return float(spearmanr(np.asarray(scores), true_query_degree(split, view.candidates, bundle))[0])
