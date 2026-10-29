"""H hub-first and H-static (DESIGN v3.1 §6.2): rank candidates by the predicted query-side SL degree q̂
(hsl.strategies.degree: unfitted prior for n < h_switch_n, ridge model from then on); randomised greedy (top-k
uniform) with ≤ cluster_cap_per_batch genes per k-means cluster (hsl.features clusters via view.clusters).
H-static uses the prior computed once at its first call (round 0) and never reads screened rows."""
from __future__ import annotations

import numpy as np

from hsl.strategies import base, degree


class Hub(base.Selector):
    name = "H"
    reads = frozenset({"rows", "smf"})

    def __init__(self, static: bool = False, **kw):
        super().__init__(static=bool(static), **kw)
        self.static = bool(static)
        if self.static:
            self.name = "H-static"
            self.reads = frozenset({"smf"})
        self._prior = None

    def reset(self, candidates, rng) -> None:
        super().reset(candidates, rng)
        self._prior = None

    # ---- scores
    def scores(self, view) -> tuple[np.ndarray, str]:
        """Score per view.candidates (higher = more hub-like) and the phase used."""
        cand = np.asarray(view.candidates, dtype=str)
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        if self.static:
            if self._prior is None:                     # round 0: prior over the round-0 candidates, frozen
                self._prior = dict(zip(cand, degree.prior_scores(view, cand, feats, pool)))
            return np.array([self._prior[g] for g in cand]), "prior_static"
        if view.n_screened < int(base.protocol_selector()["h_switch_n"]):
            return degree.prior_scores(view, cand, feats, pool), "prior"
        q, info = degree.fit_predict(view, cand, feats, pool)
        self._diag(**{f"qhat_{k}": v for k, v in info.items()})
        return q, "model"

    def gains(self, view) -> np.ndarray:
        return self.scores(view)[0]

    def _choose(self, view, k, rng) -> list:
        g, phase = self.scores(view)
        self._diag(phase=phase)
        return self._greedy(g, k, rng, clusters=view.clusters(view.candidates))
