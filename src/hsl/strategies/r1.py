"""R1 Billmann-like (status quo, static; DESIGN v3.1 §6.2): round-robin over SMF × abundance tertile bins, spread
across k-means clusters; trajectory r = seeded permutation of the within-bin order. Computed once, at the first call
(round 0), from external features and array SMF only; never reads screened rows.

  SMF bin: tertile of the candidate's array SMF; edges = tertiles of all non-NaN array SMF values of the view (the
           library). Off-array candidates have no SMF: each is assigned an SMF tertile uniformly at random from the
           trajectory's selection stream (review R1-1: a ridge-imputed SMF regresses to the mean and put 347/410
           off-array genes in one tertile, halving their share). No missingness stratum.
  abundance bin: tertile among the round-0 candidates (missing abundance -> median).
  clusters: k-means (k = protocol kmeans_k, cluster_seed) on the external-kernel embedding of the VISIBLE universe
           (candidates ∪ screened ∪ visible arrays: non-test genes; dev: also excluding D_val) —
           hsl.strategies.kernels.visible_clusters.
  order:   within each bin a seeded permutation (selection stream), then stably re-ordered by the occurrence index of
           each gene's cluster within the bin (cluster round-robin); bins visited round-robin in fixed order.
"""
from __future__ import annotations

import numpy as np

from hsl.strategies import base, kernels


class Billmann(base.Selector):
    name = "R1"
    static = True
    reads = frozenset({"smf"})

    def reset(self, candidates, rng) -> None:
        super().reset(candidates, rng)
        self._order, self._pos, self._rng = None, 0, rng

    def bins(self, view, rng=None) -> tuple[np.ndarray, np.ndarray]:
        """(bin id per candidate, cluster per candidate). Off-array SMF tertiles are drawn from `rng` (default: the
        selection stream given to reset())."""
        sp = base.knobs()["spec"]
        ns, na = int(sp["n_smf_bins"]), int(sp["n_abund_bins"])
        cand = np.asarray(view.candidates, dtype=str)
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        smf = base.candidate_smf(view)
        off = ~np.isfinite(smf)
        rng = self._rng if rng is None else rng
        ref = np.asarray(view.array_smf, dtype=np.float64)
        e_s = np.nanquantile(ref, np.arange(1, ns) / ns)
        ab = base.abundance(feats)[np.searchsorted(pool, cand)]
        ab = np.where(np.isfinite(ab), ab, np.nanmedian(ab))
        e_a = np.quantile(ab, np.arange(1, na) / na)
        b_s = np.searchsorted(e_s, np.where(off, 0.0, smf), side="right")
        if off.any():
            b_s[off] = rng.integers(0, ns, int(off.sum()))
        b_a = np.searchsorted(e_a, ab, side="right")
        U, lab = kernels.visible_clusters(view, lambda v, g: self._feats.get(v, g))
        cl = lab[np.searchsorted(U, cand)]
        return b_s * na + b_a, cl

    def order(self, view, rng) -> list:
        cand = np.asarray(view.candidates, dtype=str)
        b, cl = self.bins(view, rng)
        queues = []
        for h in range(int(b.max()) + 1 if len(b) else 0):
            idx = np.flatnonzero(b == h)
            if not len(idx):
                continue
            idx = idx[rng.permutation(len(idx))]
            seen: dict = {}
            occ = np.empty(len(idx), dtype=np.int64)
            for j, i in enumerate(idx):
                occ[j] = seen.get(cl[i], 0)
                seen[cl[i]] = occ[j] + 1
            queues.append(list(idx[np.argsort(occ, kind="stable")]))
        out = []
        while any(queues):
            for q in queues:
                if q:
                    out.append(str(cand[q.pop(0)]))
        return out

    def _choose(self, view, k, rng) -> list:
        if self._order is None:
            self._order = self.order(view, rng)
            self._diag(phase="order_built")
        pos = {g: i for i, g in enumerate(np.asarray(view.candidates, dtype=str))}
        out = []
        while len(out) < k and self._pos < len(self._order):
            g = self._order[self._pos]
            self._pos += 1
            if g in pos:
                out.append(pos[g])
        return out

    def gains(self, view) -> np.ndarray:
        raise NotImplementedError("R1 is an ordering, not a gain-based selector")
