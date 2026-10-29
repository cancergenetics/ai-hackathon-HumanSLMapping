"""M max posterior variance (DESIGN v3.1 §6.2; review r2 item 15): pick argmax σ²(g) = k(g, g) − k_gS(K_SS + λI)⁻¹k_Sg
under the GI-only kernel (the mk GI 'sym' kernel; hsl.strategies.kernels.gi_inputs), λ = lam·w_gi. Within a batch the
picked genes are added as pseudo-observations (σ²(x | S ∪ g) = σ²(x | S) − cov(x, g | S)²/(σ²(g | S) + λ)).
Round 0 (S = ∅) = the first batch of the R0 permutation of the same trajectory: reset() draws the R0 permutation
exactly as harness.RandomSelector.reset does from the same selection stream. Off-array candidates (no GI column)
are never picked by the criterion (gain −inf; count reported in diagnostics). If the GI kernel is degenerate (no
screened on-array gene), M continues along the R0 permutation (on-array genes only)."""
from __future__ import annotations

import numpy as np

from hsl import harness
from hsl.strategies import base, kernels


class MaxVar(base.Selector):
    name = "M"
    reads = frozenset({"rows"})

    def __init__(self, lam: float | None = None, **kw):
        self.lam = float(base.knobs()["M"]["lam"] if lam is None else lam)
        super().__init__(lam=self.lam, **kw)
        self._r0 = None

    def reset(self, candidates, rng) -> None:
        super().reset(candidates, rng)
        self._r0 = harness.RandomSelector()
        self._r0.reset(candidates, rng)                 # same draw as R0 on this selection stream

    def r0_order(self) -> list:
        return list(self._r0._perm)

    def _state(self, view):
        w_gi = float(base.knobs()["kernels"]["w_gi"])
        lam = self.lam * w_gi
        Zs, Zc = kernels.gi_inputs(view, w_gi)
        on = np.asarray(view.candidate_col) >= 0
        if Zs.shape[1] == 0 or not np.any(Zs):
            return None
        Kss = Zs @ Zs.T
        A = np.linalg.inv(Kss + lam * np.eye(len(Kss)))
        Kcs = Zc @ Zs.T
        KA = Kcs @ A
        var = np.einsum("ij,ij->i", Zc, Zc) - np.einsum("ij,ij->i", KA, Kcs)
        st = {"var": np.maximum(var, 0.0), "C": [], "d": []}

        def gains():
            return np.where(on, st["var"], -np.inf)

        def update(g, i):
            c = Zc @ Zc[i] - KA @ Kcs[i]                     # cov(·, i | S)
            for cj, dj in zip(st["C"], st["d"]):             # condition on earlier picks of this batch
                c = c - cj * cj[i] / dj
            d = st["var"][i] + lam
            st["var"] = np.maximum(st["var"] - c ** 2 / d, 0.0)
            st["C"].append(c)
            st["d"].append(d)
            return gains()
        return gains, update

    def gains(self, view) -> np.ndarray:
        st = self._state(view)
        return st[0]() if st is not None else np.full(len(view.candidates), np.nan)

    def _choose(self, view, k, rng) -> list:
        cand = np.asarray(view.candidates, dtype=str)
        on = np.asarray(view.candidate_col) >= 0
        self._diag(n_offarray_candidates=int((~on).sum()))
        pos = {g: i for i, g in enumerate(cand)}
        if view.n_screened == 0:
            self._diag(phase="r0_first_batch")
            return [pos[g] for g in self._r0.select(None, k, rng)]
        st = self._state(view)
        if st is None:
            self._diag(phase="r0_continuation")
            return [pos[g] for g in self.r0_order() if g in pos and on[pos[g]]][:k]
        gains, update = st
        self._diag(phase="maxvar")
        return self._greedy(gains(), k, rng, inner=update)
