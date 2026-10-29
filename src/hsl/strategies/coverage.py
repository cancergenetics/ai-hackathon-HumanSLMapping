"""C coverage (DESIGN v3.1 §6.2): facility location over the predictor's `mk` kernel K = K_gi + Σ w_k K_ext,k
(hsl.strategies.kernels: GI 'sym' kernel over the screened on-array coordinates + tier-1 external kernels), target
pool = unscreened candidates. Objective F(S) = Σ_{x∈cand} max(0, max_{s∈S} K(x, s)); gain of c =
Σ_x max(0, K(x, c) − cov_x). Exact greedy gains are recomputed after every pick (vectorised, O(|cand|²) per pick —
cheaper than lazy-greedy bookkeeping at this size; identical result), with the randomised top-k rule.
Round 0 (S = ∅): the GI part is empty, so C runs on the external kernel alone (the reported "static part").
C-static (not an arm; diagnostic) uses the external kernel only at every round and never reads screened rows."""
from __future__ import annotations

import numpy as np

from hsl.strategies import base, kernels


class Coverage(base.Selector):
    name = "C"
    reads = frozenset({"rows"})

    def __init__(self, static: bool = False, normalize: bool | None = None, **kw):
        nm = bool(base.knobs()["C"]["normalize"] if normalize is None else normalize)
        super().__init__(static=bool(static), normalize=nm, **kw)
        self.static, self.normalize = bool(static), nm
        if self.static:
            self.name, self.reads = "C-static", frozenset()

    def kernel_blocks(self, view) -> tuple[np.ndarray, np.ndarray]:
        """(K_cc [n_cand x n_cand], K_cs [n_cand x n_screened])."""
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        cand, S = np.asarray(view.candidates, dtype=str), np.asarray(view.screened_orfs, dtype=str)
        K = kernels.external_kernel(view, pool, feats)
        pc, ps = np.searchsorted(pool, cand), np.searchsorted(pool, S)
        Kcc = np.array(K[np.ix_(pc, pc)])
        Kcs = np.array(K[np.ix_(pc, ps)])
        dS = np.array(K[ps, ps])
        if not self.static:
            Zs, Zc = kernels.gi_inputs(view)
            if Zs.shape[1]:
                Kcc += Zc @ Zc.T
                Kcs += Zc @ Zs.T
                dS += np.einsum("ij,ij->i", Zs, Zs)
        if self.normalize:
            dc = np.sqrt(np.clip(np.diag(Kcc), 1e-12, None))
            ds = np.sqrt(np.clip(dS, 1e-12, None))
            Kcc = Kcc / np.outer(dc, dc)
            Kcs = Kcs / np.outer(dc, ds) if len(ds) else Kcs
        return Kcc, Kcs

    def _state(self, view):
        Kcc, Kcs = self.kernel_blocks(view)
        cov = np.maximum(Kcs.max(1), 0.0) if Kcs.shape[1] else np.zeros(len(Kcc))
        st = {"cov": cov}

        def gains():
            return np.maximum(Kcc - st["cov"][:, None], 0.0).sum(0)

        def update(g, i):
            st["cov"] = np.maximum(st["cov"], Kcc[:, i])
            return gains()
        return gains, update

    def gains(self, view) -> np.ndarray:
        return self._state(view)[0]()

    def _choose(self, view, k, rng) -> list:
        gains, update = self._state(view)
        self._diag(phase="static_external" if self.static else ("external_only" if view.n_screened == 0 else "mk"))
        return self._greedy(gains(), k, rng, inner=update)
