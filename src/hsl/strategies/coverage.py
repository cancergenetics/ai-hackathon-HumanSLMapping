"""C coverage (DESIGN v3.1 §6.2): facility location over the predictor's `mk` kernel K = K_gi + Σ w_k K_ext,k
(hsl.strategies.kernels: GI 'sym' kernel over the screened on-array coordinates + tier-1 external kernels), target
pool = unscreened candidates. Objective F(S) = Σ_{x∈cand} max(0, max_{s∈S} K(x, s)); gain of c =
Σ_x max(0, K(x, c) − cov_x). Exact greedy gains are recomputed after every pick (vectorised, O(|cand|²) per pick —
cheaper than lazy-greedy bookkeeping at this size; identical result), with the randomised top-k rule.
Round 0 (S = ∅): the GI part is empty, so C runs on the external kernel alone (the reported "static part").
C-static (not an arm; diagnostic) uses the external kernel only at every round and never reads screened rows.
Off-array route calibration (review C-1; knob offarray_calibration): off-array candidates have no GI column, so their
GI rows are 0 and their coverage gains are systematically smaller (0% off-array picks after round 0 on a null toy;
`normalize` over-corrects to 68–86%). With 'quantile', each off-array gene's gain is replaced by the quantile of the
ON-array candidates' gains that matches its percentile among the remaining OFF-array candidates' gains (route-
stratified within the joint greedy; recomputed after every pick). A common-route reference (external-only coverage,
as F/H) was tried and drifts low along C's own trajectory: off-array genes are picked on the external route only, so
the remaining ones are depleted on it while on-array genes (picked on mk gains, ρ ≈ 0.1 with the external route) are
not (null toy: 0.11 vs base 0.17–0.19). On-array gains are untouched. Round 0 (S = ∅) needs no calibration (one
route)."""
from __future__ import annotations

import numpy as np

from hsl.strategies import base, kernels


class Coverage(base.Selector):
    name = "C"
    reads = frozenset({"rows"})

    def __init__(self, static: bool = False, normalize: bool | None = None, **kw):
        nm = bool(base.knobs()["C"]["normalize"] if normalize is None else normalize)
        self.calibrate = base.offarray_calibration(kw.pop("offarray_calibration", None))
        super().__init__(static=bool(static), normalize=nm,
                         offarray_calibration="quantile" if self.calibrate else "none", **kw)
        self.static, self.normalize = bool(static), nm
        if self.static:
            self.name, self.reads = "C-static", frozenset()

    def kernel_blocks(self, view, gi: bool | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(K_cc [n_cand x n_cand], K_cs [n_cand x n_screened]); gi=False: external part only."""
        gi = (not self.static) if gi is None else bool(gi)
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        cand, S = np.asarray(view.candidates, dtype=str), np.asarray(view.screened_orfs, dtype=str)
        K = kernels.external_kernel(view, pool, feats)
        pc, ps = np.searchsorted(pool, cand), np.searchsorted(pool, S)
        Kcc = np.array(K[np.ix_(pc, pc)])
        Kcs = np.array(K[np.ix_(pc, ps)])
        dS = np.array(K[ps, ps])
        if gi:
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
        on = np.asarray(view.candidate_col) >= 0
        Kcc, Kcs = self.kernel_blocks(view)
        calib = (self.calibrate and not self.static and view.n_screened > 0 and on.any() and (~on).any())
        st = {"cov": np.maximum(Kcs.max(1), 0.0) if Kcs.shape[1] else np.zeros(len(Kcc)),
              "picked": np.zeros(len(Kcc), dtype=bool)}
        self._diag(offarray_route_calibrated=bool(calib))

        def gains():
            g = np.maximum(Kcc - st["cov"][:, None], 0.0).sum(0)
            return base.calibrate_routes(g, None, on, exclude=st["picked"]) if calib else g

        def update(g, i):
            st["cov"] = np.maximum(st["cov"], Kcc[:, i])
            st["picked"][i] = True
            return gains()
        return gains, update

    def gains(self, view) -> np.ndarray:
        return self._state(view)[0]()

    def _choose(self, view, k, rng) -> list:
        gains, update = self._state(view)
        self._diag(phase="static_external" if self.static else ("external_only" if view.n_screened == 0 else "mk"))
        return self._greedy(gains(), k, rng, inner=update)
