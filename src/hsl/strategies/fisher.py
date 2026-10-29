"""F Fisher-D (main hybrid) and F-static (DESIGN v3.1 §6.2; review r1 item 1, r2 item 13).

Objective: Σ_{x∈P} log det M_x, M_x = δI_r + Σ_{s∈S} w_sx u_s u_sᵀ, w_sx = p̂_sx(1 − p̂_sx),
  P = on-array candidates ∪ visible array-only columns (= visible columns that are neither candidates nor screened).
Fitted phase (n ≥ h_switch_n):
  u_g = (1, ũ_g); ũ = r−1 truncated-SVD factors of the visible screened block R (row-centred clipped −ε, the mk GI
  signal; R ≈ U Σ Vᵀ): screened s -> its row factor U_s; on-array candidate g -> û_g = B·v_g, v_g = Uᵀ R[:, col_g]
  (its column factor) with B the (r−1)×(r−1) least-squares map fitted on the screened on-array genes (v_s -> U_s);
  off-array candidates -> KRR from the external kernel (fitted on all screened genes) into the same space. ũ whitened
  to unit mean square per component (over screened ∪ candidates).
  p̂_sx = clip(q̂_s·ĉ_x/π̄_S, 1e-4, 0.5); q̂_s = observed smoothed row rate for screened s, sigmoid(degree model) for
  candidates; δ = mean p̂(1 − p̂) over S × P.
  Gain of g (determinant lemma): Σ_x log(1 + w_gx·u_gᵀ M_x⁻¹ u_g); after each pick M_x⁻¹ is downdated by
  Sherman–Morrison (within the batch only; M_x is rebuilt each round).
Cold start (n < h_switch_n, incl. S = ∅) and F-static (always): p̂ ≡ const, so the objective is |P|·log det(I + Σ_s
  u_s u_sᵀ) with ũ = top r−1 kernel-PCA components of the external kernel over the pool (whitened); F-static never
  reads screened rows.
r = 1 (u ≡ 1): gains are monotone in q̂_g, i.e. exactly H's ordering (built-in ablation); in the cold start the r = 1
  objective is constant, so F(r = 1) uses H's prior there (and F-static(r = 1) = H-static).
Off-array route calibration (review F-1; knob offarray_calibration): the KRR route systematically yields smaller
  u_gᵀM_x⁻¹u_g than the B·v route (a route-swap control moves on-array genes from the base share to ~1%), i.e. it
  acts as an off-array score term (DESIGN §1). With 'quantile' (route-stratified, as H and C), each off-array gene's
  gain is replaced by the quantile of the ON-array candidates' gains matching its percentile among the remaining
  OFF-array candidates' gains; recomputed after every pick of the batch (picked genes leave both references). A
  common-route variant (KRR route for every candidate) was tried and drifts low along F's trajectory (0.09–0.12 vs
  base 0.19 on dev1). On-array gains and every M_x update are unchanged.
Batch cap as H (hsl.strategies.base.batch_cap).
"""
from __future__ import annotations

import numpy as np

from hsl.strategies import base, degree, kernels


class Fisher(base.Selector):
    name = "F"
    reads = frozenset({"rows", "smf"})

    def __init__(self, static: bool = False, r: int | None = None, **kw):
        cfg = base.knobs()["F"]
        self.r = int(cfg["r"] if r is None else r)
        self.krr_lam = float(kw.pop("krr_lam", cfg["krr_lam"]))
        self.b_ridge = float(kw.pop("b_ridge", cfg["b_ridge"]))
        self.calibrate = base.offarray_calibration(kw.pop("offarray_calibration", None))
        super().__init__(static=bool(static), r=self.r, krr_lam=self.krr_lam, b_ridge=self.b_ridge,
                         offarray_calibration="quantile" if self.calibrate else "none", **kw)
        self.static = bool(static)
        if self.static:
            self.name = "F-static"
            self.reads = frozenset() if self.r > 1 else frozenset({"smf"})
        self._hub = None

    def reset(self, candidates, rng) -> None:
        super().reset(candidates, rng)
        self._hub = None

    def _hub_sel(self):
        from hsl.strategies.hub import Hub
        if self._hub is None:
            self._hub = Hub(static=self.static, deterministic=self.deterministic, top_k=self.top_k)
            self._hub._feats = self._feats
        return self._hub

    # ------------------------------------------------------------------ phase dispatch
    def _phase(self, view) -> str:
        if self.static or view.n_screened < int(base.protocol_selector()["h_switch_n"]):
            return "cold"
        return "fitted"

    def gains(self, view) -> np.ndarray:
        st = self._state(view)
        return st["gains"]()

    def _choose(self, view, k, rng) -> list:
        st = self._state(view)
        self._diag(phase=st["phase"])
        return self._greedy(st["gains"](), k, rng, clusters=view.clusters(view.candidates), inner=st["update"])

    def _state(self, view) -> dict:
        phase = self._phase(view)
        if self.r <= 1 and phase == "cold":
            hub = self._hub_sel()
            g = hub.scores(view)[0]
            return {"phase": "cold_r1_prior", "gains": lambda: g, "update": None}
        return self._cold_state(view) if phase == "cold" else self._fitted_state(view)

    # ------------------------------------------------------------------ cold start / static
    def _cold_state(self, view) -> dict:
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        cand, S = np.asarray(view.candidates, dtype=str), np.asarray(view.screened_orfs, dtype=str)
        E = kernels.external_kpca(view, pool, feats, self.r - 1)
        E = E / np.sqrt(np.maximum((E ** 2).mean(0), 1e-300))          # whiten (unit mean square) over the pool
        u = np.column_stack([np.ones(len(pool)), E])
        uc, us = u[np.searchsorted(pool, cand)], u[np.searchsorted(pool, S)]
        Minv = np.linalg.inv(np.eye(self.r) + us.T @ us)
        st = {"Minv": Minv}

        def gains():
            q = np.einsum("gi,ij,gj->g", uc, st["Minv"], uc)
            return np.log1p(np.maximum(q, 0))

        def update(g, i):
            b = st["Minv"] @ uc[i]
            st["Minv"] = st["Minv"] - np.outer(b, b) / (1.0 + uc[i] @ b)
            return gains()
        return {"phase": "cold", "gains": gains, "update": update}

    # ------------------------------------------------------------------ fitted phase
    def embeddings(self, view, feats, pool, krr_all: bool = False):
        """(u_S [n x r], u_C [n_cand x r]) with u = (1, ũ); with krr_all also u_K [n_cand x r]: the external-KRR route
        for EVERY candidate (same whitening), used as the common route of the off-array calibration."""
        r1 = self.r - 1
        S, cand = np.asarray(view.screened_orfs, dtype=str), np.asarray(view.candidates, dtype=str)
        n, nc = len(S), len(cand)
        if r1 <= 0:
            return (np.ones((n, 1)), np.ones((nc, 1))) + ((np.ones((nc, 1)),) if krr_all else ())
        vis = np.asarray(view.visible_cols)
        import warnings
        from ysl.models import baselines as B, completion as CM
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # upstream mk GI signal: −clip(ε), row-centred
            R = CM.centre(-np.clip(np.asarray(view.E, dtype=np.float64)[:, vis], *B.EPS_CLIP), "row")
        G = R @ R.T
        w, V = np.linalg.eigh((G + G.T) / 2)
        top = np.argsort(w)[::-1][: min(r1, n)]
        U = V[:, top]                                            # row factors of the screened genes [n x r1']
        sign = np.sign(U[np.abs(U).argmax(0), np.arange(U.shape[1])])
        sign[sign == 0] = 1.0
        U = U * sign
        r1e = U.shape[1]
        col_of = {c: j for j, c in enumerate(vis)}
        # column factors v_x = Uᵀ R[:, x] (coefficients of column x on the row-factor basis)
        ccol, scol = np.asarray(view.candidate_col), np.asarray(view.screened_col)
        cpos = np.array([col_of.get(c, -1) for c in ccol])
        spos = np.array([col_of.get(c, -1) for c in scol])
        uS = U.copy()
        uC = np.zeros((nc, r1e))
        on_s, on_c = spos >= 0, cpos >= 0
        if on_s.sum() >= 2 and on_c.any():
            Vs = (R[:, spos[on_s]].T @ U)                         # [n_on x r1e]
            A = Vs.T @ Vs
            reg = self.b_ridge * (np.trace(A) / max(r1e, 1) or 1.0)
            Bt = np.linalg.solve(A + reg * np.eye(r1e), Vs.T @ U[on_s])   # U_on ≈ Vs Bᵀ
            uC[on_c] = (R[:, cpos[on_c]].T @ U) @ Bt
        off = ~on_c
        uK = None
        if off.any() or krr_all:                                 # external-kernel KRR into the same space
            K = kernels.external_kernel(view, pool, feats)
            ps, pc = np.searchsorted(pool, S), np.searchsorted(pool, cand)
            Kss = K[np.ix_(ps, ps)]
            lam_k = self.krr_lam * (np.trace(Kss) / max(n, 1) or 1.0)
            alpha = np.linalg.solve(Kss + lam_k * np.eye(n), U)
            uK = K[np.ix_(pc, ps)] @ alpha
            uC[off] = uK[off]
        allu = np.vstack([uS, uC])
        scale = np.sqrt(np.maximum((allu ** 2).mean(0), 1e-300))
        uS, uC = uS / scale, uC / scale
        uK = uK / scale if uK is not None else None

        def pad(a, m):
            if a is None:
                return None
            if r1e < r1:                                        # n < r−1 (cannot happen for n ≥ 30, r ≤ 31)
                a = np.hstack([a, np.zeros((m, r1 - r1e))])
            return np.column_stack([np.ones(m), a])
        out = (pad(uS, n), pad(uC, nc))
        return out + ((pad(uK, nc),) if krr_all else ())

    def target_cols(self, view) -> np.ndarray:
        """P: on-array candidates' columns ∪ visible columns that are not pool genes (array-only)."""
        vis = np.asarray(view.visible_cols)
        pool_cols = np.concatenate([np.asarray(view.candidate_col), np.asarray(view.screened_col)])
        ao = np.setdiff1d(vis, pool_cols[pool_cols >= 0])
        cc = np.asarray(view.candidate_col)
        return np.union1d(cc[cc >= 0], ao)

    def _fitted_state(self, view) -> dict:
        lo, hi = base.knobs()["spec"]["p_clip"]
        pool = base.pool_of(view)
        feats = self._feats.get(view, pool)
        cand = np.asarray(view.candidates, dtype=str)
        qlogit, info = degree.fit_predict(view, cand, feats, pool)
        self._diag(**{f"qhat_{k}": v for k, v in info.items()})
        q_c = degree.sigmoid(qlogit)
        q_s = degree.screened_rate(view)
        c_hat, pi = degree.shrunk_col_rate(view)
        P = self.target_cols(view)
        ratio = c_hat[P] / pi
        p_s = np.clip(q_s[:, None] * ratio[None, :], lo, hi)
        W_s = p_s * (1 - p_s)                                   # [n x |P|]
        delta = float(W_s.mean())
        p_c = np.clip(q_c[:, None] * ratio[None, :], lo, hi)
        W_c = p_c * (1 - p_c)                                   # [n_cand x |P|]
        on = np.asarray(view.candidate_col) >= 0
        calib = self.calibrate and on.any() and (~on).any() and self.r > 1
        uS, uC = self.embeddings(view, feats, pool)
        r = uS.shape[1]
        outer = np.einsum("si,sj->sij", uS, uS).reshape(len(uS), r * r)
        M = (W_s.T @ outer).reshape(len(P), r, r) + delta * np.eye(r)[None]
        Minv = np.linalg.inv(M)                                  # [|P| x r x r]
        uu = np.einsum("gi,gj->gij", uC, uC).reshape(len(uC), r * r)
        Q = uu @ Minv.reshape(len(P), r * r).T                   # q_gx = u_gᵀ M_x⁻¹ u_g  [n_cand x |P|]
        st = {"Minv": Minv, "Q": Q, "picked": np.zeros(len(uC), dtype=bool)}
        self._diag(n_P=int(len(P)), delta=delta, r=int(r), offarray_route_calibrated=bool(calib))
        self._qhat = qlogit

        def gains():
            g = np.log1p(W_c * np.maximum(st["Q"], 0)).sum(1)
            return base.calibrate_routes(g, None, on, exclude=st["picked"]) if calib else g

        def update(g, i):
            Minv, Q = st["Minv"], st["Q"]
            b = np.einsum("xij,j->xi", Minv, uC[i])              # M_x⁻¹ u_g  [|P| x r]
            f = W_c[i] / (1.0 + W_c[i] * Q[i])                   # [|P|]
            cross = uC @ b.T                                     # u_hᵀ M_x⁻¹ u_g  [n_cand x |P|]
            st["Q"] = Q - f[None, :] * cross ** 2
            st["Minv"] = Minv - f[:, None, None] * np.einsum("xi,xj->xij", b, b)
            st["picked"][i] = True
            return gains()
        return {"phase": "fitted", "gains": gains, "update": update, "M": M, "W_c": W_c, "uC": uC}
