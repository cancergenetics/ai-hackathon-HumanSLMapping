"""FROZEN evaluation module (DESIGN v3 §3 metric, §6.4 statistics, §11 verification). Pure numpy/scipy/sklearn; no
project data access. Changing anything here after P0 requires a PI decision recorded in RESEARCH_TREE.md.

Metric
  average_precision(y, s, w=None)   pooled AP, tied scores grouped at ONE threshold (= sklearn
                                    average_precision_score semantics, incl. sample weights; no stable-sort
                                    tie-breaking)
  macro_ap(y, s, q)                 mean per-query AP over queries with >= 1 positive (+ count of zero-positive ones)
  auroc(y, s)                       sklearn roc_auc_score (ties count 1/2)
  precision_at_k(y, s, q, k)        per-query precision of the top k cells; cells tied at the cut-off counted in
                                    expectation (order-invariant)
  lift(ap, prevalence); summarise(y, s, q)
AP engine (exact weighted AP, DESIGN v3 §6.4)
  WeightedAP(y, s, q, n_queries)    scores sorted ONCE (descending), tie groups kept. For query weights w:
                                    AP = Σ_g (ΔTP_g / TP_tot) · TP_g / (TP_g + FP_g) over tie groups g, cumulative
                                    weighted TP/FP at tie-group ends. Implementation: the cells between two consecutive
                                    positive-containing tie groups are pooled per query into a sparse
                                    [segment x query] count matrix, so a block of replicates costs one sparse product;
                                    this is algebraically identical (FP only enters at positive groups).
                                    `ap_gather(w)` is the literal per-cell gather/cumsum reference.
  SortedAP(y, s)                    one sort; AP for arbitrary per-cell weights (two-way query x array bootstrap)
  bootstrap_weights(strata, B, rng) stratified query bootstrap (multinomial within strata)
  ReplicateEngine                   replicate array A[arm][b, r, n] (shared weights w^b for all arms) + point
                                    estimates est[arm][r, n]
Inference (all driven by the replicate array)
  welch(...)                        Δ̂ = mean_r A − mean_r B; SE² = V_Q + a + b (a = σ̂²_A/R_A, b = σ̂²_B/R_B),
                                    ν = SE⁴ / [a²/(R_A−1) + b²/(R_B−1)] (V_Q treated as known); p = 1 − T_ν(Δ̂/SE)
  paired(...)                       same-S predictor-change contrast: Var = Var(d_r)/R + V_Q(d̄)
  joint_resample(...)               joint bootstrap: query weights + whole trajectories resampled within arm
  gap_closure, g_star, n_req (PAVA + log-n interpolation, left endpoint of the reached flat segment, censored
  replicates coded 0 / +inf), pava, holm, bonferroni_bands, percentile_ci
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["average_precision", "macro_ap", "auroc", "precision_at_k", "lift", "summarise", "WeightedAP",
           "SortedAP", "bootstrap_weights", "ReplicateEngine", "welch", "paired", "joint_resample", "gap_closure",
           "g_star", "pava", "n_req", "NReq", "holm", "bonferroni_bands", "percentile_ci"]


def _protocol() -> dict:
    from hsl.config import protocol
    return protocol()


def _yv(y, s):
    y = np.asarray(y)
    s = np.asarray(s, dtype=np.float64)
    if y.shape != s.shape or y.ndim != 1:
        raise ValueError(f"y and s must be 1-D of equal length, got {y.shape} and {s.shape}")
    if not np.isfinite(s).all():
        raise ValueError("scores must be finite")
    if not np.isin(y, (0, 1)).all():
        raise ValueError("labels must be 0/1")
    return y.astype(bool), s


def _tie_ends(ss: np.ndarray) -> np.ndarray:
    """Last index of each tie group of a descending-sorted score vector."""
    return np.r_[np.flatnonzero(ss[1:] != ss[:-1]), len(ss) - 1] if len(ss) else np.zeros(0, dtype=np.int64)


def _ap_from_group_cums(tp: np.ndarray, fp: np.ndarray) -> float:
    if len(tp) == 0 or tp[-1] <= 0:
        return float("nan")
    tot = tp + fp
    prec = np.divide(tp, tot, out=np.zeros_like(tp), where=tot > 0)
    return float(np.sum(np.diff(np.r_[0.0, tp]) * prec) / tp[-1])


# ============================================================================ metric
def average_precision(y, s, w=None) -> float:
    """Σ over distinct score thresholds (descending) of ΔRecall × Precision; all cells with the same score enter at
    one threshold. Optional non-negative cell weights. NaN if there is no positive weight."""
    y, s = _yv(y, s)
    w = np.ones(len(s)) if w is None else np.asarray(w, dtype=np.float64)
    if w.shape != s.shape or (w < 0).any():
        raise ValueError("weights must be non-negative and aligned with s")
    order = np.argsort(-s, kind="mergesort")
    ss, yy, ww = s[order], y[order], w[order]
    ends = _tie_ends(ss)
    return _ap_from_group_cums(np.cumsum(ww * yy)[ends], np.cumsum(ww * ~yy)[ends])


def _groups(q):
    q = np.asarray(q)
    order = np.argsort(q, kind="mergesort")
    qs = q[order]
    starts = np.r_[0, np.flatnonzero(qs[1:] != qs[:-1]) + 1] if len(qs) else np.zeros(0, dtype=int)
    ends = np.r_[starts[1:], len(qs)]
    return order, qs[starts] if len(qs) else qs, starts, ends


def _key(k):
    return k.item() if hasattr(k, "item") else k


def macro_ap(y, s, q) -> dict:
    """Mean per-query AP over queries with >= 1 positive. {macro, n_queries, n_used, n_zero_pos, per_query}."""
    y, s = _yv(y, s)
    order, keys, st, en = _groups(q)
    per = {}
    for k, a, b in zip(keys, st, en):
        idx = order[a:b]
        per[_key(k)] = average_precision(y[idx], s[idx]) if y[idx].any() else np.nan
    vals = np.array(list(per.values()), dtype=float)
    used = ~np.isnan(vals)
    return {"macro": float(vals[used].mean()) if used.any() else float("nan"), "n_queries": int(len(vals)),
            "n_used": int(used.sum()), "n_zero_pos": int((~used).sum()), "per_query": per}


def auroc(y, s) -> float:
    from sklearn.metrics import roc_auc_score
    y, s = _yv(y, s)
    if y.all() or not y.any():
        return float("nan")
    return float(roc_auc_score(y, s))


def precision_at_k(y, s, q, k: int) -> dict:
    """{query: precision among its top min(k, n_q) cells}, ties at the cut-off counted by their expected share."""
    y, s = _yv(y, s)
    order, keys, st, en = _groups(q)
    out = {}
    for key, a, b in zip(keys, st, en):
        idx = order[a:b]
        sq, yq = s[idx], y[idx]
        kk = min(int(k), len(sq))
        v = np.sort(sq)[::-1][kk - 1]
        gt, eq = sq > v, sq == v
        exp_pos = yq[gt].sum() + (kk - gt.sum()) * yq[eq].sum() / eq.sum()
        out[_key(key)] = float(exp_pos / kk)
    return out


def lift(ap: float, prevalence: float) -> float:
    return float(ap / prevalence) if prevalence > 0 else float("nan")


def summarise(y, s, q, ks=None) -> dict:
    """Standard metric bundle for one score vector on one evaluation set."""
    ks = _protocol()["eval"]["precision_at_k"] if ks is None else ks
    y, s = _yv(y, s)
    ap = average_precision(y, s)
    prev = float(y.mean()) if len(y) else float("nan")
    m = macro_ap(y, s, q)
    out = {"n_cells": int(len(y)), "n_pos": int(y.sum()), "prevalence": prev, "auprc": ap, "lift": lift(ap, prev),
           "macro_auprc": m["macro"], "macro_n_used": m["n_used"], "macro_n_zero_pos": m["n_zero_pos"],
           "auroc": auroc(y, s)}
    for k in ks:
        p = np.array(list(precision_at_k(y, s, q, k).values()))
        out[f"p_at_{k}"] = float(p.mean()) if len(p) else float("nan")
    return out


# ============================================================================ AP engine
class WeightedAP:
    """Exact AP of one score vector under arbitrary query weights (see module docstring)."""

    def __init__(self, y, s, q, n_queries: int):
        from scipy import sparse
        y, s = _yv(y, s)
        q = np.asarray(q, dtype=np.int64)
        if q.shape != s.shape or (len(q) and ((q < 0).any() or (q >= n_queries).any())):
            raise ValueError("q must hold query indices 0..n_queries-1 aligned with s")
        self.n_queries = int(n_queries)
        order = np.argsort(-s, kind="mergesort")
        ss = s[order]
        self.q = q[order].astype(np.int32)
        self.y = y[order]
        self.ends = _tie_ends(ss)
        # tie group of each cell; positive-containing groups; segment = index of the next positive group at/after
        gid = np.zeros(len(ss), dtype=np.int64)
        if len(ss):
            gid[1:] = np.cumsum(ss[1:] != ss[:-1])
        pos_groups = np.unique(gid[self.y])
        seg = np.searchsorted(pos_groups, gid, side="left")       # cells after the last positive group -> n_seg
        keep = seg < len(pos_groups)
        n_seg = len(pos_groups)
        self.n_seg = n_seg
        shape = (max(n_seg, 1), self.n_queries)
        yk, sk, qk = self.y[keep], seg[keep], self.q[keep]
        self._pos = sparse.csr_matrix((np.ones(int(yk.sum())), (sk[yk], qk[yk])), shape=shape)
        self._neg = sparse.csr_matrix((np.ones(int((~yk).sum())), (sk[~yk], qk[~yk])), shape=shape)

    def ap(self, W) -> np.ndarray:
        """AP for each row of query weights W [b x n_queries] (or one vector) -> [b]."""
        W = np.atleast_2d(np.asarray(W, dtype=np.float64))
        if self.n_seg == 0:
            return np.full(len(W), np.nan)
        P = np.asarray((self._pos @ W.T).T)        # [b x n_seg]: weighted positives of positive group k
        N = np.asarray((self._neg @ W.T).T)        # weighted negatives entering at or before group k (after k-1)
        tp = np.cumsum(P, axis=1)
        tot = tp + np.cumsum(N, axis=1)
        prec = np.divide(tp, tot, out=np.zeros_like(tp), where=tot > 0)
        Ptot = tp[:, -1]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(Ptot > 0, (P * prec).sum(1) / Ptot, np.nan)

    def ap_gather(self, w) -> float:
        """Literal reference: gather w[q] per sorted cell, cumulate at tie-group ends."""
        wc = np.asarray(w, dtype=np.float64)[self.q]
        return _ap_from_group_cums(np.cumsum(wc * self.y)[self.ends], np.cumsum(wc * ~self.y)[self.ends])


class SortedAP:
    """Sort once; exact tie-grouped AP for arbitrary non-negative per-cell weights (two-way bootstrap w_t · v_a)."""

    def __init__(self, y, s):
        y, s = _yv(y, s)
        self.order = np.argsort(-s, kind="mergesort")
        self.y = y[self.order]
        self.ends = _tie_ends(s[self.order])

    def ap(self, w) -> float:
        w = np.asarray(w, dtype=np.float64)[self.order]
        return _ap_from_group_cums(np.cumsum(w * self.y)[self.ends], np.cumsum(w * ~self.y)[self.ends])


def bootstrap_weights(strata, B: int, rng) -> np.ndarray:
    """Stratified query bootstrap weights [B x n_queries] (int32 counts): within each stratum of n_h queries, n_h
    draws with replacement. `rng` = numpy Generator, SeedSequence or int."""
    rng = rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)
    strata = np.asarray(strata)
    W = np.zeros((int(B), len(strata)), dtype=np.int32)
    for h in np.unique(strata):
        idx = np.flatnonzero(strata == h)
        W[:, idx] = rng.multinomial(len(idx), np.full(len(idx), 1.0 / len(idx)), size=int(B))
    return W


@dataclass
class ReplicateEngine:
    """One replicate array per arm: A[arm][b, r, j] = AP of trajectory r at budget grid[j] under query weights W[b]
    (the SAME W for all arms); est[arm][r, j] = AP with unit weights. Missing (r, j) stay NaN."""
    W: np.ndarray
    grid: tuple
    R: dict                               # arm -> number of trajectories
    chunk: int | None = None
    A: dict = field(default_factory=dict)
    est: dict = field(default_factory=dict)

    def __post_init__(self):
        self.W = np.asarray(self.W)
        self.grid = tuple(int(n) for n in self.grid)
        self.chunk = int(_protocol()["bootstrap"]["chunk"] if self.chunk is None else self.chunk)
        for arm, r in self.R.items():
            self.A[arm] = np.full((len(self.W), int(r), len(self.grid)), np.nan)
            self.est[arm] = np.full((int(r), len(self.grid)), np.nan)

    def add(self, arm, r: int, n: int, y, s, q) -> None:
        j = self.grid.index(int(n))
        e = WeightedAP(y, s, q, self.W.shape[1])
        self.A[arm][:, r, j] = np.concatenate([e.ap(self.W[i:i + self.chunk])
                                               for i in range(0, len(self.W), self.chunk)])
        self.est[arm][r, j] = float(e.ap(np.ones(self.W.shape[1]))[0])

    def at(self, arm, n: int):
        """(replicates [B x R], estimates [R]) of `arm` at budget n."""
        j = self.grid.index(int(n))
        return self.A[arm][:, :, j], self.est[arm][:, j]


# ============================================================================ inference
def _tdist_p(t, df):
    from scipy import stats
    return float(stats.t.sf(t, df)) if np.isfinite(df) else float(stats.norm.sf(t))


def welch(boot_a, est_a, boot_b, est_b, margin: float = 0.0) -> dict:
    """Unpaired arm contrast Δ̂ = mean_r est_a − mean_r est_b (H1: Δ > margin).
    boot_* [B x R_*]: replicate APs under shared query weights; V_Q = Var_b(mean_r boot_a − mean_r boot_b).
    SE² = V_Q + a + b, a = σ̂²_a/R_a, b = σ̂²_b/R_b; ν = SE⁴/[a²/(R_a−1) + b²/(R_b−1)] (V_Q known; ∞ if a = b = 0)."""
    from scipy import stats
    ea, eb = np.atleast_1d(np.asarray(est_a, float)), np.atleast_1d(np.asarray(est_b, float))
    ba = np.asarray(boot_a, float).reshape(len(boot_a), -1)
    bb = np.asarray(boot_b, float).reshape(len(boot_b), -1)
    d_hat = float(ea.mean() - eb.mean())
    VQ = float(np.var(ba.mean(1) - bb.mean(1), ddof=1))
    Ra, Rb = len(ea), len(eb)
    a = float(np.var(ea, ddof=1) / Ra) if Ra > 1 else 0.0
    b = float(np.var(eb, ddof=1) / Rb) if Rb > 1 else 0.0
    se2 = VQ + a + b
    den = (a ** 2 / (Ra - 1) if Ra > 1 else 0.0) + (b ** 2 / (Rb - 1) if Rb > 1 else 0.0)
    nu = float(se2 ** 2 / den) if den > 0 else float("inf")
    se = float(np.sqrt(se2))
    t = (d_hat - margin) / se if se > 0 else float(np.sign(d_hat - margin) * np.inf)
    crit = float(stats.t.ppf(0.975, nu)) if np.isfinite(nu) else float(stats.norm.ppf(0.975))
    return {"delta": d_hat, "margin": float(margin), "V_Q": VQ, "a": a, "b": b, "se": se, "df": nu, "t": float(t),
            "p_one_sided": _tdist_p(t, nu), "ci95": (d_hat - crit * se, d_hat + crit * se), "R_a": Ra, "R_b": Rb,
            "traj_var_missing": bool(Ra < 2 or Rb < 2)}


def paired(boot_c, est_c, boot_p, est_p) -> dict:
    """Paired predictor-change contrast on the SAME trajectories (child c vs parent p): d_r = est_c,r − est_p,r;
    Var = Var(d_r)/R + V_Q(d̄), V_Q(d̄) = Var_b(mean_r (boot_c − boot_p)); ν = Var²/[(Var(d)/R)²/(R−1)]."""
    from scipy import stats
    d = np.asarray(est_c, float) - np.asarray(est_p, float)
    R = len(d)
    dd = np.asarray(boot_c, float).reshape(len(boot_c), -1) - np.asarray(boot_p, float).reshape(len(boot_p), -1)
    VQ = float(np.var(dd.mean(1), ddof=1))
    tr = float(np.var(d, ddof=1) / R) if R > 1 else 0.0
    var = VQ + tr
    nu = float(var ** 2 / (tr ** 2 / (R - 1))) if R > 1 and tr > 0 else float("inf")
    se = float(np.sqrt(var))
    dbar = float(d.mean())
    t = dbar / se if se > 0 else float(np.sign(dbar) * np.inf)
    crit = float(stats.t.ppf(0.975, nu)) if np.isfinite(nu) else float(stats.norm.ppf(0.975))
    return {"delta": dbar, "rel_gain": float(dbar / np.mean(est_p)), "V_Q": VQ, "traj": tr, "se": se, "df": nu,
            "p_one_sided": _tdist_p(t, nu), "ci95": (dbar - crit * se, dbar + crit * se), "R": R}


def joint_resample(A_arm: np.ndarray, rng) -> np.ndarray:
    """Joint bootstrap of one arm: for replicate b (query weights w^b already in A_arm[b]) resample R whole
    trajectories with replacement -> trajectory-mean curve M[b, j] (NaN-aware over budgets)."""
    rng = rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)
    B, R, _ = A_arm.shape
    idx = rng.integers(0, R, size=(B, R))
    return np.nanmean(A_arm[np.arange(B)[:, None], idx, :], axis=1)


def gap_closure(a_s, a_r_lo, a_r_hi):
    """G = (A_s(n) − A_r(n)) / (A_r(n_hi) − A_r(n)) (vectorised over replicates)."""
    a_s, lo, hi = (np.asarray(x, float) for x in (a_s, a_r_lo, a_r_hi))
    with np.errstate(divide="ignore", invalid="ignore"):
        return (a_s - lo) / (hi - lo)


def g_star(r_lo, r_hi, margin_frac: float | None = None):
    """G* = 1 − m·r_hi/(r_hi − r_lo): the gap closure a 'matches within m·r_hi' claim implies (discriminating if ≥ 0.5)."""
    m = float(_protocol()["stats"]["match_margin_frac"] if margin_frac is None else margin_frac)
    lo, hi = np.asarray(r_lo, float), np.asarray(r_hi, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 - m * hi / (hi - lo)


def pava(y, w=None) -> np.ndarray:
    """Weighted least-squares non-decreasing fit (pool adjacent violators)."""
    y = np.asarray(y, dtype=float)
    w = np.ones_like(y) if w is None else np.asarray(w, dtype=float)
    vals, wts, cnt = [], [], []
    for yi, wi in zip(y, w):
        vals.append(yi)
        wts.append(wi)
        cnt.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, w2, c2 = vals.pop(), wts.pop(), cnt.pop()
            v1, w1, c1 = vals.pop(), wts.pop(), cnt.pop()
            vals.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wts.append(w1 + w2)
            cnt.append(c1 + c2)
    return np.repeat(vals, cnt)


@dataclass(frozen=True)
class NReq:
    n: float          # random-equivalent budget; 0 if left-censored, +inf if right-censored
    censored: str     # "" | "left" (target < fit(n_min)) | "right" (target > fit(n_max))


def n_req(target: float, grid_n, curve, w=None) -> NReq:
    """Budget at which the PAVA fit f of the random curve (trajectory means on the dense grid), linearly interpolated
    in log n, first reaches `target`: if f[k−1] < target ≤ f[k] the value is interpolated on (n[k−1], n[k]]; reaching
    a flat segment gives its LEFT endpoint. Outside [f(n_min), f(n_max)]: censored, coded 0 / +inf (never
    extrapolated). Also used for n_80 / n_90 (target = frac · fit at the pool maximum)."""
    n = np.asarray(grid_n, dtype=float)
    if not (np.diff(n) > 0).all():
        raise ValueError("grid_n must be strictly increasing")
    curve = np.asarray(curve, dtype=float)
    ok = np.isfinite(curve)
    if not ok.all():
        n, curve = n[ok], curve[ok]
        w = None if w is None else np.asarray(w, float)[ok]
    f = pava(curve, w)
    t = float(target)
    if not np.isfinite(t) or len(f) == 0:
        return NReq(float("nan"), "")
    if t < f[0]:
        return NReq(0.0, "left")
    if t > f[-1]:
        return NReq(float("inf"), "right")
    k = int(np.argmax(f >= t))
    if k == 0 or f[k - 1] >= t or f[k] == t:          # target reached exactly at a grid point / flat segment start
        return NReq(float(n[k]), "")
    x0, x1 = np.log(n[k - 1]), np.log(n[k])
    frac = (t - f[k - 1]) / (f[k] - f[k - 1])
    return NReq(float(np.exp(x0 + frac * (x1 - x0))), "")


def holm(p, alpha: float | None = None) -> dict:
    """Holm step-down: adjusted p_(i) = max_{j<=i} min(1, (m − j + 1) p_(j)); reject where adjusted <= alpha."""
    alpha = float(_protocol()["stats"]["alpha_familywise"] if alpha is None else alpha)
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    adj_sorted = np.maximum.accumulate(np.minimum(1.0, (m - np.arange(m)) * p[order]))
    adj = np.empty(m)
    adj[order] = adj_sorted
    return {"p_adj": adj, "reject": adj <= alpha, "alpha": alpha}


def bonferroni_bands(est, se, df=None, level: float | None = None) -> dict:
    """Simultaneous bands over K checkpoints: est_k ± t_{1−(1−level)/(2K), ν_k} · se_k."""
    from scipy import stats
    level = float(_protocol()["stats"]["ci_level"] if level is None else level)
    est, se = np.asarray(est, float), np.asarray(se, float)
    K = len(est)
    q = 1 - (1 - level) / (2 * K)
    df = np.full(K, np.inf) if df is None else np.asarray(df, float)
    crit = np.array([stats.t.ppf(q, d) if np.isfinite(d) else stats.norm.ppf(q) for d in df])
    return {"lo": est - crit * se, "hi": est + crit * se, "crit": crit}


def percentile_ci(x, level: float | None = None) -> tuple[float, float]:
    """Percentile interval of finite values (±inf censored codes are kept: they bound the interval)."""
    level = float(_protocol()["stats"]["ci_level"] if level is None else level)
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    a = (1 - level) / 2
    return float(np.quantile(x, a)), float(np.quantile(x, 1 - a))
