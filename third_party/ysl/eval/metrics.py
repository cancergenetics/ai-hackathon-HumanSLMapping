"""Evaluation metrics (DESIGN.md §1, §8, §9). PI-owned; frozen by config/frozen.yaml (AGENTS.md rule 1).

Inputs are flat arrays over the scored (observed) test cells.

  average_precision(y, s, w)  : step-wise AP exactly as sklearn.metrics.average_precision_score, incl. ties (all cells
                                with equal score form one threshold) and sample weights. The pooled AUPRC (headline)
                                is sklearn's own function; tests assert both agree to 1e-12.
  pooled_auprc / pooled_auroc : sklearn average_precision_score / roc_auc_score over all cells.
  macro_auprc(y, s, groups)   : mean per-group AP over groups with >= 1 positive (vectorised: one sort by
                                (group, -score), tie blocks within groups, segment sums), plus its exact random
                                baseline: mean over the same groups of E[AP | random order] =
                                (H_n + (m-1)(n-H_n)/(n-1)) / n for a group with n cells and m positives.
  precision_at_k(y, s, ks)    : tie-averaged precision@k = (positives strictly above the k-th score's tie block +
                                (k - n_above) * positive fraction of that block) / k. Independent of row order.
  SortedAP                    : one sort of (y, s, cluster) supporting AP under arbitrary per-cluster weights,
                                batched (sparse segment matrices) — used by the query-ORF cluster bootstrap.
  cluster_bootstrap / paired_bootstrap / joint_bootstrap : query-ORF cluster bootstrap (resample ORFs with
                                replacement, all their cells together; arrays fixed), paired across methods (same
                                weights), and jointly across splits (one ORF weight vector applied in all splits).
  stratified                  : AUPRC, prevalence, lift (= AUPRC / prevalence), n, n_pos per boolean stratum.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from sklearn.metrics import average_precision_score, roc_auc_score


def check_scores(y, s):
    y = np.asarray(y).astype(bool).ravel()
    s = np.asarray(s, dtype=np.float64).ravel()
    if y.shape != s.shape:
        raise ValueError("labels and scores differ in shape")
    if not np.isfinite(s).all():
        raise ValueError("scores contain NaN/inf; every scored cell needs a finite score (no silent fallbacks)")
    return y, s


_check = check_scores  # backwards-compatible name


# ----------------------------------------------------------------------------- AP
def average_precision(y, s, w=None) -> float:
    """sklearn-identical AP (ties = one threshold), optional non-negative sample weights."""
    y, s = check_scores(y, s)
    w = np.ones(len(y)) if w is None else np.asarray(w, dtype=np.float64).ravel()
    order = np.argsort(-s, kind="mergesort")
    s_o, y_o, w_o = s[order], y[order], w[order]
    ends = np.flatnonzero(np.r_[s_o[1:] != s_o[:-1], True])
    tp = np.cumsum(w_o * y_o)[ends]
    fp = np.cumsum(w_o * ~y_o)[ends]
    if tp[-1] <= 0:
        return float("nan")
    ps = tp + fp
    prec = np.divide(tp, ps, out=np.zeros_like(tp), where=ps > 0)
    rec = tp / tp[-1]
    return float(np.sum(np.diff(np.r_[0.0, rec]) * prec))


def pooled_auprc(y, s) -> float:
    y, s = check_scores(y, s)
    return float(average_precision_score(y, s)) if y.any() else float("nan")


def pooled_auroc(y, s) -> float:
    y, s = check_scores(y, s)
    return float(roc_auc_score(y, s)) if 0 < y.sum() < len(y) else float("nan")


# ----------------------------------------------------------------------------- macro AP
def expected_random_ap(n, m) -> np.ndarray:
    """Exact E[AP] of a uniformly random ranking (no ties) of n items with m >= 1 positives."""
    n = np.asarray(n, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    h = _harmonic(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = (h + np.where(n > 1, (m - 1) * (n - h) / (n - 1), 0.0)) / n
    return e


def _harmonic(n: np.ndarray) -> np.ndarray:
    from scipy.special import digamma
    return digamma(np.asarray(n, dtype=np.float64) + 1.0) + np.euler_gamma


def macro_auprc(y, s, groups) -> dict:
    """Mean per-group AP over groups with >= 1 positive (sklearn tie handling within each group) + exact random
    baseline over the same groups."""
    y, s = check_scores(y, s)
    g = np.asarray(groups).ravel()
    if g.shape != y.shape:
        raise ValueError("groups differ in shape")
    _, gi = np.unique(g, return_inverse=True)
    order = np.lexsort((-s, gi))
    gi_o, s_o, y_o = gi[order], s[order], y[order].astype(np.float64)
    n = len(y)
    new_group = np.r_[True, gi_o[1:] != gi_o[:-1]]
    block_start = new_group | np.r_[True, s_o[1:] != s_o[:-1]]
    bstart = np.flatnonzero(block_start)
    bend = np.r_[bstart[1:], n] - 1
    bgroup = gi_o[bstart]
    tp_blk = np.add.reduceat(y_o, bstart)
    cnt_blk = np.diff(np.r_[bstart, n]).astype(np.float64)
    # cumulative within group: global cumsum minus the value before the group start
    gstart_cell = np.flatnonzero(new_group)
    cum_tp = np.cumsum(tp_blk)
    cum_n = np.cumsum(cnt_blk)
    first_blk = np.searchsorted(bstart, gstart_cell)           # first block of each group
    off_tp = np.r_[0.0, cum_tp][first_blk]
    off_n = np.r_[0.0, cum_n][first_blk]
    grp_of_blk = np.searchsorted(gstart_cell, bstart, side="right") - 1
    TP = cum_tp - off_tp[grp_of_blk]
    N = cum_n - off_n[grp_of_blk]
    P_g = np.bincount(bgroup, weights=tp_blk)
    n_g = np.bincount(bgroup, weights=cnt_blk)
    contrib = tp_blk * TP / N
    ap_g = np.bincount(bgroup, weights=contrib) / np.where(P_g > 0, P_g, 1.0)
    has = P_g > 0
    base = expected_random_ap(n_g[has], P_g[has])
    return {"macro_auprc": float(ap_g[has].mean()) if has.any() else float("nan"),
            "macro_auprc_random_baseline": float(base.mean()) if has.any() else float("nan"),
            "n_groups_with_pos": int(has.sum()), "n_groups": int(len(P_g))}


# ----------------------------------------------------------------------------- precision@k
def precision_at_k(y, s, ks=(100, 1000, 10000)) -> dict:
    """Tie-averaged precision@k (expected precision over random tie-breaking); row-order independent."""
    y, s = check_scores(y, s)
    order = np.argsort(-s, kind="mergesort")
    s_o, y_o = s[order], y[order].astype(np.float64)
    cum = np.r_[0.0, np.cumsum(y_o)]
    out = {}
    for k in ks:
        if k > len(y):
            continue
        sk = s_o[k - 1]
        lo = np.searchsorted(-s_o, -sk, side="left")   # first index of the tie block of the k-th score
        hi = np.searchsorted(-s_o, -sk, side="right")  # one past its last index
        pos_block = cum[hi] - cum[lo]
        out[str(k)] = float((cum[lo] + (k - lo) * pos_block / (hi - lo)) / k)
    return out


# ----------------------------------------------------------------------------- strata
def stratified(y, s, strata: dict) -> dict:
    """strata: name -> boolean mask over cells. AUPRC, prevalence, lift, n, n_pos per stratum."""
    y, s = check_scores(y, s)
    out = {}
    for name, m in strata.items():
        m = np.asarray(m, dtype=bool).ravel()
        yy, ss = y[m], s[m]
        prev = float(yy.mean()) if m.any() else float("nan")
        ap = pooled_auprc(yy, ss) if yy.any() else float("nan")
        out[name] = {"auprc": ap, "prevalence": prev, "lift": ap / prev if prev and prev > 0 else float("nan"),
                     "n": int(m.sum()), "n_pos": int(yy.sum())}
    return out


# ----------------------------------------------------------------------------- weighted AP on a single sort
class SortedAP:
    """Sort once; evaluate AP under per-cluster weights W [n_clusters x B] for B replicates at a time.

    Cells are grouped into segments ending at the last cell of each tie block that contains >= 1 positive; with
    sparse count matrices Cpos/Cneg [n_segments x n_clusters], TP/FP at those thresholds are cumsums of Cpos @ W and
    Cneg @ W, and AP = sum_b (dTP_b / TP_total) * TP_b / (TP_b + FP_b) — exactly sklearn's step-wise AP with ties.
    """

    def __init__(self, y, s, clusters, n_clusters: int | None = None):
        y, s = check_scores(y, s)
        c = np.asarray(clusters).ravel()
        if c.shape != y.shape or not np.issubdtype(c.dtype, np.integer):
            raise ValueError("clusters must be integer codes, one per cell")
        self.n_clusters = int(c.max()) + 1 if n_clusters is None else int(n_clusters)
        order = np.argsort(-s, kind="mergesort")
        s_o, y_o, c_o = s[order], y[order], c[order]
        ends = np.flatnonzero(np.r_[s_o[1:] != s_o[:-1], True])
        blk = np.repeat(np.arange(len(ends)), np.diff(np.r_[-1, ends]))
        pos_blocks = np.unique(blk[y_o])
        seg_end = ends[pos_blocks]                                     # last cell index of each positive block
        seg = np.searchsorted(seg_end, np.arange(len(y_o)), side="left")  # segment of each cell
        keep = seg < len(seg_end)                                      # cells after the last positive block: unused
        nseg = len(seg_end)
        self.Cpos = sp.csr_matrix((np.ones(int((keep & y_o).sum())), (seg[keep & y_o], c_o[keep & y_o])),
                                  shape=(nseg, self.n_clusters))
        self.Cneg = sp.csr_matrix((np.ones(int((keep & ~y_o).sum())), (seg[keep & ~y_o], c_o[keep & ~y_o])),
                                  shape=(nseg, self.n_clusters))

    def ap(self, W: np.ndarray) -> np.ndarray:
        """W: [n_clusters] or [n_clusters x B] non-negative weights -> AP per replicate."""
        W = np.asarray(W, dtype=np.float64)
        one = W.ndim == 1
        W = W[:, None] if one else W
        dtp = np.asarray(self.Cpos @ W)
        dfp = np.asarray(self.Cneg @ W)
        TP = np.cumsum(dtp, axis=0)
        FP = np.cumsum(dfp, axis=0)
        tot = TP[-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            prec = np.where(TP + FP > 0, TP / (TP + FP), 0.0)
            out = np.where(tot > 0, (dtp * prec).sum(0) / tot, np.nan)
        return out[0] if one else out


def bootstrap_weights(pool: np.ndarray, n_clusters: int, n_boot: int, rng: np.random.Generator,
                      batch: int = 100):
    """Yield weight batches W [n_clusters x b]: counts of each cluster in a resample (with replacement, size
    len(pool)) of the cluster codes in `pool`; clusters outside the pool get weight 0."""
    pool = np.asarray(pool)
    done = 0
    while done < n_boot:
        b = min(batch, n_boot - done)
        draws = pool[rng.integers(0, len(pool), size=(b, len(pool)))]
        W = np.zeros((n_clusters, b))
        for j in range(b):
            W[:, j] = np.bincount(draws[j], minlength=n_clusters)
        yield W
        done += b


def bootstrap_p_two_sided(d: np.ndarray) -> float:
    """Two-sided bootstrap p-value for H0: Δ = 0 from replicate differences: min(1, 2 min(P(d <= 0), P(d >= 0)))."""
    d = np.asarray(d, dtype=np.float64)
    d = d[np.isfinite(d)]
    return float(min(1.0, 2.0 * min((d <= 0).mean(), (d >= 0).mean())))


def holm(pvals: dict) -> dict:
    """Holm step-down adjusted p-values: adj_(i) = max_{j<=i} min(1, (m - j + 1) p_(j)) over ascending p."""
    keys = sorted(pvals, key=lambda k: pvals[k])
    m = len(keys)
    out, running = {}, 0.0
    for j, k in enumerate(keys):
        running = max(running, min(1.0, (m - j) * pvals[k]))
        out[k] = running
    return out


def _summ(x: np.ndarray, point: float) -> dict:
    x = x[np.isfinite(x)]
    return {"point": float(point), "boot_mean": float(x.mean()), "boot_sd": float(x.std(ddof=1)),
            "ci95": [float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))], "n_boot": int(len(x))}


def cluster_bootstrap(y, s, clusters, n_boot: int = 1000, seed: int = 0) -> dict:
    """Query-ORF cluster bootstrap of the pooled AUPRC (clusters = integer ORF codes per cell)."""
    sa = SortedAP(y, s, clusters)
    pool = np.unique(np.asarray(clusters))
    rng = np.random.default_rng(seed)
    reps = np.concatenate([sa.ap(W) for W in bootstrap_weights(pool, sa.n_clusters, n_boot, rng)])
    return _summ(reps, float(sa.ap(np.ones(sa.n_clusters))))


def paired_bootstrap(y, s_a, s_b, clusters, n_boot: int = 1000, seed: int = 0) -> dict:
    """Paired cluster bootstrap of ΔAUPRC = AP(s_a) - AP(s_b) with identical ORF weights per replicate."""
    A, B = SortedAP(y, s_a, clusters), SortedAP(y, s_b, clusters)
    n = max(A.n_clusters, B.n_clusters)
    A.n_clusters = B.n_clusters = n
    A.Cpos.resize(A.Cpos.shape[0], n); A.Cneg.resize(A.Cneg.shape[0], n)  # noqa: E702
    B.Cpos.resize(B.Cpos.shape[0], n); B.Cneg.resize(B.Cneg.shape[0], n)  # noqa: E702
    pool = np.unique(np.asarray(clusters))
    rng = np.random.default_rng(seed)
    d, ra, rb = [], [], []
    for W in bootstrap_weights(pool, n, n_boot, rng):
        a, b = A.ap(W), B.ap(W)
        ra.append(a), rb.append(b), d.append(a - b)
    d = np.concatenate(d)
    one = np.ones(n)
    pa, pb = float(A.ap(one)), float(B.ap(one))
    out = _summ(d, pa - pb)
    out.update({"auprc_a": pa, "auprc_b": pb, "relative_gain": (pa - pb) / pb if pb else float("nan"),
                "frac_boot_gt0": float((d > 0).mean()), "ci_excludes_0": bool(out["ci95"][0] > 0 or out["ci95"][1] < 0),
                "p_two_sided": bootstrap_p_two_sided(d)})
    return out


def joint_bootstrap(splits: list[tuple], n_clusters: int, n_boot: int = 1000, seed: int = 0) -> dict:
    """Joint multi-split bootstrap of the mean pooled AUPRC over splits. `splits` = [(y, s, cluster_codes), ...] with
    cluster codes GLOBAL (same ORF -> same code in every split). Each replicate resamples the union of test ORFs
    once and applies the same ORF weights in every split (DESIGN §4/§9)."""
    sas = [SortedAP(y, s, c, n_clusters=n_clusters) for y, s, c in splits]
    pool = np.unique(np.concatenate([np.unique(np.asarray(c)) for _, _, c in splits]))
    rng = np.random.default_rng(seed)
    reps, n_missing = [], 0
    for W in bootstrap_weights(pool, n_clusters, n_boot, rng):
        A = np.stack([sa.ap(W) for sa in sas])
        n_missing += int(np.isnan(A).any(axis=0).sum())   # replicates where a split had no positive weight
        reps.append(np.nanmean(A, axis=0))
    one = np.ones(n_clusters)
    point = float(np.mean([sa.ap(one) for sa in sas]))
    return _summ(np.concatenate(reps), point) | {"n_splits": len(sas), "n_pool_orfs": int(len(pool)),
                                                 "n_reps_with_a_split_undefined": n_missing}


def joint_paired_bootstrap(splits: list[tuple], n_clusters: int, n_boot: int = 1000, seed: int = 0) -> dict:
    """Joint multi-split paired bootstrap of Δ(mean AUPRC): splits = [(y, s_a, s_b, cluster_codes), ...]."""
    A = [SortedAP(y, a, c, n_clusters=n_clusters) for y, a, _, c in splits]
    B = [SortedAP(y, b, c, n_clusters=n_clusters) for y, _, b, c in splits]
    pool = np.unique(np.concatenate([np.unique(np.asarray(c)) for *_, c in splits]))
    rng = np.random.default_rng(seed)
    d = []
    for W in bootstrap_weights(pool, n_clusters, n_boot, rng):
        d.append(np.nanmean(np.stack([a.ap(W) - b.ap(W) for a, b in zip(A, B)]), axis=0))
    one = np.ones(n_clusters)
    pa = float(np.mean([a.ap(one) for a in A]))
    pb = float(np.mean([b.ap(one) for b in B]))
    d = np.concatenate(d)
    out = _summ(d, pa - pb)
    out.update({"auprc_a": pa, "auprc_b": pb, "relative_gain": (pa - pb) / pb if pb else float("nan"),
                "ci_excludes_0": bool(out["ci95"][0] > 0 or out["ci95"][1] < 0), "p_two_sided": bootstrap_p_two_sided(d),
                "n_splits": len(A), "n_pool_orfs": int(len(pool))})
    return out


# ----------------------------------------------------------------------------- convenience
def evaluate(y, s, groups=None, strata: dict | None = None) -> dict:
    """Flat-array metrics (kept for toy use). The harness entry point is ysl.eval.evaluator.Evaluator."""
    y, s = check_scores(y, s)
    out = {"auprc": pooled_auprc(y, s), "auroc": pooled_auroc(y, s), "prevalence": float(y.mean()),
           "n": int(len(y)), "n_pos": int(y.sum()), "precision_at_k": precision_at_k(y, s)}
    out["lift"] = out["auprc"] / out["prevalence"] if out["prevalence"] > 0 else float("nan")
    if groups is not None:
        out.update(macro_auprc(y, s, groups))
    if strata:
        out["strata"] = stratified(y, s, strata)
    return out
