"""Q3 analysis 4 — what prediction can contribute (DESIGN v3.1 §7.4).

    .venv/bin/python -m analysis.q3.predcontrib --source dev|final|dryrun [--out DIR]

(a) Replicate-like reference: on the on-array target stratum (dev E_VV / final E_TT: target x on-array target
    columns) each cell (t, a) has a reciprocal measurement (a, t) of the same pair in the other direction (a
    different strain pair, the closest within-yeast replicate). Reference 'predictor' = the reciprocal label:
    precision = P(AB SL | BA SL), recall = P(BA SL | AB SL), on cells where both directions are measured.
    P* on the same cells: precision at the reference's recall, recall at the reference's precision, AUPRC.
    HAP1 context: Billmann per-screen precision 0.61 / recall 0.37 against replicate screens.
(b) Score thresholds chosen by leave-batch-out CV on SCREENED rows (analysis.q3.lbo cells; never targets):
    threshold where CV recall first reaches r (r in recall_targets) or the max-recall threshold with CV precision
    >= p; then applied to the same trajectory/budget's P* scores on the targets: realised precision / recall.
(c) Expected true SLs among the top-k predictions per unscreened gene (E_fixed and on-array columns), vs the
    random expectation (row prevalence x k).
Outputs: pred_reciprocal.csv, pred_vs_reference.csv, pred_topk.csv, pred_thresholds.csv, pred_summary.json,
fig_pred_pr.png
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from analysis.q3 import common as C

BUDGETS = (160, 330, 660)
ARMS = ("R0", "H", "C", "F_r16")
KS = (5, 10, 20, 50, 100)
RECALL_TARGETS = (0.05, 0.10, 0.20)
PRECISION_TARGETS = (0.10, 0.20)
MIN_PRED_POS = 20


# ============================================================================ PR helpers (tie-grouped)
def pr_points(y, s):
    """(threshold, precision, recall, n_pred) at every tie-group end, scores descending."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=np.float64)
    o = np.argsort(-s, kind="mergesort")
    ss, yy = s[o], y[o]
    ends = np.r_[np.flatnonzero(ss[1:] != ss[:-1]), len(ss) - 1]
    tp = np.cumsum(yy)[ends]
    npred = ends + 1
    P = max(int(y.sum()), 1)
    return ss[ends], tp / npred, tp / P, npred


def precision_at_recall(y, s, r):
    th, p, rc, _ = pr_points(y, s)
    k = np.flatnonzero(rc >= r)
    return float(p[k[0]]) if len(k) else float("nan")


def recall_at_precision(y, s, p_min, min_pred=1):
    th, p, rc, npred = pr_points(y, s)
    ok = (p >= p_min) & (npred >= min_pred)
    return float(rc[ok].max()) if ok.any() else 0.0


def threshold_for_recall(y, s, r):
    th, p, rc, _ = pr_points(y, s)
    k = np.flatnonzero(rc >= r)
    return float(th[k[0]]) if len(k) else float("nan")


def threshold_for_precision(y, s, p_min, min_pred=MIN_PRED_POS):
    th, p, rc, npred = pr_points(y, s)
    ok = np.flatnonzero((p >= p_min) & (npred >= min_pred))
    return float(th[ok[-1]]) if len(ok) else float("nan")


def apply_threshold(y, s, tau):
    if not np.isfinite(tau):
        return float("nan"), float("nan"), 0
    m = np.asarray(s) >= tau
    tp = int(np.asarray(y)[m].sum())
    return (tp / m.sum() if m.sum() else float("nan")), tp / max(int(np.sum(y)), 1), int(m.sum())


# ============================================================================ (a) reciprocal pairs
def reciprocal_index(L: dict) -> tuple[np.ndarray, np.ndarray]:
    """For cells of the on-array target stratum: indices i (cell (t, a)) and j (its reciprocal (a, t)) where both
    are evaluation cells."""
    targets = np.asarray(L["targets"], dtype=str)
    arrays = np.asarray(L["array_orfs"], dtype=str)
    ncol = len(arrays)
    key = L["rows"].astype(np.int64) * ncol + L["cols"]
    order = np.argsort(key)
    skey = key[order]
    t_row = {g: i for i, g in enumerate(targets)}
    a_col = {g: j for j, g in enumerate(arrays)}
    idx = np.flatnonzero(L["OA"])
    t_names = targets[L["rows"][idx]]
    a_names = arrays[L["cols"][idx]]
    rr = np.array([t_row.get(a, -1) for a in a_names])
    cc = np.array([a_col.get(t, -1) for t in t_names])
    ok = (rr >= 0) & (cc >= 0)
    rk = rr.astype(np.int64) * ncol + cc
    pos = np.searchsorted(skey, rk)
    pos = np.clip(pos, 0, len(skey) - 1)
    ok &= skey[pos] == rk
    return idx[ok], order[pos[ok]]


def reference(L) -> dict:
    i, j = reciprocal_index(L)
    y, yr = L["y"][i].astype(bool), L["y"][j].astype(bool)
    tp = int((y & yr).sum())
    po = float((y == yr).mean())
    pe = float(y.mean() * yr.mean() + (1 - y.mean()) * (1 - yr.mean()))
    return {"n_cells": int(len(i)), "n_pos": int(y.sum()), "prevalence": float(y.mean()),
            "ref_precision": tp / max(int(yr.sum()), 1), "ref_recall": tp / max(int(y.sum()), 1),
            "kappa": (po - pe) / (1 - pe) if pe < 1 else float("nan"), "_i": i, "_j": j}


# ============================================================================ (c) top-k
def topk(y, s, rows, ks=KS) -> dict:
    o = np.lexsort((-np.asarray(s, np.float64), rows))
    r_sorted = rows[o]
    first = np.r_[0, np.flatnonzero(np.diff(r_sorted)) + 1]
    start = np.repeat(first, np.diff(np.r_[first, len(o)]))
    rank = np.arange(len(o)) - start
    yy = y[o].astype(float)
    n_rows = len(first)
    cnt = np.bincount(r_sorted, minlength=rows.max() + 1)
    pos = np.bincount(r_sorted, weights=yy, minlength=rows.max() + 1)
    present = cnt > 0
    prev = pos[present] / cnt[present]
    out = {}
    for k in ks:
        m = rank < k
        hits = np.bincount(r_sorted[m], weights=yy[m], minlength=rows.max() + 1)[present]
        kk = np.minimum(k, cnt[present])
        out[k] = {"true_sl_per_gene": float(hits.mean()), "random_expect": float((prev * kk).mean()),
                  "precision": float(hits.sum() / kk.sum()), "frac_genes_ge1": float((hits >= 1).mean()),
                  "n_genes": int(n_rows)}
    return out


# ============================================================================ main
def analyse(source: str, out=None) -> dict:
    from hsl import eval as E
    t0 = time.perf_counter()
    O = C.out_dir(source, out)
    K = C.hap1_constants()["billmann"]
    rows = C.rows(source)
    rows = rows[rows.arm.isin(ARMS) & rows.n.isin(BUDGETS)]
    refs, vs, tk = [], [], []
    for sp, gsp in rows.groupby("split"):
        L = C.eval_labels(source, sp)
        ref = reference(L)
        i, j = ref.pop("_i"), ref.pop("_j")
        refs.append({"split": sp, **ref})
        oa = L["OA"]
        for _, x in gsp.iterrows():
            s = C.load_pred(x.pred_file)
            if len(s) != len(L["y"]):
                raise AssertionError(f"{x.pred_file}: cells != labels")
            yi, si = L["y"][i], s[i]
            vs.append({"split": sp, "arm": x.arm, "traj": x.traj, "n": x.n,
                       "auprc_paired": E.average_precision(yi, si),
                       "prec_at_ref_recall": precision_at_recall(yi, si, ref["ref_recall"]),
                       "recall_at_ref_precision": recall_at_precision(yi, si, ref["ref_precision"], MIN_PRED_POS),
                       "prec_at_hap1_recall": precision_at_recall(yi, si, K["per_screen_recall"]),
                       "recall_at_hap1_precision": recall_at_precision(yi, si, K["per_screen_precision"],
                                                                       MIN_PRED_POS),
                       "ref_precision": ref["ref_precision"], "ref_recall": ref["ref_recall"],
                       "prevalence_paired": ref["prevalence"]})
            for setn, m in (("E_fixed", np.ones(len(s), bool)), ("on_array_cols", oa)):
                for k, v in topk(L["y"][m], s[m], L["rows"][m]).items():
                    tk.append({"split": sp, "arm": x.arm, "traj": x.traj, "n": x.n, "set": setn, "k": k, **v})
    refs, vs, tk = pd.DataFrame(refs), pd.DataFrame(vs), pd.DataFrame(tk)
    refs.to_csv(O / "pred_reciprocal.csv", index=False)
    vsm = vs.groupby(["arm", "n", "split"]).mean(numeric_only=True).groupby(["arm", "n"]).mean().reset_index()
    vsm.to_csv(O / "pred_vs_reference.csv", index=False)
    tkm = (tk.groupby(["arm", "n", "set", "k", "split"]).mean(numeric_only=True)
           .groupby(["arm", "n", "set", "k"]).mean().reset_index())
    tkm.to_csv(O / "pred_topk.csv", index=False)
    th = thresholds(source, O, rows)
    summ = {"source": source, "reference": refs.drop(columns=[], errors="ignore").to_dict("records"),
            "hap1_context": {"precision": K["per_screen_precision"], "recall": K["per_screen_recall"]},
            "n_threshold_rows": int(len(th)), "runtime_s": time.perf_counter() - t0}
    C.write_json(O / "pred_summary.json", summ)
    fig_pr(source, O, rows, refs)
    return summ


def thresholds(source: str, O, rows) -> pd.DataFrame:
    """(b): thresholds from LBO-CV cells of screened rows, applied to the same trajectory/budget's target scores."""
    p = O / "lbo" / "lbo_rows.jsonl"
    if not p.exists():
        pd.DataFrame().to_csv(O / "pred_thresholds.csv", index=False)
        return pd.DataFrame()
    lbo = pd.DataFrame([json.loads(x) for x in p.read_text().splitlines() if x.strip()])
    recs = []
    labels = {}
    for _, z in lbo.iterrows():
        tgt = rows[(rows.arm == z.arm) & (rows.split == z.split) & (rows.traj == z.traj) & (rows.n == z.n)]
        if tgt.empty:
            continue
        with np.load(C.ROOT / z.cells_file) as cz:
            ycv, scv = cz["y"], cz["s"]
        L = labels.setdefault(z.split, C.eval_labels(source, z.split))
        s = C.load_pred(tgt.pred_file.iloc[0])
        for kind, lvl, tau in ([("recall", r, threshold_for_recall(ycv, scv, r)) for r in RECALL_TARGETS] +
                               [("precision", q, threshold_for_precision(ycv, scv, q)) for q in PRECISION_TARGETS]):
            pc, rc, nc = apply_threshold(ycv, scv, tau)
            pt, rt, nt = apply_threshold(L["y"], s, tau)
            po, ro, no = apply_threshold(L["y"][L["OA"]], s[L["OA"]], tau)
            recs.append({"arm": z.arm, "split": z.split, "traj": z.traj, "n": z.n, "target": kind, "level": lvl,
                         "tau": tau, "cv_precision": pc, "cv_recall": rc, "target_precision": pt,
                         "target_recall": rt, "target_pred_per_gene": nt / len(L["targets"]),
                         "oa_precision": po, "oa_recall": ro, "cv_prevalence": float(ycv.mean()),
                         "target_prevalence": float(L["y"].mean())})
    T = pd.DataFrame(recs)
    T.to_csv(O / "pred_thresholds.csv", index=False)
    return T


def fig_pr(source, O, rows, refs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})
    K = C.hap1_constants()["billmann"]
    sp = sorted(rows.split.unique())[0]
    L = C.eval_labels(source, sp)
    ref = reference(L)
    i = ref["_i"]
    f, ax = plt.subplots(figsize=(5.5, 4.2), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color="#e8e7e2", lw=0.6)
    style = {("R0", 160): ("#9a9992", "random, n=160"), ("R0", 660): ("#52514e", "random, n=660"),
             ("H", 160): ("#86b6ef", "H, n=160"), ("H", 660): ("#2a78d6", "H, n=660")}
    for (a, n), (col, lab) in style.items():
        x = rows[(rows.arm == a) & (rows.n == n) & (rows.split == sp)]
        if x.empty:
            continue
        s = C.load_pred(x.pred_file.iloc[0])
        _, p, r, _ = pr_points(L["y"][i], s[i])
        ax.plot(r, p, color=col, lw=2, label=lab)
    ax.plot([ref["ref_recall"]], [ref["ref_precision"]], "s", color="#eb6834", ms=8,
            label="reciprocal measurement (BA→AB)")
    ax.plot([K["per_screen_recall"]], [K["per_screen_precision"]], "D", color="#1baf7a", ms=8,
            label="HAP1 single screen vs replicates")
    ax.axhline(ref["prevalence"], color="#c3c2b7", lw=0.8, ls=":")
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_ylim(0, 0.7)
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title(f"P* on on-array target pairs vs replicate-like references ({source}, {sp}, trajectory 0)",
                 fontsize=8)
    f.tight_layout()
    f.savefig(O / "fig_pred_pr.png", dpi=150, facecolor="#fcfcfb")
    plt.close(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    s = analyse(a.source, a.out)
    print(json.dumps({k: v for k, v in s.items() if k != "reference"}, default=str))
    print(pd.DataFrame(s["reference"]).round(4).to_string())


if __name__ == "__main__":
    main()
