"""Dev vs final side by side for the Q3 headline numbers (read from results/q3_dev and results/q3_final outputs).

    .venv/bin/python -m analysis.q3.compare        -> results/q3_final/dev_vs_final.csv and .md
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from analysis.q3 import common as C


def _rd(p):
    try:
        d = pd.read_csv(p)
        return d if len(d) else None
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def metrics(O) -> dict:
    m = {}
    t = _rd(O / "hap1_translation.csv")
    if t is not None:
        for st in ("OA", "all"):
            x = t[(t.stratum == st) & t.N_yeast_label.str.contains("primary") & (t.N_HAP1_label == "mid")]
            for tg in ("n80", "n90"):
                r = x[(x.target == tg) & (x.b_label == "measured")].iloc[0]
                m[f"random {tg} yeast queries ({st})"] = r.n_yeast
                m[f"random {tg} HAP1 lines ({st}, N_HAP1 8,500, b measured)"] = (
                    f"{r.n_HAP1:,.0f} [{r.lo:,.0f}, {r.hi:,.0f}]")
                m[f"random {tg} HAP1 lines ({st}, b = 0 bound)"] = x[(x.target == tg) &
                                                                     (x.b_label == "b=0 (count)")].n_HAP1.iloc[0]
                m[f"random {tg} HAP1 lines ({st}, b = 1 bound)"] = x[(x.target == tg) &
                                                                     (x.b_label == "b=1 (fraction)")].n_HAP1.iloc[0]
    c = _rd(O / "hap1_costs.csv")
    if c is not None:
        c = c[(c.stratum == "OA") & (c.N_HAP1 == 8500) & (c.b_label == "measured")]
        for arm in ("H", "C", "F_r16", "H__b110", "C__b110", "H-static", "R1"):
            for tg in ("n80", "n90"):
                r = c[(c.arm == arm) & (c.target == tg)]
                if len(r):
                    r = r.iloc[0]
                    v = (f"{r.lines:,.0f} lines, {r.rounds:.0f} rounds, saving {100 * r.saving_vs_random:.0f}%"
                         if np.isfinite(r.lines) else f"not reached by 660 ({r.censored})")
                    m[f"{arm} budget for random {tg} AUPRC (OA)"] = v
    q = _rd(O / "hap1_nreq.csv")
    if q is not None:
        q = q[q.N_HAP1 == 8500]
        for arm in C.FINALISTS:
            for n in (160, 660):
                r = q[(q.arm == arm) & (q.n_yeast == n)]
                if len(r):
                    r = r.iloc[0]
                    m[f"{arm} n_req at n = {n} (E_fixed; ratio)"] = (
                        f"{r.n_req_yeast:,.0f} ({r.ratio:.2f}×)" if np.isfinite(r.n_req_yeast) else
                        "> random pool max (right-censored)")
    b = _rd(O / "hap1_billmann222.csv")
    if b is not None:
        b = b[(b.stratum == "OA") & (b.N_HAP1 == 8500) & (b.b_label == "measured")].set_index("arm")
        for arm in ("R0", "H", "R1"):
            if arm in b.index:
                m[f"Billmann 222 lines: {arm} AUPRC % of pool-max random (OA)"] = 100 * b.loc[arm, "auprc_frac_of_max"]
    f = _rd(O / "frontier.csv")
    if f is not None:
        for arm in ("H", "C", "F_r16", "R1"):
            for n in (160, 660):
                r = f[(f.arm == arm) & (f.n == n)]
                if len(r):
                    r = r.iloc[0]
                    m[f"{arm} n = {n}: SL cells detected ×random / AUPRC ×random"] = (
                        f"{r.cells_detected_rel_random:.2f} / {r.auprc_rel_random:.2f}")
    p = _rd(O / "choice_profile.csv")
    if p is not None:
        for feat in ("smf", "conservation", "expr_sd", "abundance", "apms_degree"):
            r = p[(p.arm == "H") & (p.feature == feat)]
            if len(r):
                m[f"H early picks: SMD {feat}"] = r.smd.iloc[0]
    g = _rd(O / "choice_rule_gain.csv")
    if g is not None:
        for rule in ("forward_k1", "forward_k3"):
            for n in (160, 330, 660):
                r = g[(g.rule == rule) & (g.n == n)]
                if len(r):
                    r = r.iloc[0]
                    m[f"static rule {rule} n = {n}: gain vs random / recovered of H's gain"] = (
                        f"{100 * r.auprc_rel_gain_rule:+.1f}% / {r.auprc_recovered_vs_H:.2f} "
                        f"(splits {int(r.n_splits)}, fits {int(r.n_fits)})")
    rc = _rd(O / "pred_reciprocal.csv")
    if rc is not None:
        m["reciprocal-reference precision = recall (on-array stratum)"] = rc.ref_precision.mean()
    v = _rd(O / "pred_vs_reference.csv")
    if v is not None:
        for arm in ("R0", "H"):
            for n in (160, 660):
                r = v[(v.arm == arm) & (v.n == n)]
                if len(r):
                    r = r.iloc[0]
                    m[f"P* {arm} n = {n}: precision at reference recall"] = r.prec_at_ref_recall
                    m[f"P* {arm} n = {n}: precision at recall 0.37 / recall at precision 0.61"] = (
                        f"{r.prec_at_hap1_recall:.3f} / {r.recall_at_hap1_precision:.3f}")
    k = _rd(O / "pred_topk.csv")
    if k is not None:
        for arm in ("R0", "H"):
            r = k[(k.arm == arm) & (k.n == 660) & (k.k == 10) & (k.set == "E_fixed")]
            if len(r):
                m[f"true SLs in top-10 per unscreened gene, {arm} n = 660"] = r.true_sl_per_gene.iloc[0]
    th = _rd(O / "pred_thresholds.csv")
    if th is not None:
        for arm in ("R0", "H"):
            r = th[(th.arm == arm) & (th.n == 660) & (th.target == "recall") & (th.level == 0.1)]
            if len(r):
                m[f"CV threshold at CV recall 0.10, {arm} n = 660: CV precision → target precision"] = (
                    f"{r.cv_precision.mean():.3f} → {r.target_precision.mean():.3f}")
    e = _rd(O / "stop_estimator.csv")
    if e is not None:
        for arm in ("R0", "H"):
            r = e[e.arm == arm]
            if len(r):
                m[f"LBO-CV relative bias range, {arm}"] = f"{100 * r.rel_bias.min():+.0f}% to {100 * r.rel_bias.max():+.0f}%"
    a = _rd(O / "stop_rule_applied.csv")
    if a is not None:
        m["frozen stopping rule: mean |log2(n_stop / knee)|"] = a.log2_err.abs().mean()
        m["frozen stopping rule: triggered share / AUPRC at stop ÷ AUPRC(660)"] = (
            f"{a.triggered.mean():.2f} / {a.true_at_stop_frac.mean():.3f}")
        m["frozen stopping rule: trajectories"] = len(a)
        kn = _rd(O / "stop_knees.csv")
        if kn is not None:
            m["true knee range (n)"] = f"{kn.knee.min():.0f}–{kn.knee.max():.0f}"
            m["fixed n = 490 comparator: mean |log2(490 / knee)|"] = np.abs(np.log2(490 / kn.knee)).mean()
    return m


def main():
    D = metrics(C.ROOT / "results" / "q3_dev")
    F = metrics(C.ROOT / "results" / "q3_final")
    keys = list(dict.fromkeys(list(D) + list(F)))

    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, (float, np.floating)):
            return f"{v:,.0f}" if abs(v) >= 100 else f"{v:.3f}"
        if isinstance(v, (int, np.integer)):
            return f"{v:,}"
        return str(v)
    rows = [{"metric": k, "dev": fmt(D.get(k)), "final_T": fmt(F.get(k))} for k in keys]
    df = pd.DataFrame(rows)
    out = C.out_dir("final")
    df.to_csv(out / "dev_vs_final.csv", index=False)
    md = ["| metric | dev (D_val, dev1–3) | final (T) |", "|---|---|---|"] + [
        f"| {r['metric']} | {r['dev']} | {r['final_T']} |" for r in rows]
    (out / "dev_vs_final.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
