"""Compose <out>/summary.md from the Q3 analysis outputs (any source). Numbers only from the CSV/JSON files.

    .venv/bin/python -m analysis.q3.report --source dev|final|dryrun [--out DIR]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
import yaml

from analysis.q3 import common as C


def f(x, k=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(x):
        return "∞" if x > 0 else "—"
    return f"{x:.{k}f}"


def i0(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "∞" if x == np.inf else ("—" if not np.isfinite(x) else f"{x:,.0f}")


def table(df: pd.DataFrame, cols, heads=None, fmt=None) -> str:
    heads = heads or cols
    fmt = fmt or {}
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(fmt.get(c, lambda v: f(v) if isinstance(v, (float, np.floating)) else str(v))(r[c])
                                     for c in cols) + " |")
    return "\n".join(out)


def section_hap1(O) -> list[str]:
    L = ["## 1. HAP1 translation (§7.3–7.4)", ""]
    if not (O / "hap1_translation.csv").exists():
        return L + ["(not run)", ""]
    s = json.loads((O / "hap1_summary.json").read_text())
    t = pd.read_csv(O / "hap1_translation.csv")
    c = pd.read_csv(O / "hap1_costs.csv")
    b = pd.read_csv(O / "hap1_billmann222.csv")
    m = pd.read_csv(O / "hap1_marginal.csv")
    L += [f"Yeast P* random: n_80 = {i0(s['n_yeast']['n80'])}, n_90 = {i0(s['n_yeast']['n90'])} queries "
          f"(candidate pool ≈ {i0(s['pool_n'])}). Slope b (universe subsampling, re-bootstrapped): "
          f"n_80 {f(s['b']['n80']['b'], 2)} [{f(s['b']['n80']['ci95'][0], 2)}, {f(s['b']['n80']['ci95'][1], 2)}], "
          f"n_90 {f(s['b']['n90']['b'], 2)} [{f(s['b']['n90']['ci95'][0], 2)}, {f(s['b']['n90']['ci95'][1], 2)}].", "",
          "**Random selection, HAP1 lines needed** (N_yeast = 3,283; 95% CI from the joint bootstrap of n_yeast and "
          "b; screens = lines × 1.34):", ""]
    tt = t[t.N_yeast_label.str.contains("primary")].copy()
    tt["lines"] = [f"{i0(a)} [{i0(lo)}, {i0(hi)}]" if bl == "measured" else i0(a)
                   for a, lo, hi, bl in zip(tt.n_HAP1, tt.lo, tt.hi, tt.b_label)]
    tt["screens_"] = tt.screens.map(i0)
    L += [table(tt, ["target", "N_HAP1", "b_label", "lines", "screens_"],
                ["target", "N_HAP1", "scaling", "HAP1 lines", "screens"], {"N_HAP1": i0}), ""]
    tn = t[(t.N_HAP1_label == "mid") & (t.b_label == "measured")].copy()
    tn["lines"] = [f"{i0(a)} [{i0(lo)}, {i0(hi)}]" for a, lo, hi in zip(tn.n_HAP1, tn.lo, tn.hi)]
    L += ["Denominator sensitivity (N_HAP1 = 8,500, measured b):", "",
          table(tn, ["target", "N_yeast_label", "lines"], ["target", "N_yeast", "HAP1 lines"]), ""]
    cc = c[(c.N_HAP1 == 8500) & (c.b_label == "measured")].copy()
    cc["lines_ci"] = [f"{i0(a)} [{i0(lo)}, {i0(hi)}]" for a, lo, hi in zip(cc.lines, cc.lines_lo, cc.lines_hi)]
    cc["saving"] = cc.saving_vs_random.map(lambda v: f"{100 * v:.0f}%" if np.isfinite(v) else "—")
    cc["censored"] = cc.censored.fillna("")
    cc["rounds_"] = [("1 (one-shot)" if a == "random" else i0(r)) for a, r in zip(cc.arm, cc.rounds)]
    L += ["**Costs, random vs strategies** (N_HAP1 = 8,500, measured b). A strategy's budget is the n at which its "
          "own mean curve reaches the random target AUPRC, i.e. its random-equivalent budget n_req equals n_X. "
          "Rounds = sequential screening rounds (batch of 10/50/110 yeast queries = batch × (N_HAP1/N_yeast)^b HAP1 "
          "lines per round). 'right' = target not reached by n = 660; 'left' = reached before the first checkpoint "
          "(value is an upper bound).", "",
          table(cc, ["target", "arm", "n_yeast", "lines_ci", "screens", "rounds_", "lines_per_round", "saving",
                     "censored"],
                ["target", "arm", "yeast n", "HAP1 lines [95% CI]", "screens", "rounds", "lines/round",
                 "saving vs random", "censored"],
                {"n_yeast": i0, "screens": i0, "lines_per_round": i0}), ""]
    bb = b[(b.N_HAP1 == 8500) & b.arm.isin(["R0", "R1", "H", "C", "F_r16", "H-static"])].copy()
    bb["frac"] = [f"{100 * a:.0f}%" + (f" [{100 * lo:.0f}, {100 * hi:.0f}]" if np.isfinite(lo) and lo != hi else "")
                  for a, lo, hi in zip(bb.auprc_frac_of_max, bb.lo, bb.hi)]
    bb["sl"] = bb.sl_detected_frac.map(lambda v: f"{100 * v:.1f}%" if np.isfinite(v) else "—")
    L += ["**Billmann's 222 lines on the curves** (N_HAP1 = 8,500): equivalent yeast budget, AUPRC as % of the "
          "random pool-max AUPRC, strong SLs detected (% of the dev network, cell level).", "",
          table(bb, ["b_label", "arm", "n_yeast_equiv", "frac", "sl"],
                ["scaling", "arm", "yeast n", "AUPRC % of max", "SLs detected"], {"n_yeast_equiv": lambda v: f(v, 0)}),
          ""]
    mm = m[(m.N_HAP1 == 8500) & (m.b_label == "measured") & m.arm.isin(["R0", "H"])].copy()
    mm = mm.pivot_table(index="hap1_lines", columns="arm",
                        values=["auprc_frac_of_max", "d_auprc_frac_per100", "d_sl_detected_frac_per100"]).reset_index()
    mm.columns = ["_".join(x for x in c if x) for c in mm.columns]
    cols = [c for c in mm.columns]
    L += ["**Marginal returns per +100 HAP1 lines** (N_HAP1 = 8,500, measured b; Δ in percentage points of the "
          "random pool-max AUPRC and of the network's strong SLs; H only within its measured range ≤ 660 yeast):", "",
          table(mm, cols, cols, {c: (lambda v: f"{100 * v:.1f}" if np.isfinite(v) else "—") for c in cols
                                 if c != "hap1_lines"} | {"hap1_lines": i0}), ""]
    L += ["![HAP1 translation](fig_hap1_curves.png)", ""]
    return L


def section_frontier(O) -> list[str]:
    L = ["## 2. Discovery frontier (§7.2)", ""]
    if not (O / "frontier.csv").exists():
        return L + ["(not run)", ""]
    a = pd.read_csv(O / "frontier.csv")
    a = a[a.n.isin([160, 330, 660]) & a.arm.isin(["R0", "R1", "H", "C", "F_r16", "H-static", "H__b110"])]
    L += ["Counts recomputed from the stored selections on the dev network (T genes excluded on both sides). "
          "Cell level = strong-SL cells in screened rows; pair level 'detected' = SL called in a screened row, "
          "'touched' = SL pair (either direction) involving a screened gene.", "",
          table(a.sort_values(["n", "arm"]), ["n", "arm", "auprc", "cells_detected_frac", "pairs_detected_frac",
                                               "pairs_touched_frac", "detected_over_touched",
                                               "cells_detected_rel_random", "auprc_rel_random"],
                ["n", "arm", "AUPRC", "cells det.", "pairs det.", "pairs touched", "det./touched",
                 "cells det. ×random", "AUPRC ×random"], {"n": i0}), "", "![frontier](fig_frontier.png)", ""]
    return L


def section_choice(O) -> list[str]:
    L = ["## 3. Choice of query lines (§7.4)", ""]
    if not (O / "choice_profile.csv").exists():
        return L + ["(not run)", ""]
    p = pd.read_csv(O / "choice_profile.csv")
    pv = p.pivot_table(index="feature", columns="arm", values="smd")
    keep = [a for a in ("R0", "H", "C", "F_r16", "H-static", "R1") if a in pv.columns]
    pv = pv[keep].reindex(pv[keep[1] if len(keep) > 1 else keep[0]].abs().sort_values(ascending=False).index)
    pv = pv.reset_index()
    L += ["Standardised mean difference (first 160 picks vs the candidate pool; mean over trajectories; R0 = null "
          "control). SMF is the array strain's single-mutant fitness (HAP1 single-KO fitness analogue); true SL "
          "degree is evaluation-only context.", "", table(pv, ["feature"] + keep), ""]
    if (O / "choice_rules_heldout.csv").exists():
        r = pd.read_csv(O / "choice_rules_heldout.csv")
        L += ["Leave-split-out rank rules (fitted on the finalists' pick frequencies in the other splits, minus "
              "this split's D_val genes):", "",
              table(r, ["eval_split", "rule", "features", "auroc_vs_picked_any", "spearman_vs_pickfreq",
                        "spearman_vs_true_degree"],
                    ["split", "rule", "features", "AUROC vs picked", "ρ vs pick freq.", "ρ vs true degree"]), ""]
        R = json.loads((O / "choice_rules.json").read_text())
        if "all" in R:
            L += ["Rules fitted on all three dev splits (standardised features; the one to carry to HAP1): " +
                  "; ".join(f"k={p['k']}: " + ", ".join(f"{a} {v:+.2f}" for a, v in p["coef"].items())
                            for p in R["all"]["forward"]), ""]
    if (O / "choice_rule_gain.csv").exists():
        g = pd.read_csv(O / "choice_rule_gain.csv")
        L += ["**Static rank rule vs the adaptive strategies** (randomised top-20 static picks from the leave-split-out "
              "rule, P* refit; recovered = (rule − random)/(strategy − random)):", "",
              table(g, ["rule", "n", "auprc_rule", "auprc_random", "auprc_H", "auprc_H_static", "auprc_rel_gain_rule",
                        "auprc_rel_gain_H", "auprc_recovered_vs_H", "auprc_recovered_vs_finalists",
                        "auprc_recovered_H_static_vs_H", "R_seeds"],
                    ["rule", "n", "rule", "random", "H", "H-static", "rule gain", "H gain", "recovered vs H",
                     "recovered vs finalists", "H-static recovered", "seeds"], {"n": i0}), ""]
    L += ["![profile](fig_choice_profile.png)", ""]
    return L


def section_pred(O) -> list[str]:
    L = ["## 4. What prediction can contribute (§7.4)", ""]
    if not (O / "pred_vs_reference.csv").exists():
        return L + ["(not run)", ""]
    r = pd.read_csv(O / "pred_reciprocal.csv")
    v = pd.read_csv(O / "pred_vs_reference.csv")
    k = pd.read_csv(O / "pred_topk.csv")
    L += ["Replicate-like reference on the on-array target stratum: the reciprocal measurement of the same pair "
          "(BA → AB). Both directions are evaluation cells, so its precision equals its recall.", "",
          table(r, ["split", "n_cells", "n_pos", "prevalence", "ref_precision", "ref_recall", "kappa"],
                fmt={"n_cells": i0, "n_pos": i0}), "",
          "P* on the same cells (mean over trajectories and splits). HAP1 context: Billmann single screen vs "
          "replicates, precision 0.61 at recall 0.37.", "",
          table(v, ["arm", "n", "auprc_paired", "prec_at_ref_recall", "recall_at_ref_precision",
                    "prec_at_hap1_recall", "recall_at_hap1_precision"],
                ["arm", "n", "AUPRC (paired cells)", "precision @ ref recall", "recall @ ref precision",
                 "precision @ recall 0.37", "recall @ precision 0.61"], {"n": i0}), ""]
    kk = k[k.k.isin([10, 50]) & k.arm.isin(["R0", "H"])]
    L += ["Expected true strong SLs among the top-k predictions per unscreened gene:", "",
          table(kk, ["set", "arm", "n", "k", "true_sl_per_gene", "random_expect", "precision", "frac_genes_ge1"],
                ["columns", "arm", "n", "k", "true SLs / gene", "random expectation", "precision",
                 "genes with ≥1"], {"n": i0, "k": i0}), ""]
    if (O / "pred_thresholds.csv").exists() and (O / "pred_thresholds.csv").stat().st_size > 5:
        t = pd.read_csv(O / "pred_thresholds.csv")
        if len(t):
            tm = t.groupby(["arm", "n", "target", "level"]).mean(numeric_only=True).reset_index()
            L += ["Thresholds chosen by leave-batch-out CV on screened rows, applied to the targets (mean over the "
                  "LBO trajectories):", "",
                  table(tm, ["arm", "n", "target", "level", "cv_precision", "cv_recall", "target_precision",
                             "target_recall", "oa_precision", "oa_recall", "target_pred_per_gene"],
                        ["arm", "n", "set by", "level", "CV precision", "CV recall", "target precision",
                         "target recall", "on-array precision", "on-array recall", "calls / gene"],
                        {"n": i0}), ""]
    L += ["![PR](fig_pred_pr.png)", ""]
    return L


def section_stop(O) -> list[str]:
    L = ["## 5. Online stopping estimator (§7.4)", ""]
    if not (O / "stop_estimator.csv").exists():
        return L + ["(not run)", ""]
    e = pd.read_csv(O / "stop_estimator.csv")
    L += ["LBO-CV AUPRC on screened rows vs the trajectory's held-out AUPRC (E_fixed), per n:", "",
          table(e[e.arm == "all"], ["n", "R", "cv_mean", "true_mean", "bias", "mae", "rel_bias", "rel_mae", "pearson",
                                    "spearman", "pearson_split_centred", "fit_s_mean"], fmt={"n": i0, "R": i0}), "",
          table(e[e.arm != "all"], ["arm", "n", "R", "cv_mean", "true_mean", "rel_bias", "rel_mae", "pearson"],
                fmt={"n": i0, "R": i0}), ""]
    if (O / "stop_rule_grid.csv").exists() and (O / "stop_rule_grid.csv").stat().st_size > 5:
        g = pd.read_csv(O / "stop_rule_grid.csv").head(8)
        L += ["Plateau-rule grid (best first; mean |log2(n_stop / n_knee)|):", "",
              table(g, ["delta", "m", "mae_log2", "bias_log2", "trig", "stop_frac"]), ""]
    if (O / "stop_knees.csv").exists() and (O / "stop_knees.csv").stat().st_size > 5:
        kn = pd.read_csv(O / "stop_knees.csv")
        L += ["True knees per trajectory (0.95 × AUPRC(660)):", "", table(kn, list(kn.columns)), ""]
    if (O / "stop_rule_lso.csv").exists() and (O / "stop_rule_lso.csv").stat().st_size > 5:
        L += ["Leave-split-out check of the rule choice:", "", table(pd.read_csv(O / "stop_rule_lso.csv"),
                                                                    ["held_out", "delta", "m", "mae_log2", "bias_log2",
                                                                     "R"]), ""]
    rf = C.Q3 / "stopping_rule.yaml"
    if rf.exists():
        r = yaml.safe_load(rf.read_text())
        L += [f"**Frozen rule** (`analysis/q3/stopping_rule.yaml`, frozen {r['frozen_at']}): stop at the first "
              f"checkpoint n ≥ {r['rule']['n_min']} whose last {r['rule']['m']} relative CV gain(s) per doubling are "
              f"< {r['rule']['delta']}; otherwise stop at 660 (flagged).", ""]
    if (O / "stop_rule_applied.csv").exists():
        a = pd.read_csv(O / "stop_rule_applied.csv")
        L += ["Frozen rule applied:", "", table(a, list(a.columns)), ""]
    L += ["![stopping](fig_stopping.png)", ""]
    return L


def compose(source: str, out=None) -> str:
    O = C.out_dir(source, out)
    head = [f"# Q3 analyses — {source.upper()} {'(DEV DATA ONLY)' if source == 'dev' else ''}", "",
            "DESIGN v3.1 §7.2–7.4. Code: `analysis/q3/` (one command per analysis, `--source dev|final|dryrun`); "
            "rerun everything with `bash analysis/q3/run_all.sh <source>`. Predictor: frozen P* "
            "(stack_lr_wr2_lb). All numbers below are generated from the CSV/JSON files in this directory.", ""]
    body = section_hap1(O) + section_frontier(O) + section_choice(O) + section_pred(O) + section_stop(O)
    txt = "\n".join(head + body)
    notes = O / "summary_notes.md"
    if notes.exists():
        txt = txt + "\n" + notes.read_text()
    (O / "summary.md").write_text(txt)
    return txt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    compose(a.source, a.out)
    print("wrote summary.md")


if __name__ == "__main__":
    main()
