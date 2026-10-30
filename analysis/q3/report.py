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
    keep = [k for k, c in enumerate(cols) if c in df.columns]          # tolerate columns a source does not have
    cols, heads = [cols[k] for k in keep], [heads[k] for k in keep]
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
    ny = s["n_yeast"]
    L += ["Yeast P* random n_80 / n_90 by AUPRC stratum (pool maximum = 100%): " +
          "; ".join(f"{st} {i0(v['n80'])} / {i0(v['n90'])}" for st, v in ny.items()) +
          f" (candidate pool ≈ {i0(s['pool_n'])}). DESIGN §7.4 uses the on-array stratum (OA = E_VV on dev, E_TT at "
          "the final) for the HAP1 translation: OA is the primary below; 'all' (E_fixed, the Q1 endpoint with its "
          "frozen joint-bootstrap CI) and on-array rows are sensitivities. OA/onrows CIs carry over the 'all' "
          "stratum's relative CI width (approximation). Slope b (universe subsampling, re-bootstrapped): "
          f"n_80 {f(s['b']['n80']['b'], 2)} [{f(s['b']['n80']['ci95'][0], 2)}, {f(s['b']['n80']['ci95'][1], 2)}], "
          f"n_90 {f(s['b']['n90']['b'], 2)} [{f(s['b']['n90']['ci95'][0], 2)}, {f(s['b']['n90']['ci95'][1], 2)}].", "",
          "**Random selection, HAP1 lines needed** (N_yeast = 3,283; 95% CI from the joint bootstrap of n_yeast and "
          "b; screens = lines × 1.34):", ""]
    sens = t[t.N_yeast_label.str.contains("primary") & (t.N_HAP1_label == "mid") & (t.b_label == "measured")].copy()
    sens["lines"] = [f"{i0(a)} [{i0(lo)}, {i0(hi)}]" for a, lo, hi in zip(sens.n_HAP1, sens.lo, sens.hi)]
    t = t[t.stratum == "OA"]
    c = c[c.stratum == "OA"]
    b = b[b.stratum == "OA"]
    m = m[m.stratum == "OA"]
    tt = t[t.N_yeast_label.str.contains("primary")].copy()
    tt["lines"] = [f"{i0(a)} [{i0(lo)}, {i0(hi)}]" if bl == "measured" else i0(a)
                   for a, lo, hi, bl in zip(tt.n_HAP1, tt.lo, tt.hi, tt.b_label)]
    tt["screens_"] = tt.screens.map(i0)
    L += [table(tt, ["target", "N_HAP1", "b_label", "lines", "screens_"],
                ["target", "N_HAP1", "scaling", "HAP1 lines", "screens"], {"N_HAP1": i0}), ""]
    L += ["Stratum sensitivity (N_yeast 3,283, N_HAP1 8,500, measured b):", "",
          table(sens, ["stratum", "target", "n_yeast", "lines"], ["AUPRC stratum", "target", "yeast n", "HAP1 lines"],
                {"n_yeast": i0}), ""]
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
    L += ["![HAP1 translation (OA stratum)](fig_hap1_curves.png) (E_fixed version: fig_hap1_curves_all.png)", ""]
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
                        "auprc_recovered_H_static_vs_H", "n_splits", "n_fits", "offarray_share"],
                    ["rule", "n", "rule", "random", "H", "H-static", "rule gain", "H gain", "recovered vs H",
                     "recovered vs finalists", "H-static recovered", "splits", "fits", "off-array share"],
                    {"n": i0, "n_splits": i0, "n_fits": i0}), "",
              "n = 660 is partial: those fits were stopped by the 3 GB RSS watchdog (P* at n = 660 on the full dev "
              "target set peaks at ≈3.0–3.2 GB), so only the splits listed were completed.", ""]
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
                                    "spearman", "pearson_split_centred", "lift_rel_bias", "pearson_lift", "fit_s_mean"],
                fmt={"n": i0, "R": i0}), "",
          table(e[e.arm != "all"], ["arm", "n", "R", "cv_mean", "true_mean", "rel_bias", "rel_mae", "pearson",
                                    "lift_rel_bias", "lift_rel_mae", "pearson_lift"],
                fmt={"n": i0, "R": i0}), ""]
    if (O / "stop_rule_grid.csv").exists() and (O / "stop_rule_grid.csv").stat().st_size > 5:
        g = pd.read_csv(O / "stop_rule_grid.csv").head(8)
        L += ["Plateau-rule grid (best first; mean |log2(n_stop / n_knee)|; stop_frac = AUPRC at stop / AUPRC(660)):",
              "", table(g, ["family", "statistic", "delta", "m", "mae_log2", "bias_log2", "trig", "stop_frac"],
                        ["family", "statistic", "δ (plateau) / τ (projection)", "m / window", "mean |log2 err|",
                         "mean log2 err", "triggered", "AUPRC at stop / AUPRC(660)"]), ""]
    if (O / "stop_fixed_n_baseline.csv").exists() and (O / "stop_fixed_n_baseline.csv").stat().st_size > 5:
        L += ["Comparator without online information: always stop at a fixed n (mean |log2(n / knee)|):", "",
              table(pd.read_csv(O / "stop_fixed_n_baseline.csv"), ["n_fixed", "mae_log2", "bias_log2"],
                    fmt={"n_fixed": i0}), ""]
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
        R_ = r["rule"]
        how = (f"stop at the first checkpoint n ≥ {R_['n_min']} whose last {R_['m_or_window']} relative gains per "
               f"doubling of {R_['statistic']} are all < {R_['param']}" if R_["family"] == "plateau" else
               f"stop at the first checkpoint n ≥ {R_['n_min']} where {R_['statistic']} ≥ {R_['param']} × its "
               f"log-linear projection (last {R_['m_or_window']} checkpoints) to the planned maximum "
               f"n = {R_['n_max_planned']}")
        L += [f"**Frozen rule** (`analysis/q3/stopping_rule.yaml`, frozen {r['frozen_at']}): {R_['family']} family; "
              f"{how}; otherwise stop at 660 (flagged).", ""]
    if (O / "stop_rule_applied.csv").exists() and (O / "stop_rule_applied.csv").stat().st_size > 5:
        a = pd.read_csv(O / "stop_rule_applied.csv")
        L += ["Frozen rule applied:", "", table(a, list(a.columns)), ""]
    L += ["![stopping](fig_stopping.png)", ""]
    return L


NOTES = """## Ambiguities in DESIGN §7.2–7.4 and how they were resolved (none silently)

1. **HAP1 constants are not in protocol.yaml.** DESIGN §5 and the brief place them there; the frozen file has no HAP1
   block. They are in `analysis/q3/hap1_constants.yaml` with their LITERATURE.md sources (N_HAP1 7,500 / 8,500 / 9,500;
   replicate factor 298/222 = 1.342; Billmann precision 0.61 / recall 0.37). Q3 seeds live there too
   (`config/seeds.yaml` is frozen).
2. **Which AUPRC stratum drives the translation.** §7.4 names the on-array stratum as the HAP1 estimate; §7.3/Q1 use
   E_fixed. Primary = on-array target columns (OA: E_VV on dev, E_TT at the final); E_fixed ('all') and on-array
   target rows ('onrows') are sensitivities.
3. **N_yeast vs the curve's anchor.** n_80 is anchored on the AUPRC at the candidate-pool maximum (dev ≈ 2,100;
   final 2,627), while §7.4 fixes N_yeast = 3,283. 3,283 is primary; 4,458 (Q ∪ A, DESIGN sensitivity) and
   N_yeast = pool size (added sensitivity) are reported. With b ≈ 1 the pool-size version is ~1.6× larger: this is
   the largest single lever on the HAP1 number.
4. **Bootstrap of n_yeast.** Only the 95% CIs of the frozen joint bootstrap are stored, not its replicates: n_yeast
   replicates are lognormal matched to each CI (per split, averaged over splits on dev). b replicates re-run P1's
   within-level draw bootstrap exactly (point b reproduces slope.json to 1e-9). OA/onrows strata carry over the
   'all' stratum's relative CI width (approximation).
5. **'Use the strategy's n_req'.** The strategy's cost for target X is the budget n_s at which its own mean (PAVA)
   curve reaches the random target AUPRC, so n_req(n_s) = n_X. Reached before the first checkpoint → first
   checkpoint reported as an upper bound ('left'); not reached by 660 → censored ('right').
6. **Sequential rounds.** Random = 1 round (all lines in parallel). Adaptive = n_s / batch rounds; batch 10 / 50 /
   110 yeast queries = batch × (N_HAP1/N_yeast)^b HAP1 lines per round. The 110 arm (≤ 6 rounds to 660, §6.3) is
   read from its own stored curve (checkpoints 160 / 330 / 660).
7. **Pair-level discovery.** An unordered pair is SL if any measured direction among pool rows is a strong SL;
   'detected' = called in a screened row; 'touched' = involves a screened gene. T genes are dropped on both sides
   (final mode needs no T label). Stored P1/P4 cell counts include T columns, so recomputed cell counts are 0–25%
   lower (tested).
8. **'Genes picked early'.** First 160 picks of each stored trajectory vs that trajectory's candidate pool;
   standardised mean difference, CI by resampling trajectories within split. SMF exists only for on-array genes.
9. **Rank rule.** Target = the finalists' (H, C, F_r16) pick frequency in their first 160 (not the true degree, which
   needs labels); fitted leave-split-out with the evaluation split's D_val genes removed from the rule's training rows.
   Genes without SMF (off-array, a yeast artefact) are scored route-stratified as in the strategies' D5 calibration
   (ranked by the rule minus its SMF term, mapped onto the on-array score quantiles); otherwise a static SMF rule never
   picks them. Static selection = randomised top-20 picks (protocol top-k), seeds 0–2; P* refit outside the harness
   (no ledger writes), bit-identical to the stored AUPRC on dev1 r0 n = 160.
10. **Replicate reference.** The reciprocal label (BA) is a binary predictor of AB on on-array target cells where both
    directions are evaluation cells; by symmetry its precision equals its recall. P* is compared at matched recall
    and matched precision on the same cells. Billmann's 0.61 / 0.37 is a same-direction replicate screen (stricter
    reference) and is context only.
11. **LBO-CV cells.** Held-out batch rows × columns not screened in S_n (no reciprocal feature, as in E_fixed), not
    test genes, not self, measured. Batches = consecutive blocks of 10 in selection order; fold = batch mod 5
    (interleaved). Dev columns include D_val genes (unscreened arrays a lab would see); the final has no D_val.
12. **True knee and the stopping-rule family.** True knee on the trajectory's PAVA curve over the Q2 checkpoints,
    log-n interpolated (grid knee also given). DESIGN names a 'plateau rule' but no family. Two families were scored
    on dev: plateau (relative gain per doubling < δ for m consecutive checkpoints) and projection (current value ≥ τ ×
    log-linear projection to the planned maximum budget 660, i.e. 'within 1 − τ of where we would end'; this matches
    the knee definition, which is itself relative to AUPRC(660)). Each family uses either raw CV AUPRC or CV lift.
    Lift is included because adaptive arms screen less SL-rich rows as n grows, so their raw CV AUPRC falls while the
    held-out AUPRC rises. All rules stop at n ≥ 60, else at 660 (flagged). The choice is made on the same dev
    trajectories (in-sample), with a leave-split-out check and a no-information comparator (always stop at a fixed n).
13. **LBO scope.** For cost, LBO runs for R0 and H, trajectories 0–1, dev1–3 only (12 trajectories, one P* job at a
    time); the CV-threshold analysis covers those trajectories only.
"""


def key_findings(O) -> list[str]:
    """Headline numbers, read from the output files (no hardcoded values)."""
    L = ["## Key findings", ""]
    try:
        t = pd.read_csv(O / "hap1_translation.csv")
        x = t[(t.stratum == "OA") & t.N_yeast_label.str.contains("primary") & (t.N_HAP1_label == "mid")]
        m = x[x.b_label == "measured"].set_index("target")
        b0 = x[x.b_label == "b=0 (count)"].set_index("target")
        c = pd.read_csv(O / "hap1_costs.csv")
        c = c[(c.stratum == "OA") & (c.N_HAP1 == 8500) & (c.b_label == "measured")].set_index(["target", "arm"])
        L.append(f"- **Random screening, HAP1 lines** (N_HAP1 8,500, measured b, on-array stratum): n_80 ≈ "
                 f"{i0(m.loc['n80', 'n_HAP1'])} [{i0(m.loc['n80', 'lo'])}, {i0(m.loc['n80', 'hi'])}], n_90 ≈ "
                 f"{i0(m.loc['n90', 'n_HAP1'])} [{i0(m.loc['n90', 'lo'])}, {i0(m.loc['n90', 'hi'])}] lines "
                 f"(count-scaling bound b = 0: {i0(b0.loc['n80', 'n_HAP1'])} / {i0(b0.loc['n90', 'n_HAP1'])}).")
        for arm in ("H", "H__b110"):
            if ("n80", arm) in c.index:
                L.append(f"- **{arm}** reaches the same AUPRC with {i0(c.loc[('n80', arm), 'lines'])} (n_80) / "
                         f"{i0(c.loc[('n90', arm), 'lines'])} (n_90) lines in {i0(c.loc[('n80', arm), 'rounds'])} / "
                         f"{i0(c.loc[('n90', arm), 'rounds'])} rounds (saving "
                         f"{100 * c.loc[('n80', arm), 'saving_vs_random']:.0f}% / "
                         f"{100 * c.loc[('n90', arm), 'saving_vs_random']:.0f}%).")
        bb = pd.read_csv(O / "hap1_billmann222.csv")
        bb = bb[(bb.stratum == "OA") & (bb.N_HAP1 == 8500) & (bb.b_label == "measured")].set_index("arm")
        L.append(f"- **Billmann's 222 lines** ≈ {bb.loc['R0', 'n_yeast_equiv']:.0f} yeast queries at measured b: random "
                 f"selection gives {100 * bb.loc['R0', 'auprc_frac_of_max']:.0f}% of the pool-max AUPRC, H "
                 f"{100 * bb.loc['H', 'auprc_frac_of_max']:.0f}%.")
    except (FileNotFoundError, KeyError, IndexError, pd.errors.EmptyDataError):
        pass
    try:
        g = pd.read_csv(O / "choice_rule_gain.csv").set_index(["rule", "n"])
        k = g.loc[("forward_k3", 160)]
        L.append(f"- **Static 3-feature rank rule** (low single-mutant fitness, low expression variability, high "
                 f"conservation) recovers {100 * k.auprc_recovered_vs_H:.0f}% of H's AUPRC gain at n = 160 "
                 f"(H-static {100 * k.auprc_recovered_H_static_vs_H:.0f}%).")
    except (FileNotFoundError, KeyError, pd.errors.EmptyDataError):
        pass
    try:
        e = pd.read_csv(O / "stop_estimator.csv")
        r0 = e[(e.arm == "R0")]
        h = e[(e.arm == "H")]
        L.append(f"- **Online estimator**: LBO-CV AUPRC on screened rows has relative bias "
                 f"{100 * r0.rel_bias.min():+.0f}% to {100 * r0.rel_bias.max():+.0f}% under random selection but "
                 f"{100 * h.rel_bias.min():+.0f}% to {100 * h.rel_bias.max():+.0f}% under H (SL-rich screened rows).")
        f_ = pd.read_csv(O / "stop_fixed_n_baseline.csv").sort_values("mae_log2").iloc[0]
        s_ = json.loads((O / "stop_summary.json").read_text())["best"]
        kn = pd.read_csv(O / "stop_knees.csv").knee / 660
        L.append(f"- **Stopping**: true knees fall at {kn.min():.2f}–{kn.max():.2f} × the 660 cap; the best online rule "
                 f"has mean "
                 f"|log2 error| {s_['mae_log2']:.2f} vs {f_.mae_log2:.2f} for always stopping at n = {f_.n_fixed:.0f}.")
    except (FileNotFoundError, KeyError, TypeError, IndexError, pd.errors.EmptyDataError):
        pass
    return L + [""]


def compose(source: str, out=None) -> str:
    O = C.out_dir(source, out)
    head = [f"# Q3 analyses — {source.upper()} {'(DEV DATA ONLY)' if source == 'dev' else ''}", "",
            "DESIGN v3.1 §7.2–7.4. Code: `analysis/q3/` (one command per analysis, `--source dev|final|dryrun`); "
            "rerun everything with `bash analysis/q3/run_all.sh <source>`. Predictor: frozen P* "
            "(stack_lr_wr2_lb). All numbers below are generated from the CSV/JSON files in this directory.", "",
            "**On the final outputs** (only after `scripts/score_final.py` has finished): "
            "(1) `bash analysis/q3/run_all.sh final` (hap1, frontier, choice profile, predcontrib; minutes); "
            "(2) `SOURCE=final SPLITS=final TRAJS=\"0 1\" ARMS=\"R0 H\" bash analysis/q3/run_lbo.sh` (LBO-CV refits on "
            "the dev bundle, one P* job at a time, ~25 min per trajectory); (3) `bash analysis/q3/run_all.sh final` "
            "again, which applies the FROZEN stopping rule (`analysis/q3/stopping_rule.yaml`, never re-tuned) and the "
            "CV thresholds. Outputs go to results/q3_final/.", ""]
    body = key_findings(O) + section_hap1(O) + section_frontier(O) + section_choice(O) + section_pred(O) + section_stop(O)
    txt = "\n".join(head + body) + "\n" + NOTES
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
