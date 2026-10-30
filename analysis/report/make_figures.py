"""REPORT.md figures (Q1/Q2). Reads committed CSVs in results/final/ only; writes PNGs to results/report_figs/.

No model fits, no label access, no hsl import. Run:  .venv/bin/python analysis/report/make_figures.py

Colours: the validated default categorical palette of the dataviz skill (references/palette.md, light mode;
slots 1-3 validate all-pairs for CVD). Finalists H, C, F_r16 take slots 1-3 in that fixed order everywhere;
random/R1 are neutral grey (comparators); static variants reuse their adaptive parent's hue with a dashed line.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "2")

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FINAL = ROOT / "results" / "final"
OUT = ROOT / "results" / "report_figs"

# ---- design tokens (dataviz reference palette, light mode) ----
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e1"
S1, S2, S3, S4, S5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
GREY_DARK = "#52514e"   # random comparator
GREY_MID = "#9a9893"    # R1 / second comparator
ARM_COLOUR = {"H": S1, "C": S2, "F_r16": S3}
ARM_LABEL = {"H": "H (hub-first)", "C": "C (coverage)", "F_r16": "F (Fisher-D, r=16)"}
POOL = 2627

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5,
    "axes.titlesize": 10.5, "axes.titleweight": "semibold", "axes.titlelocation": "left",
    "legend.frameon": False, "lines.linewidth": 2, "lines.solid_capstyle": "round",
    "text.color": INK,
})


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(FINAL / name)


def _logx(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,}" for t in ticks])
    ax.minorticks_off()


def _end_label(ax, x, y, text, dy=0.0, colour=INK2):
    ax.annotate(text, (x, y), xytext=(6, dy), textcoords="offset points", va="center", ha="left",
                fontsize=8.5, color=colour)


# ---------------------------------------------------------------- Fig 1: Q1 scaling on T
def fig_q1():
    q1 = _read("q1_curves.csv")
    nx = _read("q1_nx.csv").set_index("predictor")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    series = [("pstar", "P* (stack_lr_wr2_lb)", S1), ("m1_hgb", "m1_hgb (R = 5)", S5),
              ("stack_lr", "stack_lr", S2), ("B1_array_rate", "B1 array rate", S3),
              ("B2_product", "B2 product", S4)]
    for key, lab, col in series:
        d = q1[q1.predictor == key].sort_values("n")
        ax.fill_between(d.n, d.ci_lo, d.ci_hi, color=col, alpha=0.10, linewidth=0)
        ax.plot(d.n, d.auprc_mean, color=col, label=lab, marker="o", ms=3.5,
                markeredgecolor=SURFACE, markeredgewidth=1)
        _end_label(ax, d.n.iloc[-1], d.auprc_mean.iloc[-1], f"{lab.split(' (')[0]} {d.auprc_mean.iloc[-1]:.3f}")
    prev = q1.prevalence.iloc[0]
    ax.axhline(prev, color=INK2, linewidth=1, linestyle=(0, (1, 2)))
    ax.text(11, prev * 1.07, f"prevalence {prev:.4f}", fontsize=8, color=INK2)
    p = nx.loc["pstar"]
    for col_n, lo, hi, tag in (("n80", "n80_lo", "n80_hi", "n_80"), ("n90", "n90_lo", "n90_hi", "n_90")):
        ax.axvspan(p[lo], p[hi], color=S1, alpha=0.06, linewidth=0)
        ax.axvline(p[col_n], color=S1, linewidth=1)
        ax.text(p[col_n], 0.0035, f" {tag} = {p[col_n]:,.0f}", fontsize=8, color=INK2, rotation=90,
                va="bottom", ha="right")
    ax.axvline(160, color=INK2, linewidth=0.8, linestyle=(0, (1, 2)))
    ax.axvline(660, color=INK2, linewidth=0.8, linestyle=(0, (1, 2)))
    ax.text(160, 0.112, "160 (4.9%)", fontsize=8, color=INK2, ha="center")
    ax.text(660, 0.112, "660 (20.1%)", fontsize=8, color=INK2, ha="center")
    _logx(ax, [10, 30, 100, 300, 1000, 2627])
    ax.set_xlim(8, 9000)
    ax.set_ylim(0, 0.118)
    ax.set_xlabel("query genes screened (random order; dev pool max 2,627)")
    ax.set_ylabel("pooled AUPRC on test set T (E_fixed)")
    ax.set_title("Q1. Random screening: AUPRC keeps rising to the full pool (no plateau)")
    ax.legend(loc="upper left", fontsize=8.5)
    fig.text(0.01, 0.005, "Mean over R permutations (R = 40 at n ≤ 660, 10 above, 5 at 2,627); bands = joint-bootstrap 95% CI. "
             "\nn_80/n_90: budget reaching 80/90% of the full-pool AUPRC (shaded = CI). Source: results/final/q1_curves.csv, q1_nx.csv.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 0.97, 1))
    fig.savefig(OUT / "fig1_q1_scaling.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- Fig 2: Q2 curves and Δ(n)
def fig_q2():
    q2 = _read("q2_curves.csv")
    q1 = _read("q1_curves.csv")
    r0 = q1[(q1.predictor == "pstar") & (q1.n >= 30) & (q1.n <= 660)].sort_values("n")
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), gridspec_kw={"width_ratios": [1.15, 1, 1]})

    ax = axes[0]
    ax.fill_between(r0.n, r0.ci_lo, r0.ci_hi, color=GREY_DARK, alpha=0.10, linewidth=0)
    ax.plot(r0.n, r0.auprc_mean, color=GREY_DARK, label="random (R0)")
    r1 = q2[q2.arm == "R1"].sort_values("n")
    ax.plot(r1.n, r1.s_mean, color=GREY_MID, linestyle=(0, (4, 2)), label="R1 Billmann-like")
    for arm in ("H", "C", "F_r16"):
        d = q2[q2.arm == arm].sort_values("n")
        ax.plot(d.n, d.s_mean, color=ARM_COLOUR[arm], label=ARM_LABEL[arm], marker="o", ms=3.5,
                markeredgecolor=SURFACE, markeredgewidth=1)
    for arm, col in (("H-static", S1), ("F-static_r16", S3)):
        d = q2[q2.arm == arm].sort_values("n")
        ax.plot(d.n, d.s_mean, color=col, linestyle=(0, (2, 2)), linewidth=1.5,
                label=f"{arm.replace('_r16', '')} (static)")
    r490 = float(r0.loc[r0.n == 490, "auprc_mean"].iloc[0])
    r660 = float(r0.loc[r0.n == 660, "auprc_mean"].iloc[0])
    ax.axhline(r490, color=GRID, linewidth=1)
    ax.axhline(r660, color=GRID, linewidth=1)
    ax.text(31, r490 + 0.0008, f"random @490 = {r490:.4f}", fontsize=7.5, color=INK2)
    ax.text(31, r660 + 0.0008, f"random @660 = {r660:.4f}", fontsize=7.5, color=INK2)
    ax.axvline(160, color=INK2, linewidth=0.8, linestyle=(0, (1, 2)))
    _logx(ax, [30, 60, 100, 160, 330, 660])
    ax.set_ylim(0.025, 0.097)
    ax.set_xlabel("query genes screened (batches of 10)")
    ax.set_ylabel("pooled AUPRC on T (E_fixed)")
    ax.set_title("a  AUPRC by budget")
    ax.legend(loc="upper left", fontsize=7.5)

    def delta_panel(ax, arms, title):
        offs = {a: o for a, o in zip(arms, (-0.035, 0.0, 0.035))}
        for arm in arms:
            d = q2[q2.arm == arm].sort_values("n")
            base = arm.replace("-static", "").replace("_r16", "")
            col = ARM_COLOUR.get(arm, {"H": S1, "F": S3}.get(base, GREY_MID))
            ls = "-" if arm in ARM_COLOUR else (0, (2, 2)) if "static" in arm else (0, (4, 2))
            if arm == "R1":
                col = GREY_MID
            x = d.n * (10 ** offs[arm])
            ax.fill_between(d.n, d.band_lo, d.band_hi, color=col, alpha=0.08, linewidth=0)
            ax.plot(d.n, d.delta, color=col, linestyle=ls, label=ARM_LABEL.get(arm, arm.replace("_r16", "")))
            ax.errorbar(x, d.delta, yerr=[d.delta - d.lo95, d.hi95 - d.delta], fmt="none", ecolor=col,
                        elinewidth=1, capsize=0, alpha=0.9)
        ax.axhline(0, color=INK2, linewidth=1)
        _logx(ax, [30, 60, 100, 160, 330, 660])
        ax.set_ylim(-0.008, 0.045)
        ax.set_xlabel("query genes screened")
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=7.5)

    delta_panel(axes[1], ["H", "C", "F_r16"], "b  Δ AUPRC vs random, finalists")
    axes[1].set_ylabel("Δ pooled AUPRC (arm − random, same n)")
    delta_panel(axes[2], ["R1", "H-static", "F-static_r16"], "c  Δ vs random: R1 and static variants")
    fig.text(0.01, 0.005, "Trajectory means (R_s = 20, R_r = 40). Bars = 95% Welch CI; shaded = Bonferroni band over the 8 checkpoints. "
             "Confirmatory test at n = 160 only. Source: results/final/q2_curves.csv, q1_curves.csv.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "fig2_q2_strategies.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- Fig 3: frontier and E_unobs
def fig_frontier():
    disc = _read("discovery.csv")
    eu = _read("e_unobs.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    styles = [("R0", "random (R0)", GREY_DARK, "-"), ("R1", "R1 Billmann-like", GREY_MID, (0, (4, 2))),
              ("H-static", "H-static", S1, (0, (2, 2))), ("H", ARM_LABEL["H"], S1, "-"),
              ("C", ARM_LABEL["C"], S2, "-"), ("F_r16", ARM_LABEL["F_r16"], S3, "-")]
    ax = axes[0]
    for arm, lab, col, ls in styles:
        d = disc[(disc.arm == arm) & (disc.n >= 30) & (disc.n <= 660)].sort_values("n")
        ax.plot(d.sl_detected_frac_mean * 100, d.auprc_mean, color=col, linestyle=ls, label=lab,
                marker="o", ms=3.5, markeredgecolor=SURFACE, markeredgewidth=1)
        for n in (160, 660):
            r = d[d.n == n]
            if len(r) and arm in ("R0", "H"):
                ax.annotate(f"n={n}", (float(r.sl_detected_frac_mean.iloc[0] * 100), float(r.auprc_mean.iloc[0])),
                            xytext=(4, -11) if arm == "R0" else (-6, 9), textcoords="offset points",
                            ha="left" if arm == "R0" else "right", fontsize=8, color=INK2)
    ax.set_xlabel("strong SLs detected in screened rows (% of all dev-pool SLs)")
    ax.set_ylabel("pooled AUPRC on T (E_fixed)")
    ax.set_title("a  Discovery frontier, n = 30 → 660")
    ax.legend(loc="lower right", fontsize=7.5)

    ax = axes[1]
    for arm, lab, col, ls in styles:
        if arm in ("H-static",):
            continue
        d = eu[(eu.arm == arm) & (eu.n <= 660)].sort_values("n")
        last = d[d.n == 660]
        ax.plot(d.n, d.auprc_mean, color=col, linestyle=ls, marker="o", ms=3.5,
                markeredgecolor=SURFACE, markeredgewidth=1,
                label=f"{lab}; prevalence @660 {float(last.prevalence_mean.iloc[0]):.4f}")
    _logx(ax, [30, 60, 100, 160, 330, 660])
    ax.set_xlim(25, 800)
    ax.set_xlabel("query genes screened")
    ax.set_ylabel("AUPRC on arm-specific unobserved test pairs")
    ax.set_title("b  E_unobs (brief's literal set; unpaired, no CI)")
    ax.legend(loc="upper left", fontsize=7.5)
    fig.text(0.01, 0.005, "E_unobs = test rows × arrays not screened by that arm, so its size and prevalence differ by arm "
             "(legend: prevalence at n = 660; random 0.0109). Source: results/final/discovery.csv, e_unobs.csv.",
             fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "fig3_frontier_eunobs.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- Fig 4: warm start and batch size
def fig_scenarios():
    c = _read("contrasts.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))

    ax = axes[0]
    for i, (arm, base) in enumerate((("H__ws220", "H"), ("C__ws220", "C"))):
        d = c[(c.a == arm) & (c.b == "R1") & (c.family == "descriptive:R1_continuation")].sort_values("n_a")
        x = d.n_a + (i - 0.5) * 12
        rel, lo, hi = d.rel * 100, d.lo / d.s_b * 100, d.hi / d.s_b * 100
        ax.errorbar(x, rel, yerr=[rel - lo, hi - rel], fmt="o", color=ARM_COLOUR[base], ms=6,
                    markeredgecolor=SURFACE, markeredgewidth=1.5, elinewidth=1.5, capsize=0,
                    label=f"{base} after 220 R1 lines")
        ax.plot(x, rel, color=ARM_COLOUR[base], linewidth=1, alpha=0.6)
    ax.axhline(0, color=INK2, linewidth=1)
    ax.set_xticks([330, 440, 490, 660])
    ax.set_xlabel("total query genes screened (first 220 = R1 panel)")
    ax.set_ylabel("% AUPRC gain over continuing R1 (95% CI)")
    ax.set_title("a  Warm start: switch to adaptive after a Billmann-like panel")
    ax.legend(loc="upper left", fontsize=8)

    ax = axes[1]
    q2 = _read("q2_curves.csv")
    markers = {10: "o", 50: "s", 110: "^"}
    for j, base in enumerate(("H", "C")):
        for k, b in enumerate((10, 50, 110)):
            if b == 10:  # batch-10 finalist: Welch CI from the Δ(n) curve table
                qq = q2[(q2.arm == base) & (q2.n.isin([160, 330, 660]))].sort_values("n")
                d = pd.DataFrame({"n_a": qq.n, "rel": qq.rel, "lo": qq.lo95, "hi": qq.hi95, "s_b": qq.r_mean})
            else:
                d = c[(c.a == f"{base}__b{b}") & (c.b == "R0") & (c.family == "descriptive:random")
                      & (c.set == "all") & (c.n_a == c.n_b)].sort_values("n_a")
            pos = {160: 0, 330: 1, 660: 2}
            x = [pos[int(n)] + (j * 3 + k - 2.5) * 0.11 for n in d.n_a]
            rel, lo, hi = d.rel * 100, d.lo / d.s_b * 100, d.hi / d.s_b * 100
            ax.errorbar(x, rel, yerr=[rel - lo, hi - rel], fmt=markers[b], color=ARM_COLOUR[base], ms=6,
                        markeredgecolor=SURFACE, markeredgewidth=1.2, elinewidth=1.3, capsize=0,
                        label=f"{base}, batch {b}")
    ax.axhline(0, color=INK2, linewidth=1)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["160", "330", "660"])
    ax.set_xlabel("query genes screened")
    ax.set_ylabel("% AUPRC gain over random at same n (95% CI)")
    ax.set_title("b  Batch size 10 / 50 / 110 (110 = 6 rounds to 660)")
    ax.set_ylim(-5, 80)
    ax.legend(loc="upper center", fontsize=7.5, ncol=3)
    fig.text(0.01, 0.005, "Descriptive arms (R = 20; R1 20, random 40). Batch-10 points are the main H and C arms (Δ(n) curve table). "
             "Source: results/final/contrasts.csv, q2_curves.csv.", fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "fig4_warmstart_batch.png", dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig_q1()
    fig_q2()
    fig_frontier()
    fig_scenarios()
    print("wrote", sorted(p.name for p in OUT.glob("*.png")))


if __name__ == "__main__":
    main()
