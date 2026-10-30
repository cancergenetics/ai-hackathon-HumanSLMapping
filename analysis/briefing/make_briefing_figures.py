"""Figures for BRIEFING.md that do not exist elsewhere: (1) the benchmark design schematic, (2) the project timeline.

The timeline phases and review times are taken from the git log (commit timestamps, 2026-10-29/30); see BRIEFING.md §2.
Run: .venv/bin/python analysis/briefing/make_briefing_figures.py  -> results/briefing_figs/
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parents[2] / "results" / "briefing_figs"
INK, MUTED, LINE = "#11201c", "#4d5c58", "#d3ddd9"
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#a9b3b0"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})


def design_schematic() -> None:
    """Query x array matrix: candidate pool rows (some screened), held-out test rows, and the evaluation block."""
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    ax.set_xlim(-3.2, 14.6), ax.set_ylim(-1.4, 11.6), ax.axis("off")
    # columns: [0,6) arrays that are also candidate queries; [6,8) test genes on the array; [8,10.5) array-only genes
    W_pool, W_test, W_ao = 6.0, 2.0, 2.5
    H_pool, H_test = 8.0, 2.0
    y_test0, y_pool0 = 0.0, 2.0
    # candidate pool rows (unscreened)
    ax.add_patch(Rectangle((0, y_pool0), W_pool + W_test + W_ao, H_pool, fc="#eef2f0", ec=LINE))
    # screened rows (revealed)
    for k, y in enumerate([8.6, 7.4, 5.6, 4.1, 2.9]):
        ax.add_patch(Rectangle((0, y), W_pool + W_test + W_ao, 0.42, fc=BLUE, ec="none", alpha=0.85))
    # held-out test rows
    ax.add_patch(Rectangle((0, y_test0), W_pool + W_test + W_ao, H_test, fc="#fbe3d8", ec=LINE))
    # evaluation block E_fixed: test rows x columns of genes that can never be screened (test genes + array-only)
    ax.add_patch(Rectangle((W_pool, y_test0), W_test + W_ao, H_test, fc="none", ec=ORANGE, lw=2.4))
    ax.text(W_pool + (W_test + W_ao) / 2, y_test0 + H_test / 2, "evaluation cells\n(held-out rows ×\nnever-screenable\ncolumns)",
            ha="center", va="center", fontsize=9, color=INK)
    # column labels
    for x0, w, lab in [(0, W_pool, "candidate genes\n(can be screened)"), (W_pool, W_test, "held-out\ngenes"),
                       (W_pool + W_test, W_ao, "array-only\ngenes")]:
        ax.text(x0 + w / 2, 10.25, lab, ha="center", va="bottom", fontsize=9, color=MUTED)
        ax.plot([x0, x0 + w], [10.15, 10.15], color=MUTED, lw=1)
    ax.text(5.25, 11.3, "array side (genes whose interactions are measured)", ha="center", fontsize=10, color=INK)
    # row labels
    ax.text(-0.25, y_pool0 + H_pool / 2, "candidate\nquery pool\n(2,627 genes)", ha="right", va="center", fontsize=9, color=MUTED)
    ax.text(-0.25, y_test0 + H_test / 2, "held-out\ntest queries\n(656 genes, 20%)", ha="right", va="center", fontsize=9, color=ORANGE)
    ax.text(-2.9, 6.0, "query side", rotation=90, ha="center", va="center", fontsize=10, color=INK)
    # annotation for screened rows
    ax.annotate("screening a query\nreveals its whole\nrow of interactions", xy=(10.5, 8.8), xytext=(11.0, 8.2),
                fontsize=9, color=BLUE, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=1))
    ax.text(5.25, -0.75, "Each strategy chooses which pool rows to screen. The predictor sees only screened rows plus "
            "human-analogue gene features,\nand must rank the true strong negatives in the orange block, "
            "which is identical for every strategy.", ha="center", va="top", fontsize=8.6, color=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / "figB1_design.png", dpi=200, facecolor="white")
    plt.close(fig)


def timeline() -> None:
    """Phases from commit times (hours since 07:45 on 29 Oct 2026); reviews and human inputs as markers."""
    t0 = 7.75

    def h(day: int, hm: str) -> float:
        hh, mm = map(int, hm.split(":"))
        return (day - 29) * 24 + hh + mm / 60 - t0

    phases = [  # (label, start, end, colour)
        ("Literature, data reconnaissance", h(29, "07:45"), h(29, "08:53"), BLUE),
        ("Design + 3 review rounds", h(29, "08:53"), h(29, "09:18"), BLUE),
        ("Infrastructure (P0) + review", h(29, "09:18"), h(29, "09:52"), BLUE),
        ("Strategy code + review", h(29, "09:52"), h(29, "11:49"), BLUE),
        ("P1 random-screening curves (dev)", h(29, "09:50"), h(30, "00:20"), BLUE),
        ("P2 predictor development + review", h(29, "10:21"), h(29, "14:19"), BLUE),
        ("P3 oracle references", h(29, "14:45"), h(29, "15:29"), BLUE),
        ("P4 strategy comparison (dev) + review", h(29, "15:29"), h(29, "22:10"), BLUE),
        ("Final pipeline build, dry run, review", h(29, "22:10"), h(30, "00:21"), BLUE),
        ("One-shot test evaluation + scoring", h(30, "00:22"), h(30, "13:40"), ORANGE),
        ("Q3 HAP1 analyses", h(30, "00:30"), h(30, "16:08"), BLUE),
        ("Independent rescoring", h(30, "13:45"), h(30, "14:28"), AQUA),
        ("Report writing + reviews", h(30, "14:30"), h(30, "16:16"), BLUE),
    ]
    reviews = [h(29, t) for t in ("09:02", "09:12", "09:18", "09:49", "10:29", "14:19", "22:10", "23:38")] + \
              [h(30, t) for t in ("00:27", "14:30", "16:16")]
    freezes = [(h(29, "09:18"), "design frozen"), (h(29, "13:36"), "predictor frozen"), (h(30, "00:21"), "finalists frozen")]
    human = [(h(29, "11:10"), "disk freed"), (h(29, "19:55"), "\"only 50% CPU?\"")]

    fig, ax = plt.subplots(figsize=(10, 5.4))
    n = len(phases)
    for i, (lab, s, e, c) in enumerate(phases):
        y = n - i
        ax.barh(y, e - s, left=s, height=0.62, color=c, alpha=0.9)
        ax.text(-0.4, y, lab, ha="right", va="center", fontsize=9)
    ax.scatter(reviews, [n + 1.0] * len(reviews), marker="D", s=26, color=INK, zorder=3)
    ax.text(-0.4, n + 1.0, "independent reviews", ha="right", va="center", fontsize=9, color=INK)
    for x, lab in freezes:
        ax.axvline(x, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
        ax.text(x + 0.2, n + 0.45, lab, fontsize=7.5, color=MUTED, va="center")
    for x, lab in human:
        ax.scatter([x], [n + 1.8], marker="v", s=40, color=GREY, zorder=3)
        ax.text(x + 0.3, n + 1.8, lab, fontsize=8, color=MUTED, va="center")
    ax.text(-0.4, n + 1.8, "human inputs (selected)", ha="right", va="center", fontsize=9, color=MUTED)
    ax.set_yticks([]), ax.set_ylim(0, n + 2.5), ax.set_xlim(0, 33)
    ticks = [0, 4.25, 8.25, 12.25, 16.25, 20.25, 24.25, 28.25, 32.25]
    ax.set_xticks(ticks, ["07:45\n29 Oct", "12:00", "16:00", "20:00", "00:00\n30 Oct", "04:00", "08:00", "12:00", "16:00"])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("wall-clock time (one 8-core laptop; ≈ 32.5 h from brief to final report)", color=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / "figB2_timeline.png", dpi=200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    design_schematic()
    timeline()
    print("wrote", sorted(p.name for p in OUT.glob("*.png")))
