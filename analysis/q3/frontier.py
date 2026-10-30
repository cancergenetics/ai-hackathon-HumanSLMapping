"""Q3 analysis 2 — discovery frontier (DESIGN v3.1 §7.2): AUPRC (E_fixed) vs SLs detected, per arm.

    .venv/bin/python -m analysis.q3.frontier --source dev|final|dryrun [--out DIR]

Counts are recomputed from the stored selections on the DEV bundle (rows of pool genes only; T rows are never
present), and pairs/cells involving a test gene (T) are dropped everywhere, so the final mode needs no T label.
  cell level (primary)  strong-SL cells (x, y) in screened rows x in S (directed, own label)
  pair level            the unordered pair {x, y} is an SL pair if any measured direction among dev rows is a strong
                        SL. 'detected' = an SL call in a screened row (x in S and Y(x, y) = 1);
                        'touched' = the SL pair involves a screened gene (x in S or y in S), whatever was called.
                        touched − detected = SL pairs a screen 'saw' but did not call (AB/BA disagreement, κ 0.14,
                        or the partner not on the array).
Fractions are of the dev network's totals (all pool-row cells / pairs, T excluded). AUPRC = stored P* E_fixed(_dev)
value of the same trajectory and budget (P1/P4 rows; final: runs.csv).
Outputs: frontier_traj.csv (arm, split, traj, n, counts, auprc), frontier.csv (arm x n means),
fig_frontier.png
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from analysis.q3 import common as C


class Network:
    """SL cells/pairs among dev-bundle rows, T genes excluded on both sides."""

    def __init__(self, Y, row_orfs, col_orfs, excluded):
        excl = set(map(str, excluded))
        row_orfs = np.asarray(row_orfs, dtype=str)
        col_orfs = np.asarray(col_orfs, dtype=str)
        rk = ~np.isin(row_orfs, list(excl))
        ck = ~np.isin(col_orfs, list(excl))
        Y = np.asarray(Y)[np.ix_(rk, ck)]
        row_orfs, col_orfs = row_orfs[rk], col_orfs[ck]
        r, c = np.nonzero(Y == 1)
        keep = row_orfs[r] != col_orfs[c]
        r, c = r[keep], c[keep]
        genes = np.unique(np.concatenate([row_orfs, col_orfs]))
        gi = {g: i for i, g in enumerate(genes)}
        a = np.array([gi[g] for g in row_orfs])[r]
        bcol = np.array([gi[g] for g in col_orfs])[c]
        G = len(genes)
        self.genes, self.gi, self.G = genes, gi, G
        self.cell_row = a                         # gene index of the screened (row) side of each SL cell
        self.cell_pair = np.minimum(a, bcol).astype(np.int64) * G + np.maximum(a, bcol)
        self.pairs = np.unique(self.cell_pair)
        self.pair_x, self.pair_y = self.pairs // G, self.pairs % G
        self.n_cells, self.n_pairs = len(self.cell_pair), len(self.pairs)

    def counts(self, S) -> dict:
        idx = np.array([self.gi[g] for g in S if g in self.gi], dtype=np.int64)
        inS = np.zeros(self.G, bool)
        inS[idx] = True
        cm = inS[self.cell_row]
        det = np.unique(self.cell_pair[cm])
        touched = inS[self.pair_x] | inS[self.pair_y]
        return {"cells_detected": int(cm.sum()), "pairs_detected": int(len(det)),
                "pairs_touched": int(touched.sum()),
                "cells_detected_frac": cm.sum() / self.n_cells, "pairs_detected_frac": len(det) / self.n_pairs,
                "pairs_touched_frac": touched.sum() / self.n_pairs}


def build_network():
    from analysis.q3 import cache
    d = cache.load()
    return Network(d["Y"], d["query_orfs"], d["array_orfs"], d["test"])


def analyse(source: str, out=None) -> pd.DataFrame:
    t0 = time.perf_counter()
    O = C.out_dir(source, out)
    net = build_network()
    rows = C.rows(source)
    rows = rows[rows.arm.isin(C.MAIN_ARMS)]
    sels = C.selections(source, arms=list(C.MAIN_ARMS) + (["R0_random"] if source == "dev" else []))
    recs = []
    for (arm, sp, r), g in rows.groupby(["arm", "split", "traj"]):
        sel = sels.get((arm, sp, int(r)))
        if sel is None:
            continue
        for _, x in g.iterrows():
            n = int(x.n)
            if n > len(sel["selected"]):
                continue
            recs.append({"arm": arm, "split": sp, "traj": int(r), "n": n, "auprc": x.auprc,
                         "sl_detected_stored": x.get("sl_detected"), **net.counts(sel["selected"][:n])})
    tr = pd.DataFrame(recs)
    tr.to_csv(O / "frontier_traj.csv", index=False)
    # mean over trajectories within split, then over splits
    m1 = tr.groupby(["arm", "split", "n"]).mean(numeric_only=True).reset_index()
    agg = m1.groupby(["arm", "n"]).agg(**{c: (c, "mean") for c in
                                          ["auprc", "cells_detected", "cells_detected_frac", "pairs_detected",
                                           "pairs_detected_frac", "pairs_touched", "pairs_touched_frac"]},
                                       n_splits=("split", "nunique")).reset_index()
    agg["detected_over_touched"] = agg.pairs_detected / agg.pairs_touched
    R = agg[agg.arm == C.RANDOM].set_index("n")
    for c in ("auprc", "cells_detected", "pairs_detected", "pairs_touched"):
        agg[f"{c}_rel_random"] = [x / R[c].get(n, np.nan) for x, n in zip(agg[c], agg.n)]
    agg.to_csv(O / "frontier.csv", index=False)
    C.write_json(O / "frontier_meta.json", {"source": source, "network_cells": net.n_cells,
                                             "network_pairs": net.n_pairs, "runtime_s": time.perf_counter() - t0})
    fig(agg, O / "fig_frontier.png", source)
    return agg


STYLE = {"R0": ("#6f6e69", "random"), "H": ("#2a78d6", "H"), "C": ("#eb6834", "C"), "F_r16": ("#1baf7a", "F_r16"),
         "H-static": ("#6250d6", "H-static"), "R1": ("#e34948", "R1 (Billmann-like)")}


def fig(agg, path, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})
    f, axes = plt.subplots(1, 2, figsize=(10, 4), facecolor="#fcfcfb")
    for ax, xc, xl in ((axes[0], "cells_detected_frac", "strong-SL cells detected (cell level), % of network"),
                       (axes[1], "pairs_touched_frac", "SL pairs touched (pair level), % of network")):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color="#e8e7e2", lw=0.6)
        for a, (col, lab) in STYLE.items():
            d = agg[(agg.arm == a) & (agg.n <= 660)].sort_values("n")
            if d.empty:
                continue
            ax.plot(d[xc] * 100, d.auprc, color=col, lw=2, marker="o", ms=3.5, label=lab)
            for nn in (160, 660):
                z = d[d.n == nn]
                if len(z) and a in ("R0", "H"):
                    ax.annotate(f"n={nn}", (z[xc].iloc[0] * 100, z.auprc.iloc[0]), xytext=(4, -10),
                                textcoords="offset points", fontsize=7, color="#52514e")
        ax.set_xlabel(xl)
        ax.set_ylabel("P* AUPRC (E_fixed)")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    f.suptitle(f"Discovery frontier ({source}): prediction quality vs SLs found, budgets 30-660", fontsize=10)
    f.tight_layout()
    f.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    C.set_threads(2)
    agg = analyse(a.source, a.out)
    print(agg[agg.n.isin([160, 660])].round(4).to_string())


if __name__ == "__main__":
    main()
