"""Q3 analysis 1 — HAP1 translation of query numbers (DESIGN v3.1 §7.3-7.4).

    .venv/bin/python -m analysis.q3.hap1 --source dev|final|dryrun [--out DIR]

n_HAP1 = n_yeast · (N_HAP1 / N_yeast)^b
  n_yeast   n_80 / n_90 of the P* random curve (budget reaching 80/90% of the AUPRC at the pool maximum; frozen
            PAVA + log-n interpolation). Dev: mean over dev1-3; final: q1_nx.csv.
  b         universe-subsampling slope (9 dev draws, results/p1/universe); bootstrap = draws resampled within level
            (as P1), OLS slope of log n_X on log N recomputed per replicate. b = 0 (count) and b = 1 (fraction)
            are shown as bounds.
  joint bootstrap: n_yeast replicate ~ lognormal matched to its 95% CI (per split, then averaged over splits on
            dev) x b replicate -> percentile CI of n_HAP1. (Approximation: the stored n_X CIs are percentile CIs of
            the frozen joint query bootstrap; their replicates are not stored.)
Strategy: for each arm, the budget n_s at which the arm's mean (PAVA) curve reaches the random target AUPRC
  (= frac · AUPRC at the pool max), i.e. the strategy budget whose random-equivalent budget n_req equals n_X.
  Censored '>max' if the arm never reaches it within its checkpoints.
Costs: lines, screens (x replicate factor) and sequential rounds (random: 1 round; adaptive: n_s / batch rounds of
  batch · scale lines). Billmann's 222 lines placed on the curves; marginal returns per +100 HAP1 lines.
Outputs (<out>/): hap1_translation.csv, hap1_costs.csv, hap1_billmann222.csv, hap1_marginal.csv, hap1_b_boot.csv,
  hap1_summary.json, fig_hap1_curves.png
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from analysis.q3 import common as C

FRACS = {"n80": 0.80, "n90": 0.90}


# ============================================================================ b bootstrap
def ols_slope(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    xm = x - x.mean()
    return float((xm * (y - y.mean())).sum() / (xm ** 2).sum())


def slope_boot(draws: pd.DataFrame, col: str, B: int, rng) -> tuple[float, np.ndarray]:
    """Point slope of log(col) on log(pool_n) and B draw-bootstrap replicates (resample draws within level)."""
    d = draws[np.isfinite(draws[col])]
    b0 = ols_slope(np.log(d.pool_n), np.log(d[col]))
    groups = [g.index.to_numpy() for _, g in d.groupby("level")]
    reps = np.empty(B)
    for i in range(B):
        idx = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in groups])
        x, y = np.log(d.loc[idx, "pool_n"]), np.log(d.loc[idx, col])
        reps[i] = ols_slope(x, y) if np.ptp(x) > 0 else np.nan
    return b0, reps


def nyeast_boot(nx: pd.DataFrame, col: str, B: int, rng) -> tuple[float, np.ndarray]:
    """Point = mean over splits; replicates: per split lognormal matched to the 95% CI, averaged over splits."""
    pts = nx[col].to_numpy(float)
    lo, hi = nx[f"{col}_lo"].to_numpy(float), nx[f"{col}_hi"].to_numpy(float)
    sd = (np.log(hi) - np.log(lo)) / (2 * 1.959964)
    mu = np.log(pts)
    reps = np.exp(mu[None, :] + sd[None, :] * rng.standard_normal((B, len(pts)))).mean(1)
    return float(pts.mean()), reps


def translate(n_yeast, N_hap1, N_yeast, b):
    return np.asarray(n_yeast, float) * (float(N_hap1) / float(N_yeast)) ** np.asarray(b, float)


def pct(x, q=(2.5, 97.5)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return (float(np.percentile(x, q[0])), float(np.percentile(x, q[1]))) if len(x) else (np.nan, np.nan)


# ============================================================================ curves
def mean_curves(rows: pd.DataFrame, arm: str) -> dict:
    """{split: DataFrame(n, auprc, sl_frac, R)} trajectory means per budget."""
    out = {}
    d = rows[rows.arm == arm]
    for sp, g in d.groupby("split"):
        m = g.groupby("n").agg(auprc=("auprc", "mean"), sl_frac=("sl_detected_frac", "mean"),
                               R=("traj", "nunique")).reset_index()
        out[sp] = m.sort_values("n")
    return out


def arm_budget_for_target(curve: pd.DataFrame, target: float) -> tuple[float, str]:
    """Strategy budget reaching `target` (frozen n_req on the arm's own PAVA curve). ('', 'left', 'right')."""
    from hsl import eval as E
    r = E.n_req(float(target), curve.n.to_numpy(float), curve.auprc.to_numpy(float))
    return float(r.n), r.censored


def auprc_at(curve: pd.DataFrame, n, col="auprc"):
    """PAVA curve value at budget n, log-n interpolation; NaN outside the measured range."""
    from hsl import eval as E
    f = E.pava(curve[col].to_numpy(float))
    x = curve.n.to_numpy(float)
    n = np.asarray(n, float)
    v = np.interp(np.log(np.clip(n, 1e-9, None)), np.log(x), f)
    return np.where((n >= x[0]) & (n <= x[-1]), v, np.nan)


# ============================================================================ main analysis
def analyse(source: str, out=None) -> dict:
    t0 = time.perf_counter()
    K = C.hap1_constants()
    O = C.out_dir(source, out)
    rng = np.random.default_rng(int(K["bootstrap"]["seed"]))
    B = int(K["bootstrap"]["B"])
    rows = C.rows(source)
    nx = C.nx_table(source)
    draws = C.universe_draws()
    NY = {"N_yeast=3283 (primary)": K["N_yeast"]["primary"], "N_yeast=4458 (Q∪A)": K["N_yeast"]["sensitivity_QuA"],
          "N_yeast=pool size": float(nx.pool_n.mean())}
    NH = K["N_HAP1"]
    rep = float(K["replicate_factor"])
    arms_all = [a for a in sorted(rows.arm.unique()) if a in C.MAIN_ARMS]
    strat_arms = [a for a in arms_all if a != C.RANDOM]
    rc = mean_curves(rows, C.RANDOM)
    anchors = dict(zip(nx.split, nx.auprc_anchor))

    # ---- b bootstrap
    bb = {}
    brows = []
    for col in FRACS:
        b0, reps = slope_boot(draws, col, B, rng)
        bb[col] = (b0, reps)
        lo, hi = pct(reps)
        brows.append({"quantity": col, "b": b0, "lo": lo, "hi": hi, "B": B, "n_draws": int(len(draws))})
    pd.DataFrame(brows).to_csv(O / "hap1_b_boot.csv", index=False)

    # ---- random n_X and strategy budgets (yeast)
    tr, costs = [], []
    for col, frac in FRACS.items():
        ny0, nyr = nyeast_boot(nx, col, B, rng)
        b0, br = bb[col]
        # strategy budget per split for the same target
        sb = {}
        for a in strat_arms:
            ac = mean_curves(rows, a)
            vals, cens = [], []
            for sp, cur in ac.items():
                if sp not in anchors:
                    continue
                n_s, c = arm_budget_for_target(cur, frac * anchors[sp])
                if c == "left":                      # reached before the first checkpoint: upper bound
                    n_s = float(cur.n.min())
                vals.append(n_s)
                cens.append(c)
            if vals:
                sb[a] = (float(np.mean(vals)) if all(c == "" for c in cens) else
                         (float("inf") if any(c == "right" for c in cens) else float(np.mean(vals))),
                         ";".join(c or "ok" for c in cens), [float(v) for v in vals],
                         int(max(cur.n.max() for cur in ac.values())))
        for nyl, NYv in NY.items():
            for nhl in ("low", "mid", "high"):
                Nh = float(NH[nhl])
                for bl, bv, brep in (("measured", b0, br), ("b=0 (count)", 0.0, None), ("b=1 (fraction)", 1.0, None)):
                    pt = float(translate(ny0, Nh, NYv, bv))
                    reps = translate(nyr, Nh, NYv, brep if brep is not None else bv)
                    lo, hi = pct(reps)
                    tr.append({"target": col, "N_yeast_label": nyl, "N_yeast": NYv, "N_HAP1_label": nhl,
                               "N_HAP1": Nh, "b_label": bl, "b": bv, "n_yeast": ny0, "n_HAP1": pt, "lo": lo,
                               "hi": hi, "screens": pt * rep, "screens_lo": lo * rep, "screens_hi": hi * rep})
                    if nyl.endswith("(primary)"):
                        scale = (Nh / NYv) ** bv
                        costs.append({"target": col, "N_HAP1": Nh, "b_label": bl, "arm": "random", "n_yeast": ny0,
                                      "lines": pt, "lines_lo": lo, "lines_hi": hi, "screens": pt * rep,
                                      "rounds": 1, "batch_yeast": None, "lines_per_round": pt, "censored": "",
                                      "saving_vs_random": 0.0})
                        for a, (ns, cens, per, nmax) in sb.items():
                            batch = 110 if "__b110" in a else (50 if "__b50" in a else 10)
                            lines = ns * scale if np.isfinite(ns) else np.inf
                            reps_s = (ns * (Nh / NYv) ** (brep if brep is not None else bv)
                                      if np.isfinite(ns) else np.full(B, np.inf))
                            slo, shi = pct(reps_s)
                            costs.append({"target": col, "N_HAP1": Nh, "b_label": bl, "arm": a, "n_yeast": ns,
                                          "n_yeast_per_split": "|".join(f"{v:.0f}" for v in per),
                                          "lines": lines, "lines_lo": slo, "lines_hi": shi,
                                          "screens": lines * rep,
                                          "rounds": (int(np.ceil(ns / batch)) if np.isfinite(ns) else None),
                                          "batch_yeast": batch, "lines_per_round": batch * scale,
                                          "censored": "" if set(cens.split(";")) == {"ok"} else cens,
                                          "max_checkpoint": nmax,
                                          "saving_vs_random": 1 - ns / ny0 if np.isfinite(ns) else np.nan})
    tr = pd.DataFrame(tr)
    costs = pd.DataFrame(costs)
    tr.to_csv(O / "hap1_translation.csv", index=False)
    costs.to_csv(O / "hap1_costs.csv", index=False)

    # ---- Billmann 222 lines + marginal returns (random and strategies), N_yeast primary
    NYp = float(K["N_yeast"]["primary"])
    bl_rows, mg = [], []
    b0 = bb["n80"][0]
    br = bb["n80"][1]
    L222 = float(K["billmann"]["lines"])
    curves = {a: mean_curves(rows, a) for a in arms_all}
    for nhl in ("low", "mid", "high"):
        Nh = float(NH[nhl])
        for bl, bv, brep in (("measured", b0, br), ("b=0 (count)", 0.0, None), ("b=1 (fraction)", 1.0, None)):
            ny = L222 * (NYp / Nh) ** bv
            nyr = L222 * (NYp / Nh) ** (brep if brep is not None else np.array([bv]))
            for a, cs in curves.items():
                fr, frr, slf = [], [], []
                for sp, cur in cs.items():
                    anc = anchors.get(sp)
                    if anc is None:
                        continue
                    fr.append(float(auprc_at(cur, ny)) / anc)
                    frr.append(auprc_at(cur, nyr) / anc)
                    slf.append(float(auprc_at(cur, ny, "sl_frac")))
                if not fr:
                    continue
                v = np.nanmean(np.vstack([np.atleast_1d(x) for x in frr]), 0)
                lo, hi = pct(v)
                bl_rows.append({"N_HAP1": Nh, "b_label": bl, "arm": a, "n_yeast_equiv": ny,
                                "auprc_frac_of_max": float(np.nanmean(fr)), "lo": lo, "hi": hi,
                                "sl_detected_frac": float(np.nanmean(slf))})
            # marginal returns per +100 lines (random + strategies within their measured range)
            for L in K["hap1_line_grid"]:
                n1, n2 = L * (NYp / Nh) ** bv, (L + 100) * (NYp / Nh) ** bv
                for a, cs in curves.items():
                    da, ds, a1 = [], [], []
                    for sp, cur in cs.items():
                        anc = anchors.get(sp)
                        if anc is None:
                            continue
                        v1, v2 = auprc_at(cur, n1), auprc_at(cur, n2)
                        s1, s2 = auprc_at(cur, n1, "sl_frac"), auprc_at(cur, n2, "sl_frac")
                        da.append(float(v2 - v1) / anc)
                        ds.append(float(s2 - s1))
                        a1.append(float(v1) / anc)
                    if not da or not np.isfinite(np.nanmean(da)):
                        continue
                    mg.append({"N_HAP1": Nh, "b_label": bl, "arm": a, "hap1_lines": L, "n_yeast": n1,
                               "auprc_frac_of_max": float(np.nanmean(a1)),
                               "d_auprc_frac_per100": float(np.nanmean(da)),
                               "d_auprc_abs_per100": float(np.nanmean(da)) * float(np.mean(list(anchors.values()))),
                               "d_sl_detected_frac_per100": float(np.nanmean(ds))})
    bl_df = pd.DataFrame(bl_rows)
    mg_df = pd.DataFrame(mg)
    bl_df.to_csv(O / "hap1_billmann222.csv", index=False)
    mg_df.to_csv(O / "hap1_marginal.csv", index=False)
    summ = {"source": source, "b": {k: {"b": v[0], "ci95": pct(v[1])} for k, v in bb.items()},
            "n_yeast": {k: float(nx[k].mean()) for k in FRACS}, "pool_n": float(nx.pool_n.mean()),
            "anchors": anchors, "arms": arms_all, "runtime_s": time.perf_counter() - t0}
    C.write_json(O / "hap1_summary.json", summ)
    fig(curves, anchors, b0, float(NH["mid"]), NYp, L222, O / "fig_hap1_curves.png", source)
    return summ


# ============================================================================ figure
ARM_STYLE = {"R0": ("#6f6e69", "random"), "H": ("#2a78d6", "H (hub-first)"), "C": ("#eb6834", "C (coverage)"),
             "F_r16": ("#1baf7a", "F_r16 (Fisher-D)"), "H-static": ("#6250d6", "H-static")}


def fig(curves, anchors, b, Nh, NY, L222, path, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})
    scale = (Nh / NY) ** b
    fig_, axes = plt.subplots(1, 2, figsize=(10, 4), facecolor="#fcfcfb")
    ends = []
    for ax, col, ylab in ((axes[0], "auprc", "AUPRC, % of pool-max random AUPRC"),
                          (axes[1], "sl_frac", "strong SLs detected, % of network")):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color="#e8e7e2", lw=0.6)
        for a, (colr, lab) in ARM_STYLE.items():
            cs = curves.get(a)
            if not cs:
                continue
            grid = sorted(set().union(*[set(c.n) for c in cs.values()]))
            ys = []
            for sp, cur in cs.items():
                v = auprc_at(cur, grid, col)
                ys.append(v / anchors[sp] * 100 if col == "auprc" else v * 100)
            y = np.nanmean(np.vstack(ys), 0)
            x = np.asarray(grid, float) * scale
            ax.plot(x, y, color=colr, lw=2, marker="o", ms=3.5, label=lab)
            if col == "auprc":
                ends.append([x[np.isfinite(y)][-1], y[np.isfinite(y)][-1], lab.split(" ")[0]])
        if ends:                                   # direct labels, pushed apart vertically (min gap 3 points)
            ends.sort(key=lambda e: e[1])
            for k in range(1, len(ends)):
                ends[k][1] = max(ends[k][1], ends[k - 1][1] + 3.0)
            for ex, ey, et in ends:
                ax.annotate(et, (ex, ey), xytext=(4, 0), textcoords="offset points", color="#52514e", fontsize=8,
                            va="center")
            ends = []
        ax.axvline(L222, color="#9a9992", lw=1, ls="--")
        ax.text(L222, ax.get_ylim()[1] * 0.97, " Billmann 222", color="#52514e", fontsize=8, va="top")
        if col == "auprc":
            for f in (80, 90):
                ax.axhline(f, color="#c3c2b7", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(f"HAP1 query lines (N_HAP1 = {Nh:.0f}, b = {b:.2f})")
        ax.set_ylabel(ylab)
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    fig_.suptitle(f"Q3 HAP1 translation ({source}): yeast P* curves rescaled by (N_HAP1/N_yeast)^b", fontsize=10,
                  color="#0b0b0b")
    fig_.tight_layout()
    fig_.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig_)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    s = analyse(a.source, a.out)
    print({k: s[k] for k in ("b", "n_yeast", "runtime_s")})


if __name__ == "__main__":
    main()
