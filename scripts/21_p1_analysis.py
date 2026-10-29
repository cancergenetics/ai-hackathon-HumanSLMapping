#!/usr/bin/env python
"""P1 analysis (DESIGN v3.1 §6.4, §7.1, §11; DEV DATA ONLY).

    pilot   power pilot on dev split 1 (R0, P*): σ_r(n), V_Q,single(n) from the frozen ReplicateEngine (B from
            protocol.yaml), g = r210 − r160, V_Q(Δ) := 2·V_Q,single(160)·|D_val|/|T|, R_r = 2·R_s,
            SE(Δ̂(160)) ≤ g/3.5  ->  results/p1/power_pilot.json + .md
    full    curves.csv -> results/p1/summary.md, summary.json, curve tables, figures; n_80/n_90 with joint-bootstrap
            CIs via the frozen PAVA/n_req inversion; runtime table; upstream-protocol sanity; saturating fits.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import math
import sys
import time

import numpy as np
import pandas as pd

from hsl import config, eval as E, harness as H, p1, views

OUT = p1.P1_DIR
N_T = 656                                                # |T| (frozen test split)
KEY_N = (30, 100, 160, 220, 330, 490, 660, 1310)
PS, REF = p1.PSTAR, p1.REF_PSTAR            # frozen P* (config/predictor.yaml) and the pre-P2 P* (reference)
PRED_ORDER = (PS, REF, "stack_lr", "B1_array_rate", "B2_product", "m1_hgb")
LABEL = {PS: f"P* {PS}", REF: f"{REF} (pre-P2 P*)", "stack_lr": "stack_lr", "B1_array_rate": "B1 array rate",
         "B2_product": "B2 product", "m1_hgb": "m1_hgb (no genome)"}
COLORS = {PS: "#2a78d6", REF: "#6250d6", "stack_lr": "#eb6834", "B1_array_rate": "#1baf7a", "B2_product": "#eda100",
          "m1_hgb": "#e87ba4"}
MARK = {PS: "o", REF: "P", "stack_lr": "s", "B1_array_rate": "^", "B2_product": "v", "m1_hgb": "D"}


def curves() -> pd.DataFrame:
    return pd.read_csv(OUT / "curves.csv")


def replicate_engine(df: pd.DataFrame, split: str, predictor: str, name: str, grid=None):
    """Frozen ReplicateEngine for one (split, predictor): A[b, r, j] with stratified query weights."""
    d = df[(df.split == split) & (df.predictor == predictor)]
    grid = sorted(d.n.unique()) if grid is None else list(grid)
    perms = sorted(d.perm.unique())
    W = p1.query_weights(split, name)
    eng = E.ReplicateEngine(W, tuple(grid), {"R0": len(perms)})
    es = H.get_eset(views.dev_bundle(), views.get_split(split), "E_fixed")
    dg = es.digest()
    for _, row in d.iterrows():
        if row.n not in grid:
            continue
        if row.eset_digest != dg:
            raise AssertionError(f"eset digest mismatch for {row.run_id}")
        s = p1.load_scores(row.pred_file)
        eng.add("R0", perms.index(row.perm), int(row.n), es.y, s, es.rows)
    # unit-weight estimate must equal the stored AUPRC (same tie-grouped AP)
    est = eng.est["R0"]
    for _, row in d.iterrows():
        if row.n in grid:
            v = est[perms.index(row.perm), list(grid).index(row.n)]
            if abs(v - row.auprc) > 1e-9:
                raise AssertionError(f"engine AP {v} != stored {row.auprc} ({row.run_id})")
    return eng, perms


# ============================================================================ power pilot
def pilot(a) -> dict:
    df = curves()
    split = "dev1"
    pred = a.predictor
    d = df[(df.split == split) & (df.predictor == pred)]
    R = len(d.perm.unique())
    t0 = time.perf_counter()
    eng, perms = replicate_engine(df, split, pred, "p1_power_pilot_dev1")
    A, est = eng.A["R0"], eng.est["R0"]
    grid = list(eng.grid)
    n_val = len(views.get_split(split).targets)
    per_n = {}
    for j, n in enumerate(grid):
        ok = np.isfinite(est[:, j])
        if ok.sum() < 2:
            continue
        per_n[int(n)] = {"R": int(ok.sum()), "mean_auprc": float(est[ok, j].mean()),
                         "sigma_r": float(est[ok, j].std(ddof=1)),
                         "V_Q_single": float(np.var(A[:, ok, j].mean(1), ddof=1)),
                         "V_Q_one_traj_mean": float(np.mean(np.var(A[:, ok, j], axis=0, ddof=1)))}
    j160, j210 = grid.index(160), grid.index(210)
    g_pair = E.paired(A[:, :, j210], est[:, j210], A[:, :, j160], est[:, j160])
    g = g_pair["delta"]
    vq_single = per_n[160]["V_Q_single"]
    VQ = 2.0 * vq_single * n_val / N_T
    sig = per_n[160]["sigma_r"]
    tgt = g / 3.5
    room = tgt ** 2 - VQ
    out = {"split": split, "predictor": pred, "R_pilot": R, "B": int(A.shape[0]), "n_val": n_val, "n_T": N_T,
           "bootstrap_seed_name": "p1_power_pilot_dev1", "per_n": per_n,
           "r160": per_n[160]["mean_auprc"], "r210": per_n[210]["mean_auprc"], "g": g,
           "g_paired_se": g_pair["se"], "g_paired_ci95": g_pair["ci95"],
           "V_Q_single_160": vq_single, "V_Q_Delta": VQ, "sigma_r_160": sig, "SE_target": tgt,
           "SE_target_sq": tgt ** 2, "rule": "SE^2 = V_Q(Δ) + σ²/R_s + σ²/R_r, R_r = 2R_s, σ_s := σ_r(160)"}
    if g <= 0:
        out.update(R_s=None, R_r=None, verdict="g <= 0 in the pilot: no finite R achieves SE <= g/3.5; Q2 claims "
                   "must be phrased as bounds (DESIGN §7.1)")
    elif room <= 0:
        out.update(R_s=None, R_r=None, verdict="V_Q(Δ) alone exceeds (g/3.5)^2: Q2 claims phrased as bounds "
                   "(DESIGN §7.1)")
    else:
        rs_exact = 1.5 * sig ** 2 / room
        rs = max(2, math.ceil(rs_exact))
        out.update(R_s_exact=rs_exact, R_s=rs, R_r=2 * rs,
                   SE_at_R=math.sqrt(VQ + sig ** 2 / rs + sig ** 2 / (2 * rs)),
                   verdict=f"provisional R_s = {rs}, R_r = {2 * rs}")
    out["SE_by_R_s"] = {str(r): math.sqrt(VQ + sig ** 2 / r + sig ** 2 / (2 * r)) for r in (3, 5, 10, 20, 40)}
    out["fraction_SE2_from_VQ_at_Rs"] = (VQ / out["SE_at_R"] ** 2) if out.get("SE_at_R") else None
    # SENSITIVITY (not the frozen rule): arms share the query-bootstrap weights, so V_Q of a contrast is reduced by
    # the replicate correlation. Proxy: two random curves on the same queries (r160 vs r210, 160 vs 220).
    m160, m210 = A[:, :, j160].mean(1), A[:, :, j210].mean(1)
    rho = float(np.corrcoef(m160, m210)[0, 1])
    vq_pair = float(np.var(m210 - m160, ddof=1)) * n_val / N_T
    sens = {"rho_boot_160_210": rho, "V_Q_Delta_paired_proxy": vq_pair,
            "V_Q_Delta_rho_formula": 2 * vq_single * (1 - rho) * n_val / N_T}
    for key, vq in (("paired_proxy", vq_pair), ("rho_formula", sens["V_Q_Delta_rho_formula"])):
        rm = tgt ** 2 - vq
        sens[f"R_s_{key}"] = (max(2, math.ceil(1.5 * sig ** 2 / rm)) if rm > 0 and g > 0 else None)
    out["sensitivity_correlated_VQ"] = sens
    out["wall_s"] = time.perf_counter() - t0
    sfx = "" if pred == PS else f"_ref_{pred}"
    (OUT / f"power_pilot{sfx}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    lines = [f"# P1 power pilot (dev split 1, R0 random, predictor {LABEL.get(pred, pred)})", "",
             f"R = {R} permutations; B = {out['B']} stratified query-bootstrap replicates (frozen ReplicateEngine, "
             f"seed `bootstrap_seed` child 'p1_power_pilot_dev1'); |D_val| = {n_val}, |T| = {N_T}. DESIGN §7.1.", "",
             "| n | R | mean AUPRC | σ_r(n) | V_Q,single(n) (traj-mean) | SD_Q single traj |", "|---|---|---|---|---|---|"]
    for n, v in per_n.items():
        lines.append(f"| {n} | {v['R']} | {v['mean_auprc']:.4f} | {v['sigma_r']:.4f} | {v['V_Q_single']:.3g} "
                     f"(SD {math.sqrt(v['V_Q_single']):.4f}) | {math.sqrt(v['V_Q_one_traj_mean']):.4f} |")
    lines += ["", f"- g = r210 − r160 = {out['r210']:.4f} − {out['r160']:.4f} = **{g:.4f}** "
              f"(paired SE {g_pair['se']:.4f}, 95% CI {g_pair['ci95'][0]:.4f} to {g_pair['ci95'][1]:.4f})",
              f"- V_Q(Δ) := 2·V_Q,single(160)·|D_val|/|T| = 2·{vq_single:.3g}·{n_val}/{N_T} = **{VQ:.3g}** "
              f"(SD {math.sqrt(VQ):.4f})",
              f"- σ_r(160) = {sig:.4f}; target SE(Δ̂(160)) ≤ g/3.5 = {tgt:.4f} (SE² {tgt ** 2:.3g})",
              f"- **{out['verdict']}**" + (f" (exact R_s {out['R_s_exact']:.2f}; SE at that R = "
                                           f"{out['SE_at_R']:.4f}; V_Q share of SE² "
                                           f"{out['fraction_SE2_from_VQ_at_Rs']:.0%})" if out.get("R_s") else ""),
              "- SE(Δ̂(160)) by R_s (R_r = 2R_s): " + ", ".join(f"R_s={k}: {v:.4f}" for k, v in
                                                            out["SE_by_R_s"].items()),
              "", "Sensitivity (NOT the frozen §7.1 rule): strategy and random arms are scored on the same queries "
              "with shared bootstrap weights, so the real V_Q(Δ) is reduced by their replicate correlation. Proxy = "
              f"two random budgets on the same queries: ρ_boot(r160, r210) = {sens['rho_boot_160_210']:.3f}; "
              f"V_Q of the paired contrast (scaled to |T|) = {sens['V_Q_Delta_paired_proxy']:.3g} "
              f"(SD {math.sqrt(sens['V_Q_Delta_paired_proxy']):.4f}); 2·V_Q,single·(1−ρ)·|D_val|/|T| = "
              f"{sens['V_Q_Delta_rho_formula']:.3g} → R_s = {sens['R_s_paired_proxy']} (paired proxy) / "
              f"{sens['R_s_rho_formula']} (ρ formula), R_r = 2R_s. A strategy whose selections differ from random "
              "decorrelates more than two nested random prefixes do, so this is optimistic; P4 re-estimates V_Q(Δ) "
              "from real strategy-vs-random replicates.",
              "", "Caveats: σ_s := σ_r (strategy trajectory variance unknown until P4); V_Q(Δ) assumes ρ = 0 "
              "(conservative); g is estimated from 10 nested-prefix trajectories (its CI above) — a smaller true g "
              "raises R steeply (R ∝ σ²/((g/3.5)² − V_Q)). R is re-estimated from P4 dev strategy runs and fixed at "
              "the freeze (may only go up)."]
    (OUT / f"power_pilot{sfx}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return out


# ============================================================================ full analysis
def n_x_table(df: pd.DataFrame, fracs=(0.8, 0.9), boot=True) -> list[dict]:
    """n_80/n_90 per (split, predictor): point = n_req(frac·PAVA(curve)[max]) on the trajectory-mean curve; P*:
    joint bootstrap (query weights + trajectories) CI."""
    rows = []
    for (split, pred), d in df.groupby(["split", "predictor"]):
        m = d.groupby("n").auprc.mean().sort_index()
        grid, curve = m.index.to_numpy(), m.to_numpy()
        pool_n = len(views.get_split(split).pool)
        if grid[-1] != pool_n:
            continue
        fit = E.pava(curve)
        base = {"split": split, "predictor": pred, "R": int(d.perm.nunique()), "auprc_pool_max_fit": float(fit[-1]),
                "auprc_pool_max_mean": float(curve[-1]), "pool_n": int(pool_n), "grid": "|".join(map(str, grid))}
        for f in fracs:
            nr = E.n_req(f * fit[-1], grid, curve)
            base[f"n{int(f * 100)}"] = nr.n
            base[f"n{int(f * 100)}_censored"] = nr.censored
        if boot and pred in (PS, REF) and len(grid) >= 10:
            eng, _ = replicate_engine(df, split, pred, f"p1_q1_nx_{split}")
            M = E.joint_resample(eng.A["R0"], np.random.default_rng(p1.boot_seed(f"p1_q1_nx_resample_{split}")))
            for f in fracs:
                vals = np.array([E.n_req(f * E.pava(M[b])[-1], grid, M[b]).n for b in range(len(M))])
                lo, hi = E.percentile_ci(vals)
                base[f"n{int(f * 100)}_ci95"] = [lo, hi]
                base[f"n{int(f * 100)}_n_censored"] = int((~np.isfinite(vals) | (vals == 0)).sum())
            base["auprc_pool_max_ci95"] = list(E.percentile_ci(M[:, -1]))
        rows.append(base)
    return rows


def sat_fit(df: pd.DataFrame) -> list[dict]:
    """a − b·n^−c per permutation (descriptive), summarised by a seed-cluster (permutation) bootstrap."""
    from scipy.optimize import curve_fit
    f = lambda n, a, b, c: a - b * n ** (-c)  # noqa: E731
    out = []
    for (split, pred), d in df.groupby(["split", "predictor"]):
        if pred not in (PS, REF):
            continue
        pars = []
        for r, dr in d.groupby("perm"):
            dr = dr.sort_values("n")
            try:
                p, _ = curve_fit(f, dr.n.to_numpy(float), dr.auprc.to_numpy(), p0=(0.08, 0.2, 0.5),
                                 bounds=([0, 0, 0.01], [1, 10, 3]), maxfev=20000)
                pars.append(p)
            except RuntimeError:
                continue
        P = np.array(pars)
        rng = np.random.default_rng(p1.boot_seed(f"p1_satfit_{split}"))
        bs = np.array([P[rng.integers(0, len(P), len(P))].mean(0) for _ in range(2000)])
        out.append({"split": split, "predictor": pred, "n_fits": len(P),
                    **{k: float(P[:, i].mean()) for i, k in enumerate("abc")},
                    **{f"{k}_ci95": [float(np.quantile(bs[:, i], .025)), float(np.quantile(bs[:, i], .975))]
                       for i, k in enumerate("abc")}})
    return out


def fmt(m, s=None, d=4):
    if m is None or (isinstance(m, float) and not np.isfinite(m)):
        return "—"
    return f"{m:.{d}f}" + (f" ± {s:.{d}f}" if s is not None and np.isfinite(s) else "")


def figures(df: pd.DataFrame, sfx: str = ""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#8a8984", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": .6})
    splits = sorted(df.split.unique())
    fig, axes = plt.subplots(1, len(splits), figsize=(4.2 * len(splits), 3.8), sharey=True, squeeze=False)
    for ax, sp in zip(axes[0], splits):
        for pred in PRED_ORDER:
            d = df[(df.split == sp) & (df.predictor == pred)]
            if not len(d):
                continue
            g = d.groupby("n").auprc.agg(["mean", "std", "count"]).sort_index()
            ax.fill_between(g.index, g["mean"] - g["std"].fillna(0), g["mean"] + g["std"].fillna(0),
                            color=COLORS[pred], alpha=.15, lw=0)
            ax.plot(g.index, g["mean"], color=COLORS[pred], lw=2, marker=MARK[pred], ms=4,
                    label=f"{LABEL[pred]} (R={int(d.perm.nunique())})")
        prev = df[df.split == sp].prevalence.iloc[0]
        ax.axhline(prev, color="#8a8984", lw=1, ls=":")
        ax.text(11, prev * 1.04, "prevalence", color="#52514e", fontsize=7)
        ax.set_xscale("log")
        ax.set_title(f"{sp} (pool {len(views.get_split(sp).pool)}; E_fixed_dev {df[df.split == sp].n_cells.iloc[0]:,} "
                     "cells)", fontsize=9)
        ax.set_xlabel("screened query genes n (log)")
    axes[0][0].set_ylabel("pooled AUPRC on E_fixed_dev (mean ± SD)")
    hl = {}
    for ax in axes[0]:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            hl.setdefault(lab.split(" (R=")[0], (h, lab))
    order = [hl[LABEL[p]] for p in PRED_ORDER if LABEL[p] in hl]
    axes[0][0].legend([h for h, _ in order], [lab for _, lab in order], fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("P1 dev random curves (R0)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / f"fig_q1_dev_curves{sfx}.png", dpi=150)
    plt.close(fig)
    # per-stratum panel for P* (mean over permutations; one line per split)
    strata = [("auprc", "all E_fixed_dev"), ("auprc_AO", "E_AO (array-only columns)"),
              ("auprc_VV", "E_VV (D_val on-array columns)"), ("auprc_on_array_rows", "on-array val genes (rows)"),
              ("auprc_off_array_rows", "off-array val genes (rows)")]
    fig, axes = plt.subplots(1, len(strata), figsize=(3.3 * len(strata), 3.4), sharey=True)
    ls = {"dev1": "-", "dev2": "--", "dev3": ":"}
    for ax, (col, title) in zip(axes, strata):
        for pred in (PS, REF, "stack_lr", "B1_array_rate"):
            for sp in splits:
                d = df[(df.split == sp) & (df.predictor == pred)]
                g = d.groupby("n")[col].mean().sort_index()
                ax.plot(g.index, g.values, color=COLORS[pred], lw=1.6, ls=ls.get(sp, "-"),
                        label=f"{LABEL[pred]} {sp}" if col == "auprc" else None)
        ax.set_xscale("log")
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("n (log)")
    axes[0].set_ylabel("pooled AUPRC (mean over permutations)")
    axes[0].legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / f"fig_q1_dev_strata{sfx}.png", dpi=150)
    plt.close(fig)


def full(a) -> dict:
    df = curves()
    interim = a.mode == "interim"
    sfx = "_interim" if interim else ""
    if interim:
        df = df[df.n <= int(config.budgets()["q2_max_n"])]
    res = {"n_rows": int(len(df))}
    # ---- curve table: mean, SD per split x predictor x n, and the across-split summary
    tab = (df.groupby(["split", "predictor", "n"])
           .agg(R=("perm", "nunique"), auprc_mean=("auprc", "mean"), auprc_sd=("auprc", "std"),
                lift_mean=("lift", "mean"), macro_mean=("macro_auprc", "mean"), auroc_mean=("auroc", "mean"),
                auprc_AO=("auprc_AO", "mean"), auprc_VV=("auprc_VV", "mean"),
                auprc_on_rows=("auprc_on_array_rows", "mean"), auprc_off_rows=("auprc_off_array_rows", "mean"),
                prec_at_r10=("prec_at_recall10", "mean"), sl_detected_frac=("sl_detected_frac_dev", "mean"),
                fit_s=("fit_predict_s", "mean"), fit_s_max=("fit_predict_s", "max"),
                rss_gb=("peak_rss_fit_gb", "max"))
           .reset_index())
    tab.to_csv(OUT / f"curve_table{sfx}.csv", index=False)
    nx, sat = [], []
    if not interim:
        nx = n_x_table(df[df.predictor != "m1_hgb"])
        nx += n_x_table(df[df.predictor == "m1_hgb"], boot=False)
        pd.DataFrame(nx).to_csv(OUT / "n80_n90.csv", index=False)
        sat = sat_fit(df)
    res.update(n_x=nx, saturating_fit=sat)
    figures(df, sfx)
    L = ["# P1 — dev random curves (R0), DEV DATA ONLY", "",
         "DESIGN v3.1 §7.1/§8. Splits dev1–3 (D_val = inner fold k, candidates = D \\ fold k); evaluation set "
         "E_fixed_dev = D_val × (A \\ D_cand \\ T) minus identity exclusions (no test-gene columns). Pooled AUPRC "
         "(tie-grouped AP, frozen `hsl.eval`). Trajectories: seeded permutations of the candidate pool "
         "(`dev_traj_base_seed`, harness R0 streams), nested prefixes; model seed = trajectory evaluation stream "
         "(at the pool maximum all permutations coincide, so the SD there is model-seed variation only). "
         "Rows: `results/p1/curves.csv`; per-(split, predictor, n) table: `results/p1/curve_table.csv`.", "",
         f"Figures: `fig_q1_dev_curves{sfx}.png` (AUPRC vs n, mean ± SD band), `fig_q1_dev_strata{sfx}.png` "
         "(strata).", ""]
    has_max = all(((df.split == sp) & (df.predictor == PS) & (df.n == len(views.get_split(sp).pool))).any()
                  for sp in df.split.unique())
    if not interim and not has_max:
        L[0] = "# P1 — dev random curves (R0), DEV DATA ONLY — PROVISIONAL (pool maximum pending)"
        L += ["**Provisional:** the n ≥ 1970 phase (batch6: 1970 and the dev-pool maximum, run alone) is held by the PI, "
              "so the 'max' column below is the largest budget run so far (1640) and n_80/n_90 (anchored at the pool "
              "maximum) are not yet available.", ""]
    if interim:
        L[0] = "# P1 — INTERIM: dev random curves n ≤ 660 (R0), DEV DATA ONLY"
        L += ["**Interim (PI decision D6):** n ≥ 820, the pool maximum, m1_hgb at 1310/max and universe subsampling "
              "are on hold until P2 is done, so n_80/n_90 (which need the pool-max anchor) are not reported here.", ""]
    # key-n table (mean ± SD over permutations; per split and the mean over splits)
    L += ["## AUPRC at key budgets (mean ± SD over permutations; E_fixed_dev)", ""]
    for sp in sorted(df.split.unique()):
        pool_n = len(views.get_split(sp).pool)
        prev = df[df.split == sp].prevalence.iloc[0]
        ns = [n for n in KEY_N if n <= int(df.n.max())] if interim else [n for n in KEY_N if n < pool_n] + [pool_n]
        L += [f"**{sp}** (prevalence {prev:.4f}, pool max {pool_n})", "",
              "| predictor | R | " + " | ".join(str(n) for n in ns) + " |", "|---|---|" + "---|" * len(ns)]
        for pred in PRED_ORDER:
            t = tab[(tab.split == sp) & (tab.predictor == pred)].set_index("n")
            if not len(t):
                continue
            L.append(f"| {LABEL[pred]} | {int(t.R.max())} | " + " | ".join(
                fmt(t.auprc_mean.get(n), t.auprc_sd.get(n)) if n in t.index else "—" for n in ns) + " |")
        L.append("")
    # mean over splits (pool max column = each split's own max)
    kn = [n for n in KEY_N if n <= int(df.n.max())]
    L += ["**Mean over splits** (mean of the per-split means; ± SD of the 3 split means; pool max = each split's)",
          "", "| predictor | " + " | ".join(str(n) for n in kn) + (" |" if interim else " | max |"),
          "|---|" + "---|" * (len(kn) + (0 if interim else 1))]
    for pred in PRED_ORDER:
        t = tab[tab.predictor == pred]
        if t.split.nunique() < 2:
            continue
        cells = []
        for n in kn:
            v = t[t.n == n].auprc_mean
            cells.append(fmt(v.mean(), v.std(ddof=1) if len(v) > 1 else None))
        if not interim:
            mx = t.loc[t.groupby("split").n.idxmax()].auprc_mean
            cells.append(fmt(mx.mean(), mx.std(ddof=1)))
        L.append(f"| {LABEL[pred]} | " + " | ".join(cells) + " |")
    L.append("")
    # paired P* vs comparators (same trajectories): relative gain, split-averaged
    pw = df.pivot_table(index=["split", "perm", "n"], columns="predictor", values="auprc")
    L += ["**P* vs baselines on the same trajectories** (mean over splits of the per-split mean relative gain; "
          "range over splits):", "", "| n | P* vs pre-P2 P* | P* vs stack_lr | P* vs B1 | P* vs B2_product |",
          "|---|---|---|---|---|"]
    for n in kn:
        cells = []
        for b in (REF, "stack_lr", "B1_array_rate", "B2_product"):
            if PS not in pw or b not in pw:
                cells.append("—")
                continue
            g = ((pw[PS] - pw[b]) / pw[b]).xs(n, level="n").dropna().groupby(level="split").mean()
            if not len(g):
                cells.append("—")
                continue
            cells.append(f"{g.mean():+.1%} [{g.min():+.1%}, {g.max():+.1%}]")
        L.append(f"| {n} | " + " | ".join(cells) + " |")
    L.append("")
    # strata for P*
    L += ["## Strata (P*, mean over permutations, split-averaged)", "",
          "| n | all | E_AO | E_VV | on-array val rows | off-array val rows | lift | macro AUPRC | AUROC | "
          "prec@10% recall | SLs detected (frac of dev rows' SLs) |", "|---|" + "---|" * 10]
    tp = tab[tab.predictor == PS].groupby("n").mean(numeric_only=True)
    for n in [n for n in KEY_N if n in tp.index]:
        r = tp.loc[n]
        L.append(f"| {n} | {r.auprc_mean:.4f} | {r.auprc_AO:.4f} | {r.auprc_VV:.4f} | {r.auprc_on_rows:.4f} | "
                 f"{r.auprc_off_rows:.4f} | {r.lift_mean:.2f} | {r.macro_mean:.4f} | {r.auroc_mean:.3f} | "
                 f"{r.prec_at_r10:.3f} | {r.sl_detected_frac:.3f} |")
    L.append("")
    # n_80 / n_90
    L += [] if interim else ["## n_80 / n_90 (budget at which the PAVA, log-n interpolated mean random curve reaches 80/90% of its value "
          "at the dev pool maximum; frozen `eval.n_req`; P*: joint bootstrap 95% CI, B = 2,000)", "",
          "| split | predictor | R | AUPRC at pool max (PAVA) | n_80 [95% CI] | n_90 [95% CI] |", "|---|---|---|---|---|---|"]
    for r in sorted(nx, key=lambda r: (PRED_ORDER.index(r["predictor"]), r["split"])):
        def c(k):
            v = r.get(f"{k}_ci95")
            return f"{r[k]:.0f}" + (f" [{v[0]:.0f}, {v[1]:.0f}]" if v else "") + \
                (f" ({r[k + '_censored']}-censored)" if r.get(k + "_censored") else "")
        ci = r.get("auprc_pool_max_ci95")
        L.append(f"| {r['split']} | {LABEL[r['predictor']]} | {r['R']} | {r['auprc_pool_max_fit']:.4f}"
                 + (f" [{ci[0]:.4f}, {ci[1]:.4f}]" if ci else "") + f" | {c('n80')} | {c('n90')} |")
    ps = [r for r in nx if r["predictor"] == PS]
    if ps:
        L += ["", f"P* mean over splits: n_80 = {np.mean([r['n80'] for r in ps]):.0f}, n_90 = "
              f"{np.mean([r['n90'] for r in ps]):.0f} (pool ≈ 2,100; as fractions of the pool: "
              f"{np.mean([r['n80'] / r['pool_n'] for r in ps]):.3f} / {np.mean([r['n90'] / r['pool_n'] for r in ps]):.3f})."]
    L.append("")
    # saturating fit
    if sat:
        L += ["## Saturating fit a − b·n^(−c) per permutation (P*, descriptive; permutation-bootstrap 95% CI)", "",
              "| split | fits | a | b | c |", "|---|---|---|---|---|"]
        for s in sat:
            L.append(f"| {s['split']} | {s['n_fits']} | " + " | ".join(
                f"{s[k]:.4f} [{s[k + '_ci95'][0]:.4f}, {s[k + '_ci95'][1]:.4f}]" for k in "abc") + " |")
        L.append("")
    # runtime
    L += ["## Runtime per fit (fit+predict s, mean over splits/permutations; peak RSS sampled during the fit, max)",
          "", "| n | " + " | ".join(LABEL[p] for p in PRED_ORDER) + " |", "|---|" + "---|" * len(PRED_ORDER)]
    rt = df.groupby(["predictor", "n"]).agg(s=("fit_predict_s", "mean"), rss=("peak_rss_fit_gb", "max"),
                                            th=("n_threads", "max"))
    ns_all = sorted(df.n.unique())
    for n in [n for n in ns_all if n in (10, 30, 100, 160, 330, 660, 980, 1310, 1640, 1970) or n > 2000]:
        L.append(f"| {int(n)} | " + " | ".join(
            (f"{rt.loc[(p, n)].s:.1f} s ({rt.loc[(p, n)].rss:.1f} GB, {int(rt.loc[(p, n)].th)} thr)"
             if (p, n) in rt.index else "—") for p in PRED_ORDER) + " |")
    L.append("")
    L += ["", "Timings are wall-clock under contention (≤ 4 P1 jobs × 2 threads, n ≤ 660; plus other agents' jobs on "
          "the same 8-core M2; load average 38–120 at times), so they are upper bounds; single-job timings "
          "(3 threads): docs/recon.md §6 (stack_lr_wr2 12–25 s at n ≤ 650)."]
    tot = df.groupby("predictor").fit_predict_s.sum() / 3600
    L.append("Total fit+predict hours by predictor: " + ", ".join(f"{LABEL[p]} {tot.get(p, 0):.2f} h"
                                                                  for p in PRED_ORDER) + ".")
    L.append("")
    # upstream-protocol sanity
    up_path = OUT / "upstream_protocol.csv"
    if up_path.exists():
        up = pd.read_csv(up_path)
        up = up[(up.n == 660) & (up.split == "dev1")]      # dev2/dev3 diagnostic rows: upstream_protocol_diag.md
        upreds = [m for m in (PS, REF, "stack_lr", "B1_array_rate") if m in set(up.predictor)]
        cols = {"E_fixed_dev.auprc": "E_fixed_dev", "E_upstream_dev.auprc": "all measured D_val cells (upstream)",
                "E_upstream_dev.auprc_recip_available": "reciprocal available",
                "E_upstream_dev.auprc_x_on_array_no_recip": "on-array, no recip.",
                "E_upstream_dev.auprc_x_off_array": "off-array", "E_upstream_dev.lift": "lift (upstream set)",
                "E_upstream_dev.prevalence": "prevalence (upstream set)"}
        ref = {"stack_lr_wr2": {"E_upstream_dev.auprc": "0.0754 (lift 6.6)",
                                "E_upstream_dev.auprc_recip_available": "0.126",
                                "E_upstream_dev.auprc_x_on_array_no_recip": "0.073",
                                "E_upstream_dev.auprc_x_off_array": "0.058", "E_upstream_dev.prevalence": "0.0115"},
               "stack_lr": {"E_upstream_dev.auprc": "0.0590 (lift 5.1)", "E_upstream_dev.auprc_recip_available":
                            "0.117", "E_upstream_dev.auprc_x_on_array_no_recip": "0.051",
                            "E_upstream_dev.auprc_x_off_array": "0.047"},
               "B1_array_rate": {"E_upstream_dev.auprc": "0.0430 (lift 3.7)",
                                 "E_upstream_dev.auprc_recip_available": "0.038",
                                 "E_upstream_dev.auprc_x_on_array_no_recip": "0.043",
                                 "E_upstream_dev.auprc_x_off_array": "0.048"}}
        L += ["## Upstream-protocol sanity (dev split 1, n = 660 ≈ 20% of 3,283; R = 10 permutations; mean ± SD)", "",
              "Upstream protocol = every measured cell of the held-out rows incl. cells whose reciprocal is known "
              "(array screened); here minus test-gene columns (dev rule) and self pairs. Upstream reference = "
              "YEAST_SLs REPORT §3 (S1 dev, 20% screened, ~2,600 held-out rows, seeds 0–4).", "",
              "| set / stratum | " + " | ".join(f"{LABEL[m]} ours | upstream" for m in upreds) + " |",
              "|---|" + "---|---|" * len(upreds)]
        res["upstream_protocol"] = {}
        for k, lab in cols.items():
            cells = []
            for m in upreds:
                d = up[up.predictor == m]
                cells += [fmt(d[k].mean(), d[k].std(ddof=1)) if k in d else "—", ref.get(m, {}).get(k, "—")]
                if k in d:
                    res["upstream_protocol"].setdefault(m, {})[k] = [float(d[k].mean()), float(d[k].std(ddof=1))]
            L.append(f"| {lab} | " + " | ".join(cells) + " |")
        if {"stack_lr_wr2", "stack_lr"} <= set(upreds):
            a_ = up[up.predictor == "stack_lr_wr2"].set_index("perm")
            b_ = up[up.predictor == "stack_lr"].set_index("perm")
            rat = (a_["E_upstream_dev.auprc"] / b_["E_upstream_dev.auprc"]).dropna()
            ratf = (a_["E_fixed_dev.auprc"] / b_["E_fixed_dev.auprc"]).dropna()
            L += ["", f"Ratio P*/stack_lr (paired by permutation): upstream-protocol set {rat.mean():.3f} ± "
                  f"{rat.std(ddof=1):.3f}, E_fixed_dev {ratf.mean():.3f} ± {ratf.std(ddof=1):.3f}; upstream dev "
                  "1.28 (REPORT §4.3)."]
            res["upstream_protocol"]["ratio_pstar_stack_lr"] = [float(rat.mean()), float(rat.std(ddof=1))]
        L += ["", "Gap diagnosis (`upstream_protocol_diag.md`): our harness reproduces upstream's own dev split 0 "
              "(P* 0.0774 vs 0.0775); dev1's target rows are the hardest fold (upstream's own predictions on dev1 rows: "
              "0.0673); dev2/dev3 upstream-protocol P* = 0.0796 / 0.0757 at n = 660."]
        d0 = up[up.predictor == upreds[0]]
        cnt = {k: d0[k].iloc[0] for k in d0.columns if k.startswith("E_upstream_dev.n_cells")}
        L += ["", "Cells (perm 0): " + ", ".join(f"{k.split('n_cells_')[1]} {int(v):,}" for k, v in cnt.items()
                                                 if "n_cells_" in k) +
              f"; set {int(d0['E_upstream_dev.n_cells'].iloc[0]):,} cells.", ""]
    L += ["## End-to-end controls (cited, not re-run)", "",
          "`tests/test_harness.py::test_end_to_end_controls_dev1_n100` (slow; P0): shuffled screened labels (ε, P, "
          "labels permuted jointly over the screened block, array SMF permuted) → AUPRC within 0.7–1.3× prevalence; "
          "oracle back door → AUPRC > 0.99, both on dev1 E_fixed_dev through `harness.run_predictor`. P1 code paths "
          "use the same `run_predictor`; every P1 engine AP (unit weights) equals the stored harness AUPRC to 1e-9.",
          ""]
    if (OUT / "power_pilot.md").exists():
        pp = json.loads((OUT / "power_pilot.json").read_text())
        L += ["## Power pilot", "", f"See `power_pilot.md`: g = {pp['g']:.4f}, V_Q(Δ) = {pp['V_Q_Delta']:.3g}, "
              f"σ_r(160) = {pp['sigma_r_160']:.4f} → {pp['verdict']}.", ""]
    (OUT / f"summary{sfx}.md").write_text("\n".join(L) + "\n")
    (OUT / f"summary{sfx}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("\n".join(L))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("mode", choices=["pilot", "full", "interim"])
    ap.add_argument("--predictor", default=PS, help="pilot: predictor (default frozen P*)")
    a = ap.parse_args(argv)
    (pilot if a.mode == "pilot" else full)(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
