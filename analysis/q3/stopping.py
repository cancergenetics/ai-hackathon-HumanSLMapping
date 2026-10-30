"""Q3 analysis 5 — online stopping estimator and plateau stopping rule (DESIGN v3.1 §7.4).

    .venv/bin/python -m analysis.q3.stopping evaluate --source dev     estimator bias/MAE/correlation + rule grid
    .venv/bin/python -m analysis.q3.stopping freeze --source dev       writes analysis/q3/stopping_rule.yaml (once)
    .venv/bin/python -m analysis.q3.stopping apply --source final      applies the FROZEN rule (never re-tuned)

Estimator: leave-batch-out CV AUPRC of P* on screened rows only (analysis.q3.lbo; LBO cells never include
targets or T). Truth: the trajectory's own held-out AUPRC on E_fixed(_dev) at the same n.
  bias(n) = mean(cv − true), MAE(n) = mean |cv − true|, relative versions / true; within-n correlation of cv with
  true across trajectories (Pearson and Spearman; raw and split-centred).
True knee (per trajectory): smallest n with AUPRC ≥ knee_frac (protocol stats.knee_frac = 0.95) × AUPRC(660),
  on the trajectory's PAVA curve over the Q2 checkpoints, log-n interpolated (frozen hsl.eval.n_req); grid knee =
  first checkpoint reaching it.
Plateau rule (family): at checkpoints n_1 < ... < n_K (Q2 grid), g_k = (cv_k / cv_{k-1} − 1) / log2(n_k / n_{k-1})
  = relative CV gain per doubling of the screened set. Stop at the first n_k ≥ n_min where the last m values of g
  are all < δ; if none, stop at the last checkpoint and flag 'not triggered'. Grid δ × m; the frozen choice minimises
  the mean |log2(n_stop / n_knee)| over all dev trajectories (ties -> smaller δ, i.e. later stopping); a
  leave-split-out check (choose on two splits, score on the third) is reported.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import time

import numpy as np
import pandas as pd
import yaml

from analysis.q3 import common as C

RULE_FILE = C.Q3 / "stopping_rule.yaml"
DELTAS = (0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30)
MS = (1, 2)
# statistic the plateau rule watches: raw LBO-CV AUPRC, or LBO-CV lift (AUPRC / CV-cell prevalence; removes the drift
# in screened-row SL density that adaptive selection causes — screened rows get less SL-rich as n grows)
STATS = ("cv_auprc", "cv_lift")
N_MIN = 60


def load_lbo(O) -> pd.DataFrame:
    p = O / "lbo" / "lbo_rows.jsonl"
    d = pd.DataFrame([json.loads(x) for x in p.read_text().splitlines() if x.strip()])
    return d.drop(columns=["fits"], errors="ignore")


def merged(source: str, O) -> pd.DataFrame:
    lbo = load_lbo(O)
    rows = C.rows(source)
    keep = ["arm", "split", "traj", "n", "auprc"] + (["prevalence"] if "prevalence" in rows.columns else [])
    m = lbo.merge(rows[keep].rename(columns={"prevalence": "true_prevalence"}), on=["arm", "split", "traj", "n"],
                  how="left")
    m = m.rename(columns={"auprc": "true_auprc"})
    if "true_prevalence" in m.columns:
        m["true_lift"] = m.true_auprc / m.true_prevalence
    return m


def knee(ns, true, frac):
    from hsl import eval as E
    ns, true = np.asarray(ns, float), np.asarray(true, float)
    f = E.pava(true)
    tgt = frac * f[-1]
    r = E.n_req(tgt, ns, true)
    grid = float(ns[np.argmax(f >= tgt)])
    return (float(r.n) if r.censored == "" else float(ns[0])), grid


def plateau_stop(ns, cv, delta, m, n_min=N_MIN) -> tuple[float, bool]:
    ns, cv = np.asarray(ns, float), np.asarray(cv, float)
    g = np.full(len(ns), np.nan)
    g[1:] = (cv[1:] / cv[:-1] - 1) / np.log2(ns[1:] / ns[:-1])
    for k in range(1, len(ns)):
        if ns[k] < n_min or k < m:
            continue
        if np.all(g[k - m + 1:k + 1] < delta):
            return float(ns[k]), True
    return float(ns[-1]), False


def projection_stop(ns, x, tau, w=3, n_max=660.0, n_min=N_MIN) -> tuple[float, bool]:
    """Stop at the first checkpoint n_k >= n_min (k >= w-1) whose value x_k >= tau * x_hat(n_max), x_hat = OLS line of
    x on log n through the last w checkpoints (incl. k), extrapolated to the planned maximum budget n_max."""
    ns, x = np.asarray(ns, float), np.asarray(x, float)
    for k in range(w - 1, len(ns)):
        if ns[k] < n_min:
            continue
        if ns[k] >= n_max:
            return float(ns[k]), False
        lx, yy = np.log(ns[k - w + 1:k + 1]), x[k - w + 1:k + 1]
        sl = np.polyfit(lx, yy, 1)
        proj = max(np.polyval(sl, np.log(n_max)), x[k])
        if x[k] >= tau * proj:
            return float(ns[k]), True
    return float(ns[-1]), False


def stop_n(ns, x, family, param, m, n_min=N_MIN):
    if family == "plateau":
        return plateau_stop(ns, x, param, m, n_min)
    if family == "projection":
        return projection_stop(ns, x, param, w=m, n_min=n_min)
    raise KeyError(family)


TAUS = (0.80, 0.85, 0.90, 0.93, 0.95, 0.97)
WS = (3, 4)
FAMILIES = [("plateau", d, m) for d in DELTAS for m in MS] + [("projection", t, w) for t in TAUS for w in WS]


def per_traj(M: pd.DataFrame, frac: float) -> list[dict]:
    out = []
    for (arm, sp, r), g in M.groupby(["arm", "split", "traj"]):
        g = g.sort_values("n")
        if not np.isfinite(g.true_auprc).all() or g.n.max() < 660:
            continue
        kn, kg = knee(g.n, g.true_auprc, frac)
        out.append({"arm": arm, "split": sp, "traj": r, "ns": g.n.to_numpy(), "cv": g.cv_auprc.to_numpy(),
                    "cv_auprc": g.cv_auprc.to_numpy(), "cv_lift": g.cv_lift.to_numpy(),
                    "true": g.true_auprc.to_numpy(), "knee": kn, "knee_grid": kg})
    return out


def grid_eval(T: list[dict]) -> pd.DataFrame:
    rec = []
    for st, (fam, d, m) in [(a, f) for a in STATS for f in FAMILIES]:
            for t in T:
                ns_, trig = stop_n(t["ns"], t[st], fam, d, m)
                rec.append({"family": fam, "statistic": st, "delta": d, "m": m, "arm": t["arm"], "split": t["split"], "traj": t["traj"],
                            "n_stop": ns_, "triggered": trig, "knee": t["knee"], "knee_grid": t["knee_grid"],
                            "log2_err": np.log2(ns_ / t["knee"]),
                            "true_at_stop_frac": float(np.interp(np.log(ns_), np.log(t["ns"]), t["true"]) /
                                                       t["true"][-1])})
    return pd.DataFrame(rec)


def choose(G: pd.DataFrame) -> pd.Series:
    s = G.assign(abs_err=G.log2_err.abs()).groupby(["family", "statistic", "delta", "m"]).agg(
        mae_log2=("abs_err", "mean"), bias_log2=("log2_err", "mean"), trig=("triggered", "mean"),
        stop_frac=("true_at_stop_frac", "mean")).reset_index()
    s = s.sort_values(["mae_log2", "family", "statistic", "delta", "m"])
    return s.iloc[0], s


def evaluate(source: str, out=None) -> dict:
    t0 = time.perf_counter()
    O = C.out_dir(source, out)
    frac = float(C.protocol()["stats"]["knee_frac"])
    M = merged(source, O)
    M["err"] = M.cv_auprc - M.true_auprc
    if "true_lift" in M.columns:
        M["lift_rel_err"] = M.cv_lift / M.true_lift - 1
    M["rel_err"] = M.err / M.true_auprc
    M["cv_c"] = M.cv_auprc - M.groupby(["split", "n"]).cv_auprc.transform("mean")
    M["true_c"] = M.true_auprc - M.groupby(["split", "n"]).true_auprc.transform("mean")
    est = []
    for (n), g in M.groupby("n"):
        rec = {"n": n, "R": len(g), "cv_mean": g.cv_auprc.mean(), "true_mean": g.true_auprc.mean(),
               "bias": g.err.mean(), "mae": g.err.abs().mean(), "rel_bias": g.rel_err.mean(),
               "rel_mae": g.rel_err.abs().mean(), "cv_prevalence": g.prevalence.mean(),
               "fit_s_mean": g.fit_s_total.mean(), "lift_rel_bias": g.get("lift_rel_err", pd.Series(np.nan)).mean(),
               "lift_rel_mae": g.get("lift_rel_err", pd.Series(np.nan)).abs().mean()}
        if len(g) >= 4:
            rec.update({"pearson_lift": g.cv_lift.corr(g.true_lift) if "true_lift" in g else np.nan,
                        "pearson": g.cv_auprc.corr(g.true_auprc), "spearman": g.cv_auprc.corr(g.true_auprc, "spearman"),
                        "pearson_split_centred": g.cv_c.corr(g.true_c)})
        est.append(rec)
    for (arm, n), g in M.groupby(["arm", "n"]):
        est.append({"arm": arm, "n": n, "R": len(g), "cv_mean": g.cv_auprc.mean(), "true_mean": g.true_auprc.mean(),
                    "bias": g.err.mean(), "mae": g.err.abs().mean(), "rel_bias": g.rel_err.mean(),
                    "rel_mae": g.rel_err.abs().mean(), "cv_prevalence": g.prevalence.mean(),
                    "lift_rel_bias": g.get("lift_rel_err", pd.Series(np.nan)).mean(),
                    "lift_rel_mae": g.get("lift_rel_err", pd.Series(np.nan)).abs().mean(),
                    "pearson_lift": g.cv_lift.corr(g.true_lift) if len(g) >= 4 and "true_lift" in g else np.nan,
                    "pearson": g.cv_auprc.corr(g.true_auprc) if len(g) >= 4 else np.nan,
                    "spearman": g.cv_auprc.corr(g.true_auprc, "spearman") if len(g) >= 4 else np.nan,
                    "pearson_split_centred": g.cv_c.corr(g.true_c) if len(g) >= 4 else np.nan})
    E_ = pd.DataFrame(est)
    E_["arm"] = E_["arm"].fillna("all")
    E_.to_csv(O / "stop_estimator.csv", index=False)
    M.to_csv(O / "stop_traj_rows.csv", index=False)
    T = per_traj(M, frac)
    G = grid_eval(T)
    G.to_csv(O / "stop_rule_grid_traj.csv", index=False)
    best, S = choose(G) if len(G) else (None, pd.DataFrame())
    S.to_csv(O / "stop_rule_grid.csv", index=False)
    # comparator: a FIXED stopping budget (no online information), each checkpoint
    fx = []
    for t in T:
        for n0 in t["ns"]:
            fx.append({"n_fixed": n0, "abs_log2_err": abs(np.log2(n0 / t["knee"])), "log2_err": np.log2(n0 / t["knee"])})
    FX = pd.DataFrame(fx)
    if len(FX):
        FX = FX.groupby("n_fixed").agg(mae_log2=("abs_log2_err", "mean"), bias_log2=("log2_err", "mean")).reset_index()
    FX.to_csv(O / "stop_fixed_n_baseline.csv", index=False)
    # leave-split-out check of the selection
    lso = []
    for sp in sorted(G.split.unique()) if len(G) else []:
        if (G.split != sp).sum() == 0:
            continue
        b2, _ = choose(G[G.split != sp])
        h = G[(G.split == sp) & (G.family == b2.family) & (G.statistic == b2.statistic) & (G.delta == b2.delta) &
              (G.m == b2.m)]
        lso.append({"held_out": sp, "family": b2.family, "statistic": b2.statistic, "delta": b2.delta, "m": b2.m, "mae_log2": h.log2_err.abs().mean(),
                    "bias_log2": h.log2_err.mean(), "R": len(h)})
    pd.DataFrame(lso).to_csv(O / "stop_rule_lso.csv", index=False)
    kn = pd.DataFrame([{k: v for k, v in t.items() if k not in ("ns", "cv", "true", "cv_auprc", "cv_lift")} for t in T])
    kn.to_csv(O / "stop_knees.csv", index=False)
    res = {"source": source, "n_traj_complete": len(T), "best": None if best is None else best.to_dict(),
           "best_fixed_n": None if not len(FX) else FX.sort_values("mae_log2").iloc[0].to_dict(),
           "lso": lso, "runtime_s": time.perf_counter() - t0}
    C.write_json(O / "stop_summary.json", res)
    if len(T):
        fig(M, T, best, O / "fig_stopping.png", source)
    return res


def freeze(source: str, out=None, force: bool = False) -> dict:
    if source != "dev":
        raise PermissionError("the stopping rule is frozen on dev only")
    if RULE_FILE.exists() and not force:
        raise FileExistsError(f"{RULE_FILE} exists: the rule is frozen (never re-tuned after the freeze)")
    O = C.out_dir(source, out)
    res = evaluate(source, out)
    b = res["best"]
    head = subprocess.run(["git", "-C", str(C.ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    rule = {
        "status": "frozen",
        "frozen_at": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "frozen_before_any_final_output": True,
        "git_head_at_freeze": head,
        "estimator": {"name": "LBO-CV AUPRC of P* on screened rows", "module": "analysis.q3.lbo",
                      "batch": 10, "k_folds": 5, "fold_rule": "batch index (i // 10) mod K, K = min(5, #batches)",
                      "cells": "fold rows x columns not screened, not test genes, not self, measured",
                      "metric": "pooled tie-grouped AP over all folds (hsl.eval.average_precision)",
                      "predictor": "frozen P* (config/predictor.yaml primary), trajectory model seed"},
        "checkpoints": [30, 60, 100, 160, 220, 330, 490, 660],
        "statistic": "x = rule.statistic: cv_auprc (LBO-CV AUPRC) or cv_lift (cv_auprc / prevalence of the CV "
                     "cells); plateau family uses g_k = (x_k / x_{k-1} - 1) / log2(n_k / n_{k-1})",
        "rule": {"family": str(b["family"]), "statistic": str(b["statistic"]),
                 "param": float(b["delta"]), "m_or_window": int(b["m"]), "n_min": N_MIN, "n_max_planned": 660,
                 "stop_plateau": "first checkpoint n_k >= n_min whose last m values of g are all < param (delta)",
                 "stop_projection": "first checkpoint n_k >= n_min with x_k >= param (tau) x x_hat(n_max), x_hat = OLS "
                                    "of x on log n over the last `m_or_window` checkpoints, extrapolated to n_max",
                 "if_never": "stop at the last checkpoint (660), flagged not_triggered"},
        "target_of_tuning": "true knee = smallest n with AUPRC >= 0.95 x AUPRC(660) (PAVA, log-n interpolated)",
        "selection": {"criterion": "min mean |log2(n_stop / n_knee)| over the dev trajectories (R0 and H, "
                                   "trajectories 0-1, dev1-3); ties -> family, statistic, param, m order",
                      "grid": {"statistics": list(STATS), "plateau_delta": list(DELTAS), "plateau_m": list(MS),
                               "projection_tau": list(TAUS), "projection_window": list(WS)},
                      "dev_mae_log2": float(b["mae_log2"]), "dev_bias_log2": float(b["bias_log2"]),
                      "dev_trigger_rate": float(b["trig"]), "n_dev_trajectories": res["n_traj_complete"],
                      "leave_split_out": res["lso"],
                      "comparator_best_fixed_n": res["best_fixed_n"],
                      "note": "the knee is defined relative to AUPRC(660); a fixed budget near 0.75 x 660 is a strong "
                              "comparator, reported alongside (stop_fixed_n_baseline.csv)",
                      "evidence": "results/q3_dev/stop_rule_grid.csv, stop_rule_grid_traj.csv, stop_rule_lso.csv"},
    }
    rule = json.loads(json.dumps(rule, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    RULE_FILE.write_text("# Q3 online stopping rule — FROZEN on dev before any final output was seen (DESIGN §7.4).\n"
                         "# Written by `python -m analysis.q3.stopping freeze`; never edited or re-tuned afterwards.\n"
                         + yaml.safe_dump(rule, sort_keys=False))
    return rule


def apply(source: str, out=None) -> pd.DataFrame:
    """Apply the frozen rule to LBO rows of `source` trajectories (no re-tuning)."""
    rule = yaml.safe_load(RULE_FILE.read_text())
    if rule.get("status") != "frozen":
        raise RuntimeError("stopping rule not frozen")
    O = C.out_dir(source, out)
    frac = float(C.protocol()["stats"]["knee_frac"])
    M = merged(source, O)
    M = M[M.n.isin(rule["checkpoints"])]
    T = per_traj(M, frac)
    rec = []
    for t in T:
        R_ = rule["rule"]
        ns_, trig = stop_n(t["ns"], t[R_["statistic"]], R_["family"], R_["param"], R_["m_or_window"], R_["n_min"])
        rec.append({"arm": t["arm"], "split": t["split"], "traj": t["traj"], "n_stop": ns_, "triggered": trig,
                    "knee": t["knee"], "knee_grid": t["knee_grid"], "log2_err": np.log2(ns_ / t["knee"]),
                    "true_at_stop_frac": float(np.interp(np.log(ns_), np.log(t["ns"]), t["true"]) / t["true"][-1])})
    A = pd.DataFrame(rec)
    A.to_csv(O / "stop_rule_applied.csv", index=False)
    return A


def fig(M, T, best, path, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})
    style = {"R0": ("#6f6e69", "random"), "H": ("#2a78d6", "H")}
    f, axes = plt.subplots(1, 2, figsize=(10, 4), facecolor="#fcfcfb")
    ax = axes[0]
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color="#e8e7e2", lw=0.6)
    for t in T:
        col, lab = style.get(t["arm"], ("#eb6834", t["arm"]))
        ax.plot(t["ns"], t["true"], color=col, lw=1.6, alpha=0.9)
        ax.plot(t["ns"], t["cv"], color=col, lw=1.2, ls="--", alpha=0.9)
    for a, (col, lab) in style.items():
        ax.plot([], [], color=col, lw=1.6, label=f"{lab}: held-out AUPRC")
        ax.plot([], [], color=col, lw=1.2, ls="--", label=f"{lab}: LBO-CV on screened rows")
    ax.set_xscale("log")
    ax.set_xticks([30, 60, 100, 160, 330, 660])
    ax.set_xticklabels(["30", "60", "100", "160", "330", "660"])
    ax.minorticks_off()
    ax.set_xlabel("queries screened (n)")
    ax.set_ylabel("AUPRC")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("online estimator vs truth (each line = one trajectory)", fontsize=9)
    ax = axes[1]
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color="#e8e7e2", lw=0.6)
    if best is not None:
        for t in T:
            col, _ = style.get(t["arm"], ("#eb6834", t["arm"]))
            ns_, _ = stop_n(t["ns"], t[best["statistic"]], best["family"], best["delta"], int(best["m"]))
            ax.plot(t["knee"], ns_, "o", color=col, ms=6, alpha=0.85)
        ax.plot([250, 700], [250, 700], color="#c3c2b7", lw=1)
        ax.set_xlim(250, 700)
        ax.set_ylim(250, 700)
        for a, (col, lab) in style.items():
            ax.plot([], [], "o", color=col, label=lab)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.set_xlabel("true knee n (0.95 x AUPRC(660))")
        ax.set_ylabel("stopping n of the rule (checkpoint)")
        ax.set_title(f"{best['family']} rule on {best['statistic']}: param {best['delta']}, m/w {int(best['m'])} (dev)",
                     fontsize=9)
    f.suptitle(f"Online stopping ({source})", fontsize=10)
    f.tight_layout()
    f.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("evaluate", "freeze", "apply"))
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "evaluate":
        print(json.dumps(evaluate(a.source, a.out), indent=1, default=str))
    elif a.cmd == "freeze":
        print(yaml.safe_dump(freeze(a.source, a.out), sort_keys=False))
    else:
        print(apply(a.source, a.out).to_string())


if __name__ == "__main__":
    main()
