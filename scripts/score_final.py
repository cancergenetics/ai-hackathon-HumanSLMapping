#!/usr/bin/env python
"""Score the final evaluation (DESIGN v3.1 §3, §6.4, §7.1, §7.2) from SAVED predictions — separate step after
`scripts/run_final.py start/work` (`run_final.py score` delegates here). No bundle, no HSL_FINAL: labels come from
<out>/run/labels.npz written by the `_labels` job; every number is recomputed from the prediction files with the
frozen engine (hsl.eval) and the dev-tested helpers (hsl.q2dev: welch_at, gap_closure_joint, n_req_joint).

Outputs (<out>/ = results/final or results/final_dryrun):
  summary.md, summary.json         family verdicts, G / G* / matches rule, n_req, Q1 n_80 / n_90, strata,
                                   sensitivities, descriptive arms, headline rule
  confirmatory.csv                 the <= 6 one-sided hypotheses (c1, a per finalist), raw p, Holm-adjusted p, reject
  contrasts.csv                    every Welch contrast (arm, predictor, set, n_a, n_b, Δ, CI, p, V_Q, ...)
  q2_curves.csv                    Δ(n) per arm with Bonferroni bands over the Q2 checkpoints
  n_req.csv, q1_curves.csv, q1_nx.csv, q1_satfit.csv, e_unobs.csv, discovery.csv, runtime.csv, runs.csv,
  selections.json, lock_record.json; figures fig_final_q1.png, fig_final_q2_delta.png, fig_final_curves.png,
  fig_final_frontier.png
  --cost: cost_projection.md / .json for the REAL final (plan of mode `final`) from measured runtimes.

    .venv/bin/python scripts/score_final.py --mode final|dryrun [--partial] [--cost]
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import gzip
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from hsl import config, eval as E, final as F, p1, q2dev, views

ROOT = F.ROOT
STAT = config.protocol()["stats"]


def _f(x, k=4):
    try:
        if x is not None and np.isinf(float(x)):
            return "∞" if float(x) > 0 else "−∞"
        return "—" if x is None or not np.isfinite(float(x)) else f"{float(x):.{k}f}"
    except (TypeError, ValueError):
        return str(x)


def _json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ============================================================================ scorer
class Scorer:
    def __init__(self, cfg, P, rows):
        self.cfg, self.P = cfg, P
        self.L = F.Labels(P)
        self.eng = F.BootEngine(P, self.L)
        self.rows = rows
        self.by = {}
        for x in rows:
            self.by.setdefault((x["cid"], x["pred"]), []).append(x)
        self._M = {}
        self._pt = {}
        self.fins = list(cfg["finalists"])
        self.alpha = float(STAT["alpha_familywise"])
        cav = cfg["descriptive"].get("deterministic_caveat") or {}
        self.caveat = {cid: cav[d["base"]] for cid, d in P["plan"]["arms"].items()
                       if d["role"] == "det" and d["base"] in cav}

    def has(self, cid, pred="pstar"):
        return bool(self.by.get((cid, pred)))

    def M(self, cid, pred, grid, stratum="all", two_way=False):
        k = (cid, pred, tuple(grid), stratum, two_way)
        if k not in self._M:
            self._M[k] = F.arm_matrix(self.by.get((cid, pred), []), list(grid), self.eng, stratum, two_way=two_way)
        return self._M[k]

    def welch(self, a, pa, na, b, pb, nb, stratum="all", scale_b=1.0, two_way=False) -> dict | None:
        As, es, _ = self.M(a, pa, [na], stratum, two_way)
        Ar, er, _ = self.M(b, pb, [nb], stratum, two_way)
        oks, okr = np.isfinite(es[:, 0]), np.isfinite(er[:, 0])
        if oks.sum() < 1 or okr.sum() < 1 or (oks.sum() + okr.sum()) < 3:
            return None
        w = q2dev.welch_at(As[:, :, 0], es[:, 0], Ar[:, :, 0], er[:, 0], scale_r=scale_b)
        sa, sb = float(np.nanmean(es[:, 0])), float(np.nanmean(er[:, 0]))
        return {"a": a, "pred_a": pa, "n_a": na, "b": b, "pred_b": pb, "n_b": nb, "set": stratum,
                "bootstrap": "two_way" if two_way else "query", "scale_b": scale_b, "s_a": sa, "s_b": sb,
                "sd_a": float(np.nanstd(es[:, 0], ddof=1)) if oks.sum() > 1 else float("nan"),
                "sd_b": float(np.nanstd(er[:, 0], ddof=1)) if okr.sum() > 1 else float("nan"),
                "min_a": float(np.nanmin(es[:, 0])), "max_a": float(np.nanmax(es[:, 0])),
                "min_b": float(np.nanmin(er[:, 0])), "max_b": float(np.nanmax(er[:, 0])),
                "delta": w["delta"], "rel": w["delta"] / (scale_b * sb) if sb else float("nan"),
                "lo": w["ci95"][0], "hi": w["ci95"][1], "se": w["se"], "df": w["df"], "p": w["p_one_sided"],
                "V_Q": w["V_Q"], "R_a": w["R_a"], "R_b": w["R_b"], "traj_var_missing": w["traj_var_missing"]}

    def paired(self, cid, pc, pp, n) -> dict | None:
        Ac, ec, tc = self.M(cid, pc, [n])
        Ap, ep, tp = self.M(cid, pp, [n])
        common = sorted(set(t for t, v in zip(tc, ec[:, 0]) if np.isfinite(v)) &
                        set(t for t, v in zip(tp, ep[:, 0]) if np.isfinite(v)))
        if len(common) < 2:
            return None
        ic, ip = [tc.index(t) for t in common], [tp.index(t) for t in common]
        w = E.paired(Ac[:, ic, 0], ec[ic, 0], Ap[:, ip, 0], ep[ip, 0])
        return {"arm": cid, "child": pc, "parent": pp, "n": n, "delta": w["delta"], "rel": w["rel_gain"],
                "lo": w["ci95"][0], "hi": w["ci95"][1], "p": w["p_one_sided"], "R": w["R"]}

    SEC = ("macro_auprc", "auroc", "p_at_10", "p_at_50")

    def secondary(self, x) -> dict:
        """DESIGN §3 secondary metrics of one unit: the harness values stored at run time (computed after the
        predictions were saved), or recomputed from the saved predictions for rows that predate them."""
        if all(x.get(k) is not None for k in self.SEC):
            return {k: float(x[k]) for k in self.SEC}
        m = E.summarise(self.L.y, self.eng.scores(x["pred_file"]), self.L.rows)
        return {k: float(m[k]) for k in self.SEC}

    def point(self, pred_file, stratum):
        """Unit-weight AP on a stratum (no bootstrap) for descriptive tables."""
        k = (pred_file, stratum)
        if k not in self._pt:
            cp = self.eng.cache / f"{Path(pred_file).stem}__{stratum}__B{self.P['B']}.npy"
            if cp.exists():
                self._pt[k] = float(np.load(cp)[0])
            else:
                s = self.eng.scores(pred_file)
                m = self.L.mask(stratum)
                y, s = (self.L.y, s) if m is None else (self.L.y[m], s[m])
                self._pt[k] = E.average_precision(y, s) if y.any() else float("nan")
        return self._pt[k]


# ============================================================================ sections
def family(S: Scorer, out: dict) -> pd.DataFrame:
    """Confirmatory family (DESIGN §6.4): (c1) s(160) > r(490), (a) Δ(160) > 0 per finalist; Holm, one-sided."""
    c1, na = S.P["c1"], int(config.budgets()["contrast_a"]["n"])
    hy = []
    for fid in S.fins:
        w1 = S.welch(fid, "pstar", int(c1["strategy_n"]), "R0", "pstar", int(c1["random_n"]))
        wa = S.welch(fid, "pstar", na, "R0", "pstar", na)
        for name, w in (("c1", w1), ("a", wa)):
            hy.append({"finalist": fid, "hypothesis": name,
                       "H1": f"AUPRC_{fid}({w['n_a']}) > AUPRC_R0({w['n_b']})" if w else "",
                       **({k: w[k] for k in ("s_a", "s_b", "delta", "lo", "hi", "se", "df", "p", "V_Q", "R_a", "R_b")}
                          if w else {"p": float("nan")})})
    df = pd.DataFrame(hy)
    ok = np.isfinite(df.p.to_numpy(float))
    h = E.holm(df.p.to_numpy(float)[ok], S.alpha)
    df["p_holm"], df["reject"] = np.nan, False
    df.loc[ok, "p_holm"] = h["p_adj"]
    df.loc[ok, "reject"] = h["reject"]
    df["alpha_familywise"] = S.alpha
    df["m_hypotheses"] = int(ok.sum())
    out["family"] = df.to_dict("records")
    return df


def brief_and_g(S: Scorer, fam: pd.DataFrame, out: dict) -> list[dict]:
    """s160 − r660, gap closure G (joint bootstrap), G*, 'matches' linear contrast and its 3-condition rule."""
    b = S.P["brief"]
    ns, nr = int(b["strategy_n"]), int(b["random_n"])
    m = float(STAT["match_margin_frac"])
    rec = []
    Ar, er, _ = S.M("R0", "pstar", [ns, nr])
    for fid in S.fins:
        As, es, _ = S.M(fid, "pstar", [ns])
        br = S.welch(fid, "pstar", ns, "R0", "pstar", nr)
        mt = S.welch(fid, "pstar", ns, "R0", "pstar", nr, scale_b=1.0 - m)
        G = q2dev.gap_closure_joint([(As, Ar, es, er)], np.random.default_rng(F.boot_seed(S.P, f"G_{fid}")))
        r160, r660 = float(np.nanmean(er[:, 0])), float(np.nanmean(er[:, 1]))
        gs = float(E.g_star(r160, r660))
        c1_rej = bool(fam[(fam.finalist == fid) & (fam.hypothesis == "c1")].reject.any())
        mt_rej = bool(mt and mt["p"] <= float(STAT["alpha_familywise"]))
        allowed = mt_rej and c1_rej and gs >= float(STAT["g_star_min"])
        rec.append({"finalist": fid, "brief_delta": br and br["delta"], "brief_lo": br and br["lo"],
                    "brief_hi": br and br["hi"], "brief_p": br and br["p"],
                    "matches_delta": mt and mt["delta"], "matches_lo": mt and mt["lo"], "matches_hi": mt and mt["hi"],
                    "matches_p": mt and mt["p"], "matches_rejected": mt_rej, "c1_rejected_holm": c1_rej,
                    "G": G["G"], "G_lo": G["G_ci"][0], "G_hi": G["G_ci"][1], "G_defined": G["defined"],
                    "G_frac_nonfinite": G["frac_nonfinite"],
                    "gap_r660_r160": G["den"], "gap_ci": G["den_ci"], "G_star": gs,
                    "G_star_discriminating": gs >= float(STAT["g_star_min"]), "matches_claim_allowed": allowed})
    out["brief_G"] = rec
    return rec


def delta_curves(S: Scorer, contrasts: list, out: dict) -> pd.DataFrame:
    """Δ(n) vs random at matched n for every strategy arm (P*), Bonferroni bands over each arm's checkpoints."""
    roles = {cid: d["role"] for cid, d in S.P["plan"]["arms"].items()}
    recs = []
    for cid, role in roles.items():
        if role in ("random", "ws", "ws_random") or not S.has(cid):
            continue
        ns = sorted({int(x["n"]) for x in S.by[(cid, "pstar")]} & set(S.P["random_grid"]))
        ws = [S.welch(cid, "pstar", n, "R0", "pstar", n) for n in ns]
        ws = [(n, w) for n, w in zip(ns, ws) if w]
        if not ws:
            continue
        bands = E.bonferroni_bands([w["delta"] for _, w in ws], [w["se"] for _, w in ws], [w["df"] for _, w in ws])
        for k, (n, w) in enumerate(ws):
            contrasts.append({**w, "family": "descriptive" if role != "finalist" else "q2_curve",
                              "caveat": S.caveat.get(cid, "")})
            recs.append({"arm": cid, "role": role, "n": n, "s_mean": w["s_a"], "r_mean": w["s_b"], "s_sd": w["sd_a"],
                         "r_sd": w["sd_b"], "s_min": w["min_a"], "s_max": w["max_a"], "r_min": w["min_b"],
                         "r_max": w["max_b"], "caveat": S.caveat.get(cid, ""), "delta": w["delta"], "rel": w["rel"], "lo95": w["lo"], "hi95": w["hi"],
                         "band_lo": bands["lo"][k], "band_hi": bands["hi"][k], "band_crit": bands["crit"][k],
                         "K": len(ws), "p": w["p"], "R_s": w["R_a"], "R_r": w["R_b"]})
    df = pd.DataFrame(recs)
    out["q2_curves"] = df.to_dict("records")
    return df


def n_req_table(S: Scorer, out: dict) -> pd.DataFrame:
    grid = S.P["random_grid"]
    Ar, er, _ = S.M("R0", "pstar", grid)
    recs = []
    for cid in S.fins + [c for c, d in S.P["plan"]["arms"].items() if d["role"] in ("R1", "static", "det")]:
        if not S.has(cid):
            continue
        ns = sorted({int(x["n"]) for x in S.by[(cid, "pstar")]})
        As, es, _ = S.M(cid, "pstar", ns)
        for j, n in enumerate(ns):
            if not np.isfinite(es[:, j]).any():
                continue
            r = q2dev.n_req_joint(As, es, Ar, er, grid, j, np.random.default_rng(F.boot_seed(S.P, f"nreq_{cid}_{n}")))
            recs.append({"arm": cid, "n": n, "n_req": r["n_req"], "censored": r["censored"], "lo": r["ci"][0],
                         "hi": r["ci"][1], "frac_left_censored": r["frac_left"], "frac_right_censored": r["frac_right"],
                         "ratio": r["n_req"] / n if np.isfinite(r["n_req"]) else float("nan")})
    df = pd.DataFrame(recs)
    out["n_req"] = df.to_dict("records")
    return df


def q1(S: Scorer, out: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Random curves per predictor: mean/SD/CI, lift, precision at 10% recall, SLs detected; n_80/n_90 (PAVA,
    log-n interpolation, censored, joint bootstrap); saturating fit a − b·n^−c per permutation (P*)."""
    prev = float(S.L.y.mean())
    curves, nx = [], []
    preds = list(S.cfg["random"]["predictors"]) + ["m1_hgb"]
    pool_n = S.P["pool_n"]
    for pk in preds:
        rows = S.by.get(("R0", pk), [])
        if not rows:
            continue
        grid = sorted({int(x["n"]) for x in rows})
        A, est, tr = S.M("R0", pk, grid)
        Mb = E.joint_resample(A, np.random.default_rng(F.boot_seed(S.P, f"q1_{pk}")))
        for j, n in enumerate(grid):
            rr = [x for x in rows if int(x["n"]) == n]
            ok = np.isfinite(est[:, j])
            pr10 = [p1.precision_at_recall(S.L.y, S.eng.scores(x["pred_file"]), 0.1) for x in rr]
            st = {f"auprc_{k}": float(np.mean([S.point(x["pred_file"], k) for x in rr]))
                  for k in ("AO", "TT", "onrows", "offrows", "tert1", "tert2", "tert3")} if pk == "pstar" else {}
            sec = pd.DataFrame([S.secondary(x) for x in rr]).mean().to_dict()
            lo, hi = q2dev.ci(Mb[:, j])
            curves.append({"predictor": pk, "n": n, "R": int(ok.sum()), "auprc_mean": float(np.nanmean(est[:, j])),
                           "auprc_sd": float(np.nanstd(est[:, j], ddof=1)) if ok.sum() > 1 else float("nan"),
                           "ci_lo": lo, "ci_hi": hi, "prevalence": prev,
                           "lift": float(np.nanmean(est[:, j])) / prev, "prec_at_recall10": float(np.mean(pr10)),
                           "sl_detected_frac": float(np.mean([x["sl_detected_frac"] for x in rr])),
                           "fit_predict_s": float(np.median([x["fit_predict_s"] for x in rr])),
                           "boot_reps_dropped": int(np.isnan(Mb[:, j]).sum()), **sec, **st})
        curve = np.nanmean(est, axis=0)
        anchor_is_pool = grid[-1] == pool_n
        base = {"predictor": pk, "R": len(tr), "grid": "|".join(map(str, grid)), "anchor_n": grid[-1],
                "anchor_is_pool_max": anchor_is_pool, "auprc_anchor_pava": float(E.pava(curve)[-1]),
                "R_at_anchor": int(np.isfinite(est[:, -1]).sum()),
                # replicates whose resampled trajectories include none observed at some budget (R < R_r there):
                # their curve has a NaN point and n_X is NaN, dropped from the CI (review L; counted here)
                "boot_reps_dropped": int(np.isnan(Mb).any(axis=1).sum())}
        for fr in (float(STAT["n80_frac"]), float(STAT["n90_frac"])):
            key = f"n{int(round(fr * 100))}"
            pt = E.n_req(fr * E.pava(curve)[-1], grid, curve)
            vals = np.array([E.n_req(fr * E.pava(Mb[b])[-1], grid, Mb[b]).n for b in range(len(Mb))])
            lo, hi = E.percentile_ci(vals)
            base.update({key: pt.n, f"{key}_censored": pt.censored, f"{key}_lo": lo, f"{key}_hi": hi,
                         f"{key}_frac_censored": float(np.mean((vals == 0) | ~np.isfinite(vals))),
                         f"{key}_frac_pool": pt.n / pool_n})
        nx.append(base)
    sat = []
    rows = S.by.get(("R0", "pstar"), [])
    if rows:
        from scipy.optimize import curve_fit
        fn = lambda n, a, b, c: a - b * n ** (-c)  # noqa: E731
        df = pd.DataFrame([{"traj": x["traj"], "n": x["n"], "auprc": S.point(x["pred_file"], "all")} for x in rows])
        nt = df.traj.nunique()
        common = sorted(n for n, k in df.groupby("n").traj.nunique().items() if k == nt)   # budgets every traj has
        df = df[df.n.isin(common)]
        pars = []
        for _, d in df.groupby("traj"):
            d = d.sort_values("n")
            if len(d) < 4:
                continue
            try:
                p, _ = curve_fit(fn, d.n.to_numpy(float), d.auprc.to_numpy(), p0=(0.08, 0.2, 0.5),
                                 bounds=([0, 0, 0.01], [1, 10, 3]), maxfev=20000)
                pars.append(p)
            except RuntimeError:
                continue
        if pars:
            Pm = np.array(pars)
            rng = np.random.default_rng(F.boot_seed(S.P, "q1_satfit"))
            bs = np.array([Pm[rng.integers(0, len(Pm), len(Pm))].mean(0) for _ in range(2000)])
            sat.append({"predictor": "pstar", "n_fits": len(Pm), "grid": "|".join(map(str, common)),
                        **{k: float(Pm[:, i].mean()) for i, k in enumerate("abc")},
                        **{f"{k}_lo": float(np.quantile(bs[:, i], .025)) for i, k in enumerate("abc")},
                        **{f"{k}_hi": float(np.quantile(bs[:, i], .975)) for i, k in enumerate("abc")}})
    out["q1_curves"], out["q1_nx"], out["q1_satfit"] = curves, nx, sat
    return pd.DataFrame(curves), pd.DataFrame(nx), pd.DataFrame(sat)


def strata_sens(S: Scorer, contrasts: list, out: dict) -> dict:
    """Strata (E_AO, E_TT, on-array rows), sensitivities (lockbox-only T, no same-chromosome cells, two-way
    bootstrap), m1_hgb and the reference predictor; the DESIGN §3 headline rule."""
    sens = S.cfg["sensitivities"]
    sets = list(sens["strata"]) + (["lockbox"] if sens.get("lockbox_only") else []) + \
        (["nosamechr"] if sens.get("no_same_chromosome") else [])
    c1 = S.P["c1"]
    ns = [int(n) for n in sens["contrast_budgets"]]
    res = {"headline": {}}
    for fid in S.fins:
        for st in ["all"] + sets:
            for n in ns:
                w = S.welch(fid, "pstar", n, "R0", "pstar", n, stratum=st)
                if w:
                    contrasts.append({**w, "family": "strata" if st in ("AO", "TT", "onrows")
                                      else ("full_set" if st == "all" else "sensitivity")})
            w = S.welch(fid, "pstar", int(c1["strategy_n"]), "R0", "pstar", int(c1["random_n"]), stratum=st)
            if w and st != "all":
                contrasts.append({**w, "family": "sensitivity_c1"})
        if sens.get("two_way_bootstrap"):
            for n_a, n_b in ((int(c1["strategy_n"]), int(c1["random_n"])), (160, 160)):
                w = S.welch(fid, "pstar", n_a, "R0", "pstar", n_b, two_way=True)
                if w:
                    contrasts.append({**w, "family": "sensitivity_two_way"})
        # m1_hgb (secondary predictor; R = 5 per arm)
        m1n = sorted({int(x["n"]) for x in S.by.get((fid, "m1_hgb"), [])})
        for n in m1n:
            for st in ("all", "AO", "TT"):
                w = S.welch(fid, "m1_hgb", n, "R0", "m1_hgb", n, stratum=st)
                if w:
                    contrasts.append({**w, "family": "m1_hgb"})
        if int(c1["strategy_n"]) in m1n:
            w = S.welch(fid, "m1_hgb", int(c1["strategy_n"]), "R0", "m1_hgb", int(c1["random_n"]))
            if w:
                contrasts.append({**w, "family": "m1_hgb_c1"})
        # headline rule: Δ(160) same sign on E_fixed, E_AO, E_TT with P* and with m1_hgb
        n = int(config.budgets()["contrast_a"]["n"])
        signs = {}
        for pk in ("pstar", "m1_hgb"):
            for st in ("all", "AO", "TT"):
                w = S.welch(fid, pk, n, "R0", pk, n, stratum=st)
                signs[f"{pk}:{st}"] = None if w is None else float(np.sign(w["delta"]))
        vals = [v for v in signs.values() if v is not None]
        res["headline"][fid] = {"signs": signs, "complete": len(vals) == 6,
                                "headline": len(vals) == 6 and len(set(vals)) == 1 and vals[0] > 0}
    # reference predictor (stack_lr_wr2): strategy vs random and paired P* − ref on the same trajectories
    for fid in S.fins:
        for n in S.P["ref"]["cps"]:
            w = S.welch(fid, "ref", n, "R0", "ref", n)
            if w:
                contrasts.append({**w, "family": "reference_predictor"})
        w = S.welch(fid, "ref", int(c1["strategy_n"]), "R0", "ref", int(c1["random_n"]))
        if w:
            contrasts.append({**w, "family": "reference_predictor_c1"})
    paired = []
    for cid in S.fins + ["R0"]:
        for n in S.P["ref"]["cps"]:
            p_ = S.paired(cid, "pstar", "ref", n)
            if p_:
                paired.append(p_)
    res["paired_pstar_vs_ref"] = paired
    rng_rows = []
    allsets = ["all"] + list(sens["strata"])              # per-stratum range (sensitivities listed separately)
    for fid in S.fins:
        for n in ns:
            ws = {st: S.welch(fid, "pstar", n, "R0", "pstar", n, stratum=st) for st in allsets}
            ws = {k: v for k, v in ws.items() if v}
            if "all" not in ws:
                continue
            w0 = ws["all"]
            rels = [v["rel"] for v in ws.values()]
            rng_rows.append({"finalist": fid, "n": n, "s": w0["s_a"], "s_min": w0["min_a"], "s_max": w0["max_a"],
                             "r": w0["s_b"], "r_min": w0["min_b"], "r_max": w0["max_b"], "rel": w0["rel"],
                             "rel_min": float(min(rels)), "rel_max": float(max(rels)), "sets": list(ws)})
    out["ranges"] = rng_rows
    out["strata_sensitivities"] = res
    return res


def descriptive(S: Scorer, contrasts: list, out: dict) -> None:
    """Static vs adaptive, deterministic vs randomised, warm start, batch size (all outside the family)."""
    arms_ = S.P["plan"]["arms"]
    for cid, d in arms_.items():
        role, base = d["role"], d["base"]
        if not S.has(cid):
            continue
        ns = sorted({int(x["n"]) for x in S.by[(cid, "pstar")]})
        refs = []
        if role in ("static", "det"):
            refs = [(base, "same_arm_adaptive_randomised")]
        elif role == "ws":
            refs = [("R1", "R1_continuation"), (f"R0__ws{S.P['ws_n0']}", "random_continuation")]
        elif role == "batch":
            refs = [(base, "same_arm_batch10"), ("R0", "random")]
        for ref, lab in refs:
            for n in ns:
                w = S.welch(cid, "pstar", n, ref, "pstar", n)
                if w:
                    contrasts.append({**w, "family": f"descriptive:{lab}", "caveat": S.caveat.get(cid, "")})
    out["descriptive_done"] = True


def per_arm_tables(S: Scorer) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    eu, disc, rt = [], [], []
    for (cid, pk), rows in S.by.items():
        df = pd.DataFrame(rows)
        for n, d in df.groupby("n"):
            if pk == "pstar" and "e_unobs" in d and d.e_unobs.notna().any():
                m = pd.DataFrame([x for x in d.e_unobs if isinstance(x, dict)])
                eu.append({"arm": cid, "n": int(n), "R": len(m), "auprc_mean": m.auprc.mean(), "auprc_sd": m.auprc.std(),
                           "prevalence_mean": m.prevalence.mean(), "lift_mean": m.lift.mean(),
                           "n_cells_mean": m.n_cells.mean()})
            if pk == "pstar":
                disc.append({"arm": cid, "n": int(n), "R": len(d), "sl_detected_mean": d.sl_detected.mean(),
                             "sl_detected_frac_mean": d.sl_detected_frac.mean(),
                             "offarray_share_mean": d.offarray_share.mean(),
                             "auprc_mean": float(np.mean([S.point(f, "all") for f in d.pred_file])),
                             **pd.DataFrame([S.secondary(x) for x in S.by[(cid, pk)] if int(x["n"]) == int(n)])
                             .mean().add_suffix("_mean").to_dict()})
            rt.append({"arm": cid, "pred": pk, "n": int(n), "R": len(d), "fit_predict_s_median": d.fit_predict_s.median(),
                       "eval_total_s_median": d.eval_total_s.median(), "unit_wall_s_median": d.unit_wall_s.median(),
                       "peak_rss_fit_gb_max": d.peak_rss_fit_gb.max(), "select_s_cum_mean": d.select_s_cum.mean()})
    return pd.DataFrame(eu), pd.DataFrame(disc), pd.DataFrame(rt)


# ============================================================================ figures
def figures(S: Scorer, OUT: Path, q1c: pd.DataFrame, q2: pd.DataFrame, disc: pd.DataFrame, tag: str) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    def _ticks(ax, ns):
        ns = sorted({int(n) for n in ns})
        ax.set_xticks(ns)
        ax.set_xticklabels([str(n) for n in ns], fontsize=7, rotation=45)
        ax.xaxis.set_minor_formatter(NullFormatter())
    pal = ["#1f4e9c", "#d1495b", "#1baf7a", "#eda100", "#6b7280", "#8e44ad", "#2a78d6", "#b07c00"]
    files = []
    if len(q1c):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for i, (pk, d) in enumerate(q1c.groupby("predictor")):
            d = d.sort_values("n")
            ax.plot(d.n, d.auprc_mean, marker="o", ms=3, color=pal[i % len(pal)], label=f"{pk} (R≤{d.R.max()})")
            ax.fill_between(d.n, d.ci_lo, d.ci_hi, color=pal[i % len(pal)], alpha=0.15, lw=0)
        ax.axhline(float(S.L.y.mean()), color="k", lw=0.8, ls=":", label="prevalence")
        ax.set_xscale("log")
        _ticks(ax, q1c.n)
        ax.set_xlabel("queries screened n (random)")
        ax.set_ylabel(f"AUPRC on {S.L.man['eset']}")
        ax.set_title(f"Q1 random curves{tag}")
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "fig_final_q1.png", dpi=150)
        plt.close(fig)
        files.append("fig_final_q1.png")
    if len(q2):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for i, (cid, d) in enumerate(q2.groupby("arm")):
            d = d.sort_values("n")
            ax.plot(d.n, d.delta, marker="o", ms=3, color=pal[i % len(pal)], label=cid)
            ax.fill_between(d.n, d.band_lo, d.band_hi, color=pal[i % len(pal)], alpha=0.12, lw=0)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xscale("log")
        _ticks(ax, q2.n)
        ax.set_xlabel("n")
        ax.set_ylabel("Δ AUPRC vs random at matched n (Bonferroni bands)")
        ax.set_title(f"Q2 strategies vs random{tag}")
        ax.legend(fontsize=7, frameon=False, ncol=2)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "fig_final_q2_delta.png", dpi=150)
        plt.close(fig)
        files.append("fig_final_q2_delta.png")
    if len(disc):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for i, (cid, d) in enumerate(disc.groupby("arm")):
            if "__" in cid:
                continue
            d = d.sort_values("n")
            ax.plot(d.sl_detected_mean, d.auprc_mean, marker="o", ms=3, color=pal[i % len(pal)], label=cid,
                    lw=2 if cid == "R0" else 1.2)
        ax.set_xlabel("strong-SL cells detected in screened rows")
        ax.set_ylabel("AUPRC (P*)")
        ax.set_title(f"Discovery frontier{tag}")
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT / "fig_final_frontier.png", dpi=150)
        plt.close(fig)
        files.append("fig_final_frontier.png")
    return files


# ============================================================================ summary
def write_summary(S, OUT, out, fam, bg, q2, nreq, q1c, nx, sat, cdf, eu, disc, rt, figs, partial, prog):
    dry = S.P["mode"] == "dryrun"
    L = [f"# Final evaluation{' — DRY RUN (pseudo-final on ' + S.P['split'] + '; numbers are meaningless)' if dry else ''}",
         ""]
    if dry:
        L += [f"> DRY RUN: {S.P['split']} stands in for T, R_s = {S.P['R_s']}, R_r = {S.P['R_r']}, reduced grids, "
              f"B = {S.P['B']}. This file proves the pipeline runs end to end; no number here is a result.", ""]
    if partial:
        L += ["> PARTIAL: scored before every unit finished (--partial).", ""]
    lk = out["lock"]
    L += [f"Lock: git `{lk['git_head'][:10]}`, created {lk['created']}, plan sha `{lk['plan_sha256'][:12]}`; units "
          f"done {prog[0]}/{prog[1]}; evaluation set `{S.L.man['eset']}` ({S.L.man['counts']['n_cells']} cells, "
          f"prevalence {S.L.man['counts']['prevalence']:.5f}). Predictor P* `{S.P['preds']['pstar'][0]}`. Inference: "
          f"B = {S.P['B']} shared stratified query-bootstrap replicates; Welch SE (V_Q + σ̂²_s/R_s + σ̂²_r/R_r).", ""]
    sw = out["scored_with"]
    L += [f"Scored with git `{sw['git']['head'][:10]}` (tree dirty: {sw['git']['dirty']}); frozen check: "
          f"{'ok' if sw['frozen_check']['ok'] else 'NOT OK — ' + sw['frozen_check']['detail']}; code SHA-256 in "
          "summary.json `scored_with`. Headline rule: the shared sign must also be positive.", ""]
    L += ["## Confirmatory family (one-sided, Holm, familywise α = " + _f(S.alpha, 3) + ")", "",
          "| finalist | H | s | r | Δ | 95% CI | p | p (Holm) | reject |", "|---|---|---|---|---|---|---|---|---|"]
    for _, r in fam.iterrows():
        L.append(f"| {r.finalist} | {r.hypothesis}: {r.H1} | {_f(r.get('s_a'))} | {_f(r.get('s_b'))} | {_f(r.get('delta'))} | "
                 f"[{_f(r.get('lo'))}, {_f(r.get('hi'))}] | {_f(r.p, 4)} | {_f(r.p_holm, 4)} | {'**yes**' if r.reject else 'no'} |")
    L += ["", "## Brief example (s160 vs r660), gap closure G, G*, 'matches 20% random'", "",
          "| finalist | s160 − r660 [CI] | G [95% CI] (defined; non-finite reps) | G* | matches: s160 − 0.95·r660 [CI], p | claim allowed |",
          "|---|---|---|---|---|---|"]
    for r in bg:
        L.append(f"| {r['finalist']} | {_f(r['brief_delta'])} [{_f(r['brief_lo'])}, {_f(r['brief_hi'])}] | "
                 f"{_f(r['G'], 2)} [{_f(r['G_lo'], 2)}, {_f(r['G_hi'], 2)}] ({'yes' if r['G_defined'] else 'no'}; "
                 f"{100 * r['G_frac_nonfinite']:.1f}%) | "
                 f"{_f(r['G_star'], 2)} | {_f(r['matches_delta'])} [{_f(r['matches_lo'])}, {_f(r['matches_hi'])}], "
                 f"p = {_f(r['matches_p'], 3)} | {'yes' if r['matches_claim_allowed'] else 'no'} |")
    L += ["", "'Matches' may be stated only if the matches contrast is rejected, c1 is rejected (Holm) and G* ≥ "
          f"{STAT['g_star_min']} (DESIGN §6.4).", ""]
    hl = out["strata_sensitivities"]["headline"]
    L += ["## Headline rule (Δ(160) same sign on E_fixed, E_AO, E_TT with P* and m1_hgb)", "",
          "| finalist | signs | complete | headline |", "|---|---|---|---|"]
    for fid, h in hl.items():
        L.append(f"| {fid} | " + ", ".join(f"{k} {'+' if v and v > 0 else ('−' if v is not None else '?')}"
                                          for k, v in h["signs"].items())
                 + f" | {h['complete']} | {'**headline**' if h['headline'] else 'qualified'} |")
    if len(q2):
        L += ["", "## Δ(n) vs random at matched n (P*; Bonferroni bands over each arm's checkpoints)", "",
              "| arm | " + " | ".join(str(n) for n in sorted(q2.n.unique())) + " |",
              "|---|" + "---|" * q2.n.nunique()]
        for cid, d in q2.groupby("arm"):
            v = {int(r.n): r for _, r in d.iterrows()}
            L.append(f"| {cid}{' †' if S.caveat.get(cid) else ''} | " + " | ".join(
                f"{_f(v[n].delta)} [{_f(v[n].band_lo)}, {_f(v[n].band_hi)}]" if n in v else "—"
                for n in sorted(q2.n.unique())) + " |")
        if S.caveat:
            L += ["", "† " + "; ".join(f"{c}: {t}" for c, t in S.caveat.items()) + "."]
    rg = out.get("ranges", [])
    if rg:
        L += ["", "## Means with ranges (finalists; trajectory min–max; relative Δ range across the strata)", "", "Per-split ranges do not exist on T (one held-out split); the dev per-split ranges are in "
              "results/p4/summary.md. Sets: " + ", ".join(rg[0]["sets"]) + ".", "",
              "| finalist | n | s mean [traj min, max] | r mean [traj min, max] | rel Δ (E_fixed) | rel Δ range over strata |",
              "|---|---|---|---|---|---|"]
        for r in rg:
            L.append(f"| {r['finalist']} | {r['n']} | {_f(r['s'])} [{_f(r['s_min'])}, {_f(r['s_max'])}] | {_f(r['r'])} "
                     f"[{_f(r['r_min'])}, {_f(r['r_max'])}] | {100 * r['rel']:+.1f}% | [{100 * r['rel_min']:+.1f}%, "
                     f"{100 * r['rel_max']:+.1f}%] |")
    if len(nreq):
        L += ["", "## Random-equivalent budget n_req(n) (inverts the P* random curve; joint bootstrap; censored "
              "replicates coded 0/∞)", "", "| arm | n | n_req [95% CI] | n_req / n | censored (left/right) |",
              "|---|---|---|---|---|"]
        for _, r in nreq.iterrows():
            L.append(f"| {r.arm} | {r.n} | {_f(r.n_req, 0)}{' (' + r.censored + ')' if r.censored else ''} "
                     f"[{_f(r.lo, 0)}, {_f(r.hi, 0)}] | {_f(r.ratio, 2)} | {_f(r.frac_left_censored, 2)}/"
                     f"{_f(r.frac_right_censored, 2)} |")
    if len(nx):
        L += ["", "## Q1: n_80 / n_90 (PAVA + log-n interpolation of the random mean curve; anchor = AUPRC at the "
              "pool maximum; joint bootstrap)", "",
              "| predictor | R (at anchor) | anchor n (pool max?) | AUPRC anchor | n_80 [CI] | n_90 [CI] | censored frac 80/90 | boot reps dropped |",
              "|---|---|---|---|---|---|---|---|"]
        for _, r in nx.iterrows():
            L.append(f"| {r.predictor} | {r.R} ({r.R_at_anchor}) | {r.anchor_n} ({'yes' if r.anchor_is_pool_max else 'NO'}) | "
                     f"{_f(r.auprc_anchor_pava)} | {_f(r.n80, 0)} [{_f(r.n80_lo, 0)}, {_f(r.n80_hi, 0)}] | "
                     f"{_f(r.n90, 0)} [{_f(r.n90_lo, 0)}, {_f(r.n90_hi, 0)}] | {_f(r.n80_frac_censored, 2)}/"
                     f"{_f(r.n90_frac_censored, 2)} | {r.boot_reps_dropped} |")
    if len(q1c):
        L += ["", "## Q1 random curves (mean, SD, joint-bootstrap CI of the trajectory mean)", "",
              "| predictor | n | R | AUPRC | SD | CI | lift | prec@10% recall | SLs detected (frac) |",
              "|---|---|---|---|---|---|---|---|---|"]
        for _, r in q1c.sort_values(["predictor", "n"]).iterrows():
            L.append(f"| {r.predictor} | {r.n} | {r.R} | {_f(r.auprc_mean)} | {_f(r.auprc_sd)} | [{_f(r.ci_lo)}, "
                     f"{_f(r.ci_hi)}] | {_f(r.lift, 2)} | {_f(r.prec_at_recall10, 3)} | {_f(r.sl_detected_frac, 3)} |")
    if len(sat):
        r = sat.iloc[0]
        L += ["", f"Saturating fit a − b·n^−c (P*, per permutation on the budgets every permutation has: {r.grid}; "
              f"seed-cluster bootstrap): a = {_f(r.a)} "
              f"[{_f(r.a_lo)}, {_f(r.a_hi)}], b = {_f(r.b)}, c = {_f(r.c, 3)} [{_f(r.c_lo, 3)}, {_f(r.c_hi, 3)}] "
              f"({int(r.n_fits)} fits; descriptive)."]
    if len(cdf):
        L += ["", "## Strata, sensitivities, secondary predictors, descriptive arms (all descriptive)", "",
              "| family | a (pred, n) | b (pred, n) | set | boot | Δ | rel | 95% CI | p | R_a/R_b |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        keep = cdf[cdf.family != "q2_curve"]
        for _, r in keep.iterrows():
            cv = " †" if isinstance(r.get("caveat"), str) and r.get("caveat") else ""
            L.append(f"| {r.family} | {r.a}{cv} ({r.pred_a}, {r.n_a}) | {r.b} ({r.pred_b}, {r.n_b}){' ×' + _f(r.scale_b, 2) if r.scale_b != 1 else ''} "
                     f"| {r.set} | {r.bootstrap} | {_f(r.delta)} | {100 * r.rel:+.1f}% | [{_f(r.lo)}, {_f(r.hi)}] | "
                     f"{_f(r.p, 3)} | {r.R_a}/{r.R_b} |")
    pr = out["strata_sensitivities"].get("paired_pstar_vs_ref", [])
    if pr:
        L += ["", "Paired P* − reference (stack_lr_wr2) on the same trajectories:", "",
              "| arm | n | Δ | rel | 95% CI | R |", "|---|---|---|---|---|---|"]
        for r in pr:
            L.append(f"| {r['arm']} | {r['n']} | {_f(r['delta'])} | {100 * r['rel']:+.1f}% | [{_f(r['lo'])}, {_f(r['hi'])}] | {r['R']} |")
    if len(eu):
        L += ["", "## E_unobs per arm (descriptive, P*; not paired)", "",
              "| arm | n | R | AUPRC | prevalence | lift |", "|---|---|---|---|---|---|"]
        for _, r in eu.sort_values(["arm", "n"]).iterrows():
            L.append(f"| {r.arm} | {r.n} | {r.R} | {_f(r.auprc_mean)} | {_f(r.prevalence_mean, 5)} | {_f(r.lift_mean, 2)} |")
    if len(disc):
        L += ["", "## Discovery, off-array share and DESIGN §3 secondary metrics (P* rows; means over trajectories)", "",
              "| arm | n | SLs detected | frac of pool SLs | off-array share | AUPRC | macro AUPRC | AUROC | p@10 | p@50 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for _, r in disc.sort_values(["arm", "n"]).iterrows():
            L.append(f"| {r.arm} | {r.n} | {_f(r.sl_detected_mean, 0)} | {_f(r.sl_detected_frac_mean, 3)} | "
                     f"{_f(r.offarray_share_mean, 3)} | {_f(r.auprc_mean)} | {_f(r.macro_auprc_mean)} | "
                     f"{_f(r.auroc_mean, 3)} | {_f(r.p_at_10_mean, 3)} | {_f(r.p_at_50_mean, 3)} |")
    L += ["", "## Runtime", "", f"See `runtime.csv` (per arm × predictor × n: fit+predict, unit wall, peak RSS). "
          f"Bootstrap vectors computed this run: {S.eng.n_vectors}.", "",
          "Files: " + ", ".join(f"`{x}`" for x in ["summary.json", "confirmatory.csv", "contrasts.csv", "q2_curves.csv",
                                                   "n_req.csv", "q1_curves.csv", "q1_nx.csv", "q1_satfit.csv",
                                                   "e_unobs.csv", "discovery.csv", "runtime.csv", "runs.csv",
                                                   "selections.json.gz", "lock_record.json", "run_ledger.jsonl"] + figs), ""]
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


# ============================================================================ main score
def score(mode: str, partial: bool = False) -> dict:
    t0 = time.time()
    cfg, P = F.mode_params(mode)
    lk = F.read_lock(P["lock"])
    if lk is None:
        raise SystemExit(f"no lock at {P['lock']}: nothing was run")
    plan = F.load_plan(P)
    if plan["sha256"] != lk["plan_sha256"]:
        raise SystemExit("plan differs from the locked plan")
    P["plan"] = plan
    if not (P["run"] / "labels.npz").exists():
        raise SystemExit("labels.npz missing: the `_labels` job has not run")
    prog = F.progress(P, plan)
    done, tot = sum(d for d, _ in prog.values()), sum(t for _, t in prog.values())
    if done < tot and not partial:
        raise SystemExit(f"{tot - done} of {tot} units not done: finish (run_final.py resume) or pass --partial")
    if done < tot and mode == "final":
        raise SystemExit("score --partial is refused in final mode (no unplanned look at T; review F3)")
    # F4: the scoring code and every frozen input must be the frozen versions; record what scored
    fz = F.check_frozen(ROOT)
    if mode == "final" and not fz.ok:
        raise SystemExit(f"frozen files changed since the freeze: {fz.detail}")
    code = {f: F.sha256_file(ROOT / f) for f in ("scripts/score_final.py", "src/hsl/final.py", "src/hsl/q2dev.py",
                                                  "src/hsl/eval.py", "src/hsl/p1.py", "config/finalists.yaml",
                                                  "config/protocol.yaml")}
    rows = {}
    for x in F.read_parts(P):
        if x.get("pred_file") and (ROOT / x["pred_file"]).exists():
            rows[x["uid"]] = x                                    # last row per unit wins (re-runs are identical)
    rows = list(rows.values())
    S = Scorer(cfg, P, rows)
    bad = [x["uid"] for x in rows if x["eset_digest"] != S.L.man["eset_digest"]]
    if bad:
        raise AssertionError(f"eset digest mismatch: {bad[:3]}")
    n_chk = F.check_inline([x for x in rows], S.eng)
    OUT = P["out"]
    OUT.mkdir(parents=True, exist_ok=True)
    out = {"mode": mode, "lock": lk, "units_done": done, "units_total": tot, "partial": bool(done < tot),
           "inline_ap_checked": n_chk, "label_manifest": S.L.man, "R_s": P["R_s"], "R_r": P["R_r"], "B": P["B"],
           "scored_with": {"git": views.git_state(), "frozen_check": fz.__dict__, "code_sha256": code}}
    contrasts = []
    fam = family(S, out)
    bg = brief_and_g(S, fam, out)
    q2 = delta_curves(S, contrasts, out)
    nreq = n_req_table(S, out)
    q1c, nx, sat = q1(S, out)
    strata_sens(S, contrasts, out)
    descriptive(S, contrasts, out)
    eu, disc, rt = per_arm_tables(S)
    cdf = pd.DataFrame(contrasts)
    fam.to_csv(OUT / "confirmatory.csv", index=False)
    cdf.to_csv(OUT / "contrasts.csv", index=False)
    q2.to_csv(OUT / "q2_curves.csv", index=False)
    nreq.to_csv(OUT / "n_req.csv", index=False)
    q1c.to_csv(OUT / "q1_curves.csv", index=False)
    nx.to_csv(OUT / "q1_nx.csv", index=False)
    sat.to_csv(OUT / "q1_satfit.csv", index=False)
    eu.to_csv(OUT / "e_unobs.csv", index=False)
    disc.to_csv(OUT / "discovery.csv", index=False)
    rt.to_csv(OUT / "runtime.csv", index=False)
    pd.DataFrame([{k: v for k, v in x.items() if k not in ("e_unobs", "params")} for x in rows]).to_csv(
        OUT / "runs.csv", index=False)
    sels = {}
    for p in sorted((P["run"] / "selections").glob("*.json")):
        d = json.loads(p.read_text())
        sels[p.stem] = {"traj_id": d["traj_id"], "model_seed": d["model_seed"], "selected": d["selected"],
                        "selection_wall_s": d["selection_wall_s"], "s0_n": d["s0_n"], "batch_size": d["batch_size"]}
    with gzip.open(OUT / "selections.json.gz", "wt") as fh:
        json.dump(sels, fh)
    shutil.copyfile(P["lock"], OUT / "lock_record.json")
    if (P["run"] / "LEDGER_final.jsonl").exists():                # crash / resume history (§10), committed copy
        shutil.copyfile(P["run"] / "LEDGER_final.jsonl", OUT / "run_ledger.jsonl")
    tag = " (DRY RUN)" if mode == "dryrun" else ""
    figs = figures(S, OUT, q1c, q2[q2.role.isin(["finalist", "R1", "static"])] if len(q2) else q2, disc, tag)
    out["score_wall_s"] = time.time() - t0
    out["n_boot_vectors_computed"] = S.eng.n_vectors
    (OUT / "summary.json").write_text(json.dumps(out, indent=1, default=_json) + "\n")
    write_summary(S, OUT, out, fam, bg, q2, nreq, q1c, nx, sat, cdf, eu, disc, rt, figs, done < tot, (done, tot))
    F.ledger_event(P, "score", units_done=done, units_total=tot, wall_s=out["score_wall_s"])
    if mode == "final":
        from hsl import harness as H
        H.append_ledger({"id": "final_evaluation", "parent": "freeze", "config_hash": lk["plan_sha256"][:16],
                         "split": "final", "kind": "final_score", "outcome": "partial" if done < tot else "scored",
                         "metrics": {"family": out["family"], "brief_G": out["brief_G"]},
                         "manifest": {"lock": lk, "run_ledger": str((P["run"] / "LEDGER_final.jsonl").relative_to(ROOT))}})
    print(f"scored {len(rows)} units ({n_chk} inline APs re-checked) in {out['score_wall_s']:.0f} s -> {OUT}/summary.md")
    return out


# ============================================================================ cost projection
def cost_projection() -> int:
    """Projected cost of the REAL final from measured runtimes (dry-run rows on dev1; P1 random-curve rows incl.
    the n >= 1970 batch; P4 selection timing), scaled to |T| / |dev1 targets| for the prediction grid."""
    cfg, P = F.mode_params("final")
    plan = F.build_plan(cfg, P)
    _, Pd = F.mode_params("dryrun")
    meas = []
    if (Pd["run"] / "parts").exists():
        for x in F.read_parts(Pd):
            meas.append({"pred": x["pred"], "n": x["n"], "t": x["unit_wall_s"], "rss": x["peak_rss_fit_gb"],
                         "src": "dryrun"})
    for x in p1.read_rows([ROOT / "results" / "final_dryrun" / "rss_calibration.jsonl"]):
        if x["split"] != "resplit":                            # same-day dev1 calibration runs (526 targets)
            meas.append({"pred": x["pred"], "n": x["n"], "t": x["total_s"], "rss": x["peak_rss_gb"], "src": "dryrun"})
    p1rows = p1.read_rows(sorted((p1.P1_DIR / "parts").glob("*.jsonl")))
    pm = {v[0]: k for k, v in P["preds"].items()}
    for x in p1rows:
        if x.get("kind", "curve") != "curve" or x.get("source") == "p2_reuse" or x.get("wall_s") is None:
            continue
        k = pm.get(x.get("predictor"))
        if k is None:
            continue
        meas.append({"pred": k, "n": int(x["n"]), "t": float(x["wall_s"]), "rss": x.get("peak_rss_fit_gb"), "src": "p1"})
    md = pd.DataFrame(meas)
    # time scaling to |T| targets: from the calibration pairs (same n, re-split ~1,050 vs dev1 526 targets, same day):
    # r = fit+predict(re-split) / fit+predict(dev1); time(T) = time(dev) * (1 + (r - 1) * (|T| - |dev|) / (|rs| - |dev|))
    cal = p1.read_rows([ROOT / "results" / "final_dryrun" / "rss_calibration.jsonl"])
    ratios, n_rs = {}, None
    for x in cal:
        if x["split"] != "resplit":
            continue
        n_rs = x["n_targets"]
        d = [y for y in cal if y["split"] != "resplit" and y["pred"] == x["pred"] and y["n"] == x["n"]]
        if d:
            ratios.setdefault(x["pred"], []).append(x["fit_predict_s"] / d[-1]["fit_predict_s"])
    ratios = {k: float(np.mean(v)) for k, v in ratios.items()}
    r_def = ratios.get("pstar")

    def scale_for(pred):
        r = ratios.get(pred, r_def)
        if r is None or not n_rs:
            return P["n_targets"] / max(Pd["n_targets"], 1)          # proportional (conservative) fallback
        return max(1.0, 1 + (r - 1) * (P["n_targets"] - Pd["n_targets"]) / (n_rs - Pd["n_targets"]))
    scale_T = {k: scale_for(k) for k in P["preds"]}
    tab = md.groupby(["src", "pred", "n"]).agg(t=("t", "median"), rss=("rss", "max"), k=("t", "size")).reset_index()

    def _curve(d, n):
        x, y, r = d.n.to_numpy(float), d.t.to_numpy(float), d.rss.to_numpy(float)
        if n <= x[-1]:
            if n < x[0]:
                return float(y[0]), float(r[0]), "interp"
            return float(np.interp(np.log(n), np.log(x), y)), float(np.interp(np.log(n), np.log(x), r)), "interp"
        m = x >= x[-1] / 3                                   # power law (time) / linear (RSS) over the top of the range
        if m.sum() >= 2:
            bt = np.polyfit(np.log(x[m]), np.log(y[m]), 1)
            br = np.polyfit(x[m], r[m], 1)
            t, rr = float(np.exp(np.polyval(bt, np.log(n)))), float(np.polyval(br, n))
        else:
            t, rr = float(y[-1] * n / x[-1]), float(r[-1] * n / x[-1])
        return max(t, float(y[-1])), max(rr, float(r[-1])), "EXTRAPOLATED"

    factors = {}
    for pk in tab.pred.unique():                             # dry-run (1 job) / P1 (loaded, <= 4 jobs) speed ratio
        a_ = tab[(tab.src == "dryrun") & (tab.pred == pk)].set_index("n").t
        b_ = tab[(tab.src == "p1") & (tab.pred == pk)].set_index("n").t
        common = sorted(set(a_.index) & set(b_.index))
        factors[pk] = float(np.median([a_[n] / b_[n] for n in common])) if common else 1.0

    def est(pred, n, how="loaded"):
        """'loaded' = P1 rows (<= 4 concurrent jobs on a shared machine), dry run where P1 has none;
        'calibrated' = P1 rows x the dry-run/P1 speed ratio of this predictor (1 job, 2 threads)."""
        d1 = tab[(tab.src == "p1") & (tab.pred == pred)].sort_values("n")
        dd = tab[(tab.src == "dryrun") & (tab.pred == pred)].sort_values("n")
        if len(d1):
            t, r, h = _curve(d1, n)
            return (t * factors.get(pred, 1.0) if how == "calibrated" else t), r, h
        if len(dd):
            t, r, h = _curve(dd, n)
            return (t / factors.get(pred, 1.0) if how == "loaded" else t), r, h
        return float("nan"), float("nan"), "none"

    # selection time per trajectory by arm role (dry-run selection files; full trajectories to 660 as in the final)
    sel_t = {}
    for p in sorted((Pd["run"] / "selections").glob("*.json")) if (Pd["run"] / "selections").exists() else []:
        d = json.loads(p.read_text())
        sel_t.setdefault(d["role"], []).append(d["selection_wall_s"])
    sel_t = {k: float(np.median(v)) for k, v in sel_t.items()}
    roles = {cid: d["role"] for cid, d in plan["arms"].items()}
    over = 15.0                                              # per-job process start + bundle + E_fixed build (s)
    g = P["compute"]["gate"]
    alone_n = int(config.protocol()["compute"]["alone_min_n"])
    m1_q1_trim = {30, 1310}                                  # §12 trigger list (m1_hgb Q1-only budgets)

    def build(how, keep=lambda c, n, pk: True):
        jobs, unit_rows = [], []
        for j in plan["jobs"]:
            us = [u for u in j["units"] if keep(u[0], int(u[2]), u[3])]
            if not us:
                continue
            T = over + sel_t.get(roles[j["cid"]], 0.0)
            for c, r, n, pk, e in us:
                t, _, _ = est(pk, int(n), how)
                t *= scale_T[pk]
                T += t
                unit_rows.append({"role": roles[c], "pred": pk, "n": int(n), "t": t, "est": how})
            jobs.append({"name": j["name"], "role": roles[j["cid"]], "group": j["group"], "cls": j["cls"], "t": T,
                         "rss": F.job_rss(P, us)[0], "big": max(int(u[2]) for u in us) >= alone_n})
        return jobs, unit_rows

    def simulate(jobs, policy, budget=None):
        """Event simulation of launch.py in plan order. 'gate' = D10 memory gate (<= max_jobs, <= max_big jobs at
        n >= 1970, budget - sum(running hints) >= hint + headroom); 'old4' = 4 jobs, n >= 1970 alone, 20 GB."""
        running, t, queue = [], 0.0, list(jobs)
        while queue or running:
            while queue:
                j = queue[0]
                busy = sum(x["rss"] for x in running)
                if policy == "gate":
                    ok = (len(running) < int(g["max_jobs"])
                          and not (j["big"] and sum(x["big"] for x in running) >= int(g["max_big"]))
                          and float(budget or g["budget_gb"]) - busy >= j["rss"] + float(g["headroom_gb"]))
                else:
                    ok = (len(running) < 4 and busy + j["rss"] <= 20.0
                          and not (j["cls"] == "alone" or any(x["cls"] == "alone" for x in running)))
                if not running or ok:
                    running.append({**j, "end": t + j["t"]})
                    queue.pop(0)
                else:
                    break
            nxt = min(x["end"] for x in running)
            t = nxt
            running = [x for x in running if x["end"] > nxt]
        return t / 3600

    def summarise(label, jobs):
        return {"scenario": label, "jobs": len(jobs), "job_hours": sum(x["t"] for x in jobs) / 3600,
                "wall_h_gate": simulate(jobs, "gate"), "wall_h_gate_lo": simulate(jobs, "gate", 12.0),
                "wall_h_gate_hi": simulate(jobs, "gate", 22.0), "wall_h_old4": simulate(jobs, "old4"),
                "mean_conc_gate": sum(x["t"] for x in jobs) / 3600 / max(simulate(jobs, "gate"), 1e-9)}

    variants = [
        ("D10 plan as configured", lambda c, n, pk: True),
        ("without the reference predictor", lambda c, n, pk: pk != "ref"),
        ("§12 trims: m1_hgb Q1 {30, 1310} and batch 50 dropped",
         lambda c, n, pk: not ((pk == "m1_hgb" and c == "R0" and n in m1_q1_trim) or c.endswith("__b50"))),
        ("both of the above", lambda c, n, pk: pk != "ref" and not (
            (pk == "m1_hgb" and c == "R0" and n in m1_q1_trim) or c.endswith("__b50"))),
    ]
    scen, uds = [], []
    for how in ("calibrated", "loaded"):
        lab = {"calibrated": "1-job speed (lower bound)", "loaded": "P1 loaded speed (upper bound)"}[how]
        for vl, keep in variants:
            jobs, unit_rows = build(how, keep)
            if vl == variants[0][0]:
                uds.append(pd.DataFrame(unit_rows))
            scen.append(summarise(f"{lab}: {vl}", jobs))
    ud = pd.concat(uds)
    by_arm = ud.groupby(["est", "role", "pred"]).agg(units=("t", "size"),
                                                      job_h=("t", lambda v: v.sum() / 3600)).reset_index()
    sc = pd.DataFrame(scen)
    OUT = ROOT / "results" / "final_dryrun"
    OUT.mkdir(parents=True, exist_ok=True)
    tb = tab.copy()
    tb["t_final_scaled"] = tb.t * tb.pred.map(scale_T)
    lines = ["# Cost projection for the REAL final (plan of mode `final`)", "",
             f"Inputs: per-unit wall time (fit + predict + E_fixed scoring + save) = median of measured runs per "
             f"(predictor, n): dry-run rows on dev1 (this machine, 1 job, 2 threads, shared with a running batch) and "
             f"P1 random-curve rows (up to 4 concurrent jobs, loaded machine), interpolated in log n; above the largest "
             f"measured n: power-law fit (time) / linear fit (RSS) over n >= max/3 (flagged). Scaled to |T| = {P['n_targets']} targets per predictor "
             f"(" + ", ".join(f"{k} ×{v:.2f}" for k, v in scale_T.items()) + f"; from same-day calibration pairs, re-split "
             f"{n_rs} vs dev1 {Pd['n_targets']} targets, fit+predict ratios " + ", ".join(f"{k} {v:.2f}" for k, v in ratios.items())
             + f"; predictors without a pair use P*'s ratio) for the larger prediction grid. Selection per trajectory = dry-run "
             f"median by arm role (full trajectories to 660): " + ", ".join(f"{k} {v:.0f} s" for k, v in sel_t.items())
             + f". Per-job overhead {over:.0f} s. Job memory = the launch.py `rss_gb` hint (finalists.yaml compute: "
               "fitted on measured peak RSS + calibrated per-target term + 0.5 GB), the same number the gate uses.", "",
             f"Plan: {plan['n_units']} units in {len(plan['jobs'])} jobs (+ labels job).", "",
             "Dry-run / P1 speed ratios per predictor: " + ", ".join(f"{k} {v:.2f}" for k, v in factors.items()) + ".", "",
             "## Units and job-hours by arm role × predictor (lower / upper bound)", "",
             "| role | predictor | units | job-h (1-job speed) | job-h (loaded) |", "|---|---|---|---|---|"]
    bw = by_arm.pivot_table(index=["role", "pred"], columns="est", values=["units", "job_h"]).reset_index()
    for _, r in bw.iterrows():
        lines.append(f"| {r[('role', '')]} | {r[('pred', '')]} | {int(r[('units', 'loaded')])} | "
                     f"{r[('job_h', 'calibrated')]:.1f} | {r[('job_h', 'loaded')]:.1f} |")
    lines += ["", f"Totals in units: {uds[0].t.sum() / 3600:.1f} (1-job speed) to {uds[1].t.sum() / 3600:.1f} (loaded) "
              f"job-h, + {len(plan['jobs']) * over / 3600:.1f} h job overhead + selection "
              f"{sum(sel_t.get(roles[j['cid']], 0) for j in plan['jobs']) / 3600:.1f} h.", "",
              "## Wall-clock (event simulation of the launch.py queue in plan order)", "",
              f"`gate` = D10 memory gate: ≤ {g['max_jobs']} jobs, ≤ {g['max_big']} jobs at n ≥ {alone_n}, start only if "
              f"{g['budget_gb']:.0f} GB (memory available to us: PI check 2026-10-29, free + inactive ≈ 15 of 24 GB, "
              "other non-project use ~0.7 GB; 12 / 22 GB shown as sensitivities) minus the running jobs' "
              f"hints ≥ hint + {g['headroom_gb']:.0f} GB. `old` = previous policy for comparison (4 jobs, n ≥ "
              f"{alone_n} alone, 20 GB). CPU contention at 6 jobs is not modelled beyond the loaded-speed bound "
              "(P1 rows ran up to 4 jobs on a shared machine).", "",
              f"| scenario | jobs | job-h | wall h, gate {g['budget_gb']:.0f} GB (mean jobs running) | gate 12 GB | gate 22 GB | "
              f"old, 4 jobs | > 10 h at {g['budget_gb']:.0f} GB? |", "|---|---|---|---|---|---|---|---|"]
    for _, r in sc.iterrows():
        lines.append(f"| {r.scenario} | {r.jobs} | {r.job_hours:.1f} | {r.wall_h_gate:.1f} ({r.mean_conc_gate:.1f}) | "
                     f"{r.wall_h_gate_lo:.1f} | {r.wall_h_gate_hi:.1f} | {r.wall_h_old4:.1f} | "
                     f"{'YES' if r.wall_h_gate > 10 else 'no'} |")
    if cal:
        lines += ["", "## Calibration runs (run_final.py calibrate; DEV ONLY; one job at a time, 2 threads)", "",
                  "| predictor | split | targets | n | peak RSS GB | fit+predict s | hint for T at this n |",
                  "|---|---|---|---|---|---|---|"]
        for x in cal:
            lines.append(f"| {x['pred']} | {x['split']} | {x['n_targets']} | {x['n']} | {x['peak_rss_gb']:.2f} | "
                         f"{x['fit_predict_s']:.0f} | {F.rss_hint_gb(P, x['pred'], x['n']):.2f} |")
    lines += ["", "## Measured per-unit times used (s, before × T scaling) and peak RSS (GB)", "",
              "| source | predictor | n | median s | × T | max RSS | k |", "|---|---|---|---|---|---|---|"]
    for _, r in tb[tb.pred.isin(["pstar", "ref", "m1_hgb", "stack_lr"])].iterrows():
        lines.append(f"| {r.src} | {r.pred} | {r.n} | {r.t:.1f} | {r.t_final_scaled:.1f} | {_f(r.rss, 2)} | {r.k} |")
    ex = sorted({(u["pred"], u["n"]) for u in unit_rows if est(u["pred"], u["n"])[2] == "EXTRAPOLATED"})
    jobs0, _ = build("calibrated")
    rss_alone = max((x["rss"] for x in jobs0 if x["cls"] == "alone"), default=float("nan"))
    lines += ["", f"Largest job memory hint (pool max; RSS extrapolated beyond n = 2103): {rss_alone:.1f} GB."]
    # score step (single process): one [1 + B] vector per prediction file + strata/two-way extras; measured 0.86 s per
    # B = 2000 vector on dev1 (748k cells), scaled by cells (E_fixed 999k); two-way 0.019 s per replicate (dev1)
    n_vec = plan["n_units"] + 6 * 2 * (len(cfg["finalists"]) * P["R_s"] + P["R_r"])
    n_tw = 2 * len(cfg["finalists"]) * P["R_s"] + 2 * P["R_r"]
    t_score = (n_vec * 0.86 * 999102 / 748157 + n_tw * P["B_two_way"] * 0.019 * 999102 / 748157) / 3600
    lines += [f"Score step (separate, 1 process, no bundle): ≈ {n_vec} bootstrap vectors + {n_tw} two-way vectors ≈ "
              f"{t_score:.1f} h (measured dry-run engine speed scaled to |E_fixed|)."]
    lines += ["", "Extrapolated (no measurement at this n): " + (", ".join(f"{p}@{n}" for p, n in ex) or "none"), ""]
    (OUT / "cost_projection.md").write_text("\n".join(lines) + "\n")
    (OUT / "cost_projection.json").write_text(json.dumps({"scenarios": sc.to_dict("records"),
                                                          "by_arm": by_arm.to_dict("records"),
                                                          "selection_s": sel_t, "scale_T": scale_T,
                                                          "speed_factors": factors},
                                                         indent=1, default=_json) + "\n")
    print("\n".join(lines))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", required=True, choices=list(F.MODES))
    ap.add_argument("--partial", action="store_true")
    ap.add_argument("--cost", action="store_true")
    a = ap.parse_args(argv)
    if a.cost:
        return cost_projection()
    score(a.mode, a.partial)
    return 0


if __name__ == "__main__":
    sys.exit(main())
