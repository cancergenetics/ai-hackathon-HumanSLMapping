#!/usr/bin/env python
"""P4 strategy development on DEV splits (DESIGN v3.1 §6, §8, §9 P4, §12; DEV DATA ONLY; never the test split).

Arms (config ids, hsl.q2dev.parse): R1, H, H-static, F_r{4,8,16}, F-static_r8, C, M_lam0.3 (+ λ grid), on dev1-3,
batch 10, trajectories r = 0..R-1 (dev trajectory seed r: same model seed as random permutation r, CRN), Q2 checkpoints
{30, 60, 100, 160, 220, 330, 490, 660} (+ 440 for R1: R1-continuation of the warm start) with the frozen P*.
Random arm = P1's frozen-P* rows (never recomputed here). Scenarios for the best <= 2 (descriptive, dev):
  <X>__det        deterministic (temperature-0) trajectory, one per split
  <X>__ws220      warm start S_0 = first 220 genes of R1 trajectory r, then X to 660; R0__ws220 = random continuation;
                  R1-continuation = the R1 arm itself (static: its 221st..660th genes)
  <X>__b50/b110   batch sizes 50 / 110 (checkpoints not divisible by b use the first n − b⌊n/b⌋ genes of the next
                  batch in greedy order: harness rule)
  <X>__feas       HAP1-feasible pool (bottom array-SMF decile of on-array pool genes removed); R0__feas reference
  <X>__min        HAP1-minimal selector features (no AP-MS); predictor unchanged

    run      --config F_r8 --split dev1 --traj 0 [--checkpoints ...] [--max-n 660]
    smoke    every listed config to a small n on dev1 (pipeline proof; results/p4/smoke/, not analysed)
    jobs     --name NAME --configs H,F_r8 [--splits ...] [--trajs 0-2] [--checkpoints ...]
    degdiag  H degree-model diagnostic: Spearman(q̂, true pool degree) along H trajectories (evaluation only)
    analyze  -> results/p4/summary.md, contrasts.csv, ranking.csv, power.json, diagnostics, figures [--ledger]
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from hsl import config, harness as H, p1, q2dev, views

OUT = q2dev.RESULTS / "p4"
PARTS = OUT / "parts"
BOOT = OUT / "boot"
PY = ".venv/bin/python"
Q2 = [int(n) for n in config.budgets()["q2_checkpoints"]]
WS = [int(n) for n in config.budgets()["warm_start_checkpoints"]]
SCEN_CPS = [160, 330, 660]
SPLITS = ("dev1", "dev2", "dev3")
MAIN = ("R1", "H", "H-static", "F_r4", "F_r8", "F_r16", "F-static_r8", "C", "M_lam0.3")
N_T = 656                                                 # |T| of the frozen test split (DESIGN §2)


def _ints(s):
    out = []
    for tok in str(s).split(","):
        if "-" in tok:
            a, b = tok.split("-")
            out += list(range(int(a), int(b) + 1))
        elif tok.strip():
            out.append(int(tok))
    return out


def default_cps(cid: str) -> list[int]:
    s = q2dev.parse(cid)
    if s.warm:
        return [n for n in WS if n > s.warm]
    if s.batch_size or s.feasible or s.minimal:
        return SCEN_CPS
    if cid == "R1":
        return sorted(set(Q2) | {440})
    return Q2


# ============================================================================ run / smoke / jobs
def cmd_run(a) -> int:
    part = PARTS / f"{a.config}_{a.split}_r{a.traj}.jsonl"
    cps = _ints(a.checkpoints) if a.checkpoints else default_cps(a.config)
    if all(q2dev.done([part], a.config, a.split, a.traj, n) for n in cps):
        print("done already", part)
        return 0
    rows = q2dev.run_config(a.config, a.split, a.traj, cps, part, max_n=a.max_n or max(cps), phase="p4",
                            n_threads=int(a.threads))
    for r in rows:
        print(r["config"], r["split"], r["traj"], r["n"], f"{r['auprc']:.4f}", f"fit {r['fit_predict_s']:.1f}s",
              f"sel {r['select_s_cum']:.1f}s", f"off {r['offarray_share']:.2f}", flush=True)
    return 0


def cmd_smoke(a) -> int:
    sm = OUT / "smoke" / "smoke.jsonl"
    for cid in a.configs.split(","):
        s = q2dev.parse(cid)
        if s.warm:
            cps, mx = [s.warm + 20], s.warm + 20
        elif s.batch_size:
            cps, mx = [30, 60], 60                         # 30 = partial batch for b = 50
        else:
            cps, mx = [30], 30
        t0 = time.perf_counter()
        rows = q2dev.run_config(cid, a.split, 0, cps, sm, max_n=mx, phase="smoke", write_log=False)
        for r in rows:
            print(f"SMOKE {cid} n={r['n']} auprc={r['auprc']:.4f} off={r['offarray_share']:.2f} "
                  f"sel={r['select_s_cum']:.1f}s fit={r['fit_predict_s']:.1f}s wall={time.perf_counter() - t0:.0f}s",
                  flush=True)
    return 0


def cmd_jobs(a) -> int:
    d = OUT / "jobs"
    d.mkdir(parents=True, exist_ok=True)
    jobs = []
    for r in _ints(a.trajs):
        for split in a.splits.split(","):
            for cid in a.configs.split(","):
                cmd = [PY, "scripts/41_p4_strategies.py", "run", "--config", cid, "--split", split, "--traj", str(r)]
                if a.checkpoints:
                    cmd += ["--checkpoints", a.checkpoints]
                jobs.append({"name": f"p4_{cid}_{split}_r{r}", "n": 660, "cmd": cmd})
    p = d / f"{a.name}.json"
    p.write_text(json.dumps(jobs, indent=1) + "\n")
    print(p, len(jobs), "jobs")
    return 0


# ============================================================================ degree diagnostic (evaluation only)
def cmd_degdiag(a) -> int:
    from hsl import oracles, strategies
    b = views.dev_bundle()
    rows = load_rows()
    out = []
    for split in SPLITS:
        sp = views.get_split(split)
        tids = sorted({(x["traj"], x["traj_id"]) for x in rows if x["config"] == "H" and x["split"] == split})
        for r, tid in tids[: int(a.max_trajs)]:
            sel_list = json.loads((H.TRAJ_DIR / split / f"{tid}.json").read_text())["selected"]
            for n in _ints(a.ns):
                S = sel_list[:n]
                v = views.SelectorView.build(b, S, sp)
                h = strategies.make("H")
                q, phase = h.scores(v)
                rho = oracles.degree_diagnostic(v, q, split)
                on = np.asarray(v.candidate_col) >= 0
                deg = oracles.true_query_degree(split, v.candidates)
                from scipy.stats import spearmanr
                rho_on = float(spearmanr(q[on], deg[on])[0])
                rho_off = float(spearmanr(q[~on], deg[~on])[0]) if (~on).sum() > 2 else float("nan")
                out.append({"split": split, "traj": r, "n": n, "phase": phase, "spearman": rho,
                            "spearman_on_array": rho_on, "spearman_off_array": rho_off,
                            "n_candidates": int(len(v.candidates))})
                print(out[-1], flush=True)
    pd.DataFrame(out).to_csv(OUT / "degree_diagnostic.csv", index=False)
    return 0


# ============================================================================ analysis
def load_rows() -> list[dict]:
    return p1.read_rows(sorted(PARTS.glob("*.jsonl")))


def _f(x, k=4):
    try:
        return "—" if x is None or not np.isfinite(float(x)) else f"{float(x):.{k}f}"
    except (TypeError, ValueError):
        return str(x)


def _pool_offarray_base() -> dict:
    on = q2dev._on_array()
    return {s: float(np.mean([g not in on for g in views.get_split(s).pool])) for s in SPLITS}


def contrasts_for(cid, rows, rnd, grid) -> dict:
    per = {}
    for split in SPLITS:
        rs = [x for x in rows if x["config"] == cid and x["split"] == split]
        if not rs:
            continue
        As, es, tr = q2dev.arm_matrix(rs, split, grid, BOOT)
        Ar, er, _ = q2dev.arm_matrix(rnd, split, grid, BOOT)
        per[split] = {"cs": q2dev.contrast_suite(As, es, Ar, er, grid, f"p4_{cid}_{split}"), "R_s": len(tr),
                      "As": As, "es": es, "Ar": Ar, "er": er}
    return per


def power(per: dict, g: float, n_t: int = N_T) -> dict:
    """DESIGN §5/§7.1 rule with P4-measured inputs: V_Q(Δ) at n = 160 from the shared-weight replicates (scaled
    |D_val|/|T|), σ_s, σ_r pooled within split; smallest R_s (R_r = 2R_s) with SE(Δ̂(160)) ≤ g/3.5."""
    vq, ss, sr, nv = [], [], [], []
    for split, d in per.items():
        w = d["cs"]["delta"].get(160)
        if w is None:
            continue
        nval = len(views.get_split(split).targets)
        vq.append(w["V_Q"] * nval / n_t)
        nv.append(nval)
        if np.isfinite(w["s_sd"]):
            ss.append(w["s_sd"] ** 2)
        sr.append(w["r_sd"] ** 2)
    if not vq:
        return {}
    VQ = float(np.mean(vq))
    s2s = float(np.mean(ss)) if ss else float("nan")
    s2r = float(np.mean(sr))
    target = (g / 3.5) ** 2
    req = None
    for R in range(2, 201):
        if VQ + s2s / R + s2r / (2 * R) <= target:
            req = R
            break
    return {"V_Q_final": VQ, "sigma_s": float(np.sqrt(s2s)), "sigma_r": float(np.sqrt(s2r)), "target_se": g / 3.5,
            "VQ_exceeds_target": bool(VQ > target), "R_s_required": req,
            "se_at": {R: float(np.sqrt(VQ + s2s / R + s2r / (2 * R))) for R in (3, 5, 10, 20, 30, 40)}}


def cmd_analyze(a) -> int:
    rows = load_rows()
    rnd = q2dev.random_rows()
    for split in SPLITS:
        q2dev.check_est([x for x in rows if x["split"] == split], split, BOOT)
    cfgs = sorted({x["config"] for x in rows})
    main = [c for c in cfgs if "__" not in c]
    grid = sorted(set(Q2) | {440})
    res = {"predictor": q2dev.PSTAR, "configs": cfgs, "B": int(config.protocol()["bootstrap"]["B_query_bootstrap"])}
    # random-curve gap for power (frozen P*, dev): g = r210 − r160, mean over splits
    rdf = pd.DataFrame(rnd)
    rm = rdf.groupby(["split", "n"]).auprc.mean().unstack()
    g_by = {s: float(rm.loc[s, 210] - rm.loc[s, 160]) for s in rm.index if 210 in rm.columns and
            np.isfinite(rm.loc[s, 210])}
    g = float(np.mean(list(g_by.values()))) if g_by else float("nan")
    res["g_r210_r160"] = {"by_split": g_by, "mean": g}
    # ---- main-arm contrasts
    crow, ranking, pw = [], [], {}
    for cid in main:
        per = contrasts_for(cid, rows, rnd, grid)
        if not per:
            continue
        rel160 = []
        pairs = []
        for split, d in per.items():
            cs = d["cs"]
            for n, w in cs["delta"].items():
                crow.append({"config": cid, "split": split, "n": n, "R_s": d["R_s"], "R_r": w["R_b"],
                             "s_mean": w["s_mean"], "r_mean": w["r_mean"], "delta": w["delta"], "rel": w["rel"],
                             "lo": w["ci95"][0], "hi": w["ci95"][1], "p": w["p_one_sided"], "V_Q": w["V_Q"],
                             "s_sd": w["s_sd"], "r_sd": w["r_sd"], "df": w["df"]})
            if 160 in cs["delta"]:
                rel160.append(cs["delta"][160]["rel"])
            if "G" in cs:
                j160, j660 = grid.index(160), grid.index(660)
                pairs.append((d["As"][:, :, [j160]], d["Ar"][:, :, [j160, j660]], d["es"][:, [j160]],
                              d["er"][:, [j160, j660]]))
            for k in ("c1", "brief", "matches"):
                if k in cs:
                    crow.append({"config": cid, "split": split, "n": k, "R_s": d["R_s"], "delta": cs[k]["delta"],
                                 "lo": cs[k]["ci95"][0], "hi": cs[k]["ci95"][1], "p": cs[k]["p_one_sided"]})
            if "G" in cs:
                crow.append({"config": cid, "split": split, "n": "G", "R_s": d["R_s"], "delta": cs["G"]["G"],
                             "lo": cs["G"]["G_ci"][0], "hi": cs["G"]["G_ci"][1], "p": None,
                             "G_star": cs["G_star"], "G_defined": cs["G"]["defined"]})
        # headline-rule strata (DESIGN §3): Δ(160) and Δ(660) on E_AO and E_VV (dev analogue of E_TT)
        for split in per:
            rs = [x for x in rows if x["config"] == cid and x["split"] == split]
            for st in ("AO", "VV", "onrows"):
                g2 = [160, 660]
                As, es, _ = q2dev.arm_matrix(rs, split, g2, BOOT, stratum=st)
                Ar, er, _ = q2dev.arm_matrix(rnd, split, g2, BOOT, stratum=st)
                for n in g2:
                    j = g2.index(n)
                    if np.isfinite(es[:, j]).any() and np.isfinite(er[:, j]).sum() >= 2:
                        w = q2dev.welch_at(As[:, :, j], es[:, j], Ar[:, :, j], er[:, j])
                        crow.append({"config": cid, "split": split, "n": f"{st}{n}", "R_s": per[split]["R_s"],
                                     "s_mean": float(np.nanmean(es[:, j])), "r_mean": float(np.nanmean(er[:, j])),
                                     "delta": w["delta"], "rel": w["delta"] / float(np.nanmean(er[:, j])),
                                     "lo": w["ci95"][0], "hi": w["ci95"][1], "p": w["p_one_sided"]})
        Gp = q2dev.gap_closure_joint(pairs, np.random.default_rng(p1.boot_seed(f"p4_G_pooled_{cid}"))) \
            if len(pairs) == 3 else None
        arm = q2dev.parse(cid).arm
        ranking.append({"config": cid, "arm": arm, "family": q2dev.family(arm), "n_splits": len(rel160),
                        "mean_rel_delta160": float(np.mean(rel160)) if rel160 else float("nan"),
                        "min_rel_delta160": float(np.min(rel160)) if rel160 else float("nan"),
                        "rel_delta160_by_split": rel160,
                        "G_pooled": Gp["G"] if Gp else None, "G_pooled_lo": Gp["G_ci"][0] if Gp else None,
                        "G_pooled_hi": Gp["G_ci"][1] if Gp else None,
                        "R_min": min(d["R_s"] for d in per.values())})
        pw[cid] = power(per, g)
    cdf = pd.DataFrame(crow)
    cdf.to_csv(OUT / "contrasts.csv", index=False)
    rk = pd.DataFrame(ranking).sort_values("mean_rel_delta160", ascending=False)
    rk.to_csv(OUT / "ranking.csv", index=False)
    res["ranking"] = rk.to_dict("records")
    res["power"] = pw
    # ---- shortlist rule (§8): rank by mean over 3 splits of relative Δ(160); ≤ 1 per family (H, F, C/M); ≤ 3
    eligible = rk[(rk.family.isin(["H", "F", "C/M"])) & (rk.n_splits == 3) & (rk.R_min >= 3)]
    short = []
    for _, r in eligible.iterrows():
        if r.family not in {s["family"] for s in short} and len(short) < 3:
            short.append({"config": r.config, "family": r.family, "mean_rel_delta160": r.mean_rel_delta160})
    res["shortlist_by_rule"] = short
    # ---- diagnostics: off-array share, runtime
    base = _pool_offarray_base()
    res["offarray_base"] = base
    df = pd.DataFrame(rows)
    diag = (df.groupby(["config", "n"]).agg(offarray=("offarray_share", "mean"), offarray_new=("offarray_share_new",
                                                                                              "mean"),
                                            sl_detected=("sl_detected", "mean"), auprc=("auprc", "mean"),
                                            fit_s=("fit_predict_s", "mean"), select_s_cum=("select_s_cum", "mean"),
                                            n_rounds=("n_rounds", "mean")).reset_index())
    diag.to_csv(OUT / "diagnostics.csv", index=False)
    tw = df.groupby(["config", "split", "traj"]).agg(wall=("traj_wall_s", "first"), sel=("select_s_cum", "max"),
                                                      rounds=("n_rounds", "max"), fit=("fit_predict_s", "sum"))
    rt = tw.groupby("config").agg(traj_wall_s=("wall", "mean"), select_s_total=("sel", "mean"),
                                  rounds=("rounds", "mean"), fit_s_total=("fit", "mean"),
                                  n_traj=("wall", "size")).reset_index()
    rt["select_s_per_round"] = rt.select_s_total / rt.rounds
    rt.to_csv(OUT / "runtime.csv", index=False)
    # ---- scenarios
    res["scenarios"] = scenarios(rows, rnd, cfgs)
    if (OUT / "degree_diagnostic.csv").exists():
        res["degree_diagnostic"] = pd.read_csv(OUT / "degree_diagnostic.csv").groupby("n")[
            ["spearman", "spearman_on_array", "spearman_off_array"]].mean().reset_index().to_dict("records")
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=_json) + "\n")
    _figures(cdf, df, rnd)
    _summary(res, cdf, rk, diag, rt, base)
    if a.ledger:
        for _, r in rk.iterrows():
            sub = cdf[(cdf.config == r.config) & (cdf.n.astype(str) == "160")]
            q2dev.ledger({"id": f"p4_{r.config}", "parent": f"R0_{q2dev.PSTAR}",
                             "config_hash": H._hash(["p4", r.config, q2dev.parse(r.config).knobs]),
                             "split": "dev1-3", "kind": "p4_strategy",
                             "outcome": "shortlisted" if r.config in [s["config"] for s in short] else "not_shortlisted",
                             "metrics": {"mean_rel_delta160": r.mean_rel_delta160,
                                         "per_split_160": sub[["split", "delta", "lo", "hi", "rel"]].to_dict("records"),
                                         "G_pooled": r.G_pooled, "power": pw.get(r.config)}})
        for cid, sc in res["scenarios"].items():
            q2dev.ledger({"id": f"p4_{cid}", "parent": q2dev.parse(cid).arm, "config_hash": H._hash(["p4", cid]),
                             "split": "dev1-3", "kind": "p4_scenario", "outcome": "descriptive", "metrics": sc})
    print((OUT / "summary.md").read_text())
    return 0


def _json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def scenarios(rows, rnd, cfgs) -> dict:
    """Descriptive scenario contrasts (per split + 3-split mean): det vs randomised; warm start vs R1-continuation
    and random continuation at 440/660; batch 50/110 vs 10; feasible pool vs full; minimal vs full features."""
    out = {}
    grid = sorted(set(Q2) | {440})
    for cid in [c for c in cfgs if "__" in c]:
        s = q2dev.parse(cid)
        base_id = cid.split("__")[0]
        if s.warm:
            if base_id == "R0":
                continue
            refs = {"R1_continuation": "R1", "random_continuation": "R0__ws220"}
            ns = [440, 660]
        elif s.feasible and base_id == "R0":
            refs, ns = {"R0_full": "R0"}, SCEN_CPS
        else:
            refs = {"same_arm_batch10_full": base_id, "random": "R0"}
            ns = Q2 if s.deterministic else SCEN_CPS
        rec = {}
        for rname, rid in refs.items():
            for split in SPLITS:
                rs = [x for x in rows if x["config"] == cid and x["split"] == split]
                rr = rnd if rid == "R0" else [x for x in rows if x["config"] == rid]
                if not rs or not [x for x in rr if x["split"] == split]:
                    continue
                As, es, _ = q2dev.arm_matrix(rs, split, grid, BOOT)
                Ar, er, _ = q2dev.arm_matrix(rr, split, grid, BOOT)
                for n in ns:
                    j = grid.index(n)
                    if np.isfinite(es[:, j]).sum() >= 1 and np.isfinite(er[:, j]).sum() >= 2:
                        w = q2dev.welch_at(As[:, :, j], es[:, j], Ar[:, :, j], er[:, j])
                        rec.setdefault(rname, {}).setdefault(str(n), {})[split] = {
                            "s": float(np.nanmean(es[:, j])), "ref": float(np.nanmean(er[:, j])),
                            "delta": w["delta"], "lo": w["ci95"][0], "hi": w["ci95"][1], "p": w["p_one_sided"],
                            "R_s": int(np.isfinite(es[:, j]).sum()), "R_ref": int(np.isfinite(er[:, j]).sum())}
            for n, d in rec.get(rname, {}).items():
                if len(d) == 3:
                    d["mean3"] = {"delta": float(np.mean([v["delta"] for k, v in d.items() if k in SPLITS])),
                                  "rel": float(np.mean([v["delta"] / v["ref"] for k, v in d.items()
                                                        if k in SPLITS]))}
        out[cid] = rec
    return out


def _figures(cdf, df, rnd) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = cdf[pd.to_numeric(cdf.n, errors="coerce").notna()].copy()
    d["n"] = d.n.astype(int)
    m = d.groupby(["config", "n"]).agg(rel=("rel", "mean"), lo=("lo", "mean"), hi=("hi", "mean"),
                                       r=("r_mean", "mean")).reset_index()
    cfgs = [c for c in MAIN if c in set(m.config)] + sorted(set(m.config) - set(MAIN))
    pal = ["#6b7280", "#d1495b", "#f4a3ad", "#2a78d6", "#1f4e9c", "#7fb3ef", "#9ec5f0", "#1baf7a", "#eda100",
           "#8e44ad", "#16a085", "#c0392b", "#7f8c8d"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, c in enumerate(cfgs):
        x = m[m.config == c].sort_values("n")
        ax.plot(x.n, 100 * x.rel, marker="o", ms=3, lw=1.6, color=pal[i % len(pal)], label=c)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks(Q2)
    ax.set_xticklabels([str(n) for n in Q2])
    ax.set_xlabel("queries screened n")
    ax.set_ylabel("relative Δ vs random at matched n (%), mean of dev1–3")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_p4_rel_delta.png", dpi=150)
    plt.close(fig)
    # absolute curves per split, main arms + random
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    rd = pd.DataFrame(rnd)
    for ax, split in zip(axes, SPLITS):
        r = rd[rd.split == split].groupby("n").auprc.mean()
        r = r[r.index <= 660]
        ax.plot(r.index, r.values, color="k", lw=2, marker="s", ms=3, label="random R0 (R = 10)")
        for i, c in enumerate(cfgs):
            if "__" in c:
                continue
            x = df[(df.config == c) & (df.split == split)].groupby("n").auprc.mean()
            ax.plot(x.index, x.values, marker="o", ms=2.5, lw=1.2, color=pal[i % len(pal)], label=c)
        ax.set_xscale("log")
        ax.set_title(split)
        ax.set_xlabel("n")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(f"AUPRC E_fixed_dev ({q2dev.PSTAR})")
    axes[0].legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_p4_curves.png", dpi=150)
    plt.close(fig)
    # discovery frontier
    fig, ax = plt.subplots(figsize=(7, 4.5))
    rr = rd[rd.n.isin(Q2)].groupby("n")[["auprc", "sl_detected"]].mean()
    ax.plot(rr.sl_detected, rr.auprc, color="k", marker="s", ms=3, lw=2, label="random")
    for i, c in enumerate(cfgs):
        if "__" in c:
            continue
        x = df[df.config == c].groupby("n")[["auprc", "sl_detected"]].mean()
        ax.plot(x.sl_detected, x.auprc, marker="o", ms=2.5, lw=1.2, color=pal[i % len(pal)], label=c)
    ax.set_xlabel("strong-SL cells detected in screened rows (mean of dev1–3)")
    ax.set_ylabel("AUPRC E_fixed_dev")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_p4_frontier.png", dpi=150)
    plt.close(fig)


def _summary(res, cdf, rk, diag, rt, base) -> None:
    L = ["# P4 strategy development (DEV ONLY; dev splits 1–3)", "",
         f"Generated by `scripts/41_p4_strategies.py analyze`. Predictor: frozen P* `{q2dev.PSTAR}`; E_fixed_dev "
         f"pooled AUPRC; random arm = P1 frozen-P* rows (R = 10 per split; never recomputed). Inference: B = {res['B']} "
         f"shared stratified query-bootstrap replicates per split (weights '{q2dev.W_NAME}'), Welch SE "
         "(V_Q + σ̂²_s/R_s + σ̂²_r/R_r), G with the joint bootstrap (DESIGN §6.4). Batch 10 unless stated. "
         "All numbers are dev-split development estimates, not test-set claims.", ""]
    L += ["## Ranking (DESIGN §8: mean over 3 splits of relative Δ(160) = (s160 − r160)/r160)", "",
          "| config | family | R_min | rel Δ(160) mean | per split | min | G pooled [95% CI] |", "|---|---|---|---|---|---|---|"]
    for _, r in rk.iterrows():
        L.append(f"| {r.config} | {r.family} | {r.R_min} | {100 * r.mean_rel_delta160:+.1f}% | "
                 + ", ".join(f"{100 * v:+.1f}%" for v in r.rel_delta160_by_split)
                 + f" | {100 * r.min_rel_delta160:+.1f}% | {_f(r.G_pooled, 2)} [{_f(r.G_pooled_lo, 2)}, "
                   f"{_f(r.G_pooled_hi, 2)}] |")
    L += ["", "Shortlist by the §8 rule (≤ 1 per family H / F / C-M, ≤ 3): "
          + (", ".join(f"**{s['config']}** ({s['family']}, {100 * s['mean_rel_delta160']:+.1f}%)"
                       for s in res["shortlist_by_rule"]) or "none") + ". See `shortlist_proposal.md`.", ""]
    d = cdf[cdf.n.astype(str) == "160"]
    L += ["## Δ(160) per split (Welch, 95% CI, one-sided p)", "",
          "| config | split | R_s | s160 | r160 | Δ | 95% CI | p | V_Q | σ_s | σ_r |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in d.sort_values(["config", "split"]).iterrows():
        L.append(f"| {r.config} | {r.split} | {r.R_s} | {_f(r.s_mean)} | {_f(r.r_mean)} | {_f(r.delta)} | "
                 f"[{_f(r.lo)}, {_f(r.hi)}] | {_f(r.p, 3)} | {r.V_Q:.2e} | {_f(r.s_sd)} | {_f(r.r_sd)} |")
    for k, title in (("c1", "c1-style: s160 − r490"), ("brief", "brief example: s160 − r660"),
                     ("matches", "'matches' contrast: s160 − 0.95·r660 (H1 > 0)"), ("G", "gap closure G (per split)")):
        dd = cdf[cdf.n.astype(str) == k]
        if not len(dd):
            continue
        L += ["", f"## {title}", "", "| config | dev1 | dev2 | dev3 |", "|---|---|---|---|"]
        for c, x in dd.groupby("config"):
            cells = []
            for s in SPLITS:
                y = x[x.split == s]
                if len(y):
                    y = y.iloc[0]
                    extra = f" p={_f(y.p, 3)}" if k != "G" else f" G*={_f(y.get('G_star'), 2)}"
                    cells.append(f"{_f(y.delta, 4 if k != 'G' else 2)} [{_f(y.lo, 4 if k != 'G' else 2)}, "
                                 f"{_f(y.hi, 4 if k != 'G' else 2)}]{extra}")
                else:
                    cells.append("—")
            L.append(f"| {c} | " + " | ".join(cells) + " |")
    L += ["", "## Headline-rule strata: relative Δ at n = 160 / 660 (mean of dev1–3; per-split CIs in contrasts.csv)",
          "", "E_AO = array-only columns; E_VV = on-array D_val columns (dev analogue of E_TT); onrows = on-array target "
          "rows only (removes the off-array query artefact).", "",
          "| config | all 160 | AO 160 | VV 160 | onrows 160 | all 660 | AO 660 | VV 660 | onrows 660 | same sign at 160 (all, AO, VV) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for c, x in cdf.groupby("config"):
        v = x.groupby(x.n.astype(str)).rel.mean()
        keys = ["160", "AO160", "VV160", "onrows160", "660", "AO660", "VV660", "onrows660"]
        sg = [np.sign(v[k]) for k in ("160", "AO160", "VV160") if k in v.index]
        L.append(f"| {c} | " + " | ".join(f"{100 * v[k]:+.1f}%" if k in v.index else "—" for k in keys)
                 + f" | {'yes' if len(sg) == 3 and len(set(sg)) == 1 else 'no'} |")
    L += ["", "## Relative Δ(n) vs random at matched n (mean of dev1–3)", "",
          "| config | " + " | ".join(str(n) for n in Q2) + " |", "|---|" + "---|" * len(Q2)]
    dn = cdf[pd.to_numeric(cdf.n, errors="coerce").notna()].copy()
    dn["n"] = dn.n.astype(int)
    for c, x in dn.groupby("config"):
        v = x.groupby("n").rel.mean()
        L.append(f"| {c} | " + " | ".join(f"{100 * v[n]:+.1f}%" if n in v.index else "—" for n in Q2) + " |")
    pw = res["power"]
    L += ["", "## Power and R (DESIGN §5, §7.1; P4-measured)", "",
          f"g = r210 − r160 (frozen-P* dev random curves, mean of splits) = {_f(res['g_r210_r160']['mean'])}; target "
          f"SE(Δ̂(160)) ≤ g/3.5 = {_f(res['g_r210_r160']['mean'] / 3.5)}. V_Q scaled to the test set as "
          f"V_Q,dev·|D_val|/|T| (|T| = {N_T}).", "",
          "| config | V_Q(Δ) final-scaled | σ_s(160) | σ_r(160) | V_Q > target² | R_s required | SE at R_s = 5 / 10 / 20 / 40 |",
          "|---|---|---|---|---|---|---|"]
    for c, p in pw.items():
        if not p:
            continue
        L.append(f"| {c} | {p['V_Q_final']:.2e} | {_f(p['sigma_s'])} | {_f(p['sigma_r'])} | {p['VQ_exceeds_target']} | "
                 f"{p['R_s_required'] or '> 200 / bounds only'} | "
                 + " / ".join(_f(p["se_at"][R]) for R in (5, 10, 20, 40)) + " |")
    L += ["", "## Selector diagnostics", "", "Off-array share of the screened set (base = off-array share of the pool: "
          + ", ".join(f"{s} {v:.3f}" for s, v in base.items()) + ")", "",
          "| config | n=30 | n=160 | n=660 |", "|---|---|---|---|"]
    for c, x in diag.groupby("config"):
        v = dict(zip(x.n, x.offarray))
        L.append(f"| {c} | " + " | ".join(_f(v.get(n), 3) for n in (30, 160, 660)) + " |")
    if "degree_diagnostic" in res:
        L += ["", "H degree model: Spearman(q̂, true query degree) over the remaining pool along H trajectories "
              "(evaluation-only diagnostic, `degree_diagnostic.csv`):", "",
              "| n | all | on-array | off-array |", "|---|---|---|---|"]
        for r in res["degree_diagnostic"]:
            L.append(f"| {int(r['n'])} | {_f(r['spearman'], 3)} | {_f(r['spearman_on_array'], 3)} | "
                     f"{_f(r['spearman_off_array'], 3)} |")
    L += ["", "## Runtime (2 BLAS threads, shared 8-core M2 under load; per trajectory)", "",
          "| config | trajectories | wall s | P* fit+predict s (sum over checkpoints) | selection s total | s/round |",
          "|---|---|---|---|---|---|"]
    for _, r in rt.iterrows():
        L.append(f"| {r.config} | {r.n_traj} | {r.traj_wall_s:.0f} | {r.fit_s_total:.0f} | {r.select_s_total:.1f} | "
                 f"{r.select_s_per_round:.2f} |")
    sc = res["scenarios"]
    if sc:
        L += ["", "## Scenarios (descriptive; Δ = scenario − reference, mean of available splits)", "",
              "| scenario | reference | n | Δ per split [95% CI] | mean Δ (rel) |", "|---|---|---|---|---|"]
        for cid, rec in sc.items():
            for ref, byn in rec.items():
                for n, d in byn.items():
                    cells = "; ".join(f"{s}: {_f(v['delta'])} [{_f(v['lo'])}, {_f(v['hi'])}] (R {v['R_s']})"
                                      for s, v in d.items() if s in SPLITS)
                    m3 = d.get("mean3")
                    L.append(f"| {cid} | {ref} | {n} | {cells} | "
                             + (f"{_f(m3['delta'])} ({100 * m3['rel']:+.1f}%)" if m3 else "—") + " |")
    L += ["", "Files: `contrasts.csv` (every Welch contrast), `ranking.csv`, `diagnostics.csv` (off-array share, SLs "
          "detected, runtime by n), `runtime.csv`, `summary.json`; figures `fig_p4_rel_delta.png`, `fig_p4_curves.png`, "
          "`fig_p4_frontier.png`.", ""]
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--config", required=True)
    r.add_argument("--split", required=True, choices=list(SPLITS))
    r.add_argument("--traj", type=int, required=True)
    r.add_argument("--checkpoints", default="")
    r.add_argument("--max-n", type=int, default=None)
    r.add_argument("--threads", default=2)
    s = sub.add_parser("smoke")
    s.add_argument("--configs", default=",".join(MAIN))
    s.add_argument("--split", default="dev1")
    j = sub.add_parser("jobs")
    j.add_argument("--name", required=True)
    j.add_argument("--configs", required=True)
    j.add_argument("--splits", default=",".join(SPLITS))
    j.add_argument("--trajs", default="0-2")
    j.add_argument("--checkpoints", default="")
    dd = sub.add_parser("degdiag")
    dd.add_argument("--ns", default="0,30,160,330,630")
    dd.add_argument("--max-trajs", default=3)
    an = sub.add_parser("analyze")
    an.add_argument("--ledger", action="store_true")
    a = ap.parse_args(argv)
    return {"run": cmd_run, "smoke": cmd_smoke, "jobs": cmd_jobs, "degdiag": cmd_degdiag,
            "analyze": cmd_analyze}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
