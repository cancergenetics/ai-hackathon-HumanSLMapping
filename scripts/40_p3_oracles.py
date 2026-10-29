#!/usr/bin/env python
"""P3 ceilings (DESIGN v3.1 §9 P3; DEV DATA ONLY): evaluation-only oracles vs the frozen-P* random arm + the §11
end-to-end controls.

  oracle_hub     hsl.oracles.OracleHubFirst: candidates ranked by TRUE query-row SL count (static); dev1-3,
                 trajectories r = 0..R-1 (tie order + model seed from dev trajectory seed r), Q2 checkpoints.
  oracle_greedy  hsl.oracles.GreedyOracle: per round of 10, 200 random remaining candidates (selection stream of
                 trajectory r), each scored by the realised E_fixed_dev AUPRC of the closed-form fold-in
                 F2_mk_foldin on S ∪ {c}; top 10 by gain. dev1 only, to n = 160; P* evaluated at {30, 60, 100, 160}.
  Both are evaluated with the frozen P* (config/predictor.yaml) at the checkpoints -> results/p3/parts/*.jsonl.
  Random arm (R0): P1's frozen-P* rows (results/p1/parts; n = 30/160/330/660 from the P2 reuse, the rest P1 batch4).

Go/no-go: G = (s160 − r160)/(r660 − r160) with the joint bootstrap (shared stratified query weights, trajectories
resampled within arm). If both oracles have G < 0.5 at n = 160, C and M are de-scoped.

    run       --config oracle_hub --split dev1 --traj 0 [--checkpoints 30,60,...] [--max-n 660]
    controls  shuffled screened labels -> AUPRC ≈ prevalence; oracle back door -> ≈ 1 (P*, dev1, n = 160)
    jobs      write results/p3/jobs/*.json (launch with scripts/launch.py)
    analyze   -> results/p3/summary.md, ceilings.csv, curves.csv, fig_p3_ceilings.png, summary.json [--ledger]
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd

from hsl import config, harness as H, p1, q2dev, views

OUT = q2dev.RESULTS / "p3"
PARTS = OUT / "parts"
BOOT = OUT / "boot"
PY = ".venv/bin/python"
Q2 = [int(n) for n in config.budgets()["q2_checkpoints"]]
GREEDY_CPS = [n for n in Q2 if n <= 160]
ORACLES = ("oracle_hub", "oracle_greedy", "oracle_greedyHO")
STRATUM_PLAN = (("oracle_hub", ("AO", "VV", "onrows")), ("oracle_greedy", ("halfA", "halfB")),
                ("oracle_greedyHO", ("halfA", "halfB")))
STRATUM_LABEL = {"AO": "E_AO (array-only columns)", "VV": "E_VV (on-array D_val columns; dev E_TT analogue)",
                 "onrows": "on-array target rows only", "halfA": "query half A (held-out oracle's selection half)",
                 "halfB": "query half B (held out from the held-out oracle)"}


def _ints(s):
    return [int(x) for x in str(s).split(",") if x.strip()]


# ============================================================================ run
def cmd_run(a) -> int:
    part = PARTS / f"{a.config}_{a.split}_r{a.traj}.jsonl"
    cps = _ints(a.checkpoints) if a.checkpoints else (GREEDY_CPS if a.config.startswith("oracle_greedy") else Q2)
    if all(q2dev.done([part], a.config, a.split, a.traj, n) for n in cps):
        print("done already", part)
        return 0
    extra = OUT / "greedy_logs" / f"{a.config}_{a.split}_r{a.traj}.json" if a.config.startswith("oracle_greedy") else None
    rows = q2dev.run_config(a.config, a.split, a.traj, cps, part, max_n=a.max_n or max(cps), phase="p3",
                            n_threads=int(a.threads), extra_log=extra)
    for r in rows:
        print(r["config"], r["split"], r["traj"], r["n"], f"{r['auprc']:.4f}", f"fit {r['fit_predict_s']:.1f}s",
              f"sel {r['select_s_cum']:.1f}s", flush=True)
    return 0


# ============================================================================ end-to-end controls (§11)
def cmd_controls(a) -> int:
    from hsl import upstream as U
    b = views.dev_bundle()
    t = p1.random_traj("dev1", 0)
    S = tuple(t.perm[:int(a.n)])
    rng = np.random.default_rng(config.seed("p3_e2e_control_seed"))
    r = b.gene_rows(np.asarray(S))
    E_, SP, Y = (np.array(getattr(b, k), copy=True) for k in ("E", "SP", "Y"))
    blk = np.ix_(r, np.arange(b.n_arrays))
    e, sp_, y = E_[blk].ravel(), SP[blk].ravel(), Y[blk].ravel()
    obs = np.flatnonzero(y >= 0)
    perm = obs[rng.permutation(len(obs))]
    for M, v in ((E_, e), (SP, sp_), (Y, y)):
        vv = v.copy()
        vv[obs] = v[perm]
        M[blk] = vv.reshape(len(r), -1)
    sb = b.with_matrices(E=E_, SP=SP, Y=Y)
    sb.array_smf = np.array(b.array_smf)[rng.permutation(b.n_arrays)]
    out = {"split": "dev1", "n": len(S), "traj": 0, "predictor": q2dev.PSTAR, "params": q2dev.PSTAR_PARAMS}
    t0 = time.perf_counter()
    res = H.run_predictor(H.RunSpec("dev1", S), q2dev.PSTAR, q2dev.PSTAR_PARAMS, {"model_seed": t.model_seed},
                          bundle=sb, save=False, ledger=False, n_threads=int(a.threads), tag="p3_control_shuffled")
    m = res.metrics["E_fixed_dev"]
    out["shuffled"] = {"auprc": m["auprc"], "prevalence": m["prevalence"], "ratio": m["auprc"] / m["prevalence"],
                       "pass": bool(0.7 * m["prevalence"] < m["auprc"] < 1.3 * m["prevalence"]),
                       "wall_s": time.perf_counter() - t0}
    del sb, E_, SP, Y

    def oracle_backdoor(view, bb):
        S_, _ = U.run_method(view, q2dev.PSTAR, model_seed=t.model_seed, n_threads=int(a.threads),
                             n_folds=10, params=q2dev.PSTAR_PARAMS)
        Yt = bb.Y[bb.gene_rows(view.target_orfs)].astype(np.float32)
        return S_ + 10.0 * np.clip(Yt, 0, 1)
    t0 = time.perf_counter()
    res = H.run_predictor(H.RunSpec("dev1", S), q2dev.PSTAR, q2dev.PSTAR_PARAMS, {"model_seed": t.model_seed},
                          score_fn=oracle_backdoor, save=False, ledger=False, n_threads=int(a.threads),
                          tag="p3_control_backdoor")
    m = res.metrics["E_fixed_dev"]
    out["backdoor"] = {"auprc": m["auprc"], "pass": bool(m["auprc"] > 0.99), "wall_s": time.perf_counter() - t0}
    # reference: the unshuffled run of the same S (P1 row, same model seed)
    ref = [x for x in q2dev.random_rows() if x["split"] == "dev1" and x["traj"] == 0 and x["n"] == len(S)]
    out["reference_auprc"] = ref[0]["auprc"] if ref else None
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "controls.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


# ============================================================================ jobs
def cmd_jobs(a) -> int:
    d = OUT / "jobs"
    d.mkdir(parents=True, exist_ok=True)
    jobs = []
    for r in _ints(a.trajs):
        for split in a.splits.split(","):
            if a.which in ("hub", "all"):
                jobs.append({"name": f"p3_oracle_hub_{split}_r{r}", "n": 660,
                             "cmd": [PY, "scripts/40_p3_oracles.py", "run", "--config", "oracle_hub", "--split",
                                     split, "--traj", str(r)]})
        if a.which in ("greedy", "all"):
            jobs.append({"name": f"p3_oracle_greedy_dev1_r{r}", "n": 160, "timeout_s": 6 * 3600,
                         "cmd": [PY, "scripts/40_p3_oracles.py", "run", "--config", "oracle_greedy", "--split",
                                 "dev1", "--traj", str(r)]})
    if a.which in ("controls", "all"):
        jobs.append({"name": "p3_controls", "n": 160, "cmd": [PY, "scripts/40_p3_oracles.py", "controls"]})
    p = d / f"{a.name}.json"
    p.write_text(json.dumps(jobs, indent=1) + "\n")
    print(p, len(jobs), "jobs")
    return 0


# ============================================================================ analysis
def load_rows() -> list[dict]:
    return p1.read_rows(sorted(PARTS.glob("*.jsonl")))


def cmd_analyze(a) -> int:
    rows = load_rows()
    rnd = q2dev.random_rows()
    grid = Q2
    res = {"predictor": q2dev.PSTAR, "weights": q2dev.W_NAME, "B": int(config.protocol()["bootstrap"]
                                                                         ["B_query_bootstrap"]), "arms": {}}
    for split in ("dev1", "dev2", "dev3"):
        n_chk = q2dev.check_est([x for x in rows + rnd if x["n"] in grid], split, BOOT)
        res.setdefault("engine_checks", {})[split] = n_chk
    curve = []
    for cfg, d in [("R0", rnd)] + [(c, [x for x in rows if x["config"] == c]) for c in ORACLES]:
        for x in d:
            if x["n"] in grid:
                curve.append({k: x.get(k) for k in ("config", "split", "traj", "n", "auprc", "auprc_AO", "auprc_VV",
                                                    "auprc_on_array_rows", "auprc_off_array_rows", "macro_auprc",
                                                    "sl_detected", "offarray_share", "fit_predict_s", "select_s_cum",
                                                    "traj_wall_s")})
    cdf = pd.DataFrame(curve)
    tab = (cdf.groupby(["config", "split", "n"]).agg(R=("auprc", "size"), auprc=("auprc", "mean"),
                                                      sd=("auprc", "std"), auprc_AO=("auprc_AO", "mean"),
                                                      auprc_VV=("auprc_VV", "mean"),
                                                      sl_detected=("sl_detected", "mean"),
                                                      fit_s=("fit_predict_s", "mean"),
                                                      select_s_cum=("select_s_cum", "mean")).reset_index())
    tab.to_csv(OUT / "curves.csv", index=False)
    # ---- contrasts and G
    ceil = []
    G_pairs = {c: [] for c in ORACLES}
    for cfg in ORACLES:
        for split in ("dev1", "dev2", "dev3"):
            rs = [x for x in rows if x["config"] == cfg and x["split"] == split]
            if not rs:
                continue
            g_s = [n for n in grid if any(x["n"] == n for x in rs)]
            As, es, tr = q2dev.arm_matrix(rs, split, grid, BOOT)
            Ar, er, _ = q2dev.arm_matrix(rnd, split, grid, BOOT)
            cs = q2dev.contrast_suite(As, es, Ar, er, grid, f"p3_{cfg}_{split}")
            j160, j660 = grid.index(160), grid.index(660)
            G_pairs[cfg].append((As[:, :, [j160]], Ar[:, :, [j160, j660]], es[:, [j160]], er[:, [j160, j660]]))
            rec = {"config": cfg, "split": split, "R_s": len(tr), "R_r": int(np.isfinite(er[:, j160]).sum()),
                   "checkpoints": g_s}
            for n in grid:
                if np.isfinite(er[:, grid.index(n)]).any():
                    rec[f"r{n}"] = float(np.nanmean(er[:, grid.index(n)]))
            for n, w in cs["delta"].items():
                rec[f"s{n}"], rec[f"r{n}"], rec[f"delta{n}"] = w["s_mean"], w["r_mean"], w["delta"]
                rec[f"delta{n}_lo"], rec[f"delta{n}_hi"] = w["ci95"]
                rec[f"rel{n}"] = w["rel"]
            for k in ("c1", "brief", "matches"):
                if k in cs:
                    rec[f"{k}_est"], rec[f"{k}_lo"], rec[f"{k}_hi"], rec[f"{k}_p"] = (
                        cs[k]["delta"], cs[k]["ci95"][0], cs[k]["ci95"][1], cs[k]["p_one_sided"])
            if "G" in cs:
                rec.update(G=cs["G"]["G"], G_lo=cs["G"]["G_ci"][0], G_hi=cs["G"]["G_ci"][1],
                           G_defined=cs["G"]["defined"], G_star=cs["G_star"])
            # runtime
            sub = cdf[(cdf.config == cfg) & (cdf.split == split)]
            rec["traj_wall_s_mean"] = float(sub.groupby("traj").traj_wall_s.first().mean())
            rec["select_s_to160"] = float(sub[sub.n == 160].select_s_cum.mean()) if (sub.n == 160).any() else None
            rec["fit_s_mean"] = float(sub.fit_predict_s.mean())
            ceil.append(rec)
    pooled = {}
    for cfg, pairs in G_pairs.items():
        if pairs:
            rng = np.random.default_rng(p1.boot_seed(f"p3_G_pooled_{cfg}"))
            pooled[cfg] = q2dev.gap_closure_joint(pairs, rng)
            pooled[cfg]["n_splits"] = len(pairs)
    cel = pd.DataFrame(ceil)
    cel.to_csv(OUT / "ceilings.csv", index=False)
    res["ceilings"] = ceil
    res["G_pooled"] = pooled
    gmin = float(config.protocol()["stats"]["g_star_min"])
    Gs = {c: pooled[c]["G"] for c in pooled if c in ("oracle_hub", "oracle_greedy")}
    res["go"] = {"rule": f"de-scope C/M iff both oracles have G < {gmin} at n = 160",
                 "G": Gs, "descope_CM": bool(len(Gs) == 2 and all(g < gmin for g in Gs.values())),
                 "complete": bool(len(Gs) == 2)}
    res["strata"] = stratum_ceilings(rows, rnd, grid)
    if (OUT / "controls.json").exists():
        res["controls"] = json.loads((OUT / "controls.json").read_text())
    (OUT / "summary.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    _figure(tab)
    _summary(res, tab, cel)
    if a.ledger:
        for cfg in pooled:
            q2dev.ledger({"id": f"p3_{cfg}", "parent": f"R0_{q2dev.PSTAR}", "config_hash": H._hash(["p3", cfg]),
                             "split": "dev1-3" if cfg == "oracle_hub" else "dev1", "kind": "p3_oracle",
                             "outcome": "ceiling", "metrics": {"G_pooled": pooled[cfg],
                                                               "per_split": [r for r in ceil if r["config"] == cfg]},
                             "decision": res["go"]})
        if "controls" in res:
            c = res["controls"]
            q2dev.ledger({"id": "p3_e2e_controls", "parent": None, "config_hash": H._hash(["p3_controls"]),
                             "split": "dev1", "kind": "p3_control",
                             "outcome": "pass" if c["shuffled"]["pass"] and c["backdoor"]["pass"] else "fail",
                             "metrics": c})
    print((OUT / "summary.md").read_text())
    return 0


def stratum_ceilings(rows, rnd, grid) -> list[dict]:
    """G(160) and s160/r160/r660 per oracle on cell strata (PI request after P3): the hub ceiling on E_AO, E_VV and
    on-array target rows (separates the off-array artefact); both greedy oracles on query halves A/B (the held-out
    oracle selected on half A only, so its half-B G is the attainable-ceiling estimate)."""
    out = []
    grid = [160, 660]
    j160, j660 = 0, 1
    for cfg, strata in STRATUM_PLAN:
        for st in strata:
            pairs = []
            for split in ("dev1", "dev2", "dev3"):
                rs = [x for x in rows if x["config"] == cfg and x["split"] == split]
                if not rs:
                    continue
                As, es, tr = q2dev.arm_matrix(rs, split, grid, BOOT, stratum=st)
                Ar, er, _ = q2dev.arm_matrix(rnd, split, grid, BOOT, stratum=st)
                if not np.isfinite(es[:, j160]).any():
                    continue
                p = (As[:, :, [j160]], Ar[:, :, [j160, j660]], es[:, [j160]], er[:, [j160, j660]])
                pairs.append(p)
                g1 = q2dev.gap_closure_joint([p], np.random.default_rng(p1.boot_seed(f"p3_G_{cfg}_{st}_{split}")))
                w = q2dev.welch_at(As[:, :, j160], es[:, j160], Ar[:, :, j160], er[:, j160])
                out.append({"config": cfg, "stratum": st, "split": split, "R_s": len(tr),
                            "s160": float(np.nanmean(es[:, j160])), "r160": float(np.nanmean(er[:, j160])),
                            "r660": float(np.nanmean(er[:, j660])), "delta160": w["delta"],
                            "delta160_lo": w["ci95"][0], "delta160_hi": w["ci95"][1],
                            "G": g1["G"], "G_lo": g1["G_ci"][0], "G_hi": g1["G_ci"][1]})
            if len(pairs) > 1:
                gp = q2dev.gap_closure_joint(pairs, np.random.default_rng(p1.boot_seed(f"p3_G_{cfg}_{st}_pooled")))
                out.append({"config": cfg, "stratum": st, "split": "pooled", "G": gp["G"], "G_lo": gp["G_ci"][0],
                            "G_hi": gp["G_ci"][1]})
    pd.DataFrame(out).to_csv(OUT / "ceilings_strata.csv", index=False)
    return out


def _figure(tab: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"R0": "#6b7280", "oracle_hub": "#d1495b", "oracle_greedy": "#2a78d6", "oracle_greedyHO": "#1baf7a"}
    lab = {"R0": "random (R0, R = 10)", "oracle_hub": "oracle hub-first (true degree)",
           "oracle_greedy": "greedy oracle, in-sample (upper bound)",
           "oracle_greedyHO": "greedy oracle, held-out halves"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, split in zip(axes, ("dev1", "dev2", "dev3")):
        for cfg in ("R0",) + ORACLES:
            d = tab[(tab.config == cfg) & (tab.split == split)].sort_values("n")
            if len(d):
                ax.errorbar(d.n, d.auprc, yerr=d.sd.fillna(0), marker="o", ms=4, lw=1.5, capsize=2,
                            color=col[cfg], label=lab[cfg])
        ax.set_xscale("log")
        ax.set_xticks([30, 60, 100, 160, 330, 660])
        ax.set_xticklabels(["30", "60", "100", "160", "330", "660"])
        ax.set_title(f"{split} (E_fixed_dev)")
        ax.set_xlabel("queries screened n")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(f"AUPRC ({q2dev.PSTAR})")
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_p3_ceilings.png", dpi=150)
    plt.close(fig)


def _f(x, k=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{k}f}"


def _summary(res, tab, cel) -> None:
    L = ["# P3 ceilings: oracle hub-first and greedy oracle vs random (DEV ONLY)", "",
         f"Generated by `scripts/40_p3_oracles.py analyze`. Predictor: frozen P* `{q2dev.PSTAR}` "
         f"(params {q2dev.PSTAR_PARAMS}); E_fixed_dev pooled AUPRC; random arm = P1 frozen-P* rows (R = 10 per split). "
         f"Inference: B = {res['B']} stratified query-bootstrap replicates per split (shared by all arms, weights "
         f"'{q2dev.W_NAME}'), Welch SE (DESIGN §6.4), G with the joint bootstrap (weights + trajectories resampled "
         "within arm). Oracles are evaluation-only diagnostics (hsl.oracles; they read pool labels).", ""]
    go = res["go"]
    L += ["## Go/no-go", "", f"Rule (DESIGN §9 P3): {go['rule']}.", ""]
    for cfg, p in res["G_pooled"].items():
        L.append(f"- **{cfg}**: G(160) = **{_f(p['G'], 2)}** (95% CI {_f(p['G_ci'][0], 2)} to {_f(p['G_ci'][1], 2)}; "
                 f"{p['n_splits']} split(s) pooled as Σ(s160 − r160)/Σ(r660 − r160); denominator CI "
                 f"{_f(p['den_ci'][0])} to {_f(p['den_ci'][1])}, G defined: {p['defined']})")
    if go["complete"]:
        L += ["", f"**Verdict: {'C and M DE-SCOPED (both G < 0.5)' if go['descope_CM'] else 'GO — C and M stay in P4 (at least one oracle has G ≥ 0.5)'}.**"]
    else:
        L += ["", "**Verdict pending: not all oracle runs are complete.**"]
    L += ["", "## Ceilings per split", "",
          "| oracle | split | R | s160 | r160 | r660 | Δ(160) [95% CI] | G [95% CI] | s160 − r490 [CI] | s160 − r660 [CI] | "
          "matches p (s160 − 0.95·r660) | G* | traj wall s | select s to 160 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in cel.iterrows():
        L.append(f"| {r.config} | {r.split} | {r.R_s} | {_f(r.get('s160'))} | {_f(r.get('r160'))} | {_f(r.get('r660'))} | "
                 f"{_f(r.get('delta160'))} [{_f(r.get('delta160_lo'))}, {_f(r.get('delta160_hi'))}] | "
                 f"{_f(r.get('G'), 2)} [{_f(r.get('G_lo'), 2)}, {_f(r.get('G_hi'), 2)}] | "
                 f"{_f(r.get('c1_est'))} [{_f(r.get('c1_lo'))}, {_f(r.get('c1_hi'))}] | "
                 f"{_f(r.get('brief_est'))} [{_f(r.get('brief_lo'))}, {_f(r.get('brief_hi'))}] | "
                 f"{_f(r.get('matches_p'), 3)} | {_f(r.get('G_star'), 2)} | {_f(r.get('traj_wall_s_mean'), 0)} | "
                 f"{_f(r.get('select_s_to160'), 0)} |")
    if res.get("strata"):
        L += ["", "## Ceilings by stratum (G at n = 160; joint bootstrap)", "",
              "The hub ceiling restricted to on-array rows / on-array columns separates the yeast off-array SL excess "
              "(D5); the greedy oracles on query halves separate in-sample optimisation from an attainable ceiling: "
              "`oracle_greedy` selected on ALL D_val queries (in-sample on both halves: an upper bound, **not "
              "attainable**), `oracle_greedyHO` selected on half A only, so its **half-B G is the attainable-ceiling "
              "estimate**.", "",
              "| oracle | stratum | split | s160 | r160 | r660 | Δ(160) [95% CI] | G [95% CI] |", "|---|---|---|---|---|---|---|---|"]
        for r in res["strata"]:
            L.append(f"| {r['config']} | {STRATUM_LABEL[r['stratum']]} | {r['split']} | {_f(r.get('s160'))} | "
                     f"{_f(r.get('r160'))} | {_f(r.get('r660'))} | {_f(r.get('delta160'))} "
                     f"[{_f(r.get('delta160_lo'))}, {_f(r.get('delta160_hi'))}] | {_f(r['G'], 2)} "
                     f"[{_f(r['G_lo'], 2)}, {_f(r['G_hi'], 2)}] |")
    L += ["", "## Mean AUPRC by budget (all arms)", "", "| arm | split | " + " | ".join(str(n) for n in Q2) + " |",
          "|---|---|" + "---|" * len(Q2)]
    for (cfg, split), d in tab.groupby(["config", "split"]):
        v = {int(n): a for n, a in zip(d.n, d.auprc)}
        L.append(f"| {cfg} | {split} | " + " | ".join(_f(v.get(n)) for n in Q2) + " |")
    L += ["", "Strata (E_AO = array-only columns, E_VV = on-array D_val columns) and runtime per checkpoint: "
          "`curves.csv`. Figure: `fig_p3_ceilings.png`.", ""]
    if "controls" in res:
        c = res["controls"]
        L += ["## End-to-end controls (DESIGN §11; P*, dev1, random permutation 0, n = 160)", "",
              f"- Shuffled screened labels (ε, P, label jointly permuted over the screened block; array SMF permuted; "
              f"seed `p3_e2e_control_seed`): AUPRC {_f(c['shuffled']['auprc'])} vs prevalence "
              f"{_f(c['shuffled']['prevalence'])} (ratio {_f(c['shuffled']['ratio'], 2)}; pass = within 0.7–1.3×: "
              f"**{c['shuffled']['pass']}**). Unshuffled reference (same S, same model seed): "
              f"{_f(c.get('reference_auprc'))}.",
              f"- Oracle back door (P* scores + 10 × true target label): AUPRC {_f(c['backdoor']['auprc'])} "
              f"(pass = > 0.99: **{c['backdoor']['pass']}**).", ""]
    L += [f"Engine check: unit-weight replicate-engine AP equals the stored harness AUPRC (|diff| ≤ 1e-9) for "
          f"{sum(res['engine_checks'].values())} prediction files.", ""]
    (OUT / "summary.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--config", required=True, choices=["oracle_hub", "oracle_greedy", "oracle_greedyHO"])
    r.add_argument("--split", required=True)
    r.add_argument("--traj", type=int, required=True)
    r.add_argument("--checkpoints", default="")
    r.add_argument("--max-n", type=int, default=None)
    r.add_argument("--threads", default=2)
    c = sub.add_parser("controls")
    c.add_argument("--n", type=int, default=160)
    c.add_argument("--threads", default=2)
    j = sub.add_parser("jobs")
    j.add_argument("--name", required=True)
    j.add_argument("--which", default="all", choices=["hub", "greedy", "controls", "all"])
    j.add_argument("--splits", default="dev1,dev2,dev3")
    j.add_argument("--trajs", default="0")
    an = sub.add_parser("analyze")
    an.add_argument("--ledger", action="store_true")
    a = ap.parse_args(argv)
    return {"run": cmd_run, "controls": cmd_controls, "jobs": cmd_jobs, "analyze": cmd_analyze}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
