#!/usr/bin/env python
"""P2 predictor development for n <= 330 (DESIGN v3.1 §4, §8 predictor-change gate, §9 P2; DEV DATA ONLY).

Shared R0 random trajectories of P1 (p1.random_traj: SeedSequence(dev_traj_base_seed).spawn, nested prefixes) with
identical model seeds; parent = stack_lr_wr2 (P*), candidates = stack_lr_wr2x switches (hsl.predictors.lowbudget).
Every candidate and the parent run through the frozen harness (harness.run_predictor; ledger=False / save=False per
run — one ledger line per CONFIGURATION is written by `gate`). Per run: E_fixed_dev point metrics + strata and the
B stratified query-bootstrap replicate APs under the split's shared weights (p1.query_weights(split, 'p2_gate')),
stored as float64 .npy (results/p2/boot/) — no prediction vectors are saved here.

    diag   --perms 0-4 --budgets 30,160,660            evaluation-only oracle decomposition on dev1 (where P* loses)
    run    --split dev1 --perms 0-9 --budgets 30,160,330,660 --configs dense,prop --job NAME
    jobs   --out FILE --configs dense,prop [--splits dev1,dev2,dev3,resplit]
    gate   [--configs ...]                             paired contrasts, gate verdicts, summary.md, ledger lines
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import hsl.predictors  # noqa: F401  (registers stack_lr_wr2x)
from hsl import config, eval as E, harness as H, p1, views
from hsl.predictors import lowbudget as LB

P2 = p1.RESULTS / "p2"
PARTS = P2 / "parts"
BOOT = P2 / "boot"
PY = ".venv/bin/python"
PARENT = "stack_lr_wr2"
WEIGHTS_NAME = "p2_gate"
# PI 2026-10-29: disk no longer constrained -> pass-2 runs save E_fixed_dev predictions (results/predictions/p2/);
# pass-1 runs (before the note) saved bootstrap replicates only
SAVE = True
H.PRED_DIR = p1.RESULTS / "predictions" / "p2"

# Candidate configurations (cap protocol.yaml gate_max_configs = 10; a P2 addition, not a DESIGN §8 rule). name -> (method, params, rationale)
CONFIGS = {
    # parent run as stack_lr_wr2x with every switch off: bit-identical to stack_lr_wr2 (tests/test_p2_predictors.py;
    # dev1 run ids reproduce P1's AUPRC exactly) and it fills the feature cache the candidates share
    "parent": ("stack_lr_wr2x", {}, "P* = stack_lr_wr2, upstream params"),
    "dense": ("stack_lr_wr2x", {"extra": ["dense"]},
              "dense epsilon-based array/query degree statistics (mild-negative rate, mean -eps; LOO for training "
              "rows) next to the label rates, which are very noisy at small n"),
    "prop": ("stack_lr_wr2x", {"extra": ["prop"]},
             "upstream query-propensity sub-model (nested cross-fitted) as LR feature (stack_lr_prop columns)"),
    "gi100": ("stack_lr_wr2x", {"gi_n0": 100.0}, "n-dependent GI kernel weight w_gi = n/(n+100) in mk"),
    "gi30": ("stack_lr_wr2x", {"gi_n0": 30.0}, "n-dependent GI kernel weight w_gi = n/(n+30) in mk"),
    "gi300": ("stack_lr_wr2x", {"gi_n0": 300.0}, "n-dependent GI kernel weight w_gi = n/(n+300) in mk"),
    "dense_prop": ("stack_lr_wr2x", {"extra": ["dense", "prop"]}, "dense + prop"),
    # pass 2 (after the R = 2-3 interim look, results/p2/contrasts_interim_R2.csv): query-level columns overfit at
    # n = 30 (dense -3%, prop -6%) -> array-side-only variants; gi300 the best GI weight so far
    "dense_arr": ("stack_lr_wr2x", {"extra": ["dense_arr"]},
                  "array-side dense degree statistics only (arr_mild, arr_negeps)"),
    "dense_arr_aprior": ("stack_lr_wr2x", {"extra": ["dense_arr", "aprior"]},
                         "dense_arr + empirical-Bayes array-degree prior (Poisson GLM on SMF spline + external "
                         "gene scalars of the array gene)"),
    "dense_arr_gi300": ("stack_lr_wr2x", {"extra": ["dense_arr"], "gi_n0": 300.0}, "dense_arr + gi300"),
    # 10th (last) configuration, added after the R = 1 look (results/p2/contrasts_interim_R1.csv): aprior and gi300
    # each added to dense_arr; their combination
    "dense_arr_aprior_gi300": ("stack_lr_wr2x", {"extra": ["dense_arr", "aprior"], "gi_n0": 300.0},
                               "dense_arr + aprior + gi300"),
}
STOPPED_FUTILITY = {"prop": "negative at n = 30 and 160 at the R = 2 interim look; costly (+70-170 s/fit)",
                    "dense_prop": "below dense at every budget at the R = 2 interim look; needs prop"}


def _ints(s: str) -> list[int]:
    out = []
    for tok in str(s).split(","):
        if "-" in tok:
            a, b = tok.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif tok.strip():
            out.append(int(tok))
    return out


_W: dict = {}


def weights(split: str) -> np.ndarray:
    if split not in _W:
        _W[split] = p1.query_weights(split, WEIGHTS_NAME)
    return _W[split]


def boot_path(split, r, n, cfg) -> Path:
    return BOOT / split / f"{cfg}_r{r}_n{n}.npy"


def replicate_ap(es, s, W) -> tuple[np.ndarray, float]:
    e = E.WeightedAP(es.y, s, es.rows, W.shape[1])
    ch = int(config.protocol()["bootstrap"]["chunk"])
    reps = np.concatenate([e.ap(W[i:i + ch]) for i in range(0, len(W), ch)])
    return reps, float(e.ap(np.ones(W.shape[1]))[0])


def run_one(t, n: int, cfg: str, n_threads: int | None) -> dict:
    method, params, _ = CONFIGS[cfg]
    S = t.perm[:int(n)]
    t0 = time.perf_counter()
    res = H.run_predictor(H.RunSpec(t.split, tuple(S), ("E_fixed",), ("E_fixed",), label=f"p2_{cfg}"), method, params,
                          {"model_seed": t.model_seed, "traj_index": t.r}, n_threads=n_threads, parent=t.traj_id,
                          save=SAVE, ledger=False, tag=f"p2_{cfg}")
    wall = time.perf_counter() - t0
    es, s = res.esets["E_fixed"], res.scores["E_fixed"]
    m = res.metrics[es.name]
    t1 = time.perf_counter()
    reps, est = replicate_ap(es, s, weights(t.split))
    boot_s = time.perf_counter() - t1
    if abs(est - m["auprc"]) > 1e-9:
        raise AssertionError(f"engine AP {est} != harness AP {m['auprc']}")
    p = boot_path(t.split, t.r, n, cfg)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.save(p, reps)
    ex = p1.extra_metrics(t.split, es, s)
    rt = res.manifest["runtime_s"]
    return {"split": t.split, "perm": t.r, "n": int(n), "config": cfg, "method": method, "params": json.dumps(params),
            "model_seed": t.model_seed, "run_id": res.run_id, "config_hash": res.manifest["config_hash"],
            **{k: m[k] for k in ("auprc", "auprc_AO", "auprc_VV", "macro_auprc", "auroc", "p_at_10", "p_at_50",
                                 "n_cells", "n_pos")},
            **ex, "fit_predict_s": rt.get("fit_predict_s"), "wall_s": wall, "boot_s": boot_s,
            "lr_coef": json.dumps((res.manifest.get("method_info") or {}).get("lr_coef")),
            "git_head": res.manifest["git"].get("head"), "git_dirty": res.manifest["git"].get("dirty"),
            "n_threads": res.manifest["n_threads"], "peak_rss_gb": H.peak_rss_gb(),
            "pred_file": (res.manifest.get("predictions") or {}).get("E_fixed", {}).get("file")}


def _done(part: Path) -> set:
    return {(r["split"], r["perm"], r["n"], r["config"]) for r in p1.read_rows([part])}


def cmd_run(a) -> int:
    part = PARTS / f"{a.job}.jsonl"
    done = _done(part)
    cfgs = [c for c in a.configs.split(",") if c]
    if "parent" not in cfgs:
        cfgs = ["parent"] + cfgs
    fails = 0
    for r in _ints(a.perms):
        t = p1.random_traj(a.split, r)
        for n in _ints(a.budgets):
            todo = [c for c in cfgs if (a.split, r, n, c) not in done and not boot_path(a.split, r, n, c).exists()]
            with LB.feature_cache():                      # features shared across configurations of this cell
                for c in todo:
                    try:
                        row = run_one(t, n, c, a.threads)
                        p1.append_row(part, row)
                        print(f"[p2] {a.split} r{r} n{n} {c}: AUPRC {row['auprc']:.4f} "
                              f"fit {row['fit_predict_s']:.1f}s boot {row['boot_s']:.1f}s", flush=True)
                    except Exception as e:  # noqa: BLE001
                        fails += 1
                        print(f"[p2] FAIL {a.split} r{r} n{n} {c}: {type(e).__name__}: {e}", flush=True)
    return min(fails, 255)


def cmd_jobs(a) -> int:
    jobs = []
    for sp in a.splits.split(","):
        budgets = "30,160,330" if sp == "resplit" else "30,160,330,660"
        for r in _ints(a.perms):
            name = f"p2_{a.tag}_{sp}_r{r}"
            cmd = [PY, "scripts/30_p2_predictor.py", "run", "--split", sp, "--perms", str(r), "--budgets", budgets,
                   "--configs", a.configs, "--job", name]
            jobs.append({"name": name, "cmd": cmd, "n": 660})
    Path(a.out).write_text(json.dumps(jobs, indent=1) + "\n")
    print(f"{len(jobs)} jobs -> {a.out}")
    return 0


# ============================================================================ diagnostics (evaluation-only)
def cmd_diag(a) -> int:
    """Where does P* lose at low n? Oracle decomposition on dev1: P* + oracle array degree (label rate of every
    array over ALL candidate-pool rows) and/or + oracle query degree (true row rate of every target and training
    row). Evaluation-only diagnostics: they read labels the predictor cannot have; never a predictor."""
    from hsl import upstream as U
    from ysl.models.registry import get
    from ysl.models.view import MethodContext
    from ysl.runner import set_thread_limits
    split = "dev1"
    sp = views.get_split(split)
    b = views.dev_bundle()
    pr = config.protocol()
    th = int(a.threads or os.environ.get("HSL_THREADS") or 2)
    set_thread_limits(th)
    Yp = b.Y[b.gene_rows(sp.pool)].astype(np.float64)
    Yp[Yp < 0] = np.nan
    g = np.nanmean(Yp)
    arr_or = np.log((np.nansum(Yp, 0) + 1.0) / ((~np.isnan(Yp)).sum(0) + 1.0 / g))

    def row_rate(orfs):
        Y = b.Y[b.gene_rows(orfs)].astype(np.float64)
        Y[Y < 0] = np.nan
        k, m = np.nansum(Y, 1), (~np.isnan(Y)).sum(1)
        q = (k + 0.5) / (m + 1.0)
        return np.log(q / (1 - q))

    def inj(which):
        def fn(v, rows):
            n_t, n_a = v.n_targets, v.n_arrays
            cols = []
            if "arr" in which:
                cols.append(np.broadcast_to(arr_or[None, :], (n_t, n_a)).ravel())
            if "qry" in which:
                cols.append(np.repeat(row_rate(v.target_orfs), n_a))
            return np.stack(cols, 1).astype(np.float32)
        return fn

    es = H.get_eset(b, sp, "E_fixed")
    on = p1.target_on_array(split)[es.rows]
    out = P2 / "diag_oracle.csv"
    rows = []
    variants = {"parent": None, "oracle_arr": inj(("arr",)), "oracle_qry": inj(("qry",)),
                "oracle_both": inj(("arr", "qry"))}
    for r in _ints(a.perms):
        t = p1.random_traj(split, r)
        for n in _ints(a.budgets):
            S = np.sort(np.asarray(t.perm[:n], dtype=str))
            view = U.construct_view(b, S, sp.targets)
            with LB.feature_cache():
                for nm, fn in variants.items():
                    ctx = MethodContext(model_seed=t.model_seed, split_seed=0,
                                        params=dict(get("stack_lr_wr2").defaults),
                                        n_folds=int(min(pr["predictor"]["n_folds"], n)), n_threads=th,
                                        chunk_targets=int(pr["predictor"]["chunk_targets"]))
                    t0 = time.perf_counter()
                    Smat = LB.stack_predict(view, ctx, inject=fn)
                    s = es.scores(Smat)
                    row = {"split": split, "perm": r, "n": n, "variant": nm, "fit_s": time.perf_counter() - t0,
                           "auprc": E.average_precision(es.y, s)}
                    for st, mk in (("AO", es.strata["AO"]), ("VV", es.strata["VV"]), ("on_rows", on),
                                   ("off_rows", ~on)):
                        row[f"auprc_{st}"] = E.average_precision(es.y[mk], s[mk])
                    row["lr_coef"] = json.dumps(ctx.info.get("lr_coef"))
                    rows.append(row)
                    print(f"[diag] r{r} n{n} {nm}: {row['auprc']:.4f}", flush=True)
            pd.DataFrame(rows).to_csv(out, index=False)
    d = pd.DataFrame(rows)
    for nm in ("oracle_arr", "oracle_qry", "oracle_both"):
        H.append_ledger({"id": f"p2_diag_{nm}", "parent": PARENT, "config_hash": H._hash(["p2_diag", nm]),
                         "split": split, "method": f"diag:{nm}", "outcome": "ok", "tag": "p2_diag",
                         "metrics": d[d.variant == nm].groupby("n").auprc.mean().to_dict(),
                         "manifest": {"git": H._git(), "perms": _ints(a.perms), "evaluation_only": True}})
    return 0


# ============================================================================ gate
def load_parts() -> pd.DataFrame:
    rows = p1.read_rows(sorted(PARTS.glob("*.jsonl")))
    d = pd.DataFrame(rows)
    return d.drop_duplicates(["split", "perm", "n", "config"], keep="last") if len(d) else d


def contrast(split: str, cfg: str, n: int, perms) -> dict | None:
    bc, bp, ec, ep = [], [], [], []
    for r in perms:
        pc, pp = boot_path(split, r, n, cfg), boot_path(split, r, n, "parent")
        if not (pc.exists() and pp.exists()):
            return None
        c, p = np.load(pc), np.load(pp)
        bc.append(c)
        bp.append(p)
    d = load_parts.cache
    sub = d[(d.split == split) & (d.n == n)]
    for r in perms:
        ec.append(float(sub[(sub.perm == r) & (sub.config == cfg)].auprc.iloc[0]))
        ep.append(float(sub[(sub.perm == r) & (sub.config == "parent")].auprc.iloc[0]))
    return E.paired(np.stack(bc, 1), np.array(ec), np.stack(bp, 1), np.array(ep))


def cmd_gate(a) -> int:
    d = load_parts()
    load_parts.cache = d
    cfgs = [c for c in (a.configs.split(",") if a.configs else d.config.unique()) if c != "parent"]
    perms = _ints(a.perms)
    gb = [int(x) for x in config.protocol()["gate"]["gate_budgets"]]
    min_gain = float(config.protocol()["gate"]["gate_min_rel_gain"])
    max_loss_660 = float(config.protocol()["gate"]["gate_max_rel_loss_660"])
    out, verdicts = [], {}
    for cfg in cfgs:
        for sp in ("dev1", "dev2", "dev3", "resplit"):
            for n in (30, 160, 330, 660):
                if sp == "resplit" and n > 330:
                    continue
                c = contrast(sp, cfg, n, perms)
                if c is None:
                    continue
                sub = d[(d.split == sp) & (d.n == n) & (d.config == cfg)]
                par = d[(d.split == sp) & (d.n == n) & (d.config == "parent")]
                out.append({"config": cfg, "split": sp, "n": n, "R": c["R"], "parent_auprc": float(par.auprc.mean()),
                            "cand_auprc": float(sub.auprc.mean()), "delta": c["delta"], "rel_gain": c["rel_gain"],
                            "ci_lo": c["ci95"][0], "ci_hi": c["ci95"][1], "se": c["se"], "V_Q": c["V_Q"],
                            "traj_var": c["traj"], "p_one_sided": c["p_one_sided"],
                            "rel_AO": float(sub.auprc_AO.mean() / par.auprc_AO.mean() - 1),
                            "rel_VV": float(sub.auprc_VV.mean() / par.auprc_VV.mean() - 1),
                            "rel_on_rows": float(sub.auprc_on_array_rows.mean() / par.auprc_on_array_rows.mean() - 1),
                            "rel_off_rows": float(sub.auprc_off_array_rows.mean() /
                                                  par.auprc_off_array_rows.mean() - 1),
                            "fit_s_cand": float(sub.fit_predict_s.mean()), "fit_s_parent": float(par.fit_predict_s.mean())})
    res = pd.DataFrame(out)
    if not len(res):
        print("no complete contrasts")
        return 1
    sfx = f"_interim_R{len(perms)}" if a.interim else ""
    res.to_csv(P2 / f"contrasts{sfx}.csv", index=False)
    for cfg, g in res.groupby("config"):
        dev = g[g.split.isin(["dev1", "dev2", "dev3"]) & g.n.isin(gb)]
        complete = len(dev) == 9
        mean_gain = float(dev.rel_gain.mean()) if len(dev) else float("nan")
        n160 = dev[dev.n == 160]
        ci_ok = int((n160.ci_lo > 0).sum()) if mean_gain >= 0 else int((n160.ci_hi < 0).sum())
        rs = g[g.split == "resplit"]
        rs_sign = float(rs.delta.mean()) if len(rs) else float("nan")
        g660 = g[(g.n == 660) & g.split.isin(["dev1", "dev2", "dev3"])]
        rel660 = float(g660.rel_gain.mean()) if len(g660) else float("nan")
        passed = bool(complete and mean_gain >= min_gain and ci_ok >= 2 and (n160.ci_lo > 0).sum() >= 2 and
                      len(rs) and np.sign(rs_sign) == np.sign(mean_gain) and rel660 >= -max_loss_660)
        verdicts[cfg] = {"complete": complete, "mean_rel_gain_gate": mean_gain, "n160_ci_excl0_pos": int(
            (n160.ci_lo > 0).sum()), "resplit_mean_delta": rs_sign, "resplit_rel_gain": float(rs.rel_gain.mean())
            if len(rs) else float("nan"), "rel_gain_660": rel660, "pass": passed,
                         "fit_s_cand_330": float(g[g.n == 330].fit_s_cand.mean()) if (g.n == 330).any() else None,
                         "fit_s_parent_330": float(g[g.n == 330].fit_s_parent.mean()) if (g.n == 330).any() else None}
    pd.DataFrame(verdicts).T.to_csv(P2 / f"verdicts{sfx}.csv")
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(json.dumps(verdicts, indent=1, default=float))
    if a.ledger and not a.interim:
        for cfg, why in STOPPED_FUTILITY.items():
            method, params, rat = CONFIGS[cfg]
            H.append_ledger({"id": f"p2_cfg_{cfg}", "parent": PARENT, "config_hash": H._hash([method, params]),
                             "split": "dev1,dev2,dev3,resplit", "method": method, "params": params,
                             "outcome": "fail_futility", "tag": "p2_config",
                             "metrics": {"interim": "results/p2/contrasts_interim_R2.csv", "reason": why},
                             "manifest": {"git": H._git(), "rationale": rat, "weights": WEIGHTS_NAME}})
        for cfg, v in verdicts.items():
            method, params, why = CONFIGS[cfg]
            g = res[res.config == cfg]
            H.append_ledger({"id": f"p2_cfg_{cfg}", "parent": PARENT, "config_hash": H._hash([method, params]),
                             "split": "dev1,dev2,dev3,resplit", "method": method, "params": params,
                             "outcome": "pass" if v["pass"] else "fail", "tag": "p2_config",
                             "metrics": {"verdict": v, "contrasts": g.drop(columns=["config"]).to_dict("records")},
                             "manifest": {"git": H._git(), "perms": perms, "rationale": why,
                                          "weights": WEIGHTS_NAME}})
    return 0


load_parts.cache = None


def cmd_summary(a) -> int:
    """results/p2/summary.md from contrasts.csv, verdicts.csv, diag_oracle.csv (+ the futility-stopped interim)."""
    res = pd.read_csv(P2 / "contrasts.csv")
    ver = pd.read_csv(P2 / "verdicts.csv", index_col=0)
    L = ["# P2 predictor development (n <= 330) — summary", "",
         "Dev splits 1-3 (R = 10 shared R0 random trajectories each, P1 permutations and model seeds) + re-split check "
         "(R = 10, n <= 330). Parent P* = stack_lr_wr2 (upstream params). Paired contrasts on E_fixed_dev: "
         "Var(d_r)/R + V_Q(d-bar), stratified query bootstrap B = 2000 shared by parent and candidate "
         f"(weights '{WEIGHTS_NAME}'). Gate (DESIGN §8): mean relative gain over n in {{30,160,330}} x 3 splits >= 2%; "
         "95% CI at n = 160 excludes 0 on >= 2/3 splits; same sign on the re-split. P2 additions, stricter than DESIGN: "
         "<= 1% relative loss at n = 660 and a cap of 10 configurations (protocol.yaml gate); the frozen P* is the "
         "passing configuration with the largest mean gain (a P2 choice, not a DESIGN rule). "
         "Generated by `scripts/30_p2_predictor.py summary`.", "",
         "**Interpretation (independent review, docs/reviews/p2_predictor_review.md).** The gain is best read as "
         "+4.0-4.5%, the best of 10 adaptively designed dev configurations without independent confirmation. "
         "Configurations 7-10 were designed after interim looks that included the re-split and n = 660, so the re-split "
         "is not a blind check. The gi300 component alone fails the gate; dense_arr_aprior alone gives +4.00%. "
         "The gain is certainly positive (about 8 SE at n = 160), but its size is optimistic by perhaps 0.3-0.5 points. "
         "It is about +0.001-0.003 AUPRC in absolute terms and has been validated on random trajectories only. The "
         "final evaluation on T is unaffected by this selection.", ""]
    if (P2 / "diag_oracle.csv").exists():
        dg = pd.read_csv(P2 / "diag_oracle.csv")
        t = dg.pivot_table(index="variant", columns="n", values="auprc")
        rel = t / t.loc["parent"] - 1
        L += ["## 1. Where P* loses at low n (evaluation-only oracle decomposition, dev1, perms 0-4)", "",
              "P* plus an oracle column (LR refit): array degree = label rate of every array over the whole candidate "
              "pool; query degree = true row rate of every target/training query. Never predictors.", "",
              "| variant | " + " | ".join(f"n={n}" for n in t.columns) + " |", "|---" * (len(t.columns) + 1) + "|"]
        for v in ("parent", "oracle_arr", "oracle_qry", "oracle_both"):
            L.append(f"| {v} | " + " | ".join(f"{t.loc[v, n]:.4f} ({rel.loc[v, n]:+.0%})" for n in t.columns) + " |")
        st = dg[dg.variant == "parent"].groupby("n")[["auprc_AO", "auprc_VV", "auprc_on_rows", "auprc_off_rows"]].mean()
        L += ["", "P* strata (same runs): " + "; ".join(
            f"n={n}: AO {r.auprc_AO:.4f}, VV {r.auprc_VV:.4f}, on-array rows {r.auprc_on_rows:.4f}, off-array rows "
            f"{r.auprc_off_rows:.4f}" for n, r in st.iterrows()), ""]
    L += ["## 2. Gate verdicts", "",
          "| configuration | mean rel. gain (30/160/330 x 3 splits) | n=160 CI > 0 (splits) | re-split rel. gain | "
          "rel. gain n=660 | fit s @330 (cand / parent) | verdict |", "|---|---|---|---|---|---|---|"]
    for cfg, v in ver.iterrows():
        L.append(f"| {cfg} | {v.mean_rel_gain_gate:+.2%} | {int(v.n160_ci_excl0_pos)}/3 | {v.resplit_rel_gain:+.2%} | "
                 f"{v.rel_gain_660:+.2%} | {v.fit_s_cand_330:.1f} / {v.fit_s_parent_330:.1f} | "
                 f"{'PASS' if str(v['pass']) == 'True' else 'fail'} |")
    for cfg, why in STOPPED_FUTILITY.items():
        L.append(f"| {cfg} | stopped for futility at R = 2-3 ({why}) | | | | | fail |")
    L += ["", "Fit s = mean harness fit_predict wall-clock with the per-cell feature cache: the parent pays the "
          "shared feature cost, a candidate only its own extra features + LR refit. Stand-alone fit time of a "
          "candidate ≈ parent + candidate.", "",
          "## 3. Paired contrasts (Δ AUPRC candidate − parent, 95% CI, relative gain)", ""]
    for cfg, g in res.groupby("config", sort=False):
        L += [f"### {cfg}", "", "| split | n=30 | n=160 | n=330 | n=660 |", "|---|---|---|---|---|"]
        for sp in ("dev1", "dev2", "dev3", "resplit"):
            cells = []
            for n in (30, 160, 330, 660):
                r = g[(g.split == sp) & (g.n == n)]
                cells.append("—" if not len(r) else
                             f"{r.delta.iloc[0]:+.4f} [{r.ci_lo.iloc[0]:+.4f}, {r.ci_hi.iloc[0]:+.4f}] "
                             f"({r.rel_gain.iloc[0]:+.1%})")
            L.append(f"| {sp} | " + " | ".join(cells) + " |")
        L.append("")
    L += ["Strata (relative change in mean AUPRC, candidate/parent − 1) and the full table: results/p2/contrasts.csv. "
          "Per-run rows: results/p2/runs.csv. Interim look: results/p2/contrasts_interim_R2.csv.", ""]
    (P2 / "summary.md").write_text("\n".join(L) + "\n")
    d = load_parts()
    d.to_csv(P2 / "runs.csv", index=False)
    print(f"summary.md: {len(L)} lines; runs.csv: {len(d)} rows")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--split", required=True)
    r.add_argument("--perms", required=True)
    r.add_argument("--budgets", default="30,160,330,660")
    r.add_argument("--configs", required=True)
    r.add_argument("--job", required=True)
    r.add_argument("--threads", type=int, default=None)
    j = sub.add_parser("jobs")
    j.add_argument("--out", required=True)
    j.add_argument("--configs", required=True)
    j.add_argument("--tag", required=True)
    j.add_argument("--splits", default="dev1,dev2,dev3,resplit")
    j.add_argument("--perms", default="0-9")
    dg = sub.add_parser("diag")
    dg.add_argument("--perms", default="0-4")
    dg.add_argument("--budgets", default="30,160,660")
    dg.add_argument("--threads", type=int, default=None)
    g = sub.add_parser("gate")
    g.add_argument("--configs", default="")
    g.add_argument("--perms", default="0-9")
    g.add_argument("--ledger", action="store_true")
    g.add_argument("--interim", action="store_true", help="futility look on a subset of perms (no ledger)")
    sub.add_parser("summary")
    a = ap.parse_args(argv)
    return {"run": cmd_run, "jobs": cmd_jobs, "diag": cmd_diag, "gate": cmd_gate, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
