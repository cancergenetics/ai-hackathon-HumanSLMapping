#!/usr/bin/env python
"""Universe subsampling — scaling law fraction vs count (DESIGN v3.1 §7.3; DEV DATA ONLY).

Levels 50% and 25% of the non-test genes (3 draws each, seed universe_subsample_seed; hsl.universe), inner folds
rebuilt inside each draw (val = fold 1), R0 random permutations of the sub-pool (dev trajectory seeds r = 0..R−1,
as for the dev splits), P* only, budget grid at the same pool fractions as the full-universe grid (~11–12 budgets incl.
the sub-pool maximum). Level 100% = dev splits 1–3 (results/p1/curves.csv, P*, same fraction grid).
n_80 / n_90: frozen PAVA/n_req on the trajectory-mean curve relative to its value at the (sub-)pool maximum. Slope b of
log n_X on log N (N = candidate-pool size) by OLS over the draws; CI by a draw bootstrap stratified by level.

    run      --level 50 --draw 0 --perms 0-9 [--budgets 15,30] --job NAME
    jobs     --out jobs.json [--levels 50,25] [--perms 0-9]
    analyze  -> results/p1/universe/summary.md, universe_curves.csv, universe_nx.csv, slope.json
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from hsl import config, eval as E, p1, universe as UV, views

OUT = p1.P1_DIR / "universe"
PARTS = OUT / "parts"
PY = ".venv/bin/python"
FULL_BASE = [30, 60, 100, 160, 220, 330, 490, 660, 980, 1310, 1640]


def _ints(s):
    out = []
    for tok in str(s).split(","):
        if "-" in tok:
            a, b = tok.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif tok.strip():
            out.append(int(tok))
    return out


def _setup(level, draw):
    U = UV.draw_universe(level, draw)
    sb = UV.sub_bundle(U)
    sp = UV.sub_split(sb, level, draw)
    es = views.e_fixed(sb, sp)
    return U, sb, sp, es


def cmd_run(a) -> int:
    p1.use_p1_prediction_dir(p1.PRED_P1 / "universe")
    U, sb, sp, es = _setup(a.level, a.draw)
    grid = _ints(a.budgets) if a.budgets else UV.scaled_grid(len(sp.pool))
    part = PARTS / f"{a.job}.jsonl"
    done = {(r["perm"], r["n"], r.get("predictor", p1.REF_PSTAR)) for r in p1.read_rows([part])}
    c = es.counts()
    for r in _ints(a.perms):
        t = p1.random_traj(sp.name, r, pool=sp.pool, arm="R0_random_universe")
        for n in sorted(grid):
            if (r, n, a.method) in done:
                continue
            with p1.PeakRSS() as pk:
                res = UV.run_sub(sb, sp, es, t.perm[:n], a.method, t.model_seed, r, t.traj_id,
                                 params=p1.PREDICTORS.get(a.method, {}))
            m = res["metrics"]
            row = {"level": a.level, "draw": a.draw, "split": sp.name, "perm": r, "n": int(n), "predictor": a.method,
                   "universe_genes": int(len(U)), "n_query_rows": int(len(sb.query_orfs)), "n_arrays": int(sb.n_arrays),
                   "n_val": int(len(sp.targets)), "pool_n": int(len(sp.pool)), "eset_digest": es.digest(),
                   "E_n_cells": c["n_cells"], "E_n_pos": c["n_pos"], "E_AO_pos": c["AO"]["n_pos"],
                   "E_VV_pos": c["VV"]["n_pos"], "model_seed": t.model_seed, "run_id": res["run_id"],
                   "traj_id": t.traj_id, "pred_file": res["pred_file"], "fit_predict_s": res["fit_predict_s"],
                   "peak_rss_fit_gb": pk.peak,
                   **{k: m.get(k) for k in ("prevalence", "auprc", "lift", "macro_auprc", "auroc", "auprc_AO",
                                            "auprc_VV")}}
            p1.append_row(part, row)
            print(f"[univ] {sp.name} r{r} n{n}: AUPRC {m['auprc']:.4f} fit {res['fit_predict_s']:.1f}s "
                  f"peak {pk.peak:.2f} GB", flush=True)
    return 0


def cmd_jobs(a) -> int:
    c = config.protocol()["compute"]
    small, large = [], []
    for lev in _ints(a.levels):
        for d in range(UV.N_DRAWS):
            _, _, sp, _ = _setup(lev, d)
            grid = UV.scaled_grid(len(sp.pool))
            seg = {"small": [n for n in grid if n <= int(c["small_job_max_n"])],
                   "large": [n for n in grid if n > int(c["small_job_max_n"])]}
            for r in _ints(a.perms):
                for s, g in seg.items():
                    if not g:
                        continue
                    name = f"univ_u{lev}d{d}_r{r}_{s}"
                    job = {"name": name, "n": max(g), "cmd": [PY, "scripts/22_universe_subsample.py", "run", "--level",
                                                               str(lev), "--draw", str(d), "--perms", str(r),
                                                               "--budgets", ",".join(map(str, g)), "--job", name]}
                    (small if s == "small" else large).append(job)
    Path(a.out).write_text(json.dumps(small + large, indent=1) + "\n")
    print(f"{len(small) + len(large)} jobs -> {a.out} (small {len(small)}, large {len(large)})")
    return 0


def _nx(curve_df: pd.DataFrame, fracs=(0.8, 0.9)) -> dict:
    m = curve_df.groupby("n").auprc.mean().sort_index()
    grid, curve = m.index.to_numpy(), m.to_numpy()
    fit = E.pava(curve)
    out = {"auprc_max_fit": float(fit[-1]), "R": int(curve_df.perm.nunique()), "grid": "|".join(map(str, grid))}
    for f in fracs:
        nr = E.n_req(f * fit[-1], grid, curve)
        out[f"n{int(f * 100)}"], out[f"n{int(f * 100)}_censored"] = nr.n, nr.censored
    return out


def _ols(x, y):
    X = np.c_[np.ones(len(x)), x]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(beta[1]), float(beta[0])


def cmd_analyze(a) -> int:
    rows = p1.read_rows(sorted(PARTS.glob("*.jsonl")))
    ud = pd.DataFrame(rows).drop_duplicates(["split", "perm", "n"], keep="last")
    ud.sort_values(["level", "draw", "perm", "n"]).to_csv(OUT / "universe_curves.csv", index=False)
    cur = pd.read_csv(p1.P1_DIR / "curves.csv")
    cur = cur[cur.predictor == p1.PSTAR]
    recs = []
    for sp, d in cur.groupby("split"):
        pool_n = len(views.get_split(sp).pool)
        keep = FULL_BASE + [pool_n]
        dd = d[d.n.isin(keep)]
        if dd.n.max() != pool_n:
            continue
        recs.append({"level": 100, "draw": int(sp[-1]), "split": sp, "pool_n": pool_n,
                     "n_val": int(len(views.get_split(sp).targets)), "E_n_pos": int(d.n_pos.iloc[0]),
                     **_nx(dd), **{f"dense_{k}": v for k, v in _nx(d).items() if k.startswith("n")}})
    for (lev, dr), d in ud.groupby(["level", "draw"]):
        pool_n = int(d.pool_n.iloc[0])
        if d.n.max() != pool_n:
            continue
        recs.append({"level": int(lev), "draw": int(dr), "split": d.split.iloc[0], "pool_n": pool_n,
                     "n_val": int(d.n_val.iloc[0]), "E_n_pos": int(d.E_n_pos.iloc[0]), **_nx(d)})
    nx = pd.DataFrame(recs).sort_values(["level", "draw"], ascending=[False, True])
    nx.to_csv(OUT / "universe_nx.csv", index=False)
    out = {"n_draws": int(len(nx)), "levels": sorted(nx.level.unique().tolist())}
    rng = np.random.default_rng(p1.boot_seed("p1_universe_slope"))
    for k in ("n80", "n90"):
        ok = np.isfinite(nx[k]) & (nx[k] > 0)
        x, y = np.log(nx.pool_n[ok].to_numpy(float)), np.log(nx[k][ok].to_numpy(float))
        lv = nx.level[ok].to_numpy()
        b, a0 = _ols(x, y)
        bs = []
        for _ in range(2000):
            idx = np.concatenate([rng.choice(np.flatnonzero(lv == L), size=(lv == L).sum()) for L in np.unique(lv)])
            if len(np.unique(x[idx])) > 1:
                bs.append(_ols(x[idx], y[idx])[0])
        lo, hi = E.percentile_ci(np.array(bs))
        out[k] = {"b": b, "intercept": a0, "ci95_draw_bootstrap": [lo, hi], "n_points": int(ok.sum()),
                  "n_boot": len(bs)}
    # fraction view: n_X / N per level
    out["frac_by_level"] = {int(L): {k: float((d[k] / d.pool_n).mean()) for k in ("n80", "n90")}
                            for L, d in nx.groupby("level")}
    (OUT / "slope.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    L = ["# Universe subsampling — scaling of n_80 / n_90 with the candidate-pool size (DESIGN §7.3, dev only)", "",
         "P* = stack_lr_wr2, R0 random permutations; level 100% = dev splits 1–3 evaluated on the same pool-fraction "
         "grid (30, 60, …, 1640, max); sub-levels: grid scaled to the sub-pool (`hsl.universe.scaled_grid`). n_X "
         "relative to the curve's own (sub-)pool maximum (frozen PAVA/n_req).", "",
         "| level | draw | N (pool) | |D_val| | E_fixed_dev positives | R | AUPRC at max (PAVA) | n_80 | n_90 | "
         "n_80/N | n_90/N |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in nx.iterrows():
        L.append(f"| {r.level}% | {r.draw} | {r.pool_n} | {r.n_val} | {r.E_n_pos} | {r.R} | {r.auprc_max_fit:.4f} | "
                 f"{r.n80:.0f}{' (' + r.n80_censored + ')' if r.n80_censored else ''} | "
                 f"{r.n90:.0f}{' (' + r.n90_censored + ')' if r.n90_censored else ''} | "
                 f"{r.n80 / r.pool_n:.3f} | {r.n90 / r.pool_n:.3f} |")
    L += ["", "Slope of log n_X on log N (OLS over draws; 95% CI = draw bootstrap stratified by level, B = 2,000):", ""]
    for k in ("n80", "n90"):
        s = out[k]
        L.append(f"- **{k}: b = {s['b']:.2f}** [{s['ci95_draw_bootstrap'][0]:.2f}, {s['ci95_draw_bootstrap'][1]:.2f}] "
                 f"({s['n_points']} draws). b ≈ 0 → count scaling, b ≈ 1 → fraction scaling.")
    L += ["", "## What was and was not restricted to the sub-universe", "",
          "- Restricted (recomputed inside U): query rows and array columns (same genes removed from both sides; T "
          "columns dropped), inner folds (rebuilt, val = fold 1), candidate pool, E_fixed_dev, and every statistic the "
          "predictor derives from the view — array rates (incl. LOO), query propensity, kNN similarities, `tied` "
          "low-rank, 2-hop scores, GI-profile kernel, external-kernel normalisation over the training genes (mean "
          "diagonal, PSD projection), λ scaling, LR feature standardisation; array SMF is per gene.",
          "- Not restricted: the upstream pair-feature store used by `mk`/`ext_off` (coexpression Fisher z, AP-MS "
          "edge score and neighbourhood cosine, ESM-2 standardised cosine, paralog identity, InterPro Jaccard). The "
          "upstream builders compute these over each SOURCE's own gene set (all genes of each expression study; all "
          "AP-MS network genes, deliberately 'so degrees are not truncated'; all ESM-embedded proteins), not over the "
          "S1 universe, so rebuilding the store on U with the upstream code would return identical values. Restricting "
          "them would require re-deriving the sources. The sub-universes therefore keep full-genome external context: "
          "**b is reported as a bound** (direction not determined; external features are per pair, so the effect is "
          "expected to be small).",
          "- Caveat (DESIGN §7.3): subsampling mimics fewer genes per module, not more modules. At 25% E_fixed_dev "
          "has ~420–500 positives (noise check above).", ""]
    (OUT / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--level", type=int, required=True)
    r.add_argument("--draw", type=int, required=True)
    r.add_argument("--perms", required=True)
    r.add_argument("--budgets", default="")
    r.add_argument("--job", required=True)
    r.add_argument("--method", default=p1.REF_PSTAR)   # default = pre-P2 P* (already-queued jobs); frozen P*: pass it
    j = sub.add_parser("jobs")
    j.add_argument("--out", required=True)
    j.add_argument("--levels", default="50,25")
    j.add_argument("--perms", default="0-9")
    sub.add_parser("analyze")
    a = ap.parse_args(argv)
    return {"run": cmd_run, "jobs": cmd_jobs, "analyze": cmd_analyze}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
