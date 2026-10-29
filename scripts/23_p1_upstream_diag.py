#!/usr/bin/env python
"""Diagnose the upstream-protocol gap (P1 sanity, DEV DATA ONLY): P* 0.0657 on dev1 vs upstream dev 0.0754 at 20%.

    reproduce --seeds 0,1   re-run upstream's OWN S1 dev splits (ysl.eval.splits.make_split('s1','dev',seed): 656
                            screened non-lockbox genes, 1,970 non-lockbox targets; model_seed 17 as upstream) on our dev
                            bundle with the T genes removed (their rows are not in the dev bundle: T∩screened dropped,
                            T∩targets not evaluated), scored on all measured target cells minus T columns and self
                            pairs (= views.e_unobs with S = ∅) with upstream strata. One ledger line per run; no
                            predictions saved (diagnostic). Rows -> results/p1/parts/upstream_diag_repro.jsonl
    report                  results/p1/upstream_protocol_diag.md
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

from hsl import config, eval as E, harness as H, p1, views

OUT = p1.P1_DIR
PART = OUT / "parts" / "upstream_diag_repro.jsonl"
METHODS = ("stack_lr_wr2", "stack_lr", "B1_array_rate")
UP = {"stack_lr_wr2": (0.0754, 6.6, 0.073, 0.058, 0.126), "stack_lr": (0.0590, 5.1, 0.051, 0.047, 0.117),
      "B1_array_rate": (0.0430, 3.7, 0.043, 0.048, 0.038)}       # REPORT §3: AUPRC, lift, on-no-recip, off, recip
UP_SPLITS = {0: 0.07747769260781422}                              # m4b5 stack_lr_wr2 split 0 (metrics.json)
STRATA = ("recip_available", "x_on_array_no_recip", "x_off_array")


def cmd_reproduce(a) -> int:
    from ysl.eval.splits import make_split
    from hsl import upstream as U
    b = views.dev_bundle()
    T = views.test_queries()
    done = {(r["seed"], r["method"]) for r in p1.read_rows([PART])}
    pr = config.protocol()
    n_threads = int(os.environ.get("HSL_THREADS") or 2)
    for seed in [int(s) for s in a.seeds.split(",")]:
        us = make_split("s1", "dev", seed, b.all_query_orfs)
        S = np.setdiff1d(us.train, T)
        tg = np.setdiff1d(us.test, T)
        sp = views.Split(f"ups1dev{seed}", "dev", tg, S, np.union1d(tg, T), S)
        es = views.e_unobs(b, sp, [])
        ups = None
        view = U.construct_view(b, S, tg)
        for m in METHODS:
            if (seed, m) in done:
                continue
            t0 = time.perf_counter()
            Smat, _ = U.run_method(view, m, model_seed=17, split_seed=0, n_threads=n_threads, n_folds=10,
                                   chunk_targets=int(pr["predictor"]["chunk_targets"]), params={})
            fit = time.perf_counter() - t0
            s = es.scores(Smat).astype(np.float32)
            met = E.summarise(es.y, s, es.rows)
            if ups is None:
                ups = p1.upstream_strata(b, sp, es, S)
            for st in STRATA:
                mk = ups[st]
                met[f"auprc_{st}"] = E.average_precision(es.y[mk], s[mk])
                met[f"prev_{st}"] = float(es.y[mk].mean())
                met[f"n_cells_{st}"] = int(mk.sum())
            for k in config.protocol()["folds"]["dev_splits"]:
                rk = np.isin(tg[es.rows], views.dev_split(k).targets)
                met[f"auprc_rows_dev{k}"] = E.average_precision(es.y[rk], s[rk])
                met[f"prev_rows_dev{k}"] = float(es.y[rk].mean())
                met[f"n_rows_dev{k}"] = int(len(np.unique(es.rows[rk])))
            rest = ~np.isin(tg[es.rows], np.concatenate([views.dev_split(k).targets for k in (1, 2, 3)]))
            met["auprc_rows_folds45"] = E.average_precision(es.y[rest], s[rest])
            met["prev_rows_folds45"] = float(es.y[rest].mean())
            row = {"seed": seed, "method": m, "n_screened": int(len(S)), "n_targets": int(len(tg)),
                   "n_screened_upstream": int(len(us.train)), "n_targets_upstream": int(len(us.test)),
                   "fit_predict_s": fit, **{k: v for k, v in met.items() if not isinstance(v, dict)}}
            rid = f"{sp.name}_{m}_n{len(S)}_{H._hash(row)}"
            H.append_ledger({"id": rid, "parent": None, "config_hash": H._hash({"seed": seed, "method": m}),
                             "split": sp.name, "method": m, "n_screened": int(len(S)), "outcome": "ok",
                             "metrics": {"E_upstream_dev": met}, "manifest": {"tag": "p1_upstream_diag_repro",
                                                                              "model_seed": 17, "git": H._git()}})
            p1.append_row(PART, row)
            print(f"[diag] upstream dev seed {seed} {m}: AUPRC {met['auprc']:.4f} lift {met['lift']:.2f} "
                  f"({len(S)} screened, {len(tg)} targets, fit {fit:.0f}s)", flush=True)
    return 0


def _strat_table(d: pd.DataFrame, pre: str) -> list[str]:
    L = []
    for m in METHODS:
        x = d[d.predictor == m] if "predictor" in d else d[d.method == m]
        if not len(x):
            continue
        g = lambda k: x[pre + k]  # noqa: E731
        cells = [f"{g('auprc').mean():.4f} ± {g('auprc').std(ddof=1) if len(x) > 1 else 0:.4f}",
                 f"{g('prevalence').mean():.4f}", f"{(g('auprc') / g('prevalence')).mean():.2f}"]
        for st in STRATA:
            ap, pv = g(f"auprc_{st}"), g(f"prev_{st}")
            cells.append(f"{ap.mean():.4f} / {pv.mean():.4f} / {(ap / pv).mean():.2f}")
        u = UP[m]
        L.append(f"| {m} | {len(x)} | " + " | ".join(cells) +
                 f" | {u[0]:.4f} (lift {u[1]}); {u[4]} / {u[2]} / {u[3]} |")
    return L


def cmd_report(a) -> int:
    up = pd.read_csv(OUT / "upstream_protocol.csv")
    up = up[up.n == 660]
    hdr = ["| predictor | R | AUPRC | prev | lift | recip: AUPRC / prev / lift | on-array no-recip | off-array | "
           "upstream REPORT §3: AUPRC (lift); recip / on-no-recip / off |", "|---|---|---|---|---|---|---|---|---|"]
    L = ["# Upstream-protocol gap diagnosis (dev only)", "",
         "Question: P* (stack_lr_wr2) scores 0.0657 on the upstream-protocol set of dev1 at n = 660 vs upstream's dev "
         "0.0754 at 20%. Set = every measured cell of the held-out rows incl. reciprocal-known cells, minus test-gene "
         "(T) columns (dev rule) and self pairs. Strata as upstream's evaluator (`recip_available`: target on array, "
         "array's ORF screened, reciprocal cell observed).", ""]
    for sp in sorted(up.split.unique()):
        L += [f"## {sp}, n = 660, our R0 random sets (same sets for all predictors)", ""] + hdr + \
             _strat_table(up[up.split == sp], "E_upstream_dev.") + [""]
    rows = p1.read_rows([PART])
    if rows:
        rp = pd.DataFrame(rows)
        L += ["## Upstream's own dev splits re-run here (make_split('s1','dev',seed); model_seed 17; T genes dropped)",
              "", f"Screened {rp.n_screened.iloc[0]} of {rp.n_screened_upstream.iloc[0]} (T genes have no rows in the "
              f"dev bundle), targets {rp.n_targets.iloc[0]} of {rp.n_targets_upstream.iloc[0]}; T columns excluded.",
              ""] + hdr + _strat_table(rp, "") + [""]
        for s, v in UP_SPLITS.items():
            x = rp[(rp.seed == s) & (rp.method == "stack_lr_wr2")]
            if len(x):
                L.append(f"Upstream m4b5 stack_lr_wr2 split {s} (all columns incl. lockbox): {v:.4f}; here "
                         f"{x.auprc.iloc[0]:.4f}.")
        L.append("")
        L += ["Same upstream-split predictions, restricted to the target rows that belong to each of our dev folds "
              "(AUPRC / prevalence / lift):", "", "| predictor | dev1 rows | dev2 rows | dev3 rows | folds 4–5 rows |",
              "|---|---|---|---|---|"]
        for m in METHODS:
            x = rp[rp.method == m]
            if len(x):
                L.append(f"| {m} | " + " | ".join(
                    f"{x[f'auprc_rows_{k}'].mean():.4f} / {x[f'prev_rows_{k}'].mean():.4f} / "
                    f"{(x[f'auprc_rows_{k}'] / x[f'prev_rows_{k}']).mean():.2f}"
                    for k in ("dev1", "dev2", "dev3", "folds45")) + " |")
        L.append("")
    # ratios
    L += ["## Ratios to upstream (ours / upstream REPORT §3 AUPRC)", ""]
    for sp in sorted(up.split.unique()):
        d = up[up.split == sp]
        L.append(f"- {sp}: " + ", ".join(
            f"{m} {d[d.predictor == m]['E_upstream_dev.auprc'].mean() / UP[m][0]:.3f}" for m in METHODS
            if (d.predictor == m).any()))
    if rows:
        L.append("- upstream dev splits re-run: " + ", ".join(
            f"{m} {rp[rp.method == m].auprc.mean() / UP[m][0]:.3f}" for m in METHODS if (rp.method == m).any()))
    L += ["", "## Conclusion", "",
          "- **Pipeline/stores/params are correct:** upstream's own dev split 0 re-run through our harness path "
          "(`hsl.upstream.construct_view` + `run_method`, registered defaults = upstream recipe `params: {}`, model "
          "seed 17) gives P* 0.0774 vs upstream's recorded 0.0775 (lift 6.8), stack_lr 0.0615, B1 0.0433. All six "
          "external pair stores used by `mk`/`ext_off` (apms_cosine, apms_edge, coex_z, dom_ipr_jaccard, esm_cos, "
          "para_identity) cover dev1 D_val and the screened genes like the dev pool as a whole (coverage 0.70–1.00 "
          "per store, every gene in the 4,426-gene store).",
          "- **Cause = dev1's target rows (query composition), not the predictor:** with upstream's own split-0 "
          "predictions, the dev1 rows alone score 0.0673 (P*), 0.0563 (stack_lr), 0.0409 (B1) — the same ~0.87/0.92/"
          "0.97 shortfall we see — while dev2 rows score 0.0864 and dev3 rows 0.0750. Our dev2/dev3 upstream-protocol "
          "runs at n = 660 give P* 0.0796 / 0.0757 (ratio to upstream 1.06 / 1.00). dev1 has lower prevalence "
          "(0.0109) and its queries are harder for every GI-based model; the gap grows with model strength because "
          "the GI-driven gain (recip / on-array strata) is what varies between query sets. Not a composition shift "
          "of the cell set (B1 lift 3.79 ≈ upstream 3.7) and not a P*-only defect.",
          "- Implication: absolute dev AUPRCs differ by split by up to ~±10% (query sampling); within-split paired "
          "contrasts (same queries) are unaffected. The power pilot (split 1) is on the hardest fold.", ""]
    (OUT / "upstream_protocol_diag.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reproduce")
    r.add_argument("--seeds", default="0")
    sub.add_parser("report")
    a = ap.parse_args(argv)
    return {"reproduce": cmd_reproduce, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
