"""Shared loaders for the Q3 analyses (DESIGN v3.1 §7.2-7.4). Nothing here is frozen; hsl.* is imported read-only.

Sources (every analysis takes --source dev|final|dryrun):
  dev     P1 random curves (results/p1/curves.csv, P* rows), P1 n_80/n_90 (results/p1/n80_n90.csv), P4 strategy rows
          (results/p4/parts/*.jsonl) and stored selections (results/p1/trajectories, results/trajectories/dev*).
  final   scripts/score_final.py outputs in results/final/ (q1_curves.csv, q1_nx.csv, runs.csv, discovery.csv,
          selections.json.gz, run/labels.npz). Read ONLY after scoring has finished; never during the live run.
  dryrun  the same files from results/final_dryrun/ (dev1 stands in for T): exercises the final adapters on dev data.

Unified tables (one row per trajectory x budget), columns:
  arm, split, traj, n, auprc, auprc_VV (on-array target columns; TT at the final), sl_detected, sl_detected_frac,
  p_at_10, p_at_50, pred_file, model_seed
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
Q3 = Path(__file__).resolve().parent
SOURCES = ("dev", "final", "dryrun")
PSTAR = "stack_lr_wr2_lb"
# arm name used throughout for the random arm (P1 'R0_random', final 'R0')
RANDOM = "R0"
# arms analysed in Q3 (dev P4 names and final cids); scenario arms (det, ws220, feas, min, M) are left to P4/score_final
MAIN_ARMS = ("R0", "H", "C", "F_r16", "H-static", "F-static_r8", "F-static_r16", "R1", "H__b110", "C__b110",
             "H__b50", "C__b50")
FINALISTS = ("H", "C", "F_r16")


def set_threads(n: int = 2) -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ[k] = str(n)
    os.environ.setdefault("HSL_THREADS", str(n))


_OUT: dict = {}


def out_dir(source: str, override: str | None = None) -> Path:
    if override:
        _OUT[source] = Path(override)
    if source in _OUT:
        p = _OUT[source]
    else:
        p = ROOT / "results" / {"dev": "q3_dev", "final": "q3_final", "dryrun": "q3_dryrun"}[source]
    p.mkdir(parents=True, exist_ok=True)
    return p


def final_dir(source: str) -> Path:
    if source == "final":
        return ROOT / "results" / "final"
    if source == "dryrun":
        return ROOT / "results" / "final_dryrun"
    raise ValueError(source)


@lru_cache(maxsize=1)
def hap1_constants() -> dict:
    return yaml.safe_load((Q3 / "hap1_constants.yaml").read_text())


@lru_cache(maxsize=1)
def protocol() -> dict:
    return yaml.safe_load((ROOT / "config" / "protocol.yaml").read_text())


def finalists() -> list[str]:
    return list(yaml.safe_load((ROOT / "config" / "finalists.yaml").read_text())["finalists"])


# ============================================================================ rows
def _p4_rows() -> pd.DataFrame:
    rows = []
    for f in sorted((ROOT / "results" / "p4" / "parts").glob("*.jsonl")):
        rows += [json.loads(ln) for ln in f.read_text().splitlines() if ln.strip()]
    d = pd.DataFrame(rows)
    d = d[d.predictor == PSTAR].rename(columns={"config": "arm_cid", "sl_detected_frac_dev": "sl_detected_frac"})
    d["arm"] = d["arm_cid"]
    return d


def _p1_random_rows() -> pd.DataFrame:
    d = pd.read_csv(ROOT / "results" / "p1" / "curves.csv")
    d = d[(d.predictor == PSTAR) & (d.kind == "curve")].copy()
    d = d.rename(columns={"perm": "traj", "sl_detected_frac_dev": "sl_detected_frac"})
    d["arm"] = RANDOM
    return d


COLS = ["arm", "split", "traj", "n", "auprc", "auprc_OA", "auprc_onrows", "auprc_AO", "sl_detected",
        "sl_detected_frac", "p_at_10", "p_at_50", "pred_file", "model_seed", "prevalence"]
# auprc_OA = on-array target-column stratum (dev E_VV, final E_TT); auprc_onrows = on-array target rows


def rows(source: str) -> pd.DataFrame:
    """All P* rows: random + strategy arms, one per (arm, split, traj, n)."""
    if source == "dev":
        d = pd.concat([_p1_random_rows(), _p4_rows()], ignore_index=True)
        d = d.rename(columns={"auprc_VV": "auprc_OA", "auprc_on_array_rows": "auprc_onrows"})
        d = d[[c for c in COLS if c in d.columns]]
    else:
        d = _final_rows(source)
    d["traj"] = d["traj"].astype(int)
    d["n"] = d["n"].astype(int)
    return d.sort_values(["arm", "split", "traj", "n"]).reset_index(drop=True)


def _final_rows(source: str) -> pd.DataFrame:
    """runs.csv of score_final.py: one row per unit (cid, traj, n, pred). AUPRC = auprc_inline (harness value on
    E_fixed; identical to the engine AP at unit weights)."""
    F = final_dir(source)
    r = pd.read_csv(F / "runs.csv")
    r = r[r["pred"] == "pstar"].copy()
    lab = json.loads((F / "run" / "label_manifest.json").read_text()) if (F / "run" / "label_manifest.json").exists() else {}
    split = lab.get("split", "final" if source == "final" else "dev1")
    out = pd.DataFrame({"arm": r["cid"], "split": split, "traj": r["traj"], "n": r["n"], "auprc": r["auprc_inline"],
                        "sl_detected": r.get("sl_detected"), "sl_detected_frac": r.get("sl_detected_frac"),
                        "pred_file": r["pred_file"], "model_seed": r["model_seed"],
                        "prevalence": r["n_pos"] / r["n_cells"]})
    # strata are not in runs.csv: unit-weight AP on the stratum recomputed from the saved predictions (cached)
    return _final_strata(source, out)


def _final_strata(source: str, d: pd.DataFrame) -> pd.DataFrame:
    from hsl import eval as E
    F = final_dir(source)
    cache = out_dir(source) / "strata_cache.csv"
    have = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=["pred_file", "auprc_OA", "auprc_onrows"])
    have = have.drop_duplicates("pred_file")
    todo = sorted(set(d.pred_file) - set(have.pred_file))
    if todo:
        with np.load(F / "run" / "labels.npz") as z:
            y, rows, TT, onr = z["y"].astype(np.int8), z["rows"], z["TT"], z["on_array_rows"]
        onm = onr[rows]
        new = []
        for pf in todo:
            s = load_pred(pf)
            new.append({"pred_file": pf, "auprc_OA": E.average_precision(y[TT], s[TT]),
                        "auprc_onrows": E.average_precision(y[onm], s[onm])})
        have = pd.concat([have, pd.DataFrame(new)], ignore_index=True)
        have.to_csv(cache, index=False)
    return d.merge(have, on="pred_file", how="left")


# ============================================================================ n_80 / n_90 (yeast)
def nx_table(source: str) -> pd.DataFrame:
    """Per split (dev) or one row (final): n80, n80_lo, n80_hi, n90, ..., pool_n (curve anchor)."""
    if source == "dev":
        d = pd.read_csv(ROOT / "results" / "p1" / "n80_n90.csv")
        d = d[d.predictor == PSTAR].copy()
        out = []
        for _, x in d.iterrows():
            lo80, hi80 = json.loads(x["n80_ci95"])
            lo90, hi90 = json.loads(x["n90_ci95"])
            out.append({"split": x["split"], "pool_n": int(x["pool_n"]), "n80": x["n80"], "n80_lo": lo80,
                        "n80_hi": hi80, "n90": x["n90"], "n90_lo": lo90, "n90_hi": hi90,
                        "auprc_anchor": x["auprc_pool_max_fit"]})
        return pd.DataFrame(out)
    F = final_dir(source)
    d = pd.read_csv(F / "q1_nx.csv")
    d = d[d.predictor == "pstar"].iloc[0]
    return pd.DataFrame([{"split": "final" if source == "final" else "dev1", "pool_n": int(d["anchor_n"]),
                          "n80": d["n80"], "n80_lo": d["n80_lo"], "n80_hi": d["n80_hi"], "n90": d["n90"],
                          "n90_lo": d["n90_lo"], "n90_hi": d["n90_hi"], "auprc_anchor": d["auprc_anchor_pava"],
                          "anchor_is_pool_max": bool(d["anchor_is_pool_max"])}])


def universe_draws() -> pd.DataFrame:
    """The 9 universe-subsampling draws (dev only, DESIGN §7.3): level, pool_n, n80, n90."""
    return pd.read_csv(ROOT / "results" / "p1" / "universe" / "universe_nx.csv")


# ============================================================================ selections
def selections(source: str, arms=None) -> dict:
    """{(arm, split, traj): {"selected": [...ordered...], "model_seed": int, "batch_size": int}}."""
    out = {}
    if source == "dev":
        for f in sorted((ROOT / "results" / "p1" / "trajectories").glob("dev*_R0_random_r*.json")):
            j = json.loads(f.read_text())
            out[(RANDOM, j["split"], int(j["traj_index"]))] = {"selected": j["selected"], "model_seed": j["model_seed"],
                                                              "batch_size": int(j.get("batch_size", 10))}
        for f in sorted((ROOT / "results" / "trajectories").glob("dev*/dev*_r*.json")):
            if f.name.endswith(".selector.json"):
                continue
            j = json.loads(f.read_text())
            arm = j["arm"]
            if arm.startswith(("smoke", "oracle")) or arm == "R0_random_smoke":
                continue
            out[(arm, j["split"], int(j["traj_index"]))] = {"selected": j["selected"], "model_seed": j["model_seed"],
                                                           "batch_size": int(j.get("batch_size") or 10)}
    else:
        F = final_dir(source)
        sel = json.loads(gzip.open(F / "selections.json.gz", "rt").read())
        split = "final" if source == "final" else "dev1"
        for k, v in sel.items():
            arm, r = k.rsplit("_r", 1)
            out[(arm, split, int(r))] = {"selected": v["selected"], "model_seed": v["model_seed"],
                                         "batch_size": int(v.get("batch_size") or 10)}
        if not any(a == RANDOM for a, _, _ in out):
            # the final stores random permutations under 'R0' too; if absent, rebuild is out of scope here
            pass
    if arms is not None:
        out = {k: v for k, v in out.items() if k[0] in set(arms)}
    return out


def dev_bundle():
    """Dev bundle (T rows physically removed) through the frozen loader; caller sets threads first."""
    sys.path[:0] = [str(ROOT / "src"), str(ROOT / "third_party")] if str(ROOT / "src") not in sys.path else []
    import hsl  # noqa: F401
    from hsl import views
    return views.dev_bundle()


def split_obj(name: str):
    from hsl import views
    if name == "final":
        return views.final_split()
    return views.get_split(name)


# ============================================================================ small numerics
def loglin_interp(x_grid, y_grid, x):
    """Linear interpolation of y on log(x) (grid increasing in x)."""
    return np.interp(np.log(np.asarray(x, float)), np.log(np.asarray(x_grid, float)), np.asarray(y_grid, float))


def inverse_curve(n_grid, curve, target: float) -> float:
    """Smallest n (log-n linear interpolation of the PAVA-isotonic curve) with curve(n) >= target, by the frozen
    hsl.eval.n_req: 0 if left-censored (target below the first point), +inf if right-censored (never reached)."""
    from hsl import eval as E
    return float(E.n_req(float(target), list(map(float, n_grid)), list(map(float, curve))).n)


def write_json(p: Path, obj) -> None:
    def d(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    p.write_text(json.dumps(obj, indent=1, default=d))


# ============================================================================ evaluation labels
def eval_labels(source: str, split: str) -> dict:
    """E_fixed(_dev) cell labels of a split: targets, rows (index into targets), cols (bundle array index), y and the
    on-array target-column stratum 'OA' (dev VV / final TT). dev: analysis.q3.cache; final/dryrun: the scorer's
    <out>/run/labels.npz (final: T labels — only after scoring has finished)."""
    from analysis.q3 import cache
    if source == "dev":
        e = cache.eset(split)
        return {"targets": e["targets"], "rows": e["rows"], "cols": e["cols"], "y": e["y"], "OA": e["VV"],
                "AO": e["AO"], "array_orfs": cache.load()["array_orfs"]}
    with np.load(final_dir(source) / "run" / "labels.npz") as z:
        return {"targets": z["targets"], "rows": z["rows"], "cols": z["cols"], "y": z["y"].astype(np.int8),
                "OA": z["TT"], "AO": z["AO"], "array_orfs": cache.load()["array_orfs"]}


def load_pred(pred_file: str, key: str = "E_fixed") -> np.ndarray:
    with np.load(ROOT / pred_file) as z:
        return np.asarray(z[key], dtype=np.float32)
