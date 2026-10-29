"""P1 helpers (DESIGN v3.1 §5, §7.1, §8, §11): dev random curves, power pilot, upstream-protocol sanity.

Builds ONLY on the frozen P0 infrastructure (hsl.harness / hsl.views / hsl.eval); nothing here is frozen.

Random trajectories (R0, shared with Q2 and the P2 gate)
  random_traj(split, r)    exactly what harness.Trajectory(RandomSelector(), split, seeds[r], r) selects:
                           seeds = SeedSequence(dev_traj_base_seed).spawn(R)[r]; selection stream (child 0) draws ONE
                           permutation of the sorted candidate pool; the evaluation stream (child 1) gives the model
                           seed; checkpoint n = first n genes. traj_id is computed with the harness' cfg hash (arm
                           "R0_random", max_n = pool size), and run_predictor gets the same seeds dict, so run ids are
                           those a harness Trajectory would produce.
Grids
  dev_grid(pool)           budgets.yaml random_grid truncated below the pool size, plus the pool maximum.
  m1_grid(pool)            budgets.yaml m1_hgb_q1_grid truncated likewise, plus the pool maximum.
Runs
  run_cell(...)            harness.run_predictor on E_fixed_dev (ledger line + float32 cell predictions saved under
                           results/predictions/p1/ by pointing harness.PRED_DIR there) + extra metrics: row strata
                           (on-/off-array target genes), precision at 10% recall, SLs detected by the screened rows,
                           peak RSS sampled DURING the fit (ps poll; ru_maxrss is process-cumulative).
  run_upstream_protocol    stack_lr_wr2 on the same screened set, scored on E_fixed_dev AND on every measured D_val
                           cell except test-gene columns (= views.e_unobs with S = ∅: incl. screened arrays, i.e.
                           reciprocal-known cells — the upstream YEAST_SLs protocol, dev-legal) with upstream strata
                           (recip_available, x_on_array_no_recip, x_off_array). Own ledger line (tag).
Bootstrap
  boot_seed(name)          SeedSequence(bootstrap_seed, spawn_key=(crc32(name),)): a named child per analysis.
  query_weights(split, name)  stratified query bootstrap weights (B = protocol B_query_bootstrap, frozen strata).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import json
import os
import resource
import subprocess
import sys
import threading
import time
import zlib
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

import numpy as np

from hsl import config, eval as E, harness as H, paths, views

RESULTS = paths.output_root() / "results"
P1_DIR = RESULTS / "p1"
PRED_P1 = RESULTS / "predictions" / "p1"
PREDICTORS = {"stack_lr_wr2": {}, "stack_lr": {}, "B1_array_rate": {}, "B2_product": {},
              "m1_hgb": {"drop": ["genome"]}}
PSTAR = "stack_lr_wr2"
# PI disk policy (2026-10-29, review P0 M2): only P* (stack_lr_wr2) and m1_hgb predictions are saved; others metrics only
NO_SAVE = ("stack_lr", "B1_array_rate", "B2_product")


def use_p1_prediction_dir(root: Path | None = None) -> Path:
    """Point the harness' prediction directory at results/predictions/p1 (or `root`) for this process."""
    H.PRED_DIR = Path(root or PRED_P1)
    return H.PRED_DIR


# ============================================================================ grids / trajectories
def dev_grid(pool_size: int) -> list[int]:
    g = [int(n) for n in config.budgets()["random_grid"] if int(n) < int(pool_size)]
    return g + [int(pool_size)]


def m1_grid(pool_size: int) -> list[int]:
    g = [int(n) for n in config.budgets()["m1_hgb_q1_grid"] if int(n) < int(pool_size)]
    return g + [int(pool_size)]


@dataclass(frozen=True)
class RandomTraj:
    split: str
    r: int
    perm: tuple
    model_seed: int
    traj_id: str
    cfg: dict


def random_traj(split_name: str, r: int, pool=None, arm: str = "R0_random") -> RandomTraj:
    """Permutation + model seed of dev random trajectory r (identical to harness.Trajectory + RandomSelector)."""
    pool = views.get_split(split_name).pool if pool is None else np.asarray(pool, dtype=str)
    ts = config.trajectory_seeds("dev", int(r) + 1)[int(r)]
    sel_ss, eval_ss = H.stream_seeds(ts)
    sel = H.RandomSelector()
    sel.reset(np.setdiff1d(pool, []), np.random.default_rng(sel_ss))
    perm = tuple(str(g) for g in sel._perm)
    ms = H.model_seed_of(eval_ss)
    bud = config.budgets()
    cfg = {"arm": arm, "split": split_name, "traj_index": int(r),
           "traj_seed": {"entropy": str(ts.entropy), "spawn_key": list(ts.spawn_key)}, "model_seed": ms,
           "batch_size": int(bud["batch_size"]), "max_n": int(len(pool)), "s0_n": 0, "selector_params": {},
           "static": True}
    return RandomTraj(split_name, int(r), perm, ms, f"{split_name}_{arm}_r{int(r)}_{H._hash(cfg)}", cfg)


def write_traj_log(t: RandomTraj, d: Path | None = None) -> Path:
    d = Path(d or P1_DIR / "trajectories")
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{t.traj_id}.json"
    if not p.exists():
        p.write_text(json.dumps({"traj_id": t.traj_id, **t.cfg, "selected": list(t.perm)}, indent=0) + "\n")
    return p


# ============================================================================ RSS sampling
def _rss_gb(pid: int) -> float:
    try:
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
        return int(out) / 2**20 if out else 0.0
    except (OSError, ValueError):
        return 0.0


class PeakRSS:
    """Peak resident set size of this process while the block runs (ps polled every `interval` s)."""

    def __init__(self, interval: float = 0.5):
        self.interval, self.peak, self._stop = float(interval), 0.0, threading.Event()

    def _run(self):
        pid = os.getpid()
        while not self._stop.is_set():
            self.peak = max(self.peak, _rss_gb(pid))
            self._stop.wait(self.interval)

    def __enter__(self):
        self.peak = _rss_gb(os.getpid())
        self._t = threading.Thread(target=self._run, daemon=True)
        self._t.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._t.join()
        self.peak = max(self.peak, _rss_gb(os.getpid()))


def proc_maxrss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


# ============================================================================ metrics
def precision_at_recall(y, s, recall: float = 0.1) -> float:
    """Precision at the first (tie-grouped) score threshold whose recall reaches `recall`."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(s, dtype=np.float64)
    if not y.any():
        return float("nan")
    o = np.argsort(-s, kind="mergesort")
    ss, yy = s[o], y[o]
    ends = np.r_[np.flatnonzero(ss[1:] != ss[:-1]), len(ss) - 1]
    tp = np.cumsum(yy)[ends]
    k = int(np.argmax(tp / tp[-1] >= recall))
    return float(tp[k] / (ends[k] + 1))


@lru_cache(maxsize=8)
def target_on_array(split_name: str) -> np.ndarray:
    b = views.dev_bundle()
    sp = views.get_split(split_name)
    return np.array([o in b.orf_to_col for o in sp.targets])


@lru_cache(maxsize=1)
def _dev_sl_total() -> int:
    b = views.dev_bundle()
    return int((b.Y == 1).sum())


def extra_metrics(split_name: str, es, s) -> dict:
    on = target_on_array(split_name)[es.rows]
    out = {"prec_at_recall10": precision_at_recall(es.y, s, 0.1)}
    for nm, m in (("on_array_rows", on), ("off_array_rows", ~on)):
        out[f"auprc_{nm}"] = E.average_precision(es.y[m], s[m]) if es.y[m].any() else float("nan")
        out[f"prev_{nm}"] = float(es.y[m].mean()) if m.any() else float("nan")
    return out


def sl_detected(S) -> dict:
    b = views.dev_bundle()
    k = int(views.query_sl_counts(b, np.asarray(S, dtype=str)).sum())
    return {"sl_detected": k, "sl_detected_frac_dev": k / _dev_sl_total()}


# ============================================================================ runs
def run_cell(t: RandomTraj, n: int, method: str, n_threads: int | None = None, params: dict | None = None,
             label: str = "p1") -> dict:
    """One (split, trajectory, budget, predictor) evaluation on E_fixed_dev -> flat result row."""
    S = t.perm[:int(n)]
    prm = dict(PREDICTORS.get(method, {}) if params is None else params)
    t0 = time.perf_counter()
    with PeakRSS() as pk:
        res = H.run_predictor(H.RunSpec(t.split, tuple(S), ("E_fixed",), ("E_fixed",), label=label), method, prm,
                              {"model_seed": t.model_seed, "traj_index": t.r}, n_threads=n_threads, parent=t.traj_id,
                              save=method not in NO_SAVE)
    wall = time.perf_counter() - t0
    es = res.esets["E_fixed"]
    s = res.scores["E_fixed"]
    m = dict(res.metrics[es.name])
    rt = res.manifest["runtime_s"]
    row = {"split": t.split, "perm": t.r, "n": int(n), "predictor": method, "model_seed": t.model_seed,
           "run_id": res.run_id, "traj_id": t.traj_id, "eset": es.name, "eset_digest": es.digest(),
           **{k: m[k] for k in ("n_cells", "n_pos", "prevalence", "auprc", "lift", "macro_auprc", "macro_n_used",
                                "macro_n_zero_pos", "auroc", "p_at_10", "p_at_50", "auprc_AO", "prev_AO",
                                "auprc_VV", "prev_VV")},
           **extra_metrics(t.split, es, s), **sl_detected(S),
           "view_s": rt.get("view_s"), "fit_predict_s": rt.get("fit_predict_s"), "eset_s": rt.get("eset_s"),
           "save_s": rt.get("save_s"), "total_s": rt.get("total_s"), "wall_s": wall,
           "peak_rss_fit_gb": pk.peak, "proc_maxrss_gb": proc_maxrss_gb(),
           "n_threads": res.manifest["n_threads"], "pred_file": res.manifest.get("predictions", {}).get(
               "E_fixed", {}).get("file")}
    return row


def upstream_strata(b, split, es, S) -> dict:
    """Upstream evaluator strata (ysl.eval.evaluator._build_strata) on an ESet of a dev split."""
    cmap = b.orf_to_col
    col_x = np.array([cmap.get(o, -1) for o in split.targets], dtype=np.int64)
    S = np.sort(np.asarray(S, dtype=str))
    arr_tr = np.full(b.n_arrays, -1, dtype=np.int64)
    rows = b.gene_rows(S)
    for i, o in enumerate(S):
        j = cmap.get(o, -1)
        if j >= 0:
            arr_tr[j] = i
    obs_train = b.Y[rows] >= 0
    cx = col_x[es.rows]
    ta = arr_tr[es.cols]
    rec = (ta >= 0) & (cx >= 0)
    rec[rec] = obs_train[ta[rec], cx[rec]]
    on = cx >= 0
    return {"recip_available": rec, "x_on_array_no_recip": on & ~rec, "x_off_array": ~on,
            "array_screened": ta >= 0}


def run_upstream_protocol(t: RandomTraj, n: int, method: str = PSTAR, n_threads: int | None = None) -> dict:
    """Score `method` on E_fixed_dev and on the upstream-protocol dev set (all measured D_val cells, minus test-gene
    columns and identity self pairs; screened arrays INCLUDED). Mirrors harness.run_predictor (predictions saved
    before scoring; one ledger line, tag 'p1_upstream_protocol')."""
    from hsl import upstream as U
    pr = config.protocol()
    n_threads = int(n_threads or os.environ.get("HSL_THREADS") or pr["compute"]["threads_default"])
    split = views.get_split(t.split)
    b = views.dev_bundle()
    S = np.sort(np.asarray(t.perm[:int(n)], dtype=str))
    n_folds = int(min(int(pr["predictor"]["n_folds"]), max(len(S), 2)))
    prm = H._method_params(method, PREDICTORS.get(method, {}), len(S))
    seeds = {"model_seed": t.model_seed, "traj_index": t.r, "split_seed": 0}
    tag = "p1_upstream_protocol"
    cfg = {"split": t.split, "split_digest": split.digest(), "method": method, "params": prm, "seeds": seeds,
           "screened_sha": H._hash(S.tolist()), "n_folds": n_folds, "tag": tag, "score_fn": None}
    run_id = f"{t.split}_{method}_n{len(S)}_{H._hash(cfg)}"
    rec = {"id": run_id, "parent": t.traj_id, "config_hash": H._hash(cfg), "split": t.split, "method": method,
           "n_screened": int(len(S)), "outcome": "error", "metrics": {}}
    man = {"run_id": run_id, "tag": tag, "params": prm, "seeds": seeds, "n_folds_passed": n_folds,
           "n_threads": n_threads, "git": H._git(), "split_digest": split.digest()}
    stages = {}
    try:
        t0 = time.perf_counter()
        view = U.construct_view(b, S, split.targets)
        Smat, info = U.run_method(view, method, model_seed=t.model_seed, split_seed=0, n_threads=n_threads,
                                  n_folds=n_folds, chunk_targets=int(pr["predictor"]["chunk_targets"]), params=prm)
        stages["fit_predict_s"] = time.perf_counter() - t0
        es_fx = H.get_eset(b, split, "E_fixed")
        es_up = replace(views.e_unobs(b, split, []), name="E_upstream_dev")
        esets = {"E_fixed": es_fx, "E_upstream_dev": es_up}
        scores = {k: es.scores(Smat).astype(np.float32) for k, es in esets.items()}
        if method not in NO_SAVE:
            man["predictions"] = H.save_predictions(run_id, split, scores, esets)
        metrics = {}
        for k, es in esets.items():
            m = E.summarise(es.y, scores[k], es.rows)
            for st, mask in es.strata.items():
                m[f"auprc_{st}"] = E.average_precision(es.y[mask], scores[k][mask]) if es.y[mask].any() else None
            ups = upstream_strata(b, split, es, S)
            for st, mask in ups.items():
                m[f"auprc_{st}"] = E.average_precision(es.y[mask], scores[k][mask]) if es.y[mask].any() else None
                m[f"prev_{st}"] = float(es.y[mask].mean()) if mask.any() else None
                m[f"n_cells_{st}"] = int(mask.sum())
            metrics[es.name] = m
        rec.update(outcome="ok", metrics=metrics)
    finally:
        stages["total_s"] = time.perf_counter() - t0
        man["runtime_s"] = stages
        man["peak_rss_gb"] = proc_maxrss_gb()
        rec["manifest"] = man
        H.append_ledger(rec)
    flat = {"split": t.split, "perm": t.r, "n": int(n), "predictor": method, "run_id": run_id,
            "fit_predict_s": stages["fit_predict_s"]}
    for k, m in metrics.items():
        for kk, v in m.items():
            flat[f"{k}.{kk}"] = v
    return flat


# ============================================================================ bootstrap helpers
def boot_seed(name: str) -> np.random.SeedSequence:
    return np.random.SeedSequence(config.seed("bootstrap_seed"), spawn_key=(zlib.crc32(name.encode()),))


def query_weights(split_name: str, name: str, B: int | None = None) -> np.ndarray:
    B = int(config.protocol()["bootstrap"]["B_query_bootstrap"] if B is None else B)
    sp = views.get_split(split_name)
    strata = views.query_strata(sp.targets, views.dev_bundle())
    return E.bootstrap_weights(strata, B, np.random.default_rng(boot_seed(name)))


def load_scores(pred_file: str, eset_name: str = "E_fixed") -> np.ndarray:
    with np.load(paths.PROJECT_ROOT / pred_file) as z:
        return np.asarray(z[eset_name], dtype=np.float32)


# ============================================================================ part files
def append_row(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    H.append_ledger(row, path)          # same locked JSON-lines append


def read_rows(paths_) -> list[dict]:
    out = []
    for p in paths_:
        p = Path(p)
        if p.exists():
            out.extend(json.loads(ln) for ln in p.read_text().splitlines() if ln.strip())
    return out
