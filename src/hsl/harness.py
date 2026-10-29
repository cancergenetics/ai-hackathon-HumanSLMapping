"""Run harness (DESIGN v3 §5, §6, §8): predictor runs on arbitrary screened sets, manifests, the attempt ledger,
saved evaluation-cell predictions, and the Trajectory runner with seeded selection/evaluation streams.

  run_predictor(spec, method, params, seeds, bundle=None) -> RunResult
      spec = RunSpec(split, screened, eval_sets): builds the upstream ScreenView (dev: on views.dev_bundle() via
      ScreenView._construct; final: the caller must pass views.final_bundle — the harness never loads it), runs
      the registered method (≤ protocol threads), extracts scores on the requested evaluation cells, SAVES them
      (float32, evaluation cells only, results/predictions/, gitignored) BEFORE scoring, then scores (hsl.eval) and
      appends ONE line to results/LEDGER.jsonl (id, parent, config hash, split, outcome, metrics, manifest).
      Manifest: git head + dirty flag + diff hash, seeds (model_seed, split_seed, n_folds passed to run_method),
      method + params, input checksums (bundle SHA-256 recorded upstream, split digest), n screened, runtime by
      stage, peak RSS, threads.
  Selector protocol: an object with `needs_view` (bool), `static` (bool), `reset(candidates, rng)` and
      `select(view, k, rng) -> ordered list of k genes` (greedy order). RandomSelector = nested prefixes of a
      seeded permutation of the candidate pool.
  Trajectory(selector, split, seed, ...).run(): S_0 then rounds of `batch_size`; selection uses ONLY the
      SelectorView and the selection stream; checkpoints n evaluate the predictors on the first n genes of the
      ordered selection (a checkpoint not divisible by the batch size takes the first n − b⌊n/b⌋ genes of the next
      batch in greedy order) with a separate evaluation stream, so adding a checkpoint never changes a trajectory.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import fcntl
import hashlib
import json
import os
import resource
import sys
import time
import traceback
from dataclasses import asdict, dataclass, field
from functools import lru_cache

import numpy as np

from hsl import config, paths, views

RESULTS = paths.output_root() / "results"
LEDGER = RESULTS / "LEDGER.jsonl"
PRED_DIR = RESULTS / "predictions"
TRAJ_DIR = RESULTS / "trajectories"


# ============================================================================ helpers
def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


def _hash(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


@lru_cache(maxsize=1)
def _git() -> dict:
    return views.git_state()


def append_ledger(rec: dict, path=None) -> None:
    path = LEDGER if path is None else path
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(rec, sort_keys=True, default=_json_default) + "\n"
    with open(path, "a", encoding="utf-8") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(line)
        f.flush()
        fcntl.flock(f, fcntl.LOCK_UN)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    return str(o)


_ESETS: dict = {}


def get_eset(b, split: views.Split, name: str, screened=None) -> views.ESet:
    """Cached evaluation sets (E_fixed per split; E_unobs per screened set)."""
    key = (id(b), split.digest(), name, None if name != "E_unobs" else _hash(sorted(screened)))
    if key not in _ESETS:
        if len(_ESETS) > 16:
            _ESETS.clear()
        _ESETS[key] = views.eset_for(b, split, name, screened)
    return _ESETS[key]


def _bundle_for(split: views.Split, bundle):
    if split.kind == "final":
        if bundle is None:
            raise PermissionError("final runs need the caller (scripts/run_final.py) to pass views.final_bundle")
        return bundle
    return views.dev_bundle() if bundle is None else bundle


# ============================================================================ predictor runs
@dataclass(frozen=True)
class RunSpec:
    split: str                          # dev1 | dev2 | dev3 | resplit | final
    screened: tuple                     # screened genes (any order; stored sorted)
    eval_sets: tuple = ("E_fixed",)     # "E_fixed" (= E_fixed_dev on dev splits), "E_unobs"
    save_sets: tuple = ("E_fixed",)     # which evaluation sets' cell predictions are saved (float32)
    label: str = ""


@dataclass
class RunResult:
    run_id: str
    scores: dict                        # eval set name -> float32 scores aligned with the ESet cells
    metrics: dict
    manifest: dict
    esets: dict = field(default_factory=dict)


def _method_params(method: str, params: dict | None, n: int) -> dict:
    p = dict(params or {})
    pr = config.protocol()["predictor"]
    if method.startswith("m1_hgb") and n < int(pr["m1_hgb_small_n"]) and "prop_folds" not in p:
        p["prop_folds"] = int(pr["m1_hgb_small_n_prop_folds"])
    return p


def run_predictor(spec: RunSpec, method: str, params: dict | None = None, seeds: dict | None = None, *,
                  bundle=None, n_threads: int | None = None, parent: str | None = None, save: bool = True,
                  ledger: bool = True, score_fn=None, tag: str = "") -> RunResult:
    """See module docstring. `score_fn(view, ctx_info) -> S` replaces the upstream method (test-only controls such
    as the oracle back door); it is recorded in the manifest."""
    from hsl import eval as E
    from hsl import upstream as U
    t_all = time.perf_counter()
    seeds = dict(seeds or {})
    model_seed = int(seeds.setdefault("model_seed", 0))
    split_seed = int(seeds.setdefault("split_seed", 0))
    pr = config.protocol()
    n_threads = int(n_threads or os.environ.get("HSL_THREADS") or pr["compute"]["threads_default"])
    split = views.get_split(spec.split)
    b = _bundle_for(split, bundle)
    S = np.sort(np.asarray(list(spec.screened), dtype=str))
    if not np.isin(S, split.pool).all():
        raise ValueError("screened genes must belong to the split's candidate pool")
    n = len(S)
    n_folds = int(min(int(pr["predictor"]["n_folds"]), max(n, 2)))
    prm = _method_params(method, params, n)
    stages = {}
    config_id = {"split": spec.split, "split_digest": split.digest(), "method": method, "params": prm,
                 "seeds": seeds, "screened_sha": _hash(S.tolist()), "n_folds": n_folds, "tag": tag,
                 "score_fn": getattr(score_fn, "__name__", None)}
    config_hash = _hash(config_id)
    run_id = f"{spec.split}_{method}_n{n}_{config_hash}"
    man = {"run_id": run_id, "config_hash": config_hash, "parent": parent, "split": spec.split, "kind": split.kind,
           "split_digest": split.digest(), "method": method, "params": prm, "seeds": seeds,
           "model_seed": model_seed, "split_seed_passed": split_seed, "n_folds_passed": n_folds,
           "n_threads": n_threads, "n_screened": n, "n_targets": int(len(split.targets)), "label": spec.label,
           "tag": tag, "git": _git(), "bundle_sha256": b.meta.get("bundle_sha256"),
           "bundle_view": b.meta.get("hsl_view"), "score_fn": config_id["score_fn"],
           "yeast_sls_commit": (paths.PROJECT_ROOT / "config" / "yeast_sls_commit.txt").read_text().strip()}
    rec = {"id": run_id, "parent": parent, "config_hash": config_hash, "split": spec.split, "method": method,
           "n_screened": n, "outcome": "error", "metrics": {}}
    try:
        t = time.perf_counter()
        view = U.construct_view(b, S, split.targets)
        stages["view_s"] = time.perf_counter() - t
        t = time.perf_counter()
        if score_fn is None:
            Smat, info = U.run_method(view, method, model_seed=model_seed, split_seed=split_seed, n_threads=n_threads,
                                      n_folds=n_folds, chunk_targets=int(pr["predictor"]["chunk_targets"]),
                                      params=prm)
            man["method_info"] = info.get("method_info")
        else:
            Smat = np.asarray(score_fn(view, b), dtype=np.float32)
        stages["fit_predict_s"] = time.perf_counter() - t
        if not np.isfinite(Smat).all():
            raise FloatingPointError("non-finite scores")
        t = time.perf_counter()
        esets = {nm: get_eset(b, split, nm, S) for nm in spec.eval_sets}
        scores = {nm: es.scores(Smat).astype(np.float32) for nm, es in esets.items()}
        stages["eset_s"] = time.perf_counter() - t
        if save:                                       # predictions saved BEFORE scoring
            t = time.perf_counter()
            man["predictions"] = save_predictions(run_id, split, {k: scores[k] for k in spec.save_sets if k in scores},
                                                  {k: esets[k] for k in spec.save_sets if k in esets})
            stages["save_s"] = time.perf_counter() - t
        t = time.perf_counter()
        metrics = {}
        for nm, es in esets.items():
            m = E.summarise(es.y, scores[nm], es.rows)
            for st, mask in es.strata.items():
                m[f"auprc_{st}"] = E.average_precision(es.y[mask], scores[nm][mask]) if es.y[mask].any() else None
                m[f"prev_{st}"] = float(es.y[mask].mean()) if mask.any() else None
            metrics[es.name] = m
        stages["score_s"] = time.perf_counter() - t
        rec.update(outcome="ok", metrics=metrics)
        result = RunResult(run_id, scores, metrics, man, esets)
    except Exception as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        man["traceback"] = traceback.format_exc()
        result = None
        raise_exc = e
    else:
        raise_exc = None
    finally:
        stages["total_s"] = time.perf_counter() - t_all
        man["runtime_s"] = stages
        man["peak_rss_gb"] = peak_rss_gb()
        rec["manifest"] = man
        if ledger:
            append_ledger(rec)
    if raise_exc is not None:
        raise raise_exc
    return result


def save_predictions(run_id: str, split: views.Split, scores: dict, esets: dict) -> dict:
    """float32 scores on evaluation cells (results/predictions/<split>/<run_id>.npz); the cell index of each ESet is
    written once (results/predictions/_cells/<eset digest>.npz: targets, rows, cols)."""
    out = {}
    d = PRED_DIR / split.name
    d.mkdir(parents=True, exist_ok=True)
    cells = PRED_DIR / "_cells"
    cells.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for nm, s in scores.items():
        es = esets[nm]
        dg = es.digest()
        cp = cells / f"{dg}.npz"
        if not cp.exists():
            np.savez(cp, targets=es.targets, rows=es.rows, cols=es.cols)
        arrays[nm] = np.asarray(s, dtype=np.float32)
        out[nm] = {"cells": os.path.relpath(cp, paths.PROJECT_ROOT), "eset_digest": dg, "n": int(len(s))}
    p = d / f"{run_id}.npz"
    np.savez(p, **arrays)
    for v in out.values():
        v["file"] = os.path.relpath(p, paths.PROJECT_ROOT)
    return out


# ============================================================================ selectors
class RandomSelector:
    """R0: nested prefixes of ONE seeded permutation of the (sorted) candidate pool, drawn at reset()."""
    needs_view = False
    static = True
    name = "R0_random"

    def reset(self, candidates, rng) -> None:
        c = np.sort(np.asarray(candidates, dtype=str))
        self._perm = list(c[rng.permutation(len(c))])
        self._pos = 0

    def select(self, view, k: int, rng) -> list:
        out = self._perm[self._pos:self._pos + int(k)]
        self._pos += len(out)
        return out


def randomised_greedy(gains: np.ndarray, k: int, rng, top_k: int | None = None, deterministic: bool = False,
                      update=None) -> list:
    """Pick k indices; each pick uniform among the top_k by current gain (protocol selector.randomised_greedy_top_k),
    or the argmax when deterministic (temperature 0). `update(gains, picked_idx) -> gains` refreshes gains after a
    pick (greedy within a batch); picked indices are never re-picked. Returns indices in greedy order."""
    top_k = int(config.protocol()["selector"]["randomised_greedy_top_k"] if top_k is None else top_k)
    g = np.array(gains, dtype=np.float64, copy=True)
    picked = []
    for _ in range(int(min(k, len(g)))):
        avail = np.flatnonzero(np.isfinite(g) & ~np.isin(np.arange(len(g)), picked))
        if not len(avail):
            break
        order = avail[np.lexsort((avail, -g[avail]))]      # by gain desc, index asc (deterministic ties)
        i = int(order[0]) if deterministic else int(order[rng.integers(0, min(top_k, len(order)))])
        picked.append(i)
        if update is not None:
            g = np.array(update(g, i), dtype=np.float64)
        g[i] = -np.inf
    return picked


def stream_seeds(traj_seed: np.random.SeedSequence) -> tuple[np.random.SeedSequence, np.random.SeedSequence]:
    """(selection stream, evaluation stream) children of a trajectory seed — derived without mutating it."""
    mk = lambda i: np.random.SeedSequence(traj_seed.entropy, spawn_key=tuple(traj_seed.spawn_key) + (i,))  # noqa
    return mk(0), mk(1)


def model_seed_of(eval_ss: np.random.SeedSequence) -> int:
    return int(eval_ss.generate_state(1, dtype=np.uint32)[0])


# ============================================================================ trajectories
@dataclass
class Trajectory:
    selector: object
    split: str
    traj_seed: np.random.SeedSequence
    traj_index: int = 0
    arm: str = ""
    batch_size: int | None = None
    checkpoints: tuple = ()
    max_n: int | None = None
    predictors: tuple = ()             # ((method, params), ...) evaluated at every checkpoint
    s0: tuple = ()                     # warm start (ordered)
    eval_sets: tuple = ("E_fixed",)
    save: bool = True
    ledger: bool = True
    write_log: bool = True
    bundle: object = None              # final: passed by scripts/run_final.py

    def run(self) -> dict:
        bud = config.budgets()
        b_size = int(self.batch_size or bud["batch_size"])
        split = views.get_split(self.split)
        b = _bundle_for(split, self.bundle)
        sel_ss, eval_ss = stream_seeds(self.traj_seed)
        sel_rng = np.random.default_rng(sel_ss)
        model_seed = model_seed_of(eval_ss)
        pool = split.pool
        S = [str(g) for g in self.s0]
        if not np.isin(S, pool).all() or len(set(S)) != len(S):
            raise ValueError("S_0 must be distinct pool genes")
        max_n = int(min(self.max_n or bud["q2_max_n"], len(pool)))
        cps = sorted({int(c) for c in self.checkpoints if len(S) <= int(c) <= max_n})
        sel = self.selector
        sel.reset(np.setdiff1d(pool, S), sel_rng)
        cfg = {"arm": self.arm or getattr(sel, "name", type(sel).__name__), "split": self.split,
               "traj_index": int(self.traj_index), "traj_seed": {"entropy": str(self.traj_seed.entropy),
                                                                 "spawn_key": list(self.traj_seed.spawn_key)},
               "model_seed": model_seed, "batch_size": b_size, "max_n": max_n, "s0_n": len(S),
               "selector_params": getattr(sel, "params", {}), "static": bool(getattr(sel, "static", False))}
        traj_id = f"{self.split}_{cfg['arm']}_r{self.traj_index}_{_hash(cfg)}"
        log = {"traj_id": traj_id, **cfg, "rounds": [], "checkpoints": {}, "git": _git()}
        done_cp = set()
        rnd = 0
        while len(S) < max_n:
            cand = np.setdiff1d(pool, S)
            if not len(cand):
                break
            k = int(min(b_size, len(cand)))
            t = time.perf_counter()
            view = views.SelectorView.build(b, S, split) if getattr(sel, "needs_view", True) else None
            t_view = time.perf_counter() - t
            t = time.perf_counter()
            picks = [str(g) for g in sel.select(view, k, sel_rng)]
            t_sel = time.perf_counter() - t
            if len(picks) != k or len(set(picks)) != k or not np.isin(picks, cand).all():
                raise ValueError(f"selector returned an invalid batch in round {rnd}: {picks[:5]}...")
            S.extend(picks)
            log["rounds"].append({"round": rnd, "n_before": len(S) - k, "picks": picks, "view_s": t_view,
                                  "select_s": t_sel})
            rnd += 1
            for c in cps:                          # evaluate prefixes reached by this batch (greedy order)
                if c not in done_cp and len(S) >= c:
                    log["checkpoints"][str(c)] = self._evaluate(S[:c], c, model_seed, traj_id, b)
                    done_cp.add(c)
        for c in cps:
            if c not in done_cp and len(S) >= c:
                log["checkpoints"][str(c)] = self._evaluate(S[:c], c, model_seed, traj_id, b)
        log["selected"] = S
        log["peak_rss_gb"] = peak_rss_gb()
        if self.write_log:
            d = TRAJ_DIR / self.split
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{traj_id}.json").write_text(json.dumps(log, indent=1, default=_json_default))
        return log

    def _evaluate(self, S, n, model_seed, traj_id, b) -> dict:
        out = {}
        for method, params in self.predictors:
            r = run_predictor(RunSpec(self.split, tuple(S), self.eval_sets, ("E_fixed",)), method, params,
                              {"model_seed": model_seed, "traj_index": int(self.traj_index)}, bundle=self.bundle,
                              parent=traj_id, save=self.save, ledger=self.ledger)
            out[method] = {"run_id": r.run_id, "metrics": r.metrics, "runtime_s": r.manifest["runtime_s"]}
        return out
