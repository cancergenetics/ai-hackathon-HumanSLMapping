"""Experiment runner (DESIGN.md §8, §11; r1 agent-failure R2/R7/R8). The only code that writes metrics.json,
manifest.json and the ledger.

    from ysl.runner import run
    run("stack_lr", dataset="s1", split_kind="dev", seeds=[0, 1, 2], config={"oracle": True})

For each split seed: split (committed lockbox) -> ScreenView.from_split -> method.fit_predict(view, ctx) ->
Evaluator.evaluate (+ optional oracle diagnostic, poisoning self-check, keyed prediction file). Writes
    results/<exp_id>/metrics.json   per split + mean/SD over splits
    results/<exp_id>/manifest.json  git commit, dirty-tree hash of src/scripts/config, config hash, data checksums,
                                    split fingerprints, seeds, package versions, hardware, wall-clock by stage,
                                    pytest status, peak RSS, metrics SHA-256
    results/LEDGER.jsonl            one appended line per run (failed runs included, status != complete)
    data/processed/predictions/<exp_id>/<dataset>_<kind>_<seed>.parquet  [query_strain_id, array_strain_id, score]
                                    (observed test cells only; optional) — the input of scripts/rescore.py
Final splits run ONLY inside `run_final_batch` (scripts/32_final.py): one shot over all pre-declared methods, each
from a committed recipe in config/recipes/; it refuses if results/final/LOCK exists, holds results/final/LOCK.pending
(with a batch token) while running, computes the pre-registered confirmatory analysis and writes LOCK once at the end.
The pre-registered SECONDARY finals (DESIGN §4: TheCellMap release; train fractions 10%/30%) run ONLY inside
`run_secondary_batch` (scripts/33_secondary.py), after the primary LOCK exists, with exactly the recipes, seeds,
candidates and comparator recorded in it; each (dataset, train fraction) is one-shot under
results/final/secondary_<dataset>_f<frac>/{LOCK,LOCK.pending,confirmatory.json} (docs/harness_secondary.md).
While a method's fit_predict runs, `data_guard()` blocks opening any Costanzo raw/processed/bundle/prediction file
(external features under data/processed/features/ stay readable) and hides the in-process bundle cache.
"""
from __future__ import annotations

import builtins  # noqa: F401  (documented: the guard relies on the 'open' audit event, which builtins.open raises)
import contextlib
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from ysl._paths import OUTPUT_ROOT, UPSTREAM_ROOT  # humanSLmapping PATCH: read-only upstream root, writable output
ROOT = UPSTREAM_ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[2]); data guard + recipes (read-only)
RESULTS = OUTPUT_ROOT / "results" / "ysl_runner"  # humanSLmapping PATCH (was ROOT / "results")
LEDGER = RESULTS / "LEDGER.jsonl"
LOCK = RESULTS / "final" / "LOCK"
PRED_DIR = OUTPUT_ROOT / "data" / "ysl_predictions"  # humanSLmapping PATCH (was ROOT / "data/processed/predictions")
UPSTREAM_PRED_DIR = UPSTREAM_ROOT / "data" / "processed" / "predictions"  # humanSLmapping PATCH (still protected)
TREE_DIRS = ("src", "scripts", "config")
SNAPSHOT_DIRS = ("src", "scripts", "config", "tests")
RECIPE_DIR = ROOT / "config" / "recipes"
SUMMARY_KEYS = ("auprc", "auroc", "lift", "macro_auprc", "macro_auprc_random_baseline", "prevalence")


# ----------------------------------------------------------------------------- provenance helpers
def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _file_hashes(dirs) -> dict:
    """{relative path: SHA-256} of every tracked-type file under `dirs` (ROOT-relative), in tree_hash order."""
    out = {}
    for d in dirs:
        for p in sorted((ROOT / d).rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts and p.suffix not in (".pyc",) and p.name != ".DS_Store":
                out[str(p.relative_to(ROOT))] = _sha256(p.read_bytes())
    return out


def tree_hash(dirs=TREE_DIRS) -> dict:
    """SHA-256 over (path, file SHA-256) of every tracked-type file under src/, scripts/, config/."""
    h = hashlib.sha256()
    files = _file_hashes(dirs)
    for rel, sha in files.items():
        h.update(rel.encode())
        h.update(sha.encode())
    return {"sha256": h.hexdigest(), "n_files": len(files)}


CODE_PIN_DIRS = ("src/ysl", "config")


def code_state(dirs=CODE_PIN_DIRS) -> dict:
    """Per-file SHA-256 of the code and configuration a final batch depends on (src/ysl/, config/). The primary final
    records it in LOCK; a secondary final refuses to run on different code unless given a recorded reason."""
    files = _file_hashes(dirs)
    return {"dirs": list(dirs), "sha256": _sha256(json.dumps(files, sort_keys=True).encode()), "files": files}


def git_info() -> dict:
    def g(*a):
        r = subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    status = g("status", "--porcelain", "--", *TREE_DIRS, "tests") or ""
    return {"head": g("rev-parse", "HEAD"), "dirty": bool(status), "dirty_files": status.splitlines()[:50]}


def package_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version
    out = {"python": sys.version.split()[0]}
    for p in ("numpy", "scipy", "scikit-learn", "pandas", "pyarrow", "torch", "threadpoolctl"):
        try:
            out[p] = version(p)
        except PackageNotFoundError:
            out[p] = None
    return out


def hardware() -> dict:
    mem = None
    try:
        mem = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout) / 2**30
    except Exception:  # noqa: BLE001
        pass
    return {"machine": platform.machine(), "system": f"{platform.system()} {platform.release()}",
            "processor": platform.processor(), "cpu_count": os.cpu_count(), "mem_gb": mem}


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**30 if sys.platform == "darwin" else r / 2**20


def run_pytest(scope: str = "frozen") -> dict:
    """Run pytest in a subprocess (strict: skips counted as failures) and parse the summary."""
    if scope == "none":
        return {"scope": "none", "ran": False}
    target = {"frozen": ["tests/frozen"], "all": ["tests"]}[scope]
    t = time.time()
    env = dict(os.environ, YSL_STRICT="1")
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *target], cwd=ROOT,
                       capture_output=True, text=True, env=env)
    tail = (r.stdout.strip().splitlines() or [""])[-1]
    counts = {}
    import re
    for n, k in re.findall(r"(\d+) (passed|failed|skipped|errors?|xfailed|xpassed)", tail):
        counts[k.rstrip("s") if k.startswith("error") else k] = int(n)
    return {"scope": scope, "ran": True, "exit_code": r.returncode, "summary": tail, **counts,
            "ok": r.returncode == 0 and counts.get("skipped", 0) == 0, "duration_s": round(time.time() - t, 1)}


def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def set_thread_limits(n: int) -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(k, str(n))
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(n)
    except Exception:  # noqa: BLE001
        pass


# ----------------------------------------------------------------------------- data guard (method sandbox)
_GUARD = {"active": 0, "hook": False}


def _protected(path) -> bool:
    """True for Costanzo GI data: data/raw/costanzo2016/**, data/processed/costanzo2016_*, saved predictions."""
    try:
        p = os.path.realpath(os.fsdecode(path))
    except (TypeError, ValueError):
        return False
    raw = os.path.realpath(ROOT / "data" / "raw" / "costanzo2016")
    proc = os.path.realpath(ROOT / "data" / "processed")
    preds = os.path.realpath(PRED_DIR)
    upreds = os.path.realpath(UPSTREAM_PRED_DIR)  # humanSLmapping PATCH: keep upstream predictions protected too
    if p == raw or p.startswith(raw + os.sep) or p == preds or p.startswith(preds + os.sep) \
            or p == upreds or p.startswith(upreds + os.sep):
        return True
    d, name = os.path.split(p)
    return d == proc and name.startswith("costanzo2016_")


def _audit(event, args):
    if _GUARD["active"] and event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)) \
            and _protected(args[0]):
        raise PermissionError(f"data guard: methods may not open Costanzo data files ({os.fsdecode(args[0])})")


def _guarded(fn):
    def wrapper(*a, **k):
        src = a[0] if a else next((k[x] for x in ("source", "path", "where", "path_or_paths", "filepath_or_buffer")
                                   if x in k), None)
        cands = src if isinstance(src, (list, tuple)) else [src]
        for c in cands:
            if isinstance(c, (str, bytes, os.PathLike)) and _protected(c):
                raise PermissionError(f"data guard: methods may not read Costanzo data files ({os.fsdecode(c)})")
        return fn(*a, **k)
    wrapper.__wrapped__ = fn
    return wrapper


@contextlib.contextmanager
def data_guard():
    """While active: opening Costanzo raw/processed/bundle/prediction files raises PermissionError (Python `open`
    audit event — covers open/io.open/os.open/np.load/np.memmap/pandas text readers — plus wrappers around the
    pyarrow/pandas parquet readers, which open files in C++), and ysl.data.dense's bundle cache is hidden.
    External feature files (data/processed/features/**) and data/raw/<other sources> remain readable."""
    import pandas as pd
    import pyarrow as pa
    import pyarrow.dataset as pds
    import pyarrow.parquet as pq

    import ysl.data.dense as D
    if not _GUARD["hook"]:
        sys.addaudithook(_audit)
        _GUARD["hook"] = True
    targets = [(pq, "read_table"), (pq, "ParquetFile"), (pq, "ParquetDataset"), (pq, "read_metadata"),
               (pq, "read_schema"), (pds, "dataset"), (pa, "memory_map"), (pa, "OSFile"), (pa, "input_stream"),
               (pd, "read_parquet"), (pd, "read_feather")]
    saved = [(mod, name, getattr(mod, name)) for mod, name in targets if hasattr(mod, name)]
    cache = dict(D._CACHE)
    D._CACHE.clear()
    for mod, name, fn in saved:
        setattr(mod, name, _guarded(fn))
    _GUARD["active"] += 1
    try:
        yield
    finally:
        _GUARD["active"] -= 1
        for mod, name, fn in saved:
            setattr(mod, name, fn)
        D._CACHE.clear()
        D._CACHE.update(cache)


def guarded_call(fn, view, ctx):
    """Run a method's fit_predict inside the data guard."""
    with data_guard():
        return fn(view, ctx)


# ----------------------------------------------------------------------------- code snapshot
def code_snapshot(dest: Path) -> dict:
    """tar.gz of src/, scripts/, config/, tests/ (no data, no caches) so dirty-tree runs stay provable."""
    import tarfile
    files = []
    for d in SNAPSHOT_DIRS:
        for p in sorted((ROOT / d).rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc" and p.name != ".DS_Store" \
                    and p.stat().st_size < 5_000_000:
                files.append(p)
    with tarfile.open(dest, "w:gz") as tar:
        for p in files:
            ti = tar.gettarinfo(p, arcname=str(p.relative_to(ROOT)))
            ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname = 0, 0, 0, "", ""
            with open(p, "rb") as fh:
                tar.addfile(ti, fh)
    return {"path": os.path.relpath(dest, ROOT), "sha256": _sha256(dest.read_bytes()), "n_files": len(files),
            "bytes": dest.stat().st_size}


# ----------------------------------------------------------------------------- poisoning self-check
def poisoned_bundle(bundle, split, seed: int = 999):
    """Replace every held-out (non-training) query row by finite garbage and permute its observed pattern."""
    rng = np.random.default_rng(seed)
    bad = ~np.isin(bundle.strain_orf, split.train)
    mats = {}
    for k in ("E", "SP", "DMF", "DMF_sd"):
        M = np.array(bundle.__getattribute__(k), copy=True)
        M[bad] = rng.normal(-5, 50, size=M[bad].shape).astype(np.float32)
        mats[k] = M
    Y = np.array(bundle.Y, copy=True)
    Yb = rng.integers(-1, 2, size=Y[bad].shape).astype(np.int8)
    Y[bad] = Yb
    mats["Y"] = Y
    return bundle.with_matrices(**mats)


# ----------------------------------------------------------------------------- run
DEFAULT_CONFIG = {
    "model_seed": None,        # None -> config/seeds.yaml model_seeds[0]
    "params": {},              # overrides of the method's registered defaults
    "n_folds": 10,
    "n_boot": 1000,
    "boot_seed": 0,
    "n_threads": 4,
    "chunk_targets": 250,
    "save_predictions": True,
    "oracle": False,           # dev only: oracle diagnostic (true held-out query rate) on this method's scores
    "poison_check": False,     # re-run the first seed on a bundle with poisoned held-out rows; must be bit-identical
    "pytest": "frozen",        # "frozen" | "all" | "none" | a precomputed result dict
    "headline": False,         # requires a clean git tree and a passing pytest with 0 skips
    "notes": "",
    "code_snapshot": True,     # tar.gz of src/scripts/config/tests in results/<exp_id>/ (SHA-256 in the manifest)
    "pred_format": "snappy",   # "snappy" (pandas default) | "compact" (lossless zstd + byte-stream-split; ~40% smaller)
}
PRED_FORMATS = ("snappy", "compact")


def write_predictions(df, path: Path, fmt: str = "snappy") -> None:
    """Keyed predictions [query_strain_id, array_strain_id, score] -> parquet. Both formats are lossless and are read
    identically by pandas/pyarrow (Evaluator.scores_from_keyed, scripts/rescore.py); "compact" only changes encoding."""
    if fmt == "snappy":
        df.to_parquet(path, index=False)
    elif fmt == "compact":
        df.to_parquet(path, index=False, compression="zstd", compression_level=9,
                      use_dictionary=["query_strain_id", "array_strain_id"],
                      column_encoding={"score": "BYTE_STREAM_SPLIT"})
    else:
        raise ValueError(f"pred_format must be one of {PRED_FORMATS}")


def run(method_name: str, dataset: str = "s1", split_kind: str = "dev", seeds=None, config: dict | None = None,
        exp_id: str | None = None, allow_test_only: bool = False, train_fraction: float | None = None,
        _bundle=None, _splits: dict | None = None, _final_token: str | None = None, _recipe: dict | None = None,
        _final_scope: str | None = None) -> dict:
    """`_bundle` / `_splits` ({seed: Split}) inject toy data — tests only. Final splits require `_final_token` from
    `run_final_batch` (the token in results/final/LOCK.pending) or, with `_final_scope`, from `run_secondary_batch`
    (the token in results/final/<scope>/LOCK.pending). `train_fraction` (None = config/seeds.yaml) is passed to
    make_split; on final splits only a secondary batch may set it."""
    from ysl.data.costanzo import dataset_of
    from ysl.data.dense import load_bundle
    from ysl.eval.evaluator import Evaluator
    from ysl.eval.splits import load_seeds, make_split, n_train_orfs
    from ysl.models import registry
    from ysl.models.view import MethodContext, ScreenView

    t_all = time.time()
    cfg = dict(DEFAULT_CONFIG, **(config or {}))
    set_thread_limits(int(cfg["n_threads"]))
    ds = dataset_of(dataset) if _bundle is None else _bundle.dataset
    spec = registry.get(method_name)
    if spec.test_only and not allow_test_only:
        raise ValueError(f"{method_name} is a test-only mutant")
    seeds_cfg = load_seeds()
    if seeds is None:
        seeds = seeds_cfg["dev_split_seeds"][:3] if split_kind == "dev" else seeds_cfg["final_split_seeds"]
    seeds = [int(s) for s in seeds]
    model_seed = int(cfg["model_seed"] if cfg["model_seed"] is not None else seeds_cfg["model_seeds"][0])
    if model_seed not in seeds_cfg["model_seeds"] and _bundle is None:
        raise ValueError(f"model seed {model_seed} is not in config/seeds.yaml")
    params = dict(spec.defaults, **cfg["params"])
    if split_kind not in ("dev", "final"):
        raise ValueError(split_kind)
    if split_kind == "final":
        lock = _lock_path(_final_scope)
        if lock.exists():
            raise PermissionError(f"{lock} exists: the final evaluation is one-shot")
        if _final_scope is not None and not LOCK.exists():
            raise PermissionError(f"secondary final runs need the primary {LOCK}")
        if _final_scope is None and train_fraction is not None:
            raise PermissionError("the primary final uses the train fraction in config/seeds.yaml")
        if _final_scope is not None and train_fraction is None:
            raise PermissionError("a secondary final scope needs an explicit train_fraction")
        pending = _pending_path(_final_scope)
        tok = json.loads(pending.read_text()).get("token") if pending.exists() else None
        if _final_token is None or tok != _final_token:
            raise PermissionError("final splits run only inside run_final_batch (scripts/32_final.py) or "
                                  "run_secondary_batch (scripts/33_secondary.py)")
        if cfg["oracle"]:
            raise PermissionError("oracle diagnostics are dev-only")
    if cfg["pred_format"] not in PRED_FORMATS:
        raise ValueError(f"pred_format must be one of {PRED_FORMATS}")
    resolved = {"method": method_name, "dataset": ds, "split_kind": split_kind, "seeds": seeds,
                "model_seed": model_seed, "params": params,
                **{k: cfg[k] for k in ("n_folds", "n_boot", "boot_seed", "oracle", "poison_check")}}
    if train_fraction is not None:   # only when set, so default runs keep their config hash
        resolved["train_fraction"] = float(train_fraction)
    config_sha = _sha256(json.dumps(_jsonable(resolved), sort_keys=True).encode())
    exp_id = exp_id or f"{method_name}__{ds}_{split_kind}__{time.strftime('%Y%m%d-%H%M%S')}"
    out_dir = RESULTS / exp_id
    if out_dir.exists():
        raise FileExistsError(f"{out_dir} exists; results are append-only")
    out_dir.mkdir(parents=True)

    if isinstance(cfg["pytest"], dict):
        pyt = dict(cfg["pytest"], precomputed=True)
    else:
        pyt = run_pytest(cfg["pytest"])
    git = git_info()
    if cfg["headline"] and (git["dirty"] or not pyt.get("ok")):
        raise PermissionError("headline runs need a clean tree and a passing pytest with no skips")

    snapshot = code_snapshot(out_dir / "code_snapshot.tar.gz") if cfg["code_snapshot"] else None
    lockbox_ref = None
    if _bundle is None:
        from ysl.eval.splits import read_lockbox
        lockbox_ref = read_lockbox(ds)   # independent of split.lockbox: the committed list
    stages = {}
    status, error = "running", None
    test_fps, n_lb_test = {}, {}
    per_split, fingerprints, pred_files = {}, {}, {}
    t = time.time()
    bundle = load_bundle(ds, unblind=(split_kind == "final")) if _bundle is None else _bundle
    stages["load_bundle_s"] = round(time.time() - t, 2)
    try:
        for seed in seeds:
            st = {}
            t = time.time()
            split = make_split(ds, split_kind, seed, bundle.all_query_orfs, train_fraction=train_fraction) \
                if _splits is None else _splits[seed]
            if split.kind != split_kind:
                raise ValueError("split kind mismatch")
            if train_fraction is not None and len(split.train) != n_train_orfs(len(bundle.all_query_orfs),
                                                                               train_fraction):
                raise AssertionError(f"split has {len(split.train)} training ORFs, expected "
                                     f"round({train_fraction} x {len(bundle.all_query_orfs)})")
            view = ScreenView.from_split(bundle, split)
            if not np.array_equal(view.target_orfs, split.test):
                raise AssertionError("view target order != split.test order (prediction rows would misalign)")
            ev = Evaluator(bundle, split, n_boot=cfg["n_boot"], boot_seed=cfg["boot_seed"], strata=_bundle is None)
            st["split_view_evaluator_s"] = round(time.time() - t, 2)
            ctx = MethodContext(model_seed=model_seed, split_seed=seed, params=dict(params), n_folds=cfg["n_folds"],
                                n_threads=cfg["n_threads"], chunk_targets=cfg["chunk_targets"])
            t = time.time()
            S = guarded_call(spec.fn, view, ctx)
            st["fit_predict_s"] = round(time.time() - t, 2)
            if np.shape(S) != (view.n_targets, view.n_arrays):
                raise AssertionError(f"{method_name} returned shape {np.shape(S)}, expected "
                                     f"{(view.n_targets, view.n_arrays)}")
            t = time.time()
            m = ev.evaluate(S)
            st["evaluate_s"] = round(time.time() - t, 2)
            m["method_info"] = ctx.info
            m["split"] = {"fingerprint": split.fingerprint, "n_train_orfs": int(len(split.train)),
                          "n_test_orfs": int(len(split.test)), "n_lockbox_orfs": int(len(split.lockbox)),
                          "n_lockbox_in_test": split.n_lockbox_in_test}
            fingerprints[str(seed)] = split.fingerprint
            test_fps[str(seed)] = _sha256(",".join(split.test).encode())
            lb = split.lockbox if lockbox_ref is None else lockbox_ref
            n_lb_test[str(seed)] = int(np.isin(split.test, lb).sum())
            m["split"]["test_orfs_sha256"] = test_fps[str(seed)]
            m["split"]["n_lockbox_in_test_vs_committed_list"] = n_lb_test[str(seed)]
            if train_fraction is not None:
                m["split"]["train_fraction"] = float(train_fraction)
            if cfg["oracle"]:
                from ysl.models.baselines import array_rate
                t = time.time()
                m["oracle"] = ev.oracle_diagnostic(S, array_rate(view, params.get("rate_prior", 1.0)))
                m["oracle"]["method_logit_scale"] = spec.logit_scale
                st["oracle_s"] = round(time.time() - t, 2)
            if cfg["poison_check"] and seed == seeds[0]:
                t = time.time()
                pv = ScreenView.from_split(poisoned_bundle(bundle, split), split)
                Sp = guarded_call(spec.fn, pv, MethodContext(model_seed=model_seed, split_seed=seed, params=dict(params),
                                                             n_folds=cfg["n_folds"], n_threads=cfg["n_threads"],
                                                             chunk_targets=cfg["chunk_targets"]))
                m["poison_check"] = {"view_digest_identical": pv.digest() == view.digest(),
                                     "predictions_bit_identical": bool(np.array_equal(S, Sp))}
                st["poison_check_s"] = round(time.time() - t, 2)
                del pv, Sp
                if not m["poison_check"]["predictions_bit_identical"]:
                    raise AssertionError(f"POISON CHECK FAILED for {method_name}: held-out rows influence predictions")
            if cfg["save_predictions"]:
                t = time.time()
                PRED_DIR.joinpath(exp_id).mkdir(parents=True, exist_ok=True)
                pf = PRED_DIR / exp_id / f"{ds}_{split_kind}_{seed}.parquet"
                write_predictions(ev.keyed_predictions(S), pf, cfg["pred_format"])
                pred_files[str(seed)] = {"path": os.path.relpath(pf, ROOT), "sha256": _sha256(pf.read_bytes()),
                                         "bytes": pf.stat().st_size, "format": cfg["pred_format"]}
                st["save_predictions_s"] = round(time.time() - t, 2)
            per_split[str(seed)] = m
            stages[f"seed_{seed}"] = st
            print(f"[{exp_id}] seed {seed}: AUPRC {m['auprc']:.4f} (prev {m['prevalence']:.4f}) "
                  f"fit {st['fit_predict_s']}s eval {st['evaluate_s']}s", flush=True)
            del S, view, ev
        status = "complete"
    except Exception as e:  # noqa: BLE001 — recorded in the ledger, then re-raised
        status, error = "failed", f"{type(e).__name__}: {e}"
        raise
    finally:
        summary = summarise(per_split)
        metrics = {"exp_id": exp_id, "method": method_name, "dataset": ds, "split_kind": split_kind, "seeds": seeds,
                   "status": status, "summary": summary, "splits": per_split}
        mtext = json.dumps(_jsonable(metrics), indent=1)
        if status == "complete":
            _atomic_write(out_dir / "metrics.json", mtext)
        stages["total_s"] = round(time.time() - t_all, 2)
        manifest = {
            "exp_id": exp_id, "status": status, "error": error, "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "cmd": " ".join(sys.argv), "resolved_config": resolved, "config_sha256": config_sha, "notes": cfg["notes"],
            "git": git, "tree": tree_hash(), "code_snapshot": snapshot, "recipe": _recipe,
            "data": bundle.meta, "split_fingerprints": fingerprints, "test_orfs_sha256": test_fps,
            "predictions": pred_files, "packages": package_versions(), "hardware": hardware(),
            "threads": {"n_threads": cfg["n_threads"], "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
            "wall_clock_s": stages, "peak_rss_gb": round(peak_rss_gb(), 2),
            "peak_rss_note": "process-wide peak (includes earlier methods run in the same process)", "pytest": pyt,
            "metrics_sha256": _sha256(mtext.encode()) if status == "complete" else None,
        }
        _atomic_write(out_dir / "manifest.json", json.dumps(_jsonable(manifest), indent=1))
        ledger_line = {
            "timestamp": manifest["timestamp"], "exp_id": exp_id, "kind": "evaluation", "status": status,
            "method": method_name, "dataset": ds, "split_kind": split_kind, "seeds": seeds, "model_seed": model_seed,
            "config_sha256": config_sha, "tree_sha256": manifest["tree"]["sha256"], "git_head": git["head"],
            "dirty": git["dirty"], "split_fingerprints": fingerprints, "test_orfs_sha256": test_fps,
            "n_lockbox_in_test": int(sum(n_lb_test.values())), "n_lockbox_in_test_by_seed": n_lb_test,
            "code_snapshot_sha256": snapshot["sha256"] if snapshot else None,
            "auprc": {k: per_split[k]["auprc"] for k in per_split},
            "auprc_mean": summary.get("auprc", {}).get("mean"), "auprc_sd": summary.get("auprc", {}).get("sd"),
            "pytest_ok": pyt.get("ok"), "error": error,
        }
        append_ledger(ledger_line)
    return metrics


def append_ledger(line: dict) -> None:
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with open(LEDGER, "a", encoding="utf-8") as f:
        f.write(json.dumps(_jsonable(line), sort_keys=True) + "\n")


def summarise(per_split: dict) -> dict:
    """Mean and SD (ddof=1 when >1 split) over splits of the headline metrics, P@k and stratum AUPRC/lift."""
    if not per_split:
        return {}
    ms = list(per_split.values())

    def agg(vals):
        v = np.array([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
        return {"mean": float(v.mean()) if len(v) else None,
                "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0 if len(v) else None, "n": int(len(v))}
    out = {k: agg([m.get(k) for m in ms]) for k in SUMMARY_KEYS}
    out["precision_at_k"] = {k: agg([m["precision_at_k"].get(k) for m in ms]) for k in ms[0]["precision_at_k"]}
    if "strata" in ms[0]:
        out["strata"] = {s: {f: agg([m["strata"][s][f] for m in ms if s in m["strata"]])
                             for f in ("auprc", "prevalence", "lift")} for s in ms[0]["strata"]}
    if "oracle" in ms[0]:
        out["oracle"] = {k: agg([m["oracle"][k] for m in ms]) for k, v in ms[0]["oracle"].items()
                         if isinstance(v, float)}
    out["runtime_note"] = "per-stage wall clock in manifest.json"
    return out


# ----------------------------------------------------------------------------- one-shot final evaluation
def _final_dir(scope: str | None = None) -> Path:
    """results/final/ for the primary final; results/final/<scope>/ for a secondary one (derived from LOCK at call time
    so tests can redirect LOCK)."""
    return LOCK.parent if scope is None else LOCK.parent / scope


def _lock_path(scope: str | None = None) -> Path:
    return LOCK if scope is None else _final_dir(scope) / "LOCK"


def _pending_path(scope: str | None = None) -> Path:
    return _final_dir(scope) / "LOCK.pending"


def load_recipe(path, check_committed: bool = True) -> dict:
    """A frozen recipe: config/recipes/<name>.yaml with keys method, params (dict), model_seed (optional),
    dev_exp_id (optional, the dev run it was selected on), notes. Must be committed and unmodified."""
    import yaml
    p = Path(path).resolve()
    if p.parent != RECIPE_DIR.resolve():
        raise PermissionError(f"recipes must live in {os.path.relpath(RECIPE_DIR, ROOT)}/ (got {path})")
    if not p.exists():
        raise FileNotFoundError(p)
    if check_committed:
        rel = os.path.relpath(p, ROOT)
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=ROOT, capture_output=True)
        clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT)
        if tracked.returncode != 0 or clean.returncode != 0:
            raise PermissionError(f"recipe {rel} must be committed and unmodified before a final run")
    raw = p.read_bytes()
    r = yaml.safe_load(raw) or {}
    unknown = set(r) - {"method", "params", "model_seed", "dev_exp_id", "notes"}
    if "method" not in r or unknown:
        raise ValueError(f"bad recipe {p.name}: needs 'method'; unknown keys {sorted(unknown)}")
    return {"file": os.path.relpath(p, ROOT), "sha256": _sha256(raw), "method": r["method"],
            "params": dict(r.get("params") or {}), "model_seed": r.get("model_seed"),
            "dev_exp_id": r.get("dev_exp_id"), "notes": r.get("notes", "")}




# shared by the primary and the secondary final batches ------------------------------------------------------------
def _claim_attempt(pending: Path, resume_reason: str | None) -> int:
    """Attempt number of a final batch; an unfinished batch (LOCK.pending present) is resumed only with a reason."""
    if not pending.exists():
        return 1
    if not resume_reason:
        raise PermissionError(f"{pending} exists (an earlier final batch did not finish); rerun only after an "
                              "infrastructure failure, with resume_reason")
    return int(json.loads(pending.read_text()).get("attempt", 1)) + 1


def _refuse_existing_output(out: Path) -> None:
    """A confirmatory output, once written, is never overwritten: its numbers may already have been seen. If a batch
    failed after writing it (ledger/LOCK step), LOCK.pending remains and a resume is refused here; the PI resolves it
    by hand from the saved file (docs/harness_secondary.md)."""
    if out.exists():
        raise PermissionError(f"{out} already exists: this batch's confirmatory analysis was already computed and is "
                              "never overwritten (no resume after the analysis; PI resolves by hand)")


def _load_recipes(recipe_paths: list, check_committed: bool) -> tuple[list, list]:
    recipes = [load_recipe(p, check_committed=check_committed) for p in recipe_paths]
    methods = [r["method"] for r in recipes]
    if len(set(methods)) != len(methods):
        raise ValueError("one recipe per method")
    return recipes, methods


def _provenance(pytest_scope: str, skip: bool) -> tuple[dict, dict]:
    """pytest gate (skips = failures) + git state; a final batch refuses a dirty tree or a failing/skipping pytest."""
    pyt = {"scope": "none", "ran": False} if skip else run_pytest(pytest_scope)
    git = git_info()
    if not skip and (git["dirty"] or not pyt.get("ok")):
        raise PermissionError("the final batch needs a clean committed tree and a passing pytest with no skips")
    return pyt, git


def _run_recipes(recipes: list, ds: str, seeds: list, n_boot: int, pyt: dict, token: str, attempt: int,
                 exp_id_of, label: str, train_fraction: float | None = None, scope: str | None = None,
                 pred_format: str = "snappy", _bundle=None, _splits: dict | None = None) -> dict:
    """Run every frozen recipe on the final splits (inside the batch identified by `token`)."""
    exp_ids = {}
    for r in recipes:
        cfg = {"params": r["params"], "pytest": pyt, "n_boot": n_boot, "save_predictions": True,
               "notes": f"{label} batch {token[:8]} attempt {attempt}; recipe {r['file']}"}
        if pred_format != "snappy":
            cfg["pred_format"] = pred_format
        if r["model_seed"] is not None:
            cfg["model_seed"] = r["model_seed"]
        eid = exp_id_of(r["method"])
        run(r["method"], ds, "final", seeds=seeds, config=cfg, exp_id=eid, train_fraction=train_fraction,
            _bundle=_bundle, _splits=_splits, _final_token=token, _recipe=r, _final_scope=scope)
        exp_ids[r["method"]] = eid
    return exp_ids


def _confirm_from_saved(exp_ids: dict, ds: str, seeds: list, candidates: list, comparator: str, n_boot: int,
                        train_fraction: float | None = None, _bundle=None, _splits: dict | None = None):
    """The pre-registered confirmatory analysis from the SAVED keyed predictions (the files the independent rescorer
    reads; SHA-256 checked against each manifest). Returns (confirmatory dict, bundle, {seed: Split})."""
    import pandas as pd

    from ysl.data.dense import load_bundle
    from ysl.eval.evaluator import Evaluator, confirmatory_analysis
    from ysl.eval.splits import make_split
    bundle = load_bundle(ds, unblind=True) if _bundle is None else _bundle
    evs, scores, splits = [], {m: [] for m in exp_ids}, {}
    for seed in seeds:
        split = make_split(ds, "final", seed, bundle.all_query_orfs, train_fraction=train_fraction) \
            if _splits is None else _splits[seed]
        splits[seed] = split
        ev = Evaluator(bundle, split, n_boot=0, strata=False)
        evs.append(ev)
        for m, eid in exp_ids.items():
            man = json.loads((RESULTS / eid / "manifest.json").read_text())
            pf = ROOT / man["predictions"][str(seed)]["path"]
            if _sha256(pf.read_bytes()) != man["predictions"][str(seed)]["sha256"]:
                raise AssertionError(f"prediction file changed since it was written: {pf}")
            scores[m].append(ev.scores_from_keyed(pd.read_parquet(pf)))
    return confirmatory_analysis(evs, scores, candidates, comparator, n_boot=n_boot), bundle, splits


def run_final_batch(recipe_paths: list, candidates: list, comparator: str, dataset: str = "s1", seeds=None,
                    n_boot: int = 1000, pytest_scope: str = "all", resume_reason: str | None = None,
                    _bundle=None, _splits: dict | None = None, _skip_provenance_checks: bool = False) -> dict:
    """THE final evaluation (DESIGN §4, §9): every pre-declared method (frozen recipes) on the final splits, then the
    confirmatory analysis; LOCK is written once, after everything has run. LOCK records the recipes (SHA-256),
    candidates, comparator, seeds, dataset and n_boot that the secondary finals (run_secondary_batch) must reuse."""
    import uuid

    from ysl.data.costanzo import dataset_of
    from ysl.eval.splits import load_seeds

    if LOCK.exists():
        raise PermissionError(f"{LOCK} exists: the final evaluation is one-shot")
    pending = _pending_path()
    attempt = _claim_attempt(pending, resume_reason)
    recipes, methods = _load_recipes(recipe_paths, check_committed=not _skip_provenance_checks)
    if not 1 <= len(candidates) <= 3 or not set(candidates) <= set(methods) or comparator not in methods \
            or comparator in candidates:
        raise ValueError("1-3 candidates and a distinct comparator, all among the recipes")
    ds = dataset_of(dataset) if _bundle is None else _bundle.dataset
    out = LOCK.parent / f"confirmatory_{ds}.json"
    _refuse_existing_output(out)
    seeds = [int(s) for s in (seeds if seeds is not None else load_seeds()["final_split_seeds"])]
    pyt, git = _provenance(pytest_scope, _skip_provenance_checks)
    code = code_state()
    token = uuid.uuid4().hex
    pending.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    pending.write_text(json.dumps({"token": token, "attempt": attempt, "resume_reason": resume_reason,
                                   "methods": methods, "candidates": candidates, "comparator": comparator,
                                   "recipes": recipes, "seeds": seeds, "dataset": ds, "started": stamp}, indent=1))
    exp_ids = _run_recipes(recipes, ds, seeds, n_boot, pyt, token, attempt,
                           lambda m: f"final_{m}__{ds}__{stamp}_a{attempt}", "FINAL",
                           _bundle=_bundle, _splits=_splits)
    t = time.time()
    conf, _, _ = _confirm_from_saved(exp_ids, ds, seeds, candidates, comparator, n_boot,
                                     _bundle=_bundle, _splits=_splits)
    conf.update({"dataset": ds, "token": token, "attempt": attempt, "exp_ids": exp_ids, "recipes": recipes,
                 "git": git, "tree": tree_hash(), "pytest": pyt, "runtime_s": round(time.time() - t, 1)})
    _refuse_existing_output(out)
    _atomic_write(out, json.dumps(_jsonable(conf), indent=1))
    conf_sha = _sha256(out.read_bytes())
    prim = conf["primary_lockbox_rows"]
    append_ledger({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "kind": "final_confirmatory", "dataset": ds,
                   "split_kind": "final", "exp_ids": exp_ids, "candidates": candidates, "comparator": comparator,
                   "primary_auprc_mean": {m: prim["methods"][m]["auprc_mean"] for m in methods},
                   "primary_p_holm": {c: prim["vs_comparator"][c]["p_holm"] for c in candidates},
                   "confirmatory_sha256": conf_sha, "token": token})
    LOCK.write_text(json.dumps(_jsonable({
        "token": token, "exp_ids": exp_ids, "confirmatory": os.path.relpath(out, ROOT),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "confirmatory_sha256": conf_sha, "dataset": ds,
        "seeds": seeds, "candidates": list(candidates), "comparator": comparator, "n_boot": int(n_boot),
        "recipes": recipes, "code_state": code}), indent=1))
    pending.unlink()
    return conf


# ----------------------------------------------------------------------------- secondary finals (DESIGN §4)
SECONDARY_COMBOS = (("cellmap", 0.20), ("s1", 0.10), ("s1", 0.30))   # pre-registered (dataset, train fraction)
REAL_DATASETS = ("s1", "cellmap")
KEEP_PREDICTIONS = ("all", "lockbox")
PRIMARY_LOCK_KEYS = ("recipes", "candidates", "comparator", "seeds", "dataset", "n_boot")


def canonical_fraction(train_fraction) -> float:
    f = float(train_fraction)
    r = round(f, 2)
    if not 0 < r < 1 or abs(f - r) > 1e-9:
        raise ValueError(f"train_fraction must lie in (0, 1) with at most two decimals (got {train_fraction})")
    return r


def secondary_scope(dataset: str, train_fraction) -> str:
    return f"secondary_{dataset}_f{canonical_fraction(train_fraction):.2f}"


def read_primary_lock() -> dict:
    """The protocol fixed by the primary final: recipes (file + SHA-256), candidates, comparator, seeds, dataset,
    n_boot. Taken from results/final/LOCK; for each key also present in the primary confirmatory JSON the two must
    agree (keys absent from LOCK fall back to that JSON). The JSON's SHA-256 is checked when LOCK records it."""
    if not LOCK.exists():
        raise PermissionError(f"{LOCK} is missing: secondary finals run only after the primary final "
                              "(scripts/32_final.py)")
    lock = json.loads(LOCK.read_text())
    cpath = (ROOT / lock["confirmatory"]).resolve()
    raw = cpath.read_bytes()
    if lock.get("confirmatory_sha256") and _sha256(raw) != lock["confirmatory_sha256"]:
        raise AssertionError(f"{cpath} does not match the SHA-256 recorded in {LOCK}")
    conf = json.loads(raw)
    if conf.get("token") != lock.get("token"):
        raise AssertionError(f"{cpath} is not the confirmatory output of the batch recorded in {LOCK}")
    out = {"token": lock["token"], "lock": lock, "confirmatory": os.path.relpath(cpath, ROOT),
           "confirmatory_sha256": _sha256(raw), "code_state": lock.get("code_state"), "git": conf.get("git")}
    for k in PRIMARY_LOCK_KEYS:
        a, b = lock.get(k), conf.get(k)
        if a is not None and b is not None and a != b:
            raise AssertionError(f"primary LOCK and confirmatory JSON disagree on {k!r}")
        if a is None and b is None:
            raise ValueError(f"primary LOCK records no {k!r}")
        out[k] = a if a is not None else b
    return out


def _prune_to_lockbox(exp_ids: dict, seeds: list, splits: dict, bundle, fmt: str) -> dict:
    """keep_predictions="lockbox": after the confirmatory analysis has been computed from the full files, rewrite each
    saved prediction file to its lockbox-ORF rows only (what `scripts/rescore.py --orf-list` needs to re-verify the
    primary-style endpoint). The pre-prune SHA-256 stays in each manifest; both hashes go into the secondary LOCK."""
    import pandas as pd
    rec: dict = {}
    for m, eid in exp_ids.items():
        man = json.loads((RESULTS / eid / "manifest.json").read_text())
        for seed in seeds:
            meta = man["predictions"][str(seed)]
            pf = ROOT / meta["path"]
            raw = pf.read_bytes()
            if _sha256(raw) != meta["sha256"]:
                raise AssertionError(f"prediction file changed since it was written: {pf}")
            df = pd.read_parquet(pf)
            lb_strains = bundle.strain_ids[np.isin(bundle.strain_orf, splits[seed].lockbox)]
            keep = df["query_strain_id"].astype(str).isin(set(lb_strains.tolist())).to_numpy()
            tmp = pf.with_name(pf.name + ".tmp")
            write_predictions(df[keep].reset_index(drop=True), tmp, fmt)
            os.replace(tmp, pf)
            rec.setdefault(eid, {})[str(seed)] = {
                "path": meta["path"], "sha256_full": meta["sha256"], "n_rows_full": int(len(df)),
                "bytes_full": len(raw), "sha256_lockbox_only": _sha256(pf.read_bytes()),
                "n_rows_lockbox_only": int(keep.sum()), "bytes_lockbox_only": pf.stat().st_size}
    return rec


def _compare_code(primary: dict | None, current: dict, reason: str | None) -> dict:
    """Secondary finals must run the primary's code (src/ysl/, config/). Any difference (or a primary LOCK that recorded
    no code state) is refused unless a reason is given; the reason and the changed files are returned for LOCK."""
    if primary is None:
        changed, note = None, "primary LOCK records no code_state; code identity cannot be verified"
    else:
        a, b = primary.get("files", {}), current["files"]
        changed = sorted(f for f in set(a) | set(b) if a.get(f) != b.get(f))
        note = None
    identical = changed == []
    if not identical and not (reason and reason.strip()):
        what = note or f"{len(changed)} file(s) differ from the primary final: {changed[:20]}"
        raise PermissionError(f"code differs from the primary final ({what}); rerun with code_change_reason "
                              "(--code-change-reason) to proceed, which is recorded in LOCK")
    return {"identical": identical, "primary_sha256": (primary or {}).get("sha256"), "current_sha256":
            current["sha256"], "changed_files": changed, "note": note,
            "code_change_reason": None if identical else reason}


def run_secondary_batch(recipe_paths: list, dataset: str, train_fraction: float, pytest_scope: str = "all",
                        resume_reason: str | None = None, predictions: str = "compact",
                        keep_predictions: str = "all", code_change_reason: str | None = None, _bundle=None,
                        _splits: dict | None = None, _skip_provenance_checks: bool = False) -> dict:
    """A pre-registered SECONDARY final (DESIGN §4): TheCellMap release at 20% or S1 at 10%/30% training queries.

    Requires the primary results/final/LOCK and reuses exactly its protocol: the same recipe files (SHA-256 must equal
    LOCK's; no re-tuning), seeds, candidates, comparator and n_boot, all read from LOCK. Same provenance checks as the
    primary (committed, unmodified recipes; clean tree; pytest with no skips). Same confirmatory analysis: this
    dataset's lockbox rows as the primary-style endpoint, all test rows as secondary. One-shot per (dataset, fraction):
    results/final/secondary_<ds>_f<frac>/{LOCK.pending, confirmatory.json, LOCK}; (s1, 0.20) is the primary and is
    refused. `predictions`: "compact" (lossless, ~40% smaller) or "snappy"; `keep_predictions="lockbox"` prunes the
    saved prediction files to lockbox rows after the analysis and after LOCK (see docs/harness_secondary.md).
    The code (src/ysl/, config/; per-file SHA-256) must equal the primary's `code_state` in LOCK, unless
    `code_change_reason` is given; the changed files and the reason are then recorded in this batch's LOCK."""
    import uuid

    from ysl.data.costanzo import dataset_of
    from ysl.eval.splits import load_seeds, n_train_orfs

    t_all = time.time()
    frac = canonical_fraction(train_fraction)
    ds = dataset_of(dataset) if _bundle is None else _bundle.dataset
    if _bundle is not None and dataset != ds:
        raise ValueError(f"dataset {dataset!r} != injected bundle dataset {ds!r}")
    if ds == "s1" and abs(frac - 0.20) < 1e-9:
        raise PermissionError("(s1, 0.20) is the primary final (scripts/32_final.py), not a secondary")
    if ds in REAL_DATASETS and (ds, frac) not in SECONDARY_COMBOS:
        raise PermissionError(f"({ds}, {frac:.2f}) is not a pre-registered secondary final {SECONDARY_COMBOS}")
    if predictions not in PRED_FORMATS or keep_predictions not in KEEP_PREDICTIONS:
        raise ValueError(f"predictions in {PRED_FORMATS}, keep_predictions in {KEEP_PREDICTIONS}")
    prim = read_primary_lock()
    if ds == prim["dataset"] and abs(frac - float(load_seeds()["train_query_fraction"])) < 1e-9:
        raise PermissionError(f"({ds}, {frac:.2f}) is the primary final")
    scope = secondary_scope(ds, frac)
    lock, pending = _lock_path(scope), _pending_path(scope)
    if lock.exists():
        raise PermissionError(f"{lock} exists: each secondary final is one-shot")
    out = _final_dir(scope) / "confirmatory.json"
    _refuse_existing_output(out)
    attempt = _claim_attempt(pending, resume_reason)
    recipes, methods = _load_recipes(recipe_paths, check_committed=not _skip_provenance_checks)
    want = {(r["file"], r["sha256"]) for r in prim["recipes"]}
    got = {(r["file"], r["sha256"]) for r in recipes}
    if got != want:
        raise PermissionError("secondary finals use exactly the primary's recipes (same files, same SHA-256 as LOCK); "
                              f"not in LOCK: {sorted(f for f, _ in got - want)}; "
                              f"missing: {sorted(f for f, _ in want - got)}")
    candidates, comparator = list(prim["candidates"]), prim["comparator"]
    seeds, n_boot = [int(s) for s in prim["seeds"]], int(prim["n_boot"])
    if ds in REAL_DATASETS and seeds != [int(s) for s in load_seeds()["final_split_seeds"]]:
        raise PermissionError("primary LOCK seeds are not the final split seeds of config/seeds.yaml")
    code = code_state()
    code_check = _compare_code(prim.get("code_state"), code, code_change_reason)
    pyt, git = _provenance(pytest_scope, _skip_provenance_checks)
    if code_state() != code:   # review N2: a change landing during the pytest gate must not slip past the code pin
        raise PermissionError("src/ysl or config changed while the pytest gate ran; rerun the secondary batch")
    token = uuid.uuid4().hex
    pending.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    pending.write_text(json.dumps({"token": token, "attempt": attempt, "resume_reason": resume_reason,
                                   "scope": scope, "dataset": ds, "train_fraction": frac,
                                   "primary_token": prim["token"], "methods": methods, "candidates": candidates,
                                   "comparator": comparator, "recipes": recipes, "seeds": seeds,
                                   "predictions": predictions, "keep_predictions": keep_predictions,
                                   "code_vs_primary": code_check, "started": stamp}, indent=1))
    exp_ids = _run_recipes(recipes, ds, seeds, n_boot, pyt, token, attempt,
                           lambda m: f"final2_{m}__{ds}_f{frac:.2f}__{stamp}_a{attempt}", f"SECONDARY {scope}",
                           train_fraction=frac, scope=scope, pred_format=predictions, _bundle=_bundle,
                           _splits=_splits)
    t = time.time()
    conf, bundle, splits = _confirm_from_saved(exp_ids, ds, seeds, candidates, comparator, n_boot,
                                               train_fraction=frac, _bundle=_bundle, _splits=_splits)
    n_all = len(bundle.all_query_orfs)
    n_train = [int(len(splits[s].train)) for s in seeds]
    if any(n != n_train_orfs(n_all, frac) for n in n_train):
        raise AssertionError("a final split does not have round(train_fraction x N) training ORFs")
    conf.update({"analysis": "secondary", "scope": scope, "dataset": ds, "train_fraction": frac,
                 "n_query_orfs": n_all, "n_train_orfs_per_split": n_train,
                 "lockbox_orfs_sha256": _sha256(",".join(splits[seeds[0]].lockbox).encode()),
                 "endpoint_note": "secondary final (DESIGN §4): 'primary_lockbox_rows' is the primary-style endpoint "
                                  "on this dataset's lockbox rows; the confirmatory claim is the primary final's",
                 "primary": {"token": prim["token"], "confirmatory": prim["confirmatory"],
                             "confirmatory_sha256": prim["confirmatory_sha256"]},
                 "token": token, "attempt": attempt, "exp_ids": exp_ids, "recipes": recipes, "git": git,
                 "tree": tree_hash(), "pytest": pyt, "predictions": predictions,
                 "keep_predictions": keep_predictions, "code_vs_primary": code_check,
                 "runtime_s": round(time.time() - t, 1),
                 "batch_wall_clock_s": round(time.time() - t_all, 1), "hardware": hardware()})
    _refuse_existing_output(out)
    _atomic_write(out, json.dumps(_jsonable(conf), indent=1))
    conf_sha = _sha256(out.read_bytes())
    prim_rows = conf["primary_lockbox_rows"]
    append_ledger({"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "kind": "secondary_confirmatory",
                   "scope": scope, "dataset": ds, "train_fraction": frac, "split_kind": "final", "exp_ids": exp_ids,
                   "candidates": candidates, "comparator": comparator,
                   "primary_auprc_mean": {m: prim_rows["methods"][m]["auprc_mean"] for m in methods},
                   "primary_p_holm": {c: prim_rows["vs_comparator"][c]["p_holm"] for c in candidates},
                   "confirmatory_sha256": conf_sha, "token": token, "primary_token": prim["token"]})
    prune_record = _final_dir(scope) / "pruned_predictions.json"
    lock.write_text(json.dumps(_jsonable({
        "token": token, "primary_token": prim["token"], "scope": scope, "dataset": ds, "train_fraction": frac,
        "exp_ids": exp_ids, "confirmatory": os.path.relpath(out, ROOT), "confirmatory_sha256": conf_sha,
        "seeds": seeds, "candidates": candidates, "comparator": comparator, "n_boot": n_boot, "recipes": recipes,
        "predictions": predictions, "keep_predictions": keep_predictions,
        "pruned_predictions": os.path.relpath(prune_record, ROOT) if keep_predictions == "lockbox" else None,
        "code_state": code, "code_vs_primary": code_check, "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}),
        indent=1))
    pending.unlink()
    # pruning only after the results are sealed by LOCK: a failure here can no longer lead to a re-run
    if keep_predictions == "lockbox":
        _atomic_write(prune_record, json.dumps({"status": "running"}))
        pruned = _prune_to_lockbox(exp_ids, seeds, splits, bundle, predictions)
        _atomic_write(prune_record, json.dumps(_jsonable({"status": "complete", "files": pruned}), indent=1))
    return conf
