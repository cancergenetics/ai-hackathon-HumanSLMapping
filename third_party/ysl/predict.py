"""User-facing prediction tool (scripts/50_predict.py; docs/TOOL.md): score unscreened query genes against an array.

Use case: "I have screened these query genes (an SGA matrix: rows = screened queries, columns = array strains, eps
and P). Predict the stringent-negative GI / synthetic-lethal score for these other query genes against the array."

Two input modes, one prediction path:
  (a) Costanzo mode  `costanzo_view(dataset, screened, targets)`: the screened rows come from the (blinded) S1 or
      CellMap bundle. The ScreenView is built EXACTLY as the runner builds it for a split: a dev-kind `Split` with
      train = screened, test = targets and the committed lockbox, then `ScreenView.from_split(bundle, split)`.
      With screened/targets = the train/test ORFs of a harness split the view (digest) and the scores are those of
      the harness run (tests/test_predict.py).
  (b) User-matrix mode `user_view(matrix, targets, array_smf)`: long-format rows (query_orf, array_orf, epsilon,
      pvalue) -> gene-level dense rows in a `_UserSource`, passed to the same private ScreenView constructor that
      `from_split` uses. No Costanzo file is read. Arrays = array ORFs of the kept (screened) rows, sorted.
Then `predict(view, method, ...)` runs the registered method's fit_predict inside the runner's data guard
(`ysl.runner.guarded_call`: Costanzo raw/processed/bundle/prediction files cannot be opened while the method runs)
and, in user mode, an extra guard that blocks re-opening the user's own input files.

Held-out rows (DESIGN.md §3, AGENTS.md rule 5): a target's own rows are never read by the method.
  * Costanzo mode: the view physically contains only the screened rows; an ORF listed both as screened and as a
    target is REJECTED (two explicit lists contradict each other).
  * User-matrix mode: matrix rows whose query_orf is a target are IGNORED - dropped when the file is read, before the
    array set, SMF or anything else is derived (counted in the report; predictions are bit-identical to a matrix
    without those rows, tested). This allows "hide some of my own screens and see how well they are predicted".
Gene IDs are systematic ORF names; retired names are mapped to current SGD names for external-feature joins only
(`ysl.data.features.to_current`, when the SGD feature table is present). Methods that use bundled external features
refuse targets (and, for array-side feature methods, arrays) outside the external-feature gene universe unless
`allow_missing_features=True`.
"""
from __future__ import annotations

import contextlib
import dataclasses
import hashlib
import json
import os
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ysl.data.dense import signed_neglog10p
from ysl.eval import labels as L
from ysl.models import registry
from ysl.models.view import MethodContext, ScreenView

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
DEFAULT_METHOD = "stack_lr_wr2"
DEFAULT_THREADS = 3
MATRIX_COLUMNS = ("query_orf", "array_orf", "epsilon", "pvalue")
ORF_RE = r"^(Y[A-P][LR]\d{3}[WC](-[A-Z])?|Q\d{4}|R\d{4}[WC]?)$"   # as ysl.data.features.ORF_RE

# Which bundled external (non-GI) features a method's module reads (model code, read-only survey 2026-10-27):
#   "query": gene x gene features between screened queries and targets (pair store / gene table / co-fitness store)
#   "array": additionally target x array pair features
FEATURE_NEEDS = {
    "ysl.models.baselines": (),
    "ysl.models.completion": (),
    "ysl.models.foldin2": ("query",),
    "ysl.models.tier2": ("query",),
    "ysl.models.propensity": ("query",),
    "ysl.models.external": ("query", "array"),
    "ysl.models.combiner": ("query", "array"),
}


class InputError(ValueError):
    """Bad or inconsistent user input (raised with a message meant for the user)."""


# ----------------------------------------------------------------------------- small helpers
def _sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(1 << 24):
            h.update(b)
    return h.hexdigest()


def read_orf_list(path) -> np.ndarray:
    """One ORF per line (first whitespace/tab/comma-separated token); '#' comments and blank lines ignored."""
    out = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.split("#", 1)[0].strip()
        if ln:
            out.append(ln.replace(",", " ").split()[0].upper())
    if not out:
        raise InputError(f"{path}: no ORFs found")
    arr = np.asarray(out, dtype=str)
    u, c = np.unique(arr, return_counts=True)
    if (c > 1).any():
        warnings.warn(f"{path}: {int((c > 1).sum())} duplicated ORFs ignored (e.g. {u[c > 1][:3].tolist()})")
    return u


def _check_orf_names(names, what: str) -> None:
    names = pd.Index(np.asarray(names, dtype=str))
    bad = names[~names.str.match(ORF_RE)]
    if len(bad):
        raise InputError(f"{len(bad)} {what} are not systematic ORF names (e.g. {bad[:5].tolist()}). Use current "
                         "SGD systematic names such as YAL002W (strip strain suffixes like '_sn273'; no gene symbols).")


def current_names(orfs) -> tuple[np.ndarray, int, bool]:
    """(current SGD names, number remapped, mapping available). Identity when the SGD feature table is absent."""
    orfs = np.asarray(orfs, dtype=str)
    try:
        from ysl.data import features as F
        cur = F.to_current(pd.Index(orfs)).to_numpy(dtype=str)
        return cur, int((cur != orfs).sum()), True
    except (FileNotFoundError, OSError):
        return orfs.copy(), 0, False


def resolve_method(method: str | None = None, recipe=None, params: dict | None = None,
                   model_seed: int | None = None) -> dict:
    """Method name, parameters (registered defaults <- recipe params <- `params`) and model seed.
    Recipe format = config/recipes/README.md (keys method, params, model_seed, dev_exp_id, notes); the recipe's method
    is used unless `method` is given and differs (then an error). Seeds must be listed in config/seeds.yaml."""
    from ysl.eval.splits import load_seeds
    rec = {}
    if recipe is not None:
        import yaml
        raw = Path(recipe).read_bytes()
        rec = yaml.safe_load(raw) or {}
        unknown = set(rec) - {"method", "params", "model_seed", "dev_exp_id", "notes"}
        if "method" not in rec or unknown:
            raise InputError(f"bad recipe {recipe}: needs 'method'; unknown keys {sorted(unknown)}")
        if method is not None and method != rec["method"]:
            raise InputError(f"--method {method} contradicts recipe method {rec['method']}")
        rec["sha256"] = hashlib.sha256(raw).hexdigest()
        rec["file"] = str(recipe)
    name = rec.get("method", method or DEFAULT_METHOD)
    spec = registry.get(name)
    if spec.test_only:
        raise InputError(f"{name} is a test-only mutant")
    p = dict(spec.defaults) | dict(rec.get("params") or {}) | dict(params or {})
    seeds = load_seeds()
    ms = model_seed if model_seed is not None else rec.get("model_seed")
    ms = int(seeds["model_seeds"][0] if ms is None else ms)
    if ms not in seeds["model_seeds"]:
        raise InputError(f"model seed {ms} is not in config/seeds.yaml model_seeds {seeds['model_seeds']}")
    return {"method": name, "params": p, "model_seed": ms, "spec": spec,
            "recipe": {k: rec[k] for k in ("file", "sha256", "dev_exp_id", "notes") if k in rec} or None}


def feature_needs(spec) -> tuple[str, ...]:
    return FEATURE_NEEDS.get(spec.fn.__module__, ("query", "array"))   # unknown module: assume the strictest


def check_external_features(view: ScreenView, spec, allow_missing: bool = False) -> dict:
    """Refuse targets (and arrays, for methods with target x array features) that are absent from the bundled
    external-feature gene universe (data/processed/features/pairs/genes.txt). Screened queries outside it only warn.
    Returns coverage counts for the report."""
    needs = feature_needs(spec)
    if not needs:
        return {"needed": [], "note": f"{spec.name} uses no external features"}
    from ysl.features import pairs as PF
    try:
        universe = PF.genes()
    except FileNotFoundError as e:
        raise InputError(f"{spec.name} needs the bundled external features ({e}); build them with "
                         "scripts/10_download_features.py, 11_process_features.py, 12_esm_embeddings.py and "
                         "42_build_pair_features.py") from None
    out = {"needed": list(needs), "universe_size": int(len(universe))}
    miss_t = view.target_orfs[universe.get_indexer(view.target_orfs_current) < 0]
    miss_s = view.train_orfs[universe.get_indexer(view.train_orfs_current) < 0]
    miss_a = view.array_orfs[universe.get_indexer(view.array_orfs_current) < 0] if "array" in needs \
        else np.array([], dtype=str)
    out |= {"targets_missing": int(len(miss_t)), "screened_missing": int(len(miss_s)),
            "arrays_missing": int(len(miss_a))}
    problems = []
    if len(miss_t):
        problems.append(f"{len(miss_t)} target(s) (e.g. {miss_t[:5].tolist()})")
    if len(miss_a):
        problems.append(f"{len(miss_a)} array gene(s) (e.g. {miss_a[:5].tolist()})")
    if problems and not allow_missing:
        raise InputError(
            f"method {spec.name} uses bundled external features, but " + " and ".join(problems) + " have none (not "
            f"in the {len(universe)}-gene external-feature universe = current names of the Costanzo S1 query and array "
            "genes). Their scores would silently lose the external part of the model. Use current systematic ORF "
            "names, drop these genes, choose a method without external features (e.g. stack_lr_wr, stack_lr, "
            "B1_array_rate), or pass --allow-missing-features to score them anyway.")
    if problems:
        warnings.warn("scoring without external features for " + " and ".join(problems))
    if len(miss_s):
        warnings.warn(f"{len(miss_s)} screened queries have no external features (they still contribute GI rows)")
    return out


# ----------------------------------------------------------------------------- (a) Costanzo mode
def costanzo_view(dataset: str, screened, targets=None, split_seed: int = 0):
    """(view, bundle, split): screened rows from the blinded bundle; targets default to every unscreened query ORF
    with rows in the bundle (lockbox ORFs excluded). Names may be the bundle's ORF names or their current SGD names."""
    from ysl.data.dense import load_bundle
    from ysl.eval.splits import Split, read_lockbox
    bundle = load_bundle(dataset)          # blinded: lockbox query rows are not even loaded
    ds = bundle.dataset
    lockbox = read_lockbox(ds)
    inv = {c: o for o, c in zip(bundle.all_query_orfs, bundle.all_query_orfs_current)}
    known = set(bundle.all_query_orfs)

    def to_bundle_names(orfs):
        return np.unique(np.array([o if o in known else inv.get(o, o) for o in np.asarray(orfs, dtype=str)],
                                  dtype=str))
    scr = to_bundle_names(screened)
    unknown = scr[~np.isin(scr, bundle.all_query_orfs)]
    if len(unknown):
        raise InputError(f"{len(unknown)} screened ORFs are not query ORFs of the {ds} dataset "
                         f"(e.g. {unknown[:5].tolist()})")
    lb = scr[np.isin(scr, lockbox)]
    if len(lb):
        raise InputError(f"{len(lb)} screened ORFs are LOCKBOX queries of {ds} (e.g. {lb[:5].tolist()}); their rows "
                         "are reserved for the final evaluation and cannot be used")
    if targets is None:
        tg = bundle.query_orfs[~np.isin(bundle.query_orfs, scr)]
    else:
        tg = to_bundle_names(targets)
        both = tg[np.isin(tg, scr)]
        if len(both):
            raise InputError(f"{len(both)} ORFs are both screened and targets (e.g. {both[:5].tolist()}); a target's "
                             "own screen must not be used to predict it")
        lbt = tg[np.isin(tg, lockbox)]
        if len(lbt):
            raise InputError(f"{len(lbt)} targets are LOCKBOX queries of {ds} (e.g. {lbt[:5].tolist()}); predicting "
                             "them is reserved for the final evaluation (scripts/32_final.py)")
    if not len(tg):
        raise InputError("no targets")
    _check_orf_names(tg[~np.isin(tg, bundle.all_query_orfs) & ~np.isin(tg, bundle.array_orfs)], "targets")
    split = Split("dev", int(split_seed), scr, tg, lockbox, ds, len(bundle.all_query_orfs),
                  meta={"source": "ysl.predict"})
    view = ScreenView.from_split(bundle, split)      # exactly as the runner
    if not np.array_equal(view.target_orfs, split.test):
        raise AssertionError("view target order != split.test order")
    return view, bundle, split


# ----------------------------------------------------------------------------- (b) user-matrix mode
def read_user_matrix(path, targets=None) -> tuple[pd.DataFrame, dict]:
    """Long-format TSV/CSV (optionally .gz) with columns query_orf, array_orf, epsilon, pvalue and optionally
    array_strain (array strain ID; one strain per array ORF; sets the column key and order); other columns are
    ignored; '#' comment lines allowed. Rows of target query ORFs are DROPPED here, before anything else is derived.
    Returns (kept rows, report)."""
    p = Path(path)
    sep = "," if p.name.lower().removesuffix(".gz").endswith(".csv") else "\t"
    df = pd.read_csv(p, sep=sep, comment="#", dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    miss = [c for c in MATRIX_COLUMNS if c not in df.columns]
    if miss:
        raise InputError(f"{path}: missing column(s) {miss}; need {list(MATRIX_COLUMNS)} (tab-separated, header row)")
    has_strain = "array_strain" in df.columns
    df = df[list(MATRIX_COLUMNS) + (["array_strain"] if has_strain else [])].copy()
    rep = {"n_rows_read": int(len(df)), "array_strain_column": has_strain}
    if has_strain:
        df["array_strain"] = df["array_strain"].str.strip()
    for c in ("query_orf", "array_orf"):
        df[c] = df[c].str.strip().str.upper()
    tg = None if targets is None else set(np.asarray(targets, dtype=str))
    if tg is not None:
        # a target is matched by its given name or its current SGD name, on either side
        tg_all = tg | set(current_names(sorted(tg))[0])
        uq = df["query_orf"].unique()
        uq_cur = dict(zip(uq, current_names(uq)[0]))
        held = df["query_orf"].isin(tg_all) | df["query_orf"].map(uq_cur).isin(tg_all)
        rep["n_rows_dropped_target_queries"] = int(held.sum())
        rep["n_target_queries_with_rows_dropped"] = int(df.loc[held, "query_orf"].nunique())
        df = df[~held]
        if held.any():
            warnings.warn(f"ignored {int(held.sum())} matrix rows of {rep['n_target_queries_with_rows_dropped']} "
                          "target query ORF(s): a target's own screen is never used to predict it")
    for c in ("epsilon", "pvalue"):
        df[c] = pd.to_numeric(df[c].replace({"": np.nan, "NA": np.nan, "NaN": np.nan, "nan": np.nan}),
                              errors="raise").astype(np.float64)
    bad = df["epsilon"].isna() | df["pvalue"].isna()
    rep["n_rows_dropped_missing_values"] = int(bad.sum())
    if bad.any():
        warnings.warn(f"dropped {int(bad.sum())} rows with missing epsilon or pvalue (treated as unmeasured)")
    df = df[~bad]
    if not len(df):
        raise InputError(f"{path}: no usable (screened) rows")
    if ((df["pvalue"] < 0) | (df["pvalue"] > 1)).any():
        raise InputError(f"{path}: pvalue must lie in [0, 1]")
    if not np.isfinite(df["epsilon"].to_numpy()).all():
        raise InputError(f"{path}: non-finite epsilon")
    _check_orf_names(df["query_orf"].unique(), "query ORFs in the matrix")
    _check_orf_names(df["array_orf"].unique(), "array ORFs in the matrix")
    dup = df.duplicated(["query_orf", "array_orf"])
    if dup.any():
        ex = df.loc[dup, ["query_orf", "array_orf"]].head(3).values.tolist()
        raise InputError(f"{path}: {int(dup.sum())} duplicated (query_orf, array_orf) pairs (e.g. {ex}); give one "
                         "gene-level measurement per pair (e.g. pick one strain/allele per gene)")
    if has_strain:
        m = df.groupby("array_orf")["array_strain"].nunique()
        m2 = df.groupby("array_strain")["array_orf"].nunique()
        if (m > 1).any() or (m2 > 1).any():
            raise InputError(f"{path}: array_strain must map 1:1 to array_orf (one array strain per array gene)")
    rep["n_rows_used"] = int(len(df))
    return df.reset_index(drop=True), rep


def read_array_smf(path) -> pd.Series:
    """TSV/CSV with columns array_orf, smf (array single-mutant fitness)."""
    p = Path(path)
    sep = "," if p.name.lower().removesuffix(".gz").endswith(".csv") else "\t"
    df = pd.read_csv(p, sep=sep, comment="#", dtype={"array_orf": str})
    df.columns = [c.strip().lower() for c in df.columns]
    if not {"array_orf", "smf"} <= set(df.columns):
        raise InputError(f"{path}: need columns array_orf, smf")
    df["array_orf"] = df["array_orf"].str.strip().str.upper()
    if df["array_orf"].duplicated().any():
        raise InputError(f"{path}: duplicated array_orf")
    return pd.Series(pd.to_numeric(df["smf"], errors="raise").to_numpy(dtype=np.float64),
                     index=pd.Index(df["array_orf"].to_numpy(dtype=str)), name="smf")


class _UserSource:
    """Gene-level rows of a user SGA matrix, in the interface ScreenView._construct expects of a row source
    (cf. view._BundleSource). Holds ONLY screened rows (target rows were dropped by read_user_matrix)."""

    def __init__(self, df: pd.DataFrame, smf: pd.Series | None = None, extra_names=()):
        self.dataset = "user"
        q = np.unique(df["query_orf"].to_numpy(dtype=str))
        # columns: sorted array ORFs, or sorted array strain IDs when given (as the Costanzo bundle: column order
        # matters only through seeded randomised components, e.g. the 'tied' low-rank feature's random start)
        key = "array_strain" if "array_strain" in df.columns else "array_orf"
        keys, first = np.unique(df[key].to_numpy(dtype=str), return_index=True)
        a = df["array_orf"].to_numpy(dtype=str)[first]
        qi = np.searchsorted(q, df["query_orf"].to_numpy(dtype=str))
        ai = np.searchsorted(keys, df[key].to_numpy(dtype=str))
        eps, pv = df["epsilon"].to_numpy(np.float64), df["pvalue"].to_numpy(np.float64)
        shape = (len(q), len(a))
        self._E = np.full(shape, np.nan, np.float32)
        self._SP = self._E.copy()
        self._Y = np.full(shape, np.nan, np.float32)
        self._E[qi, ai] = eps
        self._SP[qi, ai] = signed_neglog10p(eps, pv)
        self._Y[qi, ai] = L.stringent_negative(eps, pv).astype(np.float32)   # the project label, from float64 input
        self._q = q
        self.array_orfs = a
        self.array_strains = keys               # array strain IDs if given, else the array ORFs
        self.array_smf = (np.full(len(a), np.nan, np.float32) if smf is None
                          else smf.reindex(a).to_numpy(dtype=np.float32))
        names = np.unique(np.concatenate([q, a, np.asarray(list(extra_names), dtype=str)]))
        cur, self.n_remapped, self.mapping_available = current_names(names)
        self.current = dict(zip(names, cur))
        self.array_orfs_current = np.array([self.current[o] for o in a], dtype=str)
        self.orf_to_col = {o: j for j, o in enumerate(a)}

    @property
    def query_orfs(self) -> np.ndarray:
        return self._q

    def rows(self, orfs) -> dict:
        orfs = np.asarray(orfs, dtype=str)
        r = np.searchsorted(self._q, orfs)
        ok = r < len(self._q)
        ok[ok] = self._q[r[ok]] == orfs[ok]
        if not ok.all():
            raise KeyError(f"ORFs without rows in the user matrix: {orfs[~ok][:5].tolist()}")
        nan = np.full((len(orfs), len(self.array_orfs)), np.nan, np.float32)
        return {"E": self._E[r].copy(), "SP": self._SP[r].copy(), "Y": self._Y[r].copy(), "DMF": nan,
                "DMF_sd": nan.copy()}


def user_view(matrix, targets=None, array_smf=None):
    """(view, report): screened = every query ORF with kept rows; targets default to every array ORF that is not a
    screened query (i.e. 'screen the rest of the array as queries')."""
    tg = None if targets is None else np.unique(np.char.upper(np.asarray(targets, dtype=str)))
    if tg is not None:
        _check_orf_names(tg, "targets")
    df, rep = read_user_matrix(matrix, tg)
    smf = read_array_smf(array_smf) if array_smf is not None else None
    if smf is None:
        warnings.warn("no --array-smf: array single-mutant fitness is unknown (NaN) for every array; methods that use "
                      "it see an 'unavailable' flag")
    if tg is None:
        a, q = np.unique(df["array_orf"].to_numpy(dtype=str)), np.unique(df["query_orf"].to_numpy(dtype=str))
        tg = a[~np.isin(a, q)]
    src = _UserSource(df, smf, extra_names=tg)
    if not len(tg):
        raise InputError("no targets")
    if len(src.query_orfs) < 2:
        raise InputError("need at least 2 screened query ORFs in the matrix")
    rep |= {"n_screened": int(len(src.query_orfs)), "n_arrays": int(len(src.array_orfs)),
            "n_targets": int(len(tg)), "n_targets_on_array": int(np.isin(tg, src.array_orfs).sum()),
            "n_names_remapped_to_current": src.n_remapped, "current_name_mapping": src.mapping_available,
            "array_smf_available": int(np.isfinite(src.array_smf).sum())}
    view = ScreenView._construct(src, src.query_orfs, tg)    # the constructor ScreenView.from_split uses
    return view, rep


# ----------------------------------------------------------------------------- guard for the user's own files
_UGUARD = {"paths": set(), "hook": False}


def _uaudit(event, args):
    if _UGUARD["paths"] and event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        try:
            p = os.path.realpath(os.fsdecode(args[0]))
        except (TypeError, ValueError):
            return
        if p in _UGUARD["paths"]:
            raise PermissionError(f"data guard: methods may not open the user's input matrix ({p})")


@contextlib.contextmanager
def user_file_guard(paths):
    """While active, opening any of `paths` from Python raises PermissionError (audit hook, as runner.data_guard)."""
    if not _UGUARD["hook"]:
        sys.addaudithook(_uaudit)
        _UGUARD["hook"] = True
    add = {os.path.realpath(str(p)) for p in paths if p is not None} - _UGUARD["paths"]
    _UGUARD["paths"] |= add
    try:
        yield
    finally:
        _UGUARD["paths"] -= add


# ----------------------------------------------------------------------------- prediction
@dataclasses.dataclass
class Prediction:
    S: np.ndarray            # float32 [n_targets x n_arrays], rows = view.target_orfs, columns = view.array_orfs
    view: ScreenView
    report: dict

    def table(self, top_k: int | None = None, drop_self: bool = True) -> pd.DataFrame:
        return score_table(self.S, self.view, logit_scale=self.report.get("logit_scale", False), top_k=top_k,
                           drop_self=drop_self)


def predict(view: ScreenView, method: str | None = None, recipe=None, params: dict | None = None,
            model_seed: int | None = None, split_seed: int = 0, n_threads: int = DEFAULT_THREADS,
            n_folds: int = 10, chunk_targets: int = 250, allow_missing_features: bool = False,
            protect=()) -> Prediction:
    """Run a registered method on `view` exactly as the runner does (MethodContext + guarded_call); `protect` = extra
    files the method may not open (the user's inputs)."""
    from ysl.eval.splits import load_seeds
    from ysl.runner import guarded_call, hardware, set_thread_limits
    t0 = time.time()
    set_thread_limits(int(n_threads))
    m = resolve_method(method, recipe, params, model_seed)
    seeds = load_seeds()
    if int(split_seed) not in seeds["dev_split_seeds"] + seeds["final_split_seeds"]:
        raise InputError(f"split seed {split_seed} is not in config/seeds.yaml")
    if not np.nansum(view.Y) > 0:
        raise InputError("the screened rows contain no stringent negative interaction (P < 0.05 and eps < -0.12): "
                         "there is nothing to learn from")
    feat = check_external_features(view, m["spec"], allow_missing_features)
    ctx = MethodContext(model_seed=m["model_seed"], split_seed=int(split_seed), params=dict(m["params"]),
                        n_folds=int(n_folds), n_threads=int(n_threads), chunk_targets=int(chunk_targets))
    t1 = time.time()
    with user_file_guard(protect):
        S = guarded_call(m["spec"].fn, view, ctx)
    t2 = time.time()
    if np.shape(S) != (view.n_targets, view.n_arrays):
        raise AssertionError(f"{m['method']} returned shape {np.shape(S)}, expected {(view.n_targets, view.n_arrays)}")
    S = np.asarray(S)
    if not np.isfinite(S).all():
        raise AssertionError(f"{m['method']} returned non-finite scores")
    report = {"method": m["method"], "params": m["params"], "model_seed": m["model_seed"], "split_seed": int(split_seed),
              "recipe": m["recipe"], "logit_scale": bool(m["spec"].logit_scale), "n_folds": int(n_folds),
              "dataset": view.dataset, "n_screened": int(view.n), "n_targets": int(view.n_targets),
              "n_targets_on_array": int((view.target_col >= 0).sum()), "n_arrays": int(view.n_arrays),
              "view_digest": view.digest(), "external_features": feat, "method_info": ctx.info,
              "runtime_s": {"setup": round(t1 - t0, 2), "fit_predict": round(t2 - t1, 2)},
              "threads": int(n_threads), "hardware": hardware(), "python": platform.python_version()}
    return Prediction(S=S, view=view, report=report)


def score_table(S: np.ndarray, view: ScreenView, logit_scale: bool = False, top_k: int | None = None,
                drop_self: bool = True) -> pd.DataFrame:
    """Long table [query_orf, array_orf, score, rank_within_query (+ probability for log-odds methods)].
    rank_within_query: 1 = highest score among the target's arrays (ties broken by array order); the self pair
    (target ORF = array ORF) is dropped when `drop_self`. top_k keeps the k best arrays per target."""
    S = np.asarray(S)
    n_t, n_a = S.shape
    Sm = S.astype(np.float64, copy=True)
    selfc = view.target_col
    if drop_self:
        on = np.flatnonzero(selfc >= 0)
        Sm[on, selfc[on]] = -np.inf
    order = np.argsort(-Sm, axis=1, kind="stable")
    rank = np.empty((n_t, n_a), dtype=np.int32)
    np.put_along_axis(rank, order, np.arange(1, n_a + 1, dtype=np.int32)[None, :].repeat(n_t, 0), axis=1)
    keep = np.ones((n_t, n_a), dtype=bool)
    if drop_self:
        on = np.flatnonzero(selfc >= 0)
        keep[on, selfc[on]] = False
    if top_k is not None:
        keep &= rank <= int(top_k)
    ti, ai = np.nonzero(keep)
    out = pd.DataFrame({
        "query_orf": pd.Categorical.from_codes(ti, categories=view.target_orfs),
        "array_orf": pd.Categorical.from_codes(ai, categories=view.array_orfs),
        "score": S[ti, ai],
        "rank_within_query": rank[ti, ai]})
    if not np.array_equal(view.array_strains, view.array_orfs):
        out.insert(2, "array_strain_id", pd.Categorical.from_codes(ai, categories=view.array_strains))
    if logit_scale:
        out["probability"] = (1.0 / (1.0 + np.exp(-S[ti, ai].astype(np.float64)))).astype(np.float32)
    return out


def write_table(df: pd.DataFrame, out) -> None:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".parquet":
        df.to_parquet(out, index=False)
    elif out.name.endswith((".tsv", ".tsv.gz", ".txt")):
        df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    else:
        raise InputError(f"output must end in .parquet, .tsv or .tsv.gz (got {out.name})")


def write_report(report: dict, path) -> None:
    def js(o):
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    Path(path).write_text(json.dumps(report, indent=1, default=js))
