"""The ONLY door to Costanzo labels/data in humanSLmapping (DESIGN v3 §2, §3, §6.1, §8). FROZEN after P0.

Data access (runtime-enforced by hsl.guard's label gate: Costanzo files can only be read inside this module, or
anywhere with HSL_FINAL=1)
  dev_bundle()            upstream Bundle.blind(T): the 656 test-query rows are physically dropped from every
                          row-aligned array; names (all_query_orfs) are kept. Asserted: gene counts; Y = −1 ⇔ ε/P
                          missing; the label rule reproduces Y on every dev row; no T row present.
  final_bundle()          the unblinded bundle; requires env HSL_FINAL=1 and a clean git tree. Only
                          scripts/run_final.py may call it (grep test).
  efixed_label_summary()  label-only counts/prevalence of the FINAL E_fixed and its strata (no scores; DESIGN §3).
Splits
  inner_folds()           5 folds of D, stratified by the frozen split strata (zero + 10 SL-count bins of
                          config/split/split_manifest.json), seeded by inner_fold_seed; genes that are one current SGD
                          gene share a fold.
  dev_split(k), resplit() val = fold k (k = 1, 2, 3) / val = folds 4+5 (re-split check, candidates = folds 1–3).
  final_split()           targets = T, pool = D minus identity removals.
  query_strata(orfs, b)   frozen split stratum (0..10) of each query from its SL count in bundle b (bootstrap strata).
Gene identity (current SGD names; DESIGN v3 §3): rule (i) a query pair split across pool/targets that is one current
gene -> pool member removed; (ii) arrays whose current name equals a pool gene's -> dropped from E_fixed(_dev);
(iii) cells whose array is the row gene under current names -> dropped. `identity_report()` lists all groups.
Evaluation sets (ESet: rows index split.targets, cols index bundle arrays; measured cells only)
  e_fixed(b, split)       final: T × (A \\ D) − (ii) − (iii); dev: D_val × (A \\ D_cand \\ T) − (ii) − (iii)
                          strata: AO (array-only columns), TT (final) / VV (dev) = on-array target-set columns
  e_unobs(b, split, S)    targets × (A \\ S) [dev: \\ T as well] − identity — per arm, descriptive
SelectorView               see class docstring.
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import hashlib
import json
import os
import subprocess
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from hsl import config, guard, paths

SPLIT_DIR = paths.PROJECT_ROOT / "config" / "split"
EXPECTED_COUNTS = {"n_query": 3283, "n_array": 3827, "n_both": 2652, "n_array_only": 1175, "n_query_only": 631}


# ============================================================================ names / lists
def _read_list(name: str) -> np.ndarray:
    return np.sort(np.asarray((SPLIT_DIR / name).read_text().split(), dtype=str))


@lru_cache(maxsize=1)
def test_queries() -> np.ndarray:
    """Names of the 656 test query genes (names only: used to mask/exclude, never to read rows)."""
    return _read_list("test_queries.txt")


@lru_cache(maxsize=1)
def dev_pool() -> np.ndarray:
    return _read_list("dev_pool.txt")


def _split_manifest() -> dict:
    return json.loads((SPLIT_DIR / "split_manifest.json").read_text())


def sha256_file(p) -> str:
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


# ============================================================================ bundles
def _assert_full_names(b) -> dict:
    q, a = set(b.all_query_orfs), set(b.array_orfs)
    got = {"n_query": len(q), "n_array": len(a), "n_both": len(q & a), "n_array_only": len(a - q),
           "n_query_only": len(q - a)}
    if got != EXPECTED_COUNTS:
        raise AssertionError(f"gene counts {got} != {EXPECTED_COUNTS}")
    T, D = set(test_queries()), set(dev_pool())
    if T & D or (T | D) != q:
        raise AssertionError("split lists are not a partition of the query genes")
    return got


def _assert_rows(b, label: str) -> dict:
    """Y == −1 ⇔ ε missing ⇔ P (SP) missing (also DMF); label rule reproduces Y (float32 P round trip tolerance)."""
    from ysl.eval.labels import P_THRESHOLD, stringent_negative
    miss = np.isnan(b.E)
    if not (np.array_equal(b.Y == -1, miss) and np.array_equal(miss, np.isnan(b.SP))):
        raise AssertionError(f"{label}: Y == -1 is not equivalent to missing eps/P")
    obs = ~miss
    P = np.where(obs, 10.0 ** (-np.abs(np.nan_to_num(b.SP.astype(np.float64)))), np.nan)
    lab = stringent_negative(np.nan_to_num(b.E), np.nan_to_num(P, nan=1.0)).astype(np.int8)
    bad = obs & (lab != b.Y)
    tol = float(config.protocol()["labels"]["p_round_trip_tolerance"])
    if bad.any() and not (np.abs(P[bad] - P_THRESHOLD) < tol).all():
        raise AssertionError(f"{label}: label rule does not reproduce Y away from the P threshold")
    return {"n_rows": int(b.E.shape[0]), "n_observed": int(obs.sum()), "n_pos": int((b.Y == 1).sum()),
            "label_mismatch_at_P_threshold": int(bad.sum())}


_DEV: dict = {}


def dev_bundle():
    """Upstream S1 bundle with the T rows physically removed (Bundle.blind(T)); cached per process."""
    if "b" in _DEV:
        return _DEV["b"]
    from ysl.data.dense import load_bundle
    T = test_queries()
    with guard.label_access():
        full = load_bundle(paths.costanzo_dataset(), unblind=True, cache=False)
        counts = _assert_full_names(full)
        d = full.blind(T)
        del full
    if np.isin(d.strain_orf, T).any() or len(d.query_orfs) != len(dev_pool()):
        raise AssertionError("dev bundle still holds test rows")
    d.meta = dict(d.meta, hsl_view="dev_bundle", gene_counts=counts, row_checks=_assert_rows(d, "dev rows"))
    d.freeze()
    _DEV["b"] = d
    return d


def git_state() -> dict:
    """HEAD, dirty flag and SHA-256 of `git diff HEAD` (tracked changes) of this repository."""
    def run(*a):
        return subprocess.run(["git", "-C", str(paths.PROJECT_ROOT), *a], capture_output=True, text=True).stdout
    status = run("status", "--porcelain")
    return {"head": run("rev-parse", "HEAD").strip(), "dirty": bool(status.strip()),
            "diff_sha256": hashlib.sha256(run("diff", "HEAD").encode()).hexdigest()}


def final_bundle():
    """Unblinded bundle (T rows present). FINAL EVALUATION ONLY: requires HSL_FINAL=1 and a clean git tree."""
    if os.environ.get("HSL_FINAL") != "1":
        raise PermissionError("final_bundle() requires HSL_FINAL=1 (scripts/run_final.py only)")
    if git_state()["dirty"]:
        raise PermissionError("final_bundle() requires a clean git tree")
    from ysl.data.dense import load_bundle
    with guard.label_access():
        b = load_bundle(paths.costanzo_dataset(), unblind=True, cache=False)
    _assert_full_names(b)
    b.meta = dict(b.meta, hsl_view="final_bundle")
    return b.freeze()


# ============================================================================ gene identity
def _current_map(b) -> dict:
    d = dict(zip(b.all_query_orfs, b.all_query_orfs_current))
    d.update(zip(b.array_orfs, b.array_orfs_current))
    return d


def current_of(b, orfs) -> np.ndarray:
    m = _current_map(b)
    return np.array([m.get(o, o) for o in np.asarray(orfs, dtype=str)], dtype=str)


def cols_matching(b, genes) -> np.ndarray:
    """Bool mask over array columns: array ORF in `genes` OR array current name = current name of a gene (rule ii)."""
    genes = np.asarray(genes, dtype=str)
    return np.isin(b.array_orfs, genes) | np.isin(b.array_orfs_current, current_of(b, genes))


def identity_report(b) -> dict:
    """Every current-SGD-name group with > 1 distinct original ORF over query ∪ array ORFs."""
    groups: dict = {}
    for o in b.all_query_orfs:
        groups.setdefault(_current_map(b)[o], {"query": set(), "array": set()})["query"].add(str(o))
    for o, c in zip(b.array_orfs, b.array_orfs_current):
        groups.setdefault(str(c), {"query": set(), "array": set()})["array"].add(str(o))
    multi = {c: {k: sorted(v) for k, v in g.items()} for c, g in sorted(groups.items())
             if len(g["query"] | g["array"]) > 1}
    return {"n_groups": len(multi), "groups": multi}


def identity_pool_removals(b, pool, others) -> np.ndarray:
    """Rule (i): pool genes that are the same current gene as a gene in `others` (targets / excluded)."""
    pool = np.asarray(pool, dtype=str)
    return pool[np.isin(current_of(b, pool), current_of(b, others))]


# ============================================================================ splits
def _strata_upper_edges() -> np.ndarray:
    """Upper SL-count edge of each nonzero stratum, from the frozen split manifest labels '[a, b]' / '(a, b]'."""
    labels = [s["label"] for s in _split_manifest()["strata"]]
    return np.array([float(lab.split(",")[1].strip(" ])")) for lab in labels[1:]])


def query_sl_counts(b, orfs) -> np.ndarray:
    return (b.Y[b.gene_rows(orfs)] == 1).sum(1).astype(np.int64)


def query_strata(orfs, b=None) -> np.ndarray:
    """Frozen split stratum (0 = zero SLs, 1..10 count bins) of each query ORF (rows must be in bundle b)."""
    b = dev_bundle() if b is None else b
    c = query_sl_counts(b, orfs)
    up = _strata_upper_edges()
    k = np.searchsorted(up, c, side="left") + 1
    return np.where(c == 0, 0, np.minimum(k, len(up))).astype(np.int64)


@lru_cache(maxsize=1)
def inner_folds() -> dict:
    """{dev ORF: fold 1..K}. Per stratum (0..10 in order): ORFs sorted, shuffled with ONE Generator
    (inner_fold_seed) consumed stratum by stratum, assigned round-robin continuing across strata (balanced folds).
    Then ORFs that are one current gene are moved to the fold of their alphabetically first member."""
    b = dev_bundle()
    K = int(config.protocol()["folds"]["n_inner_folds"])
    D = dev_pool()
    st = query_strata(D, b)
    rng = np.random.default_rng(config.seed("inner_fold_seed"))
    fold, off = {}, 0
    for h in range(int(st.max()) + 1):
        orfs = np.sort(D[st == h])
        for i, o in enumerate(orfs[rng.permutation(len(orfs))]):
            fold[str(o)] = (off + i) % K + 1
        off += len(orfs)
    cur = current_of(b, D)
    for c in np.unique(cur):
        members = np.sort(D[cur == c])
        for o in members[1:]:
            fold[str(o)] = fold[str(members[0])]
    return fold


@dataclass(frozen=True)
class Split:
    name: str                 # dev1 | dev2 | dev3 | resplit | final
    kind: str                 # dev | final
    targets: np.ndarray       # sorted ORFs whose rows are predicted / evaluated (D_val or T)
    pool: np.ndarray          # candidate pool after identity removals (sorted)
    excluded: np.ndarray      # genes masked in the SelectorView (T, plus D_val in dev)
    cand_full: np.ndarray     # pool before identity removals (D_cand, or D) — defines E_fixed columns
    pool_removed: np.ndarray = field(default_factory=lambda: np.array([], dtype=str))

    def digest(self) -> str:
        h = hashlib.sha256()
        for k in ("targets", "pool", "excluded"):
            h.update(("|".join(getattr(self, k))).encode() + b"#")
        return h.hexdigest()[:16]


def _make_split(name, kind, targets, cand, excluded) -> Split:
    b = dev_bundle()
    targets, cand, excluded = (np.sort(np.asarray(x, dtype=str)) for x in (targets, cand, excluded))
    rem = identity_pool_removals(b, cand, excluded)
    pool = np.setdiff1d(cand, rem)
    if len(identity_pool_removals(b, pool, excluded)):
        raise AssertionError("identity rule (i) not empty after application")
    if np.isin(pool, excluded).any() or np.isin(targets, pool).any():
        raise AssertionError("pool overlaps targets/excluded genes")
    return Split(name, kind, targets, pool, excluded, cand, np.sort(rem))


def dev_split(k: int) -> Split:
    if int(k) not in config.protocol()["folds"]["dev_splits"]:
        raise ValueError(f"dev split {k} not in protocol folds.dev_splits")
    f = inner_folds()
    D = dev_pool()
    val = np.array([o for o in D if f[o] == int(k)])
    return _make_split(f"dev{int(k)}", "dev", val, np.setdiff1d(D, val), np.union1d(val, test_queries()))


def resplit() -> Split:
    """Re-split check (DESIGN v3 §8): val = folds 4+5, candidates = folds 1–3."""
    vf = set(config.protocol()["folds"]["confirm_val_folds"])
    f = inner_folds()
    D = dev_pool()
    val = np.array([o for o in D if f[o] in vf])
    return _make_split("resplit", "dev", val, np.setdiff1d(D, val), np.union1d(val, test_queries()))


confirm_split = resplit


def final_split() -> Split:
    T = test_queries()
    return _make_split("final", "final", T, dev_pool(), T)


def get_split(name: str) -> Split:
    if name.startswith("dev") and name[3:].isdigit():
        return dev_split(int(name[3:]))
    if name in ("resplit", "confirm"):
        return resplit()
    if name == "final":
        return final_split()
    raise KeyError(name)


# ============================================================================ evaluation sets
@dataclass(frozen=True)
class ESet:
    name: str
    split: str
    targets: np.ndarray       # rows of the score matrix (= split.targets)
    rows: np.ndarray          # int32 index into targets
    cols: np.ndarray          # int32 array column
    y: np.ndarray             # int8 0/1
    strata: dict              # stratum name -> bool mask over cells
    meta: dict

    def scores(self, S: np.ndarray) -> np.ndarray:
        return np.asarray(S)[self.rows, self.cols]

    @property
    def n(self) -> int:
        return len(self.y)

    def digest(self) -> str:
        h = hashlib.sha256("|".join(self.targets).encode())
        h.update(self.rows.tobytes())
        h.update(self.cols.tobytes())
        return h.hexdigest()[:16]

    def subset(self, mask, name: str) -> "ESet":
        m = np.asarray(mask, dtype=bool)
        return ESet(name, self.split, self.targets, self.rows[m], self.cols[m], self.y[m],
                    {k: v[m] for k, v in self.strata.items()}, dict(self.meta, parent=self.name))

    def counts(self) -> dict:
        out = {"n_cells": int(self.n), "n_pos": int(self.y.sum()),
               "prevalence": float(self.y.mean()) if self.n else float("nan"),
               "n_rows": int(len(np.unique(self.rows))), "n_cols": int(len(np.unique(self.cols)))}
        for k, m in self.strata.items():
            out[k] = {"n_cells": int(m.sum()), "n_pos": int(self.y[m].sum()),
                      "prevalence": float(self.y[m].mean()) if m.any() else float("nan"),
                      "n_cols": int(len(np.unique(self.cols[m])))}
        pos_per_row = np.bincount(self.rows, weights=self.y, minlength=len(self.targets))
        out["n_rows_zero_pos"] = int((pos_per_row[np.unique(self.rows)] == 0).sum())
        return out


def _cells(b, split: Split, col_ok: np.ndarray, name: str, extra_meta=None) -> ESet:
    rows_b = b.gene_rows(split.targets)
    Y = b.Y[rows_b]
    row_cur = current_of(b, split.targets)
    not_self = row_cur[:, None] != b.array_orfs_current[None, :]          # rule (iii) incl. plain self pairs
    M = (Y >= 0) & col_ok[None, :] & not_self
    r, c = np.nonzero(M)
    y = Y[r, c].astype(np.int8)
    q_all = set(b.all_query_orfs)
    array_only = np.array([o not in q_all for o in b.array_orfs])
    on_target = np.isin(b.array_orfs, split.targets)
    strata = {"AO": array_only[c], ("VV" if split.kind == "dev" else "TT"): on_target[c]}
    # post-application assertions (rule iii; plain self)
    if (row_cur[r] == b.array_orfs_current[c]).any():
        raise AssertionError("identity rule (iii) not empty after application")
    meta = {"split": split.name, "split_digest": split.digest(), "n_targets": int(len(split.targets)),
            "excluded_cols_by_identity": [], **(extra_meta or {})}
    return ESet(name, split.name, split.targets, r.astype(np.int32), c.astype(np.int32), y, strata, meta)


def _fixed_col_ok(b, split: Split) -> tuple[np.ndarray, list]:
    block = split.cand_full if split.kind == "final" else np.union1d(split.cand_full, test_queries())
    ok = ~cols_matching(b, block)
    by_identity = sorted(set(b.array_orfs[~ok]) - set(block))
    if cols_matching(b, split.cand_full)[ok].any():
        raise AssertionError("identity rule (ii) not empty after application")
    return ok, by_identity


def e_fixed(b, split: Split) -> ESet:
    """Primary strategy-invariant set: final E_fixed (needs a bundle holding T rows) or dev E_fixed_dev."""
    if split.kind == "final" and not np.isin(split.targets, b.query_orfs).all():
        raise PermissionError("E_fixed (final) needs a bundle with the test rows (final_bundle)")
    ok, ident = _fixed_col_ok(b, split)
    es = _cells(b, split, ok, "E_fixed" if split.kind == "final" else "E_fixed_dev")
    es.meta["identity_excluded_arrays"] = ident
    if np.isin(b.array_orfs[es.cols], dev_pool() if split.kind == "final" else split.cand_full).any():
        raise AssertionError("E_fixed contains a pool array")
    if split.kind == "dev" and np.isin(b.array_orfs[es.cols], test_queries()).any():
        raise AssertionError("E_fixed_dev contains a test-gene column")
    return es


def e_unobs(b, split: Split, screened) -> ESet:
    """Practical per-arm set: targets × arrays not screened (identity-matched), dev: minus test columns."""
    S = np.asarray(screened, dtype=str)
    ok = ~cols_matching(b, S)
    if split.kind == "dev":
        ok &= ~cols_matching(b, test_queries())
    return _cells(b, split, ok, "E_unobs", {"n_screened": int(len(S))})


def eset_for(b, split: Split, name: str, screened=None) -> ESet:
    if name in ("E_fixed", "E_fixed_dev"):
        return e_fixed(b, split)
    if name == "E_unobs":
        return e_unobs(b, split, screened)
    raise KeyError(name)


def efixed_label_summary() -> dict:
    """Label-only counts/prevalence of the FINAL E_fixed and strata (DESIGN §3: written before any scoring).
    Loads the unblinded bundle inside the gate, returns aggregate counts only, drops the bundle."""
    from ysl.data.dense import load_bundle
    with guard.label_access():
        full = load_bundle(paths.costanzo_dataset(), unblind=True, cache=False)
    try:
        es = e_fixed(full, final_split())
        out = es.counts()
        out["identity_excluded_arrays"] = es.meta["identity_excluded_arrays"]
    finally:
        del full
    return out


# ============================================================================ SelectorView
def _ro(a):
    a = np.array(a, copy=True)
    a.setflags(write=False)
    return a


class SelectorView:
    """The ONLY object a selector receives (DESIGN v3 §6.1). Built by the harness from a dev/final bundle, a
    screened set S and a Split; holds no reference to the bundle. Exposes, by construction (__slots__):
      screened_orfs, screened_col           S (sorted) and each gene's array column (−1 off-array)
      E, SP  float32 [|S| x n_arrays]       ε and signed −log10 P of SCREENED rows only; columns of excluded genes
      Y      int8    [|S| x n_arrays]       (T; + D_val in dev; identity-matched by current name) masked with the
                                            native unobserved code (NaN / −1)
      array_orfs, array_orfs_current, array_smf (NaN for excluded genes), visible_cols (int index of unmasked cols)
      candidates, candidate_col             pool \\ S (sorted) and its array column (−1 off-array)
      gene_features(orfs), similarity(kind, a, b), clusters(orfs)   whitelisted external features (hsl.features) for
                                            visible genes only (candidates, screened genes, visible array genes)
      screen_view(targets=None)             upstream ScreenView over the masked rows (for ≤ 1 stack_lr fit per round)
    Never exposed: DMF/DMF_sd, query SMF, strain metadata, per-gene measurement counts, excluded genes' features,
    unscreened rows, split files. NaN SMF of off-array genes is inherent; selectors must not turn it (or
    candidate_col < 0) into a feature (off-array indicator; DESIGN §1) — only into model routing."""

    __slots__ = ("screened_orfs", "screened_col", "E", "SP", "Y", "array_orfs", "array_orfs_current", "array_smf",
                 "visible_cols", "candidates", "candidate_col", "_visible_genes", "_dataset", "_current")

    def __init__(self, *a, **k):
        raise TypeError("built by the harness: SelectorView.build(bundle, screened, split)")

    # ---- overridable pieces (broken twins in tests override exactly one)
    @staticmethod
    def _rows(b, screened):
        r = b.gene_rows(screened)
        return {"E": np.array(b.E[r], copy=True), "SP": np.array(b.SP[r], copy=True), "Y": np.array(b.Y[r], copy=True)}

    @staticmethod
    def _masked_cols(b, split: Split) -> np.ndarray:
        return cols_matching(b, split.excluded)

    @staticmethod
    def _smf(b, masked: np.ndarray) -> np.ndarray:
        smf = np.array(b.array_smf, dtype=np.float32, copy=True)
        smf[masked] = np.nan
        return smf

    @classmethod
    def build(cls, b, screened, split: Split) -> "SelectorView":
        S = np.sort(np.asarray(list(screened), dtype=str))
        if len(np.unique(S)) != len(S):
            raise ValueError("duplicate screened genes")
        if not np.isin(S, split.pool).all():
            raise ValueError("screened genes must come from the split's candidate pool")
        v = object.__new__(cls)
        rows = cls._rows(b, S)
        masked = cls._masked_cols(b, split)
        rows["E"][:, masked] = np.nan
        rows["SP"][:, masked] = np.nan
        rows["Y"][:, masked] = -1
        for k in ("E", "SP", "Y"):
            rows[k].setflags(write=False)
            object.__setattr__(v, k, rows[k])
        cmap = {o: j for j, o in enumerate(b.array_orfs)}
        v.screened_orfs = _ro(S)
        v.screened_col = _ro(np.array([cmap.get(o, -1) for o in S], dtype=np.int64))
        v.array_orfs = _ro(b.array_orfs)
        v.array_orfs_current = _ro(b.array_orfs_current)
        v.array_smf = _ro(cls._smf(b, masked))
        v.visible_cols = _ro(np.flatnonzero(~masked))
        cand = np.setdiff1d(split.pool, S)
        v.candidates = _ro(cand)
        v.candidate_col = _ro(np.array([cmap.get(o, -1) for o in cand], dtype=np.int64))
        v._visible_genes = frozenset(cand) | frozenset(S) | frozenset(b.array_orfs[~masked])
        v._dataset = b.dataset
        v._current = {o: c for o, c in _current_map(b).items() if o in v._visible_genes}
        return v

    # ---- properties / features
    @property
    def n_screened(self) -> int:
        return len(self.screened_orfs)

    @property
    def n_arrays(self) -> int:
        return len(self.array_orfs)

    def _check_visible(self, orfs) -> np.ndarray:
        orfs = np.asarray(orfs, dtype=str)
        bad = [o for o in orfs if o not in self._visible_genes]
        if bad:
            raise KeyError(f"genes not visible to selectors: {bad[:5]}")
        return orfs

    def gene_features(self, orfs):
        from hsl import features
        return features.gene_features(self._check_visible(orfs))

    def similarity(self, kind: str, orfs_a, orfs_b) -> np.ndarray:
        from hsl import features
        return features.similarity(kind, self._check_visible(orfs_a), self._check_visible(orfs_b))

    def clusters(self, orfs) -> np.ndarray:
        from hsl import features
        return features.clusters(self._check_visible(orfs))

    def screen_view(self, targets=None):
        """Upstream ScreenView whose training rows are the MASKED screened rows (DMF/DMF_sd set to NaN) and whose
        targets are `targets` (default: candidates)."""
        from ysl.models.view import ScreenView
        tg = self.candidates if targets is None else self._check_visible(targets)
        return ScreenView._construct(_MaskedSource(self), self.screened_orfs, np.asarray(tg, dtype=str))

    def __repr__(self):
        return f"SelectorView({self.n_screened} screened, {len(self.candidates)} candidates, " \
               f"{len(self.visible_cols)}/{self.n_arrays} visible columns)"


class _MaskedSource:
    """Row source for ScreenView._construct built from a SelectorView (no bundle access)."""

    def __init__(self, v: SelectorView):
        self._v = v
        self.dataset = v._dataset
        self.array_strains = np.array(v.array_orfs, copy=True)    # strain IDs not exposed; ORF names stand in
        self.array_orfs, self.array_orfs_current = v.array_orfs, v.array_orfs_current
        self.array_smf = v.array_smf
        self.orf_to_col = {o: j for j, o in enumerate(v.array_orfs)}
        self.current = dict(v._current)

    def rows(self, orfs) -> dict:
        idx = np.searchsorted(self._v.screened_orfs, np.asarray(orfs, dtype=str))
        Y = self._v.Y[idx].astype(np.float32)
        Y[Y < 0] = np.nan
        nan = np.full(Y.shape, np.nan, dtype=np.float32)
        return {"E": np.array(self._v.E[idx]), "SP": np.array(self._v.SP[idx]), "Y": Y, "DMF": nan,
                "DMF_sd": nan.copy()}


# ============================================================================ manifest
def eset_manifest(include_final: bool = True) -> dict:
    """Sizes/prevalence of every evaluation set (dev sets from dev rows; final E_fixed label counts only),
    identity groups and removals, fold sizes. Written by scripts/04_eset_manifest.py."""
    b = dev_bundle()
    out = {"gene_counts": b.meta["gene_counts"], "dev_row_checks": b.meta["row_checks"],
           "identity": identity_report(b), "splits": {}}
    folds = inner_folds()
    out["inner_folds"] = {"seed": config.seed("inner_fold_seed"),
                          "sizes": {str(k): int(sum(1 for v in folds.values() if v == k))
                                    for k in sorted(set(folds.values()))}}
    for sp in [dev_split(k) for k in config.protocol()["folds"]["dev_splits"]] + [resplit()]:
        es = e_fixed(b, sp)
        out["splits"][sp.name] = {"n_targets": int(len(sp.targets)), "n_pool": int(len(sp.pool)),
                                  "pool_identity_removed": sp.pool_removed.tolist(),
                                  "identity_excluded_arrays": es.meta["identity_excluded_arrays"],
                                  "E_fixed_dev": es.counts(), "split_digest": sp.digest()}
    fs = final_split()
    out["final"] = {"n_targets": int(len(fs.targets)), "n_pool": int(len(fs.pool)),
                    "pool_identity_removed": fs.pool_removed.tolist()}
    if include_final:
        out["final"]["E_fixed_label_counts"] = efixed_label_summary()
    return out
