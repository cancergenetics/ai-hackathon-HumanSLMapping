"""ScreenView: the ONLY data object a method receives (DESIGN.md §3, §11; AGENTS.md rule 5).

A view physically contains only what is allowed at prediction time:
  * training (screened) query rows, gene level (canonical strain), copied:   E, SP (signed -log10 P), Y (float32
    label, NaN = unobserved), DMF, DMF_sd                                    [n_train x n_arrays]
  * array metadata: array_strains, array_orfs, array_orfs_current, array_smf (gene SMF from the array strain)
  * targets (held-out query ORFs): their NAMES (target_orfs, target_orfs_current) and array column (target_col,
    -1 if the gene has no array strain) — nothing else. x's array-side column is read from the training rows at
    column target_col; x's SMF only via its array strain (`target_smf`, NaN if off-array, DESIGN §3.4).
It never contains: held-out query rows (ε, P, DMF, sd), query-strain metadata (strain IDs, suppressor flag,
temperature), per-query observed-array counts or which cells of a held-out row are observed, query SMF.

Construction: `ScreenView.from_split(bundle, split)` and `view.without_queries(fold_orfs)` both go through the
same private constructor `_construct(source, train_orfs, target_orfs)` — the cross-fit fold view is a view whose
training rows are T \\ fold and whose targets are the fold ORFs, built exactly like a test-time view.

Method interface (registered in ysl.models.registry):
    fit_predict(view: ScreenView, ctx: MethodContext) -> float32 [len(view.target_orfs) x view.n_arrays]
    (full grid; finite everywhere; rows in the order of view.target_orfs).
"""
from __future__ import annotations

import hashlib
import zlib
from dataclasses import dataclass, field

import numpy as np

_ROW_MATS = ("E", "SP", "Y", "DMF", "DMF_sd")


# ----------------------------------------------------------------------------- method context / seeding
@dataclass
class MethodContext:
    """Everything a method may know besides the view: seeds, hyper-parameters, resources."""
    model_seed: int
    split_seed: int = 0
    params: dict = field(default_factory=dict)
    n_folds: int = 10
    n_threads: int = 4
    chunk_targets: int = 250
    info: dict = field(default_factory=dict)   # diagnostics a method may report (coefficients, sizes); not inputs

    def rng(self, tag: str) -> np.random.Generator:
        """The one seeding helper: SeedSequence(model_seed, split_seed, crc32(tag))."""
        return np.random.default_rng(np.random.SeedSequence(
            [int(self.model_seed), int(self.split_seed), zlib.crc32(tag.encode())]))

    def folds(self, view: "ScreenView", k: int | None = None) -> np.ndarray:
        """Query-fold label (0..K-1) for each training row of `view` (balanced, seeded)."""
        k = self.n_folds if k is None else k
        return (self.rng("folds").permutation(view.n) % k).astype(np.int64)

    def get(self, key, default=None):
        return self.params.get(key, default)


# ----------------------------------------------------------------------------- row sources
class _BundleSource:
    """Gene-level rows of a Bundle (canonical strain per ORF). Used only inside _construct."""

    def __init__(self, bundle):
        self._b = bundle
        self.dataset = bundle.dataset
        self.array_strains = bundle.array_strains
        self.array_orfs = bundle.array_orfs
        self.array_orfs_current = bundle.array_orfs_current
        self.array_smf = bundle.array_smf
        self.orf_to_col = bundle.orf_to_col
        self.current = bundle.current_name

    def rows(self, orfs) -> dict:
        r = self._b.gene_rows(orfs)
        Y = self._b.Y[r].astype(np.float32)
        Y[Y < 0] = np.nan
        return {"E": self._b.E[r].copy(), "SP": self._b.SP[r].copy(), "Y": Y,
                "DMF": self._b.DMF[r].copy(), "DMF_sd": self._b.DMF_sd[r].copy()}


class _ViewSource:
    """Training rows of an existing ScreenView."""

    def __init__(self, v: "ScreenView"):
        self._v = v
        self.dataset = v.dataset
        self.array_strains, self.array_orfs = v.array_strains, v.array_orfs
        self.array_orfs_current, self.array_smf = v.array_orfs_current, v.array_smf
        self.orf_to_col = v._orf_to_col
        self.current = dict(zip(v.train_orfs, v.train_orfs_current)) | dict(zip(v.target_orfs, v.target_orfs_current))

    def rows(self, orfs) -> dict:
        orfs = np.asarray(orfs, dtype=str)
        idx = np.searchsorted(self._v.train_orfs, orfs)
        ok = (idx < self._v.n)
        ok[ok] = self._v.train_orfs[idx[ok]] == orfs[ok]
        if not ok.all():
            raise KeyError(f"ORFs are not training rows of the parent view: {orfs[~ok][:5].tolist()}")
        return {k: getattr(self._v, k)[idx].copy() for k in _ROW_MATS}


# ----------------------------------------------------------------------------- the view
class ScreenView:
    __slots__ = ("dataset", "train_orfs", "train_orfs_current", "train_col", "E", "SP", "Y", "DMF", "DMF_sd",
                 "array_strains", "array_orfs", "array_orfs_current", "array_smf", "target_orfs",
                 "target_orfs_current", "target_col", "_orf_to_col", "__weakref__")

    def __init__(self, *a, **k):
        raise TypeError("build views with ScreenView.from_split(bundle, split) or view.without_queries(orfs)")

    # ---- the single constructor
    @classmethod
    def _construct(cls, source, train_orfs, target_orfs) -> "ScreenView":
        train = np.unique(np.asarray(train_orfs, dtype=str))
        targets = np.unique(np.asarray(target_orfs, dtype=str))
        if len(train) != len(train_orfs) or len(targets) != len(target_orfs):
            raise ValueError("duplicate ORFs")
        if np.isin(train, targets).any():
            raise ValueError("an ORF is both a training row and a target")
        v = object.__new__(cls)
        v.dataset = source.dataset
        rows = source.rows(train)
        for k in _ROW_MATS:
            rows[k].setflags(write=False)
            object.__setattr__(v, k, rows[k])
        v.train_orfs = train
        v.target_orfs = targets
        v._orf_to_col = dict(source.orf_to_col)
        v.train_col = np.array([v._orf_to_col.get(o, -1) for o in train], dtype=np.int64)
        v.target_col = np.array([v._orf_to_col.get(o, -1) for o in targets], dtype=np.int64)
        v.train_orfs_current = np.array([source.current.get(o, o) for o in train], dtype=str)
        v.target_orfs_current = np.array([source.current.get(o, o) for o in targets], dtype=str)
        v.array_strains = source.array_strains.copy()
        v.array_orfs = source.array_orfs.copy()
        v.array_orfs_current = source.array_orfs_current.copy()
        v.array_smf = source.array_smf.copy()
        for k in ("train_orfs", "target_orfs", "train_col", "target_col", "array_smf"):
            getattr(v, k).setflags(write=False)
        return v

    @classmethod
    def from_split(cls, bundle, split) -> "ScreenView":
        """Training rows = split.train; targets = split.test (sorted). Lockbox ORFs only with an unblinded bundle."""
        if split.dataset and split.dataset != bundle.dataset:
            raise ValueError(f"split dataset {split.dataset} != bundle dataset {bundle.dataset}")
        lb_train = np.isin(split.train, split.lockbox).any()
        lb_test = np.isin(split.test, split.lockbox).any()
        if split.kind == "dev" and (lb_train or lb_test):
            raise PermissionError("dev split touches lockbox ORFs")
        if (lb_train or lb_test) and not bundle.unblinded:
            raise PermissionError("lockbox ORFs require an unblinded bundle (final evaluation only)")
        return cls._construct(_BundleSource(bundle), split.train, split.test)

    def without_queries(self, fold_orfs) -> "ScreenView":
        """Cross-fit fold view: training rows T \\ fold, targets = fold ORFs; same constructor as test time."""
        fold = np.asarray(fold_orfs, dtype=str)
        if not np.isin(fold, self.train_orfs).all():
            raise KeyError("fold ORFs must be training rows of this view")
        keep = self.train_orfs[~np.isin(self.train_orfs, fold)]
        return ScreenView._construct(_ViewSource(self), keep, fold)

    # ---- basic properties
    @property
    def n(self) -> int:
        return len(self.train_orfs)

    @property
    def n_arrays(self) -> int:
        return len(self.array_orfs)

    @property
    def n_targets(self) -> int:
        return len(self.target_orfs)

    @property
    def orfs(self) -> np.ndarray:  # pilot-compatible alias
        return self.train_orfs

    @property
    def col(self) -> np.ndarray:  # pilot-compatible alias
        return self.train_col

    @property
    def obs(self) -> np.ndarray:
        return ~np.isnan(self.Y)

    @property
    def arr_tr(self) -> np.ndarray:
        """For each array column, the training row of that array's ORF (-1 if not a training query)."""
        out = np.full(self.n_arrays, -1, dtype=np.int64)
        on = self.train_col >= 0
        out[self.train_col[on]] = np.flatnonzero(on)
        return out

    def cols_of(self, orfs) -> np.ndarray:
        return np.array([self._orf_to_col.get(o, -1) for o in np.asarray(orfs, dtype=str)], dtype=np.int64)

    @property
    def target_smf(self) -> np.ndarray:
        """Gene SMF of each target from its ARRAY strain; NaN for targets without an array strain."""
        out = np.full(self.n_targets, np.nan, dtype=np.float32)
        on = self.target_col >= 0
        out[on] = self.array_smf[self.target_col[on]]
        return out

    # ---- additive statistics and exact leave-one-query-out versions
    def array_counts(self) -> tuple[np.ndarray, np.ndarray]:
        """(positives, observed) per array column over all training rows."""
        return np.nansum(self.Y, axis=0).astype(np.float64), self.obs.sum(0).astype(np.float64)

    def loo_array_counts(self) -> tuple[np.ndarray, np.ndarray]:
        """Exact leave-one-query-out counts [n_train x n_arrays]: row t = counts over training rows other than t."""
        k, n = self.array_counts()
        return k[None, :] - np.nan_to_num(self.Y).astype(np.float64), n[None, :] - self.obs.astype(np.float64)

    def global_counts(self) -> tuple[float, float]:
        return float(np.nansum(self.Y)), float(self.obs.sum())

    def loo_global_counts(self) -> tuple[np.ndarray, np.ndarray]:
        K, N = self.global_counts()
        return K - np.nansum(self.Y, axis=1), N - self.obs.sum(1)

    # ---- provenance
    def digest(self) -> str:
        """SHA-256 over everything the view contains (used by leakage tests)."""
        h = hashlib.sha256()
        for k in ("train_orfs", "target_orfs", "train_col", "target_col", "array_strains", "array_orfs",
                  "array_smf", "train_orfs_current", "target_orfs_current", "array_orfs_current") + _ROW_MATS:
            a = np.ascontiguousarray(getattr(self, k))
            h.update(k.encode())
            h.update(a.tobytes() if a.dtype.kind != "U" else "\x00".join(a.tolist()).encode())
        h.update(repr(sorted(self._orf_to_col.items())).encode())
        return h.hexdigest()

    def __repr__(self) -> str:
        return (f"ScreenView({self.dataset}: {self.n} training queries, {self.n_targets} targets, "
                f"{self.n_arrays} arrays)")
