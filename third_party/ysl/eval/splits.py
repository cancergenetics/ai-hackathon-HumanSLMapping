"""Query-ORF splits and the lockbox (DESIGN.md §4). PI-owned; frozen by config/frozen.yaml.

The split unit is the query ORF (original Costanzo systematic name, the prefix of the strain ID): every query strain
of an ORF falls on the same side. All functions are pure functions of (sorted ORF list, seed, committed lockbox),
so a split is reproducible from `config/seeds.yaml` + `config/lockbox_<dataset>.txt` alone.

  lockbox_orfs(orfs, seed, fraction)
        first round(fraction * N) ORFs of numpy.random.default_rng(seed).permutation(sorted ORFs). With seed 777 over
        the 3,597 TheCellMap query ORFs this is the MASTER lockbox (= the pilot's lockbox; never touched by the pilot).
  derive_s1_lockbox(master, s1_orfs, cellmap_orfs)
        S1 lockbox = (master ∩ S1 query ORFs) ∪ lockbox_orfs(S1-only ORFs, 777, 0.20)   (r2-stats N1: an independent
        S1 draw would be contaminated by the pilot). S1-only = S1 query ORFs that are not TheCellMap query ORFs.
  read_lockbox(dataset)
        the committed list config/lockbox_<dataset>.txt, verified against its SHA-256 in config/frozen.yaml.
  dev_split(orfs, seed, lockbox)
        n_train = round(train_query_fraction * N_all) training ORFs drawn with default_rng(seed).choice(...,
        replace=False) from the sorted NON-lockbox ORFs; test = all remaining non-lockbox ORFs.
  final_split(orfs, seed, lockbox)
        n_train drawn the same way from ALL sorted ORFs; test = all others (lockbox ORFs included; reported as the
        lockbox stratum = the confirmatory endpoint).
  make_split(dataset, kind, seed)
        the harness entry point: ORF universe from the bundle sidecar's dataset, committed lockbox, seed checks.

n_train is computed from the TOTAL number of query ORFs of the dataset (S1: 3,283 -> 657; CellMap: 3,597 -> 719).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
CONFIG = ROOT / "config"
SEEDS_FILE = CONFIG / "seeds.yaml"
FROZEN_FILE = CONFIG / "frozen.yaml"
LOCKBOX_FILES = {"cellmap": CONFIG / "lockbox_cellmap.txt", "s1": CONFIG / "lockbox_s1.txt"}
KINDS = ("dev", "final")


def load_seeds(path: Path | str = SEEDS_FILE) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


@dataclass(frozen=True)
class Split:
    name: str            # "dev" | "final"
    seed: int
    train: np.ndarray    # sorted query ORFs used as screened (training) queries
    test: np.ndarray     # sorted query ORFs whose rows are predicted/evaluated
    lockbox: np.ndarray  # sorted lockbox ORFs
    dataset: str = ""
    n_all: int = 0       # number of query ORFs in the dataset
    meta: dict = field(default_factory=dict, compare=False)

    @property
    def kind(self) -> str:
        return self.name

    @property
    def fingerprint(self) -> str:
        """SHA-256 over dataset, kind, seed and the sorted train and test ORF lists."""
        s = f"{self.dataset}|{self.name}|{self.seed}|" + ",".join(self.train) + "|" + ",".join(self.test)
        return sha256_bytes(s.encode())

    @property
    def n_lockbox_in_test(self) -> int:
        return int(np.isin(self.test, self.lockbox).sum())


def _sorted_unique(orfs) -> np.ndarray:
    orfs = np.asarray(list(orfs), dtype=str)
    u = np.unique(orfs)
    if len(u) != len(orfs):
        raise ValueError("query ORF list contains duplicates")
    return u


def n_train_orfs(n_orfs: int, fraction: float) -> int:
    return int(round(fraction * n_orfs))


def lockbox_orfs(query_orfs, seed: int | None = None, fraction: float | None = None) -> np.ndarray:
    """Fixed draw: first round(fraction*N) ORFs of default_rng(seed).permutation(sorted ORFs)."""
    cfg = load_seeds() if seed is None or fraction is None else {}
    seed = cfg["lockbox_seed"] if seed is None else seed
    fraction = cfg["lockbox_fraction"] if fraction is None else fraction
    orfs = _sorted_unique(query_orfs)
    perm = np.random.default_rng(seed).permutation(len(orfs))
    return np.sort(orfs[perm[: n_train_orfs(len(orfs), fraction)]])


def derive_s1_lockbox(master, s1_orfs, cellmap_orfs, seed: int | None = None,
                      fraction: float | None = None) -> tuple[np.ndarray, dict]:
    """S1 lockbox = (master ∩ S1) ∪ seed-777 20% draw over S1-only ORFs. Returns (sorted ORFs, provenance)."""
    master, s1, cm = (_sorted_unique(x) for x in (master, s1_orfs, cellmap_orfs))
    if not np.isin(master, cm).all():
        raise ValueError("master lockbox contains ORFs that are not CellMap query ORFs")
    inter = master[np.isin(master, s1)]
    s1_only = s1[~np.isin(s1, cm)]
    extra = lockbox_orfs(s1_only, seed, fraction) if len(s1_only) else np.array([], dtype=str)
    lb = np.sort(np.concatenate([inter, extra]))
    return lb, {"n_master": int(len(master)), "n_master_in_s1": int(len(inter)), "n_s1_only": int(len(s1_only)),
                "n_s1_only_drawn": int(len(extra)), "s1_only_drawn": extra.tolist(), "n_lockbox": int(len(lb))}


def write_lockbox(path: Path, orfs, header: str) -> str:
    orfs = _sorted_unique(orfs)
    text = "".join(f"# {line}\n" for line in header.strip().splitlines()) + "".join(f"{o}\n" for o in orfs)
    Path(path).write_text(text, encoding="utf-8")
    return sha256_bytes(text.encode())


def frozen_hashes(path: Path | str = FROZEN_FILE) -> dict:
    if not Path(path).exists():
        raise FileNotFoundError(f"{path} missing (scripts/31_refreeze.py writes it; PI only)")
    with open(path, encoding="utf-8") as f:
        return (yaml.safe_load(f) or {}).get("files", {})


def read_lockbox(dataset: str, verify: bool = True) -> np.ndarray:
    """Committed lockbox ORF list; with verify=True its SHA-256 must equal the entry in config/frozen.yaml."""
    from ysl.data.costanzo import dataset_of
    path = LOCKBOX_FILES[dataset_of(dataset)]
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; created once by scripts/04b_make_lockbox.py")
    raw = path.read_bytes()
    if verify:
        rel = str(path.relative_to(ROOT))
        want = frozen_hashes().get(rel)
        if want is None:
            raise ValueError(f"{rel} has no SHA-256 in config/frozen.yaml")
        if sha256_bytes(raw) != want:
            raise ValueError(f"{rel} does not match its frozen SHA-256 — the lockbox must never change")
    orfs = [ln.strip() for ln in raw.decode().splitlines() if ln.strip() and not ln.startswith("#")]
    return _sorted_unique(orfs)


def _draw(pool: np.ndarray, n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if n > len(pool):
        raise ValueError(f"cannot draw {n} training ORFs from a pool of {len(pool)}")
    idx = np.random.default_rng(seed).choice(len(pool), size=n, replace=False)
    mask = np.zeros(len(pool), dtype=bool)
    mask[idx] = True
    return pool[mask], pool[~mask]


def _check(sp: Split) -> Split:
    tr, te, lb = set(sp.train), set(sp.test), set(sp.lockbox)
    if tr & te:
        raise AssertionError("train and test ORFs overlap")
    if sp.name == "dev" and (tr & lb or te & lb):
        raise AssertionError("a dev split touches a lockbox ORF")  # final splits draw from all ORFs (DESIGN §4)
    return sp


def dev_split(query_orfs, seed: int, train_fraction: float | None = None,
              lockbox: np.ndarray | None = None, dataset: str = "") -> Split:
    """Development split: train drawn from non-lockbox ORFs, test = remaining non-lockbox ORFs.
    `lockbox=None` (legacy/tests only) draws lockbox_orfs(query_orfs); the harness always passes the committed list."""
    train_fraction = load_seeds()["train_query_fraction"] if train_fraction is None else train_fraction
    orfs = _sorted_unique(query_orfs)
    lb = lockbox_orfs(orfs) if lockbox is None else _sorted_unique(lockbox)
    pool = orfs[~np.isin(orfs, lb)]
    train, test = _draw(pool, n_train_orfs(len(orfs), train_fraction), seed)
    return _check(Split("dev", int(seed), train, test, lb, dataset, len(orfs)))


def final_split(query_orfs, seed: int, train_fraction: float | None = None,
                lockbox: np.ndarray | None = None, dataset: str = "") -> Split:
    """Final split: train drawn from all ORFs, test = all other ORFs (lockbox included)."""
    train_fraction = load_seeds()["train_query_fraction"] if train_fraction is None else train_fraction
    orfs = _sorted_unique(query_orfs)
    lb = lockbox_orfs(orfs) if lockbox is None else _sorted_unique(lockbox)
    train, test = _draw(orfs, n_train_orfs(len(orfs), train_fraction), seed)
    return _check(Split("final", int(seed), train, test, lb, dataset, len(orfs)))


def make_split(dataset: str, kind: str, seed: int, query_orfs, train_fraction: float | None = None) -> Split:
    """Harness entry point. `query_orfs` = all query ORFs of the dataset (Bundle.all_query_orfs).
    Seeds must be listed in config/seeds.yaml for the kind; the lockbox is the committed, hash-verified list."""
    from ysl.data.costanzo import dataset_of
    ds = dataset_of(dataset)
    cfg = load_seeds()
    allowed = cfg["dev_split_seeds"] if kind == "dev" else cfg["final_split_seeds"] if kind == "final" else None
    if allowed is None:
        raise ValueError(f"unknown split kind {kind!r}")
    if int(seed) not in allowed:
        raise ValueError(f"seed {seed} is not a {kind} seed in config/seeds.yaml")
    lb = read_lockbox(ds)
    if not np.isin(lb, query_orfs).all():
        raise ValueError("lockbox contains ORFs that are not query ORFs of this dataset")
    fn = dev_split if kind == "dev" else final_split
    return fn(query_orfs, seed, train_fraction=train_fraction, lockbox=lb, dataset=ds)
