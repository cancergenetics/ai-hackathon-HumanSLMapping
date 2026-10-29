"""Common selector machinery (DESIGN v3.1 §6.2): the Selector base class (harness protocol: `needs_view`, `static`,
`reset(candidates, rng)`, `select(view, k, rng) -> ordered genes`), the cluster-capped randomised greedy, knob loading
(config/strategies.yaml), and small label-free helpers (abundance, mid-rank percentiles, per-selector feature cache).

Selectors see ONLY the SelectorView (hsl.views.SelectorView) and use only its documented attributes:
  screened_orfs, screened_col, E, Y, array_smf, visible_cols, candidates, candidate_col, n_screened,
  gene_features(), similarity(), clusters().
Off-array status (candidate_col < 0, NaN SMF) is used only for ROUTING (which submodel scores a gene), never as a
feature or score term (DESIGN §1).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard)

import hashlib
import time
from functools import lru_cache

import numpy as np
import yaml

from hsl import config, harness, paths


# ============================================================================ knobs
@lru_cache(maxsize=1)
def knobs() -> dict:
    with open(paths.PROJECT_ROOT / "config" / "strategies.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def knobs_sha() -> str:
    import json
    return hashlib.sha256(json.dumps(knobs(), sort_keys=True).encode()).hexdigest()[:12]


def protocol_selector() -> dict:
    return config.protocol()["selector"]


def cache_token(view) -> str:
    """Identifies the external-feature source of a view: all real SelectorViews share the S1 tier-1 stores ("S1");
    toy views (tests) carry their own token so module-level caches never mix sources."""
    return str(getattr(view, "cache_token", "S1"))


def genes_key(view, genes, *extra) -> str:
    h = hashlib.sha256(cache_token(view).encode())
    h.update("|".join(np.asarray(genes, dtype=str)).encode())
    for e in extra:
        h.update(repr(e).encode())
    return h.hexdigest()[:24]


# ============================================================================ helpers
def midrank_pct(values: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Mid-rank empirical CDF of `values` within the finite entries of `ref` (NaN values -> NaN)."""
    r = np.sort(np.asarray(ref, dtype=np.float64)[np.isfinite(ref)])
    v = np.asarray(values, dtype=np.float64)
    if not len(r):
        return np.full(len(v), 0.5)
    lo = np.searchsorted(r, v, side="left")
    hi = np.searchsorted(r, v, side="right")
    out = (lo + hi) / (2.0 * len(r))
    return np.where(np.isfinite(v), out, np.nan)


def pool_of(view) -> np.ndarray:
    """Candidate pool of the trajectory = candidates ∪ screened (constant over a trajectory)."""
    return np.union1d(np.asarray(view.candidates, dtype=str), np.asarray(view.screened_orfs, dtype=str))


class FeatureCache:
    """Per-selector cache of whitelisted external gene features (label-free; keyed by gene list and source)."""

    def __init__(self):
        self._d: dict = {}

    def get(self, view, genes):
        genes = np.asarray(genes, dtype=str)
        key = genes_key(view, genes)
        if key not in self._d:
            if len(self._d) > 8:
                self._d.clear()
            self._d[key] = view.gene_features(genes)
        return self._d[key]


def abundance(feats, cfg: dict | None = None) -> np.ndarray:
    """Protein abundance (PaxDb Mueller 2020 log10 ppm); missing values filled from the Kulak 2014 column through a
    linear map fitted on genes measured in both (no missingness indicator is created; remaining NaN stay NaN)."""
    f = (cfg or knobs()["degree_model"])["features"]
    a = feats[f["abundance"]].to_numpy(np.float64)
    b = feats[f["abundance_fallback"]].to_numpy(np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    if both.sum() >= 10 and np.std(b[both]) > 0:
        beta = np.polyfit(b[both], a[both], 1)
        fill = ~np.isfinite(a) & np.isfinite(b)
        a = a.copy()
        a[fill] = np.polyval(beta, b[fill])
    return a


def candidate_smf(view) -> np.ndarray:
    """Array SMF of each candidate (NaN for off-array candidates: routing only)."""
    col = np.asarray(view.candidate_col)
    out = np.full(len(col), np.nan)
    on = col >= 0
    out[on] = np.asarray(view.array_smf, dtype=np.float64)[col[on]]
    return out


def capped_update(clusters: np.ndarray, cap: int, inner=None):
    """`update` callback for harness.randomised_greedy: optional inner gain refresh, then ≤ `cap` picks per cluster
    per batch (cluster −1 = unclustered, never capped). When every remaining candidate sits in a full cluster, a new
    cap cycle starts (counts reset) instead of leaving the batch to be filled at random."""
    clusters = np.asarray(clusters)
    counts: dict = {}
    st = {"g": None, "picked": [], "cycles": 0}

    def upd(g, i):
        if inner is not None:
            raw = np.array(inner(g, i), dtype=np.float64)
        else:                                   # no refresh: keep the unmasked first-pick gains ourselves
            if st["g"] is None:
                st["g"] = np.array(g, dtype=np.float64)
            raw = st["g"].copy()
        st["picked"].append(i)
        c = int(clusters[i])
        if c >= 0:
            counts[c] = counts.get(c, 0) + 1
        full = [cl for cl, n in counts.items() if n >= cap]
        if not full:
            return raw
        out = raw.copy()
        out[np.isin(clusters, full)] = -np.inf
        rest = np.ones(len(raw), dtype=bool)
        rest[st["picked"]] = False
        if not np.isfinite(out[rest]).any() and np.isfinite(raw[rest]).any():
            counts.clear()
            st["cycles"] += 1
            return raw
        return out
    upd.state = st
    return upd


def batch_cap(k: int, clusters) -> int:
    """Per-batch cluster cap: max(protocol cluster_cap_per_batch, ceil(k / #clusters among the candidates))."""
    n_cl = max(len(np.unique(np.asarray(clusters)[np.asarray(clusters) >= 0])), 1)
    return int(max(int(protocol_selector()["cluster_cap_per_batch"]), -(-int(k) // n_cl)))


def offarray_calibration(override: str | None = None) -> bool:
    """Route calibration mode (config/strategies.yaml `offarray_calibration`, DESIGN §1 / D5): 'quantile' | 'none'."""
    mode = str(override if override is not None else knobs().get("offarray_calibration", "quantile"))
    if mode not in ("quantile", "none"):
        raise ValueError(f"offarray_calibration {mode!r}")
    return mode == "quantile"


def quantile_map(main_on: np.ndarray, route_on: np.ndarray, route_off: np.ndarray) -> np.ndarray:
    """Off-array scores on the on-array scale: percentile of route_off within route_on (mid-rank; NaN -> median)
    mapped to the same quantile of main_on."""
    p = midrank_pct(route_off, route_on)
    p = np.clip(np.where(np.isfinite(p), p, 0.5), 0.0, 1.0)
    return np.quantile(np.asarray(main_on, dtype=np.float64), p)


def calibrate_routes(g_main: np.ndarray, g_alt, on: np.ndarray, exclude=None) -> np.ndarray:
    """Route calibration of gains; off-array genes get a quantile of the on-array main-route gains.
    g_alt given (D5 style): the quantile matching their percentile among ON-array genes under a COMMON alternative
      route `g_alt` (defined for every candidate).
    g_alt None ('within route'): the quantile matching their percentile among the OFF-array genes under their own
      route, i.e. off-array gains are forced onto the on-array gain distribution (route-stratified; used where the
      common-route comparison drifts along a trajectory, see coverage.py).
    `exclude` (bool) drops genes (e.g. already picked in this batch) from the references."""
    g = np.array(g_main, dtype=np.float64, copy=True)
    on = np.asarray(on, dtype=bool)
    off = ~on
    keep = np.isfinite(g) if exclude is None else np.isfinite(g) & ~np.asarray(exclude, dtype=bool)
    ref = on & keep
    if g_alt is None:
        ref_off = off & keep
        if ref.sum() >= 2 and ref_off.sum() >= 1:
            g[off] = quantile_map(g[ref], g_main[ref_off], g_main[off])
        return g
    ref &= np.isfinite(g_alt)
    if ref.sum() >= 2 and off.any():
        g[off] = quantile_map(g[ref], np.asarray(g_alt)[ref], np.asarray(g_alt)[off])
    return g


# ============================================================================ base selector
class Selector:
    """Harness selector protocol. Subclasses implement `_choose(view, k, rng) -> list[int]` (indices into
    view.candidates, greedy order). `reads` declares which SelectorView data channels the selector consumes
    ("rows" = screened-row GI data E/Y; "smf" = array SMF); the leakage tests check it both ways."""
    name = "base"
    needs_view = True
    static = False
    reads: frozenset = frozenset()

    def __init__(self, deterministic: bool = False, top_k: int | None = None, **params):
        self.deterministic = bool(deterministic)
        self.top_k = int(protocol_selector()["randomised_greedy_top_k"] if top_k is None else top_k)
        self.params = {"deterministic": self.deterministic, "top_k": self.top_k, **params,
                       "knobs_sha": knobs_sha()}               # strategies.yaml enters the trajectory id
        self.diagnostics: list = []
        self._feats = FeatureCache()
        self._g0 = self._qhat = None

    # ---- harness protocol
    def reset(self, candidates, rng) -> None:
        self.diagnostics = []
        self._round = 0
        self._g0 = self._qhat = None

    def select(self, view, k: int, rng) -> list:
        t0 = time.perf_counter()
        cand = np.asarray(view.candidates, dtype=str)
        k = int(min(k, len(cand)))
        idx = [int(i) for i in self._choose(view, k, rng)]
        if len(idx) < k:                                   # e.g. M with too few on-array candidates
            rest = np.setdiff1d(np.arange(len(cand)), idx)
            fill = rest[rng.permutation(len(rest))][: k - len(idx)]
            self._diag(n_filled_at_random=int(len(fill)))
            idx += [int(i) for i in fill]
        if len(set(idx)) != len(idx):
            raise AssertionError(f"{self.name}: duplicate picks")
        on = np.asarray(view.candidate_col) >= 0
        rec = {"picks": [str(cand[i]) for i in idx], "picks_offarray": [bool(not on[i]) for i in idx],
               "offarray_share_candidates": float((~on).mean()) if len(on) else None}
        g0 = getattr(self, "_g0", None)
        if g0 is not None and len(g0) == len(cand):         # first-pick gains / scores of this round
            g0 = np.asarray(g0, dtype=np.float64)
            order = np.lexsort((np.arange(len(g0)), -g0))
            rank = np.empty(len(g0), dtype=np.int64)
            rank[order] = np.arange(len(g0))
            rec.update(picks_gain0=[float(g0[i]) for i in idx], picks_rank0=[int(rank[i]) for i in idx])
        q = getattr(self, "_qhat", None)
        if q is not None and len(q) == len(cand):           # predicted logit query-degree q̂ of the picks
            rec["picks_qhat_logit"] = [float(q[i]) for i in idx]
        self._g0 = self._qhat = None
        self._diag(round=getattr(self, "_round", 0), n_screened=int(view.n_screened),
                   select_s=time.perf_counter() - t0, **rec)
        self._round = getattr(self, "_round", 0) + 1
        return [str(cand[i]) for i in idx]

    # ---- helpers for subclasses
    def _diag(self, **kw):
        r = getattr(self, "_round", 0)
        if not self.diagnostics or self.diagnostics[-1].get("_r") != r:
            self.diagnostics.append({"_r": r})
        self.diagnostics[-1].update(kw)

    def _greedy(self, gains, k, rng, clusters=None, inner=None) -> list:
        self._g0 = np.array(gains, dtype=np.float64, copy=True)
        upd = inner
        if clusters is not None:
            cap = batch_cap(k, clusters)
            upd = capped_update(clusters, cap, inner)
            self._diag(cluster_cap=cap)
        out = harness.randomised_greedy(gains, k, rng, top_k=self.top_k, deterministic=self.deterministic,
                                        update=upd)
        if clusters is not None:
            self._diag(cap_cycles=int(upd.state["cycles"]))
        return out

    def _choose(self, view, k, rng) -> list:          # pragma: no cover
        raise NotImplementedError

    def gains(self, view) -> np.ndarray:              # first-pick gains over view.candidates (diagnostics/tests)
        raise NotImplementedError
