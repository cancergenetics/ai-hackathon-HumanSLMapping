"""Toy worlds for strategy tests: a duck-typed SelectorView over a synthetic GI matrix (same attributes and methods
selectors are allowed to use). No real data."""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from hsl.strategies import base

_TOKENS = itertools.count()
ESM = [f"esm_pc{i}" for i in range(16)]


class ToyWorld:
    """kind = 'hub' (planted query propensity a_g, driving SLs, SMF, abundance), 'cluster' (10 planted external-
    kernel clusters; iid labels) or 'null' (no structure anywhere)."""

    def __init__(self, seed: int, kind: str = "hub", n_on: int = 200, n_off: int = 40, n_ao: int = 80,
                 n_clusters: int = 10, cluster_of_top_hubs: bool = False):
        rng = np.random.default_rng(seed)
        self.kind = kind
        on = [f"g{i:03d}" for i in range(n_on)]
        off = [f"q{i:03d}" for i in range(n_off)]
        ao = [f"a{i:03d}" for i in range(n_ao)]
        self.queries = np.array(sorted(on + off))
        self.array_orfs = np.array(sorted(on + ao))
        self.genes = np.array(sorted(set(self.queries) | set(self.array_orfs)))
        G = len(self.genes)
        gi = {g: i for i, g in enumerate(self.genes)}
        self.col = {g: j for j, g in enumerate(self.array_orfs)}
        qi = np.array([gi[g] for g in self.queries])
        ai = np.array([gi[g] for g in self.array_orfs])
        # latent structure
        a = rng.normal(0, 1, G) if kind == "hub" else np.zeros(G)
        b = 0.5 * a + 0.85 * rng.normal(0, 1, G) if kind == "hub" else np.zeros(G)
        self.a = dict(zip(self.genes, a))
        lin = -3.9 + 1.2 * a[qi][:, None] + 1.0 * b[ai][None, :]
        if kind != "hub":
            lin = np.full((len(qi), len(ai)), np.log(0.02 / 0.98))
        Y = (rng.random(lin.shape) < 1 / (1 + np.exp(-lin))).astype(np.int8)
        E = np.where(Y == 1, rng.normal(-0.3, 0.08, Y.shape), rng.normal(0, 0.05, Y.shape)).astype(np.float32)
        for r, g in enumerate(self.queries):                     # self pairs unmeasured
            if g in self.col:
                Y[r, self.col[g]] = -1
                E[r, self.col[g]] = np.nan
        self.Y, self.E = Y, E
        self.row = {g: r for r, g in enumerate(self.queries)}
        self.smf = (0.9 - 0.06 * a[ai] + rng.normal(0, 0.04, len(ai))).astype(np.float32) if kind == "hub" else \
            rng.normal(0.9, 0.06, len(ai)).astype(np.float32)
        # clusters / external kernels
        cl = rng.integers(0, n_clusters, G)
        if kind == "cluster":
            cl = np.arange(G) % n_clusters
            rng.shuffle(cl)
        if cluster_of_top_hubs:
            top = np.argsort(-a)[:30]
            cl[top] = 0
        self.cl = cl
        coex = rng.normal(0, 0.1, (G, G))
        coex = (coex + coex.T) / 2
        if kind == "cluster":
            coex = coex + 2.0 * (cl[:, None] == cl[None, :])
        np.fill_diagonal(coex, 3.0)
        self.coex = coex
        ap = (rng.random((G, G)) < 0.01).astype(np.float64)
        ap = np.maximum(ap, ap.T)
        np.fill_diagonal(ap, 0)
        self.apms = ap if kind != "cluster" else np.zeros((G, G))
        cent = rng.normal(0, 1, (n_clusters, 16))
        esm = rng.normal(0, 1, (G, 16)) if kind != "cluster" else cent[cl] + rng.normal(0, 0.3, (G, 16))
        f = base.knobs()["degree_model"]["features"]
        ab = 2 + 0.6 * a + rng.normal(0, 0.6, G) if kind == "hub" else rng.normal(2, 0.8, G)
        df = pd.DataFrame(esm, columns=ESM, index=self.genes)
        df[f["abundance"]] = ab
        df[f["abundance_fallback"]] = ab + rng.normal(0, 0.2, G)
        df[f["apms_degree"]] = np.log1p(ap.sum(1))
        df[f["coexp_degree"]] = rng.normal(0, 1, G)
        df[f["conservation"]] = rng.random(G)
        df[f["paralog_count"]] = (rng.random(G) < 0.15).astype(float)
        self.feats = df
        self.gi = gi
        self.token = f"toy{next(_TOKENS)}"

    def true_degree(self, genes) -> np.ndarray:
        return np.array([(self.Y[self.row[g]] == 1).sum() for g in genes])

    def view(self, S, pool=None) -> "ToyView":
        return ToyView(self, np.asarray(S, dtype=str), self.queries if pool is None else np.asarray(pool, dtype=str))


class ToyView:
    def __init__(self, w: ToyWorld, S, pool):
        self._w = w
        S = np.sort(np.asarray(S, dtype=str))
        self.screened_orfs = S
        self.screened_col = np.array([w.col.get(g, -1) for g in S], dtype=np.int64)
        r = [w.row[g] for g in S]
        self.E = w.E[r] if len(r) else np.zeros((0, len(w.array_orfs)), np.float32)
        self.Y = w.Y[r] if len(r) else np.zeros((0, len(w.array_orfs)), np.int8)
        self.SP = np.where(np.isnan(self.E), np.nan, 5.0 * (self.Y == 1)).astype(np.float32)
        self.array_orfs = w.array_orfs
        self.array_orfs_current = w.array_orfs
        self.array_smf = w.smf
        self.visible_cols = np.arange(len(w.array_orfs))
        self.candidates = np.setdiff1d(pool, S)
        self.candidate_col = np.array([w.col.get(g, -1) for g in self.candidates], dtype=np.int64)
        self.cache_token = w.token

    @property
    def n_screened(self):
        return len(self.screened_orfs)

    def gene_features(self, orfs):
        return self._w.feats.loc[np.asarray(orfs, dtype=str)].copy()

    def similarity(self, kind, a, b):
        ia = [self._w.gi[g] for g in np.asarray(a, dtype=str)]
        ib = [self._w.gi[g] for g in np.asarray(b, dtype=str)]
        if kind == "coexp":
            return self._w.coex[np.ix_(ia, ib)].astype(np.float32)
        if kind == "apms":
            return self._w.apms[np.ix_(ia, ib)].astype(np.float32)
        if kind == "esm":
            X = self._w.feats[ESM].to_numpy()
            X = X / np.linalg.norm(X, axis=1, keepdims=True)
            return (X[ia] @ X[ib].T).astype(np.float32)
        raise ValueError(kind)

    def clusters(self, orfs):
        return self._w.cl[[self._w.gi[g] for g in np.asarray(orfs, dtype=str)]].astype(np.int64)


def run_toy(world: ToyWorld, sel, seed: int, n_max: int, batch: int = 10) -> list:
    """Mini trajectory (same protocol as harness.Trajectory) on a toy world; validates every batch."""
    rng = np.random.default_rng(seed)
    sel.reset(world.queries, rng)
    S: list = []
    while len(S) < n_max:
        v = world.view(S)
        k = min(batch, n_max - len(S))
        picks = [str(g) for g in sel.select(v, k, rng)]
        assert len(picks) == k and len(set(picks)) == k
        assert np.isin(picks, v.candidates).all() and not np.isin(picks, S).any()
        S.extend(picks)
    return S
