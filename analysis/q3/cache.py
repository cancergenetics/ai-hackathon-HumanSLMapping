"""Small label cache of the DEV bundle for the Q3 analyses (no T rows exist in it), so analyses that only need labels
do not hold the ~1 GB bundle in memory next to a running P* job.

    .venv/bin/python -m analysis.q3.cache          -> data/q3_cache/dev_cache.npz (gitignored /data/)

Contents: dev-row label matrix Y (int8, -1 = unmeasured) with row/array names and current names, array SMF, the
test-gene names (to exclude), and per dev split the E_fixed_dev cell index (targets, rows, cols, y, AO, VV masks) and
candidate pool, all through the frozen hsl.views code.
"""
from __future__ import annotations

import time
from functools import lru_cache

import numpy as np

from analysis.q3 import common as C

PATH = C.ROOT / "data" / "q3_cache" / "dev_cache.npz"


def build() -> None:
    from hsl import views
    b = C.dev_bundle()
    d = {"Y": np.asarray(b.Y[b.gene_rows(b.query_orfs)], dtype=np.int8),   # canonical strain row / ORF
         "query_orfs": np.asarray(b.query_orfs, dtype=str),
         "array_orfs": np.asarray(b.array_orfs, dtype=str),
         "array_orfs_current": np.asarray(b.array_orfs_current, dtype=str),
         "query_orfs_current": views.current_of(b, b.query_orfs),
         "array_smf": np.asarray(b.array_smf, dtype=np.float64), "test": np.asarray(views.test_queries(), dtype=str),
         "all_query_orfs": np.asarray(b.all_query_orfs, dtype=str)}
    for k in C.protocol()["folds"]["dev_splits"]:
        sp = views.dev_split(int(k))
        es = views.e_fixed(b, sp)
        p = f"dev{k}_"
        d[p + "targets"] = np.asarray(es.targets, dtype=str)
        d[p + "rows"], d[p + "cols"], d[p + "y"] = es.rows, es.cols, es.y
        d[p + "AO"], d[p + "VV"] = es.strata["AO"], es.strata["VV"]
        d[p + "pool"] = np.asarray(sp.pool, dtype=str)
        d[p + "digest"] = np.array(es.digest())
    PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PATH, **d)


@lru_cache(maxsize=1)
def load() -> dict:
    if not PATH.exists():
        build()
    with np.load(PATH, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    if d["Y"].shape != (len(d["query_orfs"]), len(d["array_orfs"])):
        raise AssertionError("cache Y is not gene-level (one canonical row per query ORF)")
    return d


def eset(split: str) -> dict:
    d = load()
    p = f"{split}_"
    return {k: d[p + k] for k in ("targets", "rows", "cols", "y", "AO", "VV", "pool", "digest")}


if __name__ == "__main__":
    t0 = time.perf_counter()
    build()
    print(f"wrote {PATH} in {time.perf_counter() - t0:.0f} s")
