#!/usr/bin/env python
"""Fixed held-out TEST split of query genes (made once, before any modelling; PROMPT.md, RESEARCH_TREE M2).

Universe: every query ORF of the Costanzo 2016 S1 NxN (gene level as upstream: canonical strain per ORF; 3,283).
SL count per query gene: number of strong-SL cells (P < 0.05 and eps < -0.12) in its query row (upstream label).
Strata: genes with 0 SLs form their own stratum; non-zero counts are cut into deciles with pandas.qcut(q=10,
duplicates="drop") (bins (e_i, e_{i+1}], lowest bin closed; identical edges caused by ties are merged).
Per stratum: n_test = floor(0.20 * n + 0.5) (round half up).
Selection within a stratum: the stratum's ORFs (sorted) are shuffled with ONE numpy Generator
default_rng(split_seed) consumed stratum by stratum in stratum order (zeros first, then deciles ascending); ORFs in the
upstream S1 lockbox (YEAST_SLs/config/lockbox_s1.txt, SHA-verified by ysl.eval.splits.read_lockbox) are taken first in
shuffled order, then non-lockbox ORFs in shuffled order fill the stratum if it has too few lockbox members.
Rationale for preferring the lockbox: upstream hyper-parameters were tuned on non-lockbox dev splits, so test genes
drawn from the lockbox minimise upstream-tuning contamination.

Outputs (default config/split/): test_queries.txt, dev_pool.txt (one ORF per line, sorted, no header),
split_manifest.json (seed, strata, counts, SL-count distributions test vs dev, KS test, lockbox overlap, SHA-256 of
both lists and of the input data). Deterministic: re-running gives byte-identical files.
"""
from __future__ import annotations

import hsl  # noqa: F401  (YEAST_SLs write guard + paths; must be first)

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from hsl import paths
from ysl.data.dense import load_bundle
from ysl.eval.splits import read_lockbox

TEST_FRACTION = 0.20
N_DECILES = 10


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def split_seed() -> int:
    with open(paths.PROJECT_ROOT / "config" / "seeds.yaml", encoding="utf-8") as f:
        return int(yaml.safe_load(f)["split_seed"])


def sl_counts(dataset: str) -> tuple[pd.Series, dict]:
    """Strong-SL count per query ORF (all query ORFs, incl. upstream lockbox) + input provenance."""
    b = load_bundle(dataset, unblind=True)          # verifies the bundle SHA-256 against its sidecar
    rows = b.gene_rows(b.all_query_orfs)
    counts = pd.Series((b.Y[rows] == 1).sum(1).astype(int), index=pd.Index(b.all_query_orfs, name="orf"))
    prov = {"dataset": dataset,
            "bundle_file": f"YEAST_SLs/data/processed/costanzo2016_bundle_{dataset}.npz",
            "bundle_sha256_verified_on_load": b.meta["bundle_sha256"],
            "raw_nxn_sha256_recorded_by_upstream": b.meta["raw_nxn_sha256"],
            "processed_parquet_sha256_recorded_by_upstream": b.meta["processed_parquet_sha256"],
            "n_query_orfs": int(len(counts)), "label": "P < 0.05 and eps < -0.12 (ysl.eval.labels)"}
    return counts, prov


def strata(counts: pd.Series) -> tuple[pd.Series, list]:
    s = pd.Series(-1, index=counts.index, dtype=int)
    nz = counts[counts > 0]
    bins, edges = pd.qcut(nz, q=N_DECILES, labels=False, retbins=True, duplicates="drop")
    s.loc[nz.index] = bins.astype(int)
    labels = ["zero"] + [f"({edges[i]:g}, {edges[i + 1]:g}]" if i else f"[{edges[0]:g}, {edges[1]:g}]"
                         for i in range(len(edges) - 1)]
    return s + 1, labels   # 0 = zero stratum, 1..K = nonzero quantile bins


def make_split(counts: pd.Series, lockbox: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, list]:
    st, labels = strata(counts)
    rng = np.random.default_rng(seed)
    lb = set(lockbox.tolist())
    test, per = [], []
    for k in range(len(labels)):
        orfs = np.sort(st.index[st.to_numpy() == k].to_numpy(dtype=str))
        n = len(orfs)
        n_test = int(np.floor(TEST_FRACTION * n + 0.5))
        shuf = orfs[rng.permutation(n)]
        in_lb = np.array([o in lb for o in shuf], dtype=bool)
        order = np.concatenate([shuf[in_lb], shuf[~in_lb]])
        chosen = order[:n_test]
        test.extend(chosen.tolist())
        c = counts.loc[orfs]
        per.append({"stratum": k, "label": labels[k], "n": n, "n_test": n_test, "n_dev": n - n_test,
                    "sl_min": int(c.min()), "sl_max": int(c.max()),
                    "n_lockbox_in_stratum": int(in_lb.sum()),
                    "n_test_from_lockbox": int(min(in_lb.sum(), n_test)),
                    "n_test_filled_from_non_lockbox": int(max(0, n_test - in_lb.sum()))})
    test = np.sort(np.array(test, dtype=str))
    dev = np.sort(counts.index.to_numpy(dtype=str)[~np.isin(counts.index.to_numpy(dtype=str), test)])
    return test, dev, per


def describe(x: pd.Series) -> dict:
    q = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1]
    return {"n": int(len(x)), "mean": float(x.mean()), "sd": float(x.std(ddof=1)), "n_zero": int((x == 0).sum()),
            "total_SL": int(x.sum()), **{f"q{int(v * 100):02d}": float(x.quantile(v)) for v in q}}


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", type=Path, default=paths.PROJECT_ROOT / "config" / "split")
    a = ap.parse_args(argv)
    ds = paths.costanzo_dataset()
    seed = split_seed()
    counts, prov = sl_counts(ds)
    lockbox = read_lockbox(ds)                       # SHA-256 verified against upstream config/frozen.yaml
    lb_raw = paths.lockbox_file().read_bytes()
    test, dev, per = make_split(counts, lockbox, seed)

    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    t_txt = "".join(f"{o}\n" for o in test).encode()
    d_txt = "".join(f"{o}\n" for o in dev).encode()
    (out / "test_queries.txt").write_bytes(t_txt)
    (out / "dev_pool.txt").write_bytes(d_txt)
    ct, cd = counts.loc[test], counts.loc[dev]
    ks = stats.ks_2samp(ct.to_numpy(), cd.to_numpy())
    mw = stats.mannwhitneyu(ct.to_numpy(), cd.to_numpy())
    lb_in_test = np.isin(test, lockbox)
    man = {
        "description": "Fixed held-out test split of S1 query genes; stratified by #strong SLs; prefers upstream lockbox",
        "script": "scripts/01_make_split.py",
        "split_seed": seed, "seed_source": "config/seeds.yaml:split_seed",
        "test_fraction_target": TEST_FRACTION, "rounding": "per stratum n_test = floor(0.2 n + 0.5)",
        "strata_rule": "zero-SL stratum + pandas.qcut(nonzero counts, 10, duplicates='drop')",
        "universe": {"n": int(len(counts)), "definition": "all S1 NxN query ORFs (upstream gene level)"},
        "n_test": int(len(test)), "n_dev_pool": int(len(dev)),
        "strata": per,
        "sl_count_distribution": {"test": describe(ct), "dev_pool": describe(cd), "all": describe(counts)},
        "ks_test_vs_dev": {"statistic": float(ks.statistic), "pvalue": float(ks.pvalue)},
        "mannwhitney_test_vs_dev": {"statistic": float(mw.statistic), "pvalue": float(mw.pvalue)},
        "lockbox": {"file": "YEAST_SLs/config/lockbox_s1.txt", "sha256": sha256_bytes(lb_raw),
                    "n_lockbox": int(len(lockbox)), "n_test_in_lockbox": int(lb_in_test.sum()),
                    "frac_test_in_lockbox": float(lb_in_test.mean()),
                    "n_lockbox_in_dev_pool": int(np.isin(dev, lockbox).sum()),
                    "frac_lockbox_in_test": float(lb_in_test.sum() / len(lockbox))},
        "test_on_array": {"n": int(sum(o in set(load_bundle(ds, unblind=True).array_orfs) for o in test))},
        "input": prov,
        "sha256": {"test_queries.txt": sha256_bytes(t_txt), "dev_pool.txt": sha256_bytes(d_txt)},
    }
    (out / "split_manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps({k: man[k] for k in ("n_test", "n_dev_pool", "ks_test_vs_dev", "lockbox", "sha256")}, indent=1))
    return man


if __name__ == "__main__":
    main()
