"""Feature provenance (config/feature_provenance.yaml), DATA.md coverage, repository hygiene and config files."""
from __future__ import annotations

import hsl  # noqa: F401

import fnmatch
import os
import subprocess

import numpy as np
import pytest
import yaml

from hsl import config, guard, harness as H, paths, views

PROV = yaml.safe_load((paths.PROJECT_ROOT / "config" / "feature_provenance.yaml").read_text())["stores"]
ROOT = os.path.realpath(paths.yeast_sls_root())


def _uncovered(read_paths, consumer) -> list[str]:
    bad = []
    for p in read_paths:
        if not p.startswith(ROOT + os.sep):
            continue
        rel = os.path.relpath(p, ROOT)
        if not any(consumer in s["used_by"] and any(fnmatch.fnmatch(rel, pat) for pat in s["files"])
                   for s in PROV.values()):
            bad.append(rel)
    return bad


def test_provenance_yaml_fields():
    for name, s in PROV.items():
        assert {"files", "data_md", "role", "label_free", "gi_text_go_free", "human_analogue", "used_by"} <= set(s)
        if name != "costanzo_s1_bundle":
            assert s["label_free"] and s["gi_text_go_free"], name


def test_every_store_listed_in_data_md():
    text = (paths.PROJECT_ROOT / "DATA.md").read_text()
    missing = [n for n, s in PROV.items() if s["data_md"] not in text]
    assert not missing, missing


def test_selector_features_only_open_listed_stores(dev_b, pseudo):
    with guard.record_reads() as ps:
        v = views.SelectorView.build(dev_b, pseudo["S"], pseudo["split"])
        f = v.gene_features(v.candidates[:20])
        for k in ("coexp", "esm", "apms"):
            v.similarity(k, v.candidates[:5], v.candidates[:5])
        cl = v.clusters(v.candidates)
    assert not _uncovered(ps, "selector")
    assert not any(c.startswith("has_") for c in f.columns)
    assert not any("chrom" in c or c == "mid" for c in f.columns)
    assert (cl >= 0).all() and len(np.unique(cl)) > 50


def test_broken_twin_unlisted_store_is_caught(dev_b):
    from ysl.features import pairs as PF
    with guard.record_reads() as ps:
        PF.load("gen_same_chr")              # genome store: used_by = [] -> must be flagged
    assert _uncovered(ps, "selector") and _uncovered(ps, "stack_lr_wr2")


@pytest.mark.parametrize("method", ["stack_lr", "stack_lr_wr2"])
def test_predictor_only_opens_listed_stores(method):
    S = tuple(np.sort(np.random.default_rng(0).choice(views.dev_split(1).pool, 30, replace=False)))
    with guard.record_reads() as ps:
        H.run_predictor(H.RunSpec("dev1", S), method, seeds={"model_seed": 1}, save=False, ledger=False)
    assert not _uncovered(ps, method)


@pytest.mark.slow
def test_m1_hgb_only_opens_listed_stores():
    S = tuple(np.sort(np.random.default_rng(0).choice(views.dev_split(1).pool, 30, replace=False)))
    with guard.record_reads() as ps:
        H.run_predictor(H.RunSpec("dev1", S), "m1_hgb", {"drop": ["genome"]}, seeds={"model_seed": 1}, save=False,
                        ledger=False)
    assert not _uncovered(ps, "m1_hgb")


# ----------------------------------------------------------------------------- repo hygiene
def test_no_data_or_predictions_tracked():
    out = subprocess.run(["git", "-C", str(paths.PROJECT_ROOT), "ls-files"], capture_output=True, text=True).stdout
    bad = [p for p in out.splitlines() if p.startswith(("data/", "results/predictions/")) or "/predictions/" in p
           or p.endswith((".parquet", ".npz", ".npy"))]
    assert not bad, bad


def test_predictions_and_data_are_gitignored():
    for p in ("results/predictions/dev1/x.npz", "data/cache/clusters_x.npz", "results/x/y.parquet"):
        r = subprocess.run(["git", "-C", str(paths.PROJECT_ROOT), "check-ignore", "-q", p])
        assert r.returncode == 0, p


# ----------------------------------------------------------------------------- config files
def test_seeds_named_and_disjoint():
    s = config.seeds()
    for k in ("split_seed", "inner_fold_seed", "dev_traj_base_seed", "final_traj_base_seed", "bootstrap_seed",
              "universe_subsample_seed", "cluster_seed", "test_fixture_seed"):
        assert isinstance(s[k], int), k
    assert s["dev_traj_base_seed"] != s["final_traj_base_seed"]
    vals = [v for v in s.values() if isinstance(v, int)]
    assert len(vals) == len(set(vals))


def test_budgets():
    b = config.budgets()
    g = b["random_grid"]
    assert g == sorted(g) and 1640 in g and g[-1] == 2627 and 220 in g
    assert b["batch_size"] == 10 and set(b["q2_checkpoints"]) <= set(g)
    assert max(b["q2_checkpoints"]) == b["q2_max_n"] == 660
    assert b["batch_sizes_scenario"] == [10, 50, 110]
    assert b["m1_hgb_q1_grid"] == [30, 160, 490, 660, 1310, 2627]
    assert b["m1_hgb_q2_checkpoints"] == {"strategies": [160, 660], "random": [160, 490, 660]}
    assert b["warm_start_checkpoints"] == [220, 330, 440, 490, 660]


def test_protocol_keys():
    p = config.protocol()
    assert p["bootstrap"]["B_query_bootstrap"] == 2000 and p["bootstrap"]["B_two_way"] == 500
    assert p["stats"]["alpha_familywise"] == 0.025 and p["stats"]["g_star_min"] == 0.5
    assert p["selector"]["kmeans_k"] == 100 and p["selector"]["randomised_greedy_top_k"] == 20
    assert p["gate"]["gate_budgets"] == [30, 160, 330]
    assert "ap_bins" not in p["eval"]
