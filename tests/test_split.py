"""Tests of the frozen test/dev split (scripts/01_make_split.py -> config/split/)."""
from __future__ import annotations

import hsl  # noqa: F401

import hashlib
import importlib.util
import json

import numpy as np
import pytest
from scipy import stats

from hsl import paths

SPLIT = paths.PROJECT_ROOT / "config" / "split"


def _load_script():
    spec = importlib.util.spec_from_file_location("make_split", paths.PROJECT_ROOT / "scripts" / "01_make_split.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _read(p):
    return [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]


@pytest.fixture(scope="module")
def ms():
    return _load_script()


@pytest.fixture(scope="module")
def counts(ms):
    return ms.sl_counts(paths.costanzo_dataset())[0]


@pytest.fixture(scope="module")
def lists():
    return _read(SPLIT / "test_queries.txt"), _read(SPLIT / "dev_pool.txt")


def test_sorted_unique(lists):
    for lst in lists:
        assert lst == sorted(lst) and len(lst) == len(set(lst))


def test_disjoint_and_union_is_universe(lists, counts):
    test, dev = map(set, lists)
    assert not test & dev
    assert test | dev == set(counts.index)
    assert len(counts) == 3283


def test_fraction_within_one_gene_per_stratum(ms, lists, counts):
    test = set(lists[0])
    st, labels = ms.strata(counts)
    assert len(labels) >= 2 and labels[0] == "zero"
    assert (counts[st == 0] == 0).all() and (counts[st > 0] > 0).all()
    for k in range(len(labels)):
        members = st.index[st.to_numpy() == k]
        n_t = sum(o in test for o in members)
        assert abs(n_t - 0.20 * len(members)) <= 1, (k, n_t, len(members))


def test_manifest_checksums(lists):
    man = json.loads((SPLIT / "split_manifest.json").read_text())
    for name in ("test_queries.txt", "dev_pool.txt"):
        assert hashlib.sha256((SPLIT / name).read_bytes()).hexdigest() == man["sha256"][name]
    assert man["split_seed"] == 20261029
    assert man["n_test"] == len(lists[0]) and man["n_dev_pool"] == len(lists[1])


def test_deterministic_rerun(ms, tmp_path):
    ms.main(["--out-dir", str(tmp_path)])
    for name in ("test_queries.txt", "dev_pool.txt", "split_manifest.json"):
        assert (tmp_path / name).read_bytes() == (SPLIT / name).read_bytes(), name


def test_stratification_ks_not_significant(lists, counts):
    t = counts.loc[lists[0]].to_numpy()
    d = counts.loc[lists[1]].to_numpy()
    assert stats.ks_2samp(t, d).pvalue > 0.05


def test_lockbox_preferred(ms, lists, counts):
    """Every stratum takes lockbox members first: a stratum's test set contains a non-lockbox gene only if all its
    lockbox members are already in test."""
    from ysl.eval.splits import read_lockbox
    lb = set(read_lockbox(paths.costanzo_dataset()).tolist())
    test = set(lists[0])
    st, labels = ms.strata(counts)
    for k in range(len(labels)):
        members = set(st.index[st.to_numpy() == k])
        t_k, lb_k = members & test, members & lb
        if t_k - lb:
            assert lb_k <= t_k
        else:
            assert t_k <= lb_k


def test_seeds_yaml():
    import yaml
    s = yaml.safe_load((paths.PROJECT_ROOT / "config" / "seeds.yaml").read_text())
    assert s["split_seed"] == 20261029
    assert "append-only" in (paths.PROJECT_ROOT / "config" / "seeds.yaml").read_text().lower()


def test_no_yeast_paths_written():
    assert not np.any([str(paths.yeast_sls_root()) in str(p) for p in SPLIT.iterdir()])
