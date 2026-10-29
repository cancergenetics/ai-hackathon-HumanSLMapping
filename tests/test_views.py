"""hsl.views: label gate, dev bundle, folds, gene identity, evaluation-set invariants, unblinded-loader grep test.
Dev rows only; final-set checks use NAMES only (column sets), never T rows."""
from __future__ import annotations

import hsl  # noqa: F401

import hashlib
import re
from collections import Counter

import numpy as np
import pytest

from hsl import config, guard, paths, views

REQUIRED_IDENTITY = {"YFL056C", "YHR079C-B", "YIR043C", "YJL206C-A"}


# ----------------------------------------------------------------------------- runtime label gate
def test_unblinded_load_refused_outside_views(monkeypatch):
    monkeypatch.delenv("HSL_FINAL", raising=False)
    from ysl.data.dense import load_bundle
    for unblind in (True, False):
        with pytest.raises(guard.LabelGateViolation):
            load_bundle(paths.costanzo_dataset(), unblind=unblind, cache=False)


def test_costanzo_files_refused_outside_views(monkeypatch):
    monkeypatch.delenv("HSL_FINAL", raising=False)
    import pandas as pd
    from ysl.data import costanzo as C
    with pytest.raises(guard.LabelGateViolation):
        pd.read_parquet(C.PROCESSED_DIR / "costanzo2016_nxn_science_subset_queries.parquet")
    with pytest.raises(guard.LabelGateViolation):
        np.load(C.PROCESSED_DIR / "costanzo2016_bundle_s1.npz")
    with pytest.raises(guard.LabelGateViolation):
        open(C.RAW_DIR / "anything.txt", "rb")
    with pytest.raises(guard.LabelGateViolation):
        with guard.label_access():
            pass


def test_feature_stores_remain_readable():
    from ysl.data import features as F
    assert len(F.load_table("abundance")) > 1000


def test_final_bundle_requires_env(monkeypatch):
    monkeypatch.delenv("HSL_FINAL", raising=False)
    with pytest.raises(PermissionError):
        views.final_bundle()


# ----------------------------------------------------------------------------- dev bundle
def test_dev_bundle_has_no_test_rows(dev_b):
    T = views.test_queries()
    assert not np.isin(dev_b.strain_orf, T).any()
    assert len(dev_b.query_orfs) == 2627 and dev_b.E.shape == (2627, 3827)
    with pytest.raises(PermissionError):
        dev_b.gene_rows(T[:3])
    assert not dev_b.E.flags.writeable


def test_dev_bundle_assertions(dev_b):
    assert dev_b.meta["gene_counts"] == views.EXPECTED_COUNTS
    rc = dev_b.meta["row_checks"]
    assert rc["n_rows"] == 2627 and rc["label_mismatch_at_P_threshold"] <= 10
    assert set(dev_b.all_query_orfs) == set(views.test_queries()) | set(views.dev_pool())


def test_upstream_lockbox_genes_in_dev_pool_are_usable(dev_b):
    """23 upstream-lockbox genes sit in D: from_split on the blinded dev bundle refuses them; construct_view works."""
    from ysl.eval.splits import Split, read_lockbox
    from ysl.models.view import ScreenView
    from hsl import upstream as U
    lb = read_lockbox("s1")
    lbD = np.intersect1d(lb, views.dev_pool())
    assert len(lbD) == 23
    other = np.setdiff1d(views.dev_pool(), lbD)[:5]
    with pytest.raises(PermissionError):
        ScreenView.from_split(dev_b, Split("final", 0, lbD, other, lb, "s1", 3283))
    v = U.construct_view(dev_b, lbD, other)
    assert v.n == 23 and np.isfinite(v.E).any()
    v2 = U.construct_view(dev_b, other, lbD)
    assert v2.n_targets == 23


# ----------------------------------------------------------------------------- folds / splits
def test_inner_folds_balanced_stratified_deterministic(dev_b):
    f = views.inner_folds()
    D = views.dev_pool()
    assert set(f) == set(D)
    sizes = Counter(f.values())
    assert set(sizes) == {1, 2, 3, 4, 5} and max(sizes.values()) - min(sizes.values()) <= 5
    st = views.query_strata(D, dev_b)
    for h in np.unique(st):
        c = Counter(f[o] for o in D[st == h])
        assert max(c.values()) - min(c.get(k, 0) for k in range(1, 6)) <= 2, (h, c)
    views.inner_folds.cache_clear()
    assert views.inner_folds() == f


def test_identity_groups_share_a_fold(dev_b):
    f = views.inner_folds()
    D = views.dev_pool()
    cur = views.current_of(dev_b, D)
    for c in np.unique(cur):
        assert len({f[o] for o in D[cur == c]}) == 1


def test_query_strata_match_split_manifest(dev_b):
    import json
    man = json.loads((views.SPLIT_DIR / "split_manifest.json").read_text())
    st = views.query_strata(views.dev_pool(), dev_b)
    got = Counter(st.tolist())
    for s in man["strata"]:
        assert got[s["stratum"]] == s["n_dev"], s


def test_dev_splits(dev_b):
    vals = []
    for k in (1, 2, 3):
        sp = views.dev_split(k)
        vals.append(set(sp.targets))
        assert not set(sp.targets) & set(sp.pool)
        assert set(sp.targets) | set(sp.pool) | set(sp.pool_removed) == set(views.dev_pool())
        assert set(views.test_queries()) <= set(sp.excluded)
        assert not np.isin(sp.pool, views.test_queries()).any()
    assert not (vals[0] & vals[1]) and not (vals[1] & vals[2])
    rs = views.resplit()
    assert len(rs.targets) + len(rs.pool) == 2627 and not (set(rs.targets) & (vals[0] | vals[1] | vals[2]))


# ----------------------------------------------------------------------------- gene identity
def test_identity_exclusions_final_columns(dev_b):
    """Final E_fixed column rule from NAMES only: required identity arrays excluded, no dev-pool array kept."""
    fs = views.final_split()
    ok, ident = views._fixed_col_ok(dev_b, fs)
    assert REQUIRED_IDENTITY <= set(ident)
    kept = dev_b.array_orfs[ok]
    assert not np.isin(kept, views.dev_pool()).any()
    assert not np.isin(dev_b.array_orfs_current[ok], views.current_of(dev_b, views.dev_pool())).any()
    assert len(fs.pool_removed) == 0 and len(fs.pool) == 2627          # rule (i) empty for this split


def test_identity_report_lists_groups(dev_b):
    rep = views.identity_report(dev_b)
    assert rep["n_groups"] > 0
    assert "YFL056C" in rep["groups"] and "YFL057C" in rep["groups"]["YFL056C"]["query"]


def test_final_e_fixed_needs_test_rows(dev_b):
    with pytest.raises(PermissionError):
        views.e_fixed(dev_b, views.final_split())


# ----------------------------------------------------------------------------- E-set invariants (dev)
@pytest.mark.parametrize("k", [1, 2, 3])
def test_e_fixed_dev_invariants(dev_b, k):
    sp = views.dev_split(k)
    es = views.e_fixed(dev_b, sp)
    cols = dev_b.array_orfs[es.cols]
    assert not np.isin(cols, views.test_queries()).any()                 # no test-gene columns
    assert not np.isin(cols, sp.cand_full).any()                          # no candidate arrays
    assert not np.isin(dev_b.array_orfs_current[es.cols], views.current_of(dev_b, sp.cand_full)).any()
    rows = dev_b.gene_rows(es.targets[es.rows])
    assert (dev_b.Y[rows, es.cols] == es.y).all() and (es.y >= 0).all()   # measured cells, own labels
    assert not (views.current_of(dev_b, es.targets)[es.rows] == dev_b.array_orfs_current[es.cols]).any()
    c = es.counts()
    assert c["n_cells"] > 500_000 and 0.005 < c["prevalence"] < 0.02
    assert c["AO"]["n_cells"] + c["VV"]["n_cells"] == c["n_cells"]
    print(sp.name, c)


def test_e_unobs_contains_e_fixed(dev_b):
    sp = views.dev_split(1)
    S = sp.pool[:50]
    ef, eu = views.e_fixed(dev_b, sp), views.e_unobs(dev_b, sp, S)
    a = set(zip(ef.rows.tolist(), ef.cols.tolist()))
    b = set(zip(eu.rows.tolist(), eu.cols.tolist()))
    assert a <= b and len(b) > len(a)
    assert not np.isin(dev_b.array_orfs[eu.cols], S).any()


# ----------------------------------------------------------------------------- grep test: unblinded loader
UNBLINDED = re.compile(r"final_bundle\s*\(|load_s1\s*\(|unblind\s*=\s*True|label_access\s*\(|efixed_label_summary")
# Pre-P0 scripts (M0–M2: data characterisation, split creation, runtime benchmark) read the unblinded bundle before
# the gate existed; pinned by SHA-256 (they now run only with HSL_FINAL=1). Any edit fails this test.
LEGACY = {
    "scripts/00_data_characterisation.py": None,
    "scripts/01_make_split.py": None,
    "scripts/02_runtime_benchmark.py": None,
}
LEGACY_SHA_FILE = paths.PROJECT_ROOT / "tests" / "legacy_unblinded_scripts.sha256"
ALLOWED = {"scripts/run_final.py", "scripts/04_eset_manifest.py"}   # eset manifest: label counts of E_fixed only


def _legacy_pins() -> dict:
    out = {}
    for ln in LEGACY_SHA_FILE.read_text().splitlines():
        if ln.strip():
            h, p = ln.split()
            out[p] = h
    return out


def test_no_script_calls_unblinded_loader():
    pins = _legacy_pins()
    assert set(pins) == set(LEGACY)
    offenders = []
    for p in sorted((paths.PROJECT_ROOT / "scripts").rglob("*.py")):
        rel = str(p.relative_to(paths.PROJECT_ROOT))
        if UNBLINDED.search(p.read_text()):
            if rel in LEGACY:
                assert hashlib.sha256(p.read_bytes()).hexdigest() == pins[rel], f"legacy script changed: {rel}"
            elif rel not in ALLOWED:
                offenders.append(rel)
    for p in sorted((paths.PROJECT_ROOT / "src" / "hsl").rglob("*.py")):
        rel = str(p.relative_to(paths.PROJECT_ROOT))
        if rel in ("src/hsl/views.py", "src/hsl/upstream.py", "src/hsl/guard.py"):
            continue
        if UNBLINDED.search(p.read_text()):
            offenders.append(rel)
    assert not offenders, offenders


def test_manifest_script_uses_only_label_counts():
    p = paths.PROJECT_ROOT / "scripts" / "04_eset_manifest.py"
    txt = p.read_text()
    assert "final_bundle" not in txt and "load_bundle" not in txt and "label_access" not in txt
