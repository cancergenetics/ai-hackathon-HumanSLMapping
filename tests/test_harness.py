"""Harness: manifests/ledger/predictions, prediction poisoning, trajectory runner invariants, end-to-end controls.
Dev data only (dev split 1). Writes go to tmp_path (ledger, predictions, trajectory logs are redirected)."""
from __future__ import annotations

import hsl  # noqa: F401

import json

import numpy as np
import pytest

from hsl import config, harness as H, views
from hsl import eval as E


@pytest.fixture(autouse=True)
def _redirect(tmp_path, monkeypatch):
    monkeypatch.setattr(H, "LEDGER", tmp_path / "LEDGER.jsonl")
    monkeypatch.setattr(H, "PRED_DIR", tmp_path / "predictions")
    monkeypatch.setattr(H, "TRAJ_DIR", tmp_path / "trajectories")
    return tmp_path


def _S(n, seed=0, split="dev1"):
    sp = views.get_split(split)
    return tuple(np.sort(np.random.default_rng(seed).choice(sp.pool, n, replace=False)))


class AdaptiveToySelector:
    """Uses the SelectorView (array-side SL counts of candidates) + seeded randomised greedy."""
    needs_view, static, name = True, False, "toy_adaptive"

    def reset(self, candidates, rng):
        pass

    def select(self, view, k, rng):
        Y = view.Y.astype(float)
        Y[Y < 0] = 0
        col = view.candidate_col
        g = np.where(col >= 0, Y[:, np.maximum(col, 0)].sum(0), -1.0) if view.n_screened else np.zeros(len(col))
        return [str(view.candidates[i]) for i in H.randomised_greedy(g, k, rng)]


# ----------------------------------------------------------------------------- run_predictor
def test_run_predictor_manifest_ledger_predictions(_redirect):
    S = _S(30)
    r = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", seeds={"model_seed": 3})
    m = r.manifest
    for k in ("git", "seeds", "params", "n_folds_passed", "split_seed_passed", "n_threads", "runtime_s",
              "peak_rss_gb", "bundle_sha256", "split_digest", "config_hash"):
        assert k in m, k
    assert set(m["git"]) == {"head", "dirty", "diff_sha256"}
    assert m["n_screened"] == 30 and m["n_folds_passed"] == 10 and m["n_threads"] <= 3
    assert m["bundle_view"] == "dev_bundle"
    lines = (_redirect / "LEDGER.jsonl").read_text().splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert {"id", "parent", "config_hash", "split", "outcome", "metrics"} <= set(rec) and rec["outcome"] == "ok"
    pred = np.load(_redirect.parent / _redirect.name / "predictions" / "dev1" / f"{r.run_id}.npz")
    assert pred["E_fixed"].dtype == np.float32 and np.array_equal(pred["E_fixed"], r.scores["E_fixed"])
    es = r.esets["E_fixed"]
    assert abs(r.metrics["E_fixed_dev"]["auprc"] - E.average_precision(es.y, r.scores["E_fixed"])) < 1e-12
    assert r.metrics["E_fixed_dev"]["n_cells"] == es.n


def test_failed_run_is_a_ledger_line(_redirect):
    with pytest.raises(Exception):
        H.run_predictor(H.RunSpec("dev1", _S(30)), "no_such_method")
    rec = json.loads((_redirect / "LEDGER.jsonl").read_text().splitlines()[-1])
    assert rec["outcome"] == "error"


def test_final_runs_need_explicit_bundle():
    with pytest.raises(PermissionError):
        H.run_predictor(H.RunSpec("final", tuple(views.dev_pool()[:30])), "stack_lr", ledger=False, save=False)


def test_m1_hgb_small_n_params():
    assert H._method_params("m1_hgb", {}, 10)["prop_folds"] == 5
    assert "prop_folds" not in H._method_params("m1_hgb", {}, 30)
    assert H._method_params("m1_hgb", {"prop_folds": 7}, 10)["prop_folds"] == 7


# ----------------------------------------------------------------------------- prediction poisoning
def _poison_rows(b, orfs):
    from test_selector_view import poison_rows
    return poison_rows(b, orfs, seed=5)


def test_prediction_poisoning(dev_b):
    sp = views.dev_split(1)
    S = _S(30)
    kw = dict(seeds={"model_seed": 1}, save=False, ledger=False)
    base = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", bundle=dev_b, **kw).scores["E_fixed"]
    pv = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", bundle=_poison_rows(dev_b, sp.targets), **kw)
    assert np.array_equal(pv.scores["E_fixed"], base)                        # D_val rows poisoned -> identical
    unscreened = np.setdiff1d(sp.pool, S)
    pu = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", bundle=_poison_rows(dev_b, unscreened), **kw)
    assert np.array_equal(pu.scores["E_fixed"], base)                        # unscreened pool rows -> identical
    ps = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", bundle=_poison_rows(dev_b, S), **kw)
    assert not np.array_equal(ps.scores["E_fixed"], base)                    # screened rows -> change


# ----------------------------------------------------------------------------- trajectories
def _traj(sel, max_n, cps=(), preds=(), r=0, s0=(), b=None):
    seeds = config.trajectory_seeds("dev", r + 1)
    return H.Trajectory(sel, "dev1", seeds[r], r, checkpoints=cps, max_n=max_n, predictors=preds, s0=s0,
                        batch_size=b, save=False, ledger=False, write_log=False).run()


def test_random_selector_nested_prefixes():
    a = _traj(H.RandomSelector(), 30)["selected"]
    b = _traj(H.RandomSelector(), 60)["selected"]
    assert b[:30] == a and len(set(b)) == 60
    c = _traj(H.RandomSelector(), 30, r=1)["selected"]
    assert c != a                                                            # different trajectory seed


def test_adding_checkpoint_does_not_change_trajectory():
    p = (("stack_lr", {}),)
    a = _traj(AdaptiveToySelector(), 40, cps=(20,), preds=p)
    b = _traj(AdaptiveToySelector(), 40, cps=(20, 25, 40), preds=p)
    c = _traj(AdaptiveToySelector(), 40)
    assert a["selected"] == b["selected"] == c["selected"]
    assert a["checkpoints"]["20"]["stack_lr"]["metrics"] == b["checkpoints"]["20"]["stack_lr"]["metrics"]
    assert set(b["checkpoints"]) == {"20", "25", "40"}                         # 25: first 5 of the 3rd batch


def test_partial_batch_checkpoint_uses_greedy_order(monkeypatch):
    seen = []
    real = H.run_predictor

    def spy(spec, *a, **k):
        seen.append(spec.screened)
        return real(spec, *a, **k)
    monkeypatch.setattr(H, "run_predictor", spy)
    log = _traj(AdaptiveToySelector(), 30, cps=(25,), preds=(("stack_lr", {}),))
    assert list(seen[0]) == log["selected"][:25]


def test_randomised_greedy():
    g = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    assert H.randomised_greedy(g, 3, None, deterministic=True) == [0, 1, 2]
    picks = {tuple(H.randomised_greedy(g, 1, np.random.default_rng(s), top_k=2)) for s in range(50)}
    assert picks == {(0,), (1,)}
    a = H.randomised_greedy(g, 4, np.random.default_rng(1), top_k=3)
    assert a == H.randomised_greedy(g, 4, np.random.default_rng(1), top_k=3) and len(set(a)) == 4


def test_trajectory_seed_streams():
    ss = config.trajectory_seeds("dev", 3)
    fs = config.trajectory_seeds("final", 3)
    m = [H.model_seed_of(H.stream_seeds(s)[1]) for s in ss]
    assert len(set(m)) == 3 and m == [H.model_seed_of(H.stream_seeds(s)[1]) for s in config.trajectory_seeds("dev", 3)]
    assert not set(m) & {H.model_seed_of(H.stream_seeds(s)[1]) for s in fs}


# ----------------------------------------------------------------------------- slow end-to-end controls
@pytest.mark.slow
def test_end_to_end_controls_dev1_n100(dev_b):
    """Shuffled screened labels (ε, P, labels permuted jointly over the screened block; array SMF permuted) ->
    AUPRC ≈ prevalence; oracle back door (true target labels added to stack_lr scores) -> AUPRC ≈ 1."""
    from hsl import upstream as U
    S = _S(100, seed=2)
    rng = np.random.default_rng(9)
    r = dev_b.gene_rows(S)
    E_, SP, Y = (np.array(getattr(dev_b, k), copy=True) for k in ("E", "SP", "Y"))
    blk = np.ix_(r, np.arange(dev_b.n_arrays))
    e, sp_, y = E_[blk].ravel(), SP[blk].ravel(), Y[blk].ravel()
    obs = np.flatnonzero(y >= 0)
    perm = obs[rng.permutation(len(obs))]
    for M, v in ((E_, e), (SP, sp_), (Y, y)):
        vv = v.copy()
        vv[obs] = v[perm]
        M[blk] = vv.reshape(len(r), -1)
    sb = dev_b.with_matrices(E=E_, SP=SP, Y=Y)
    sb.array_smf = np.array(dev_b.array_smf)[rng.permutation(dev_b.n_arrays)]
    res = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", seeds={"model_seed": 1}, bundle=sb, save=False,
                          ledger=False)
    m = res.metrics["E_fixed_dev"]
    print("shuffled", m["auprc"], m["prevalence"])
    assert 0.7 * m["prevalence"] < m["auprc"] < 1.3 * m["prevalence"]

    def oracle_backdoor(view, b):
        S_, _ = U.run_method(view, "stack_lr", model_seed=1, n_threads=3)
        Yt = b.Y[b.gene_rows(view.target_orfs)].astype(np.float32)
        return S_ + 10.0 * np.clip(Yt, 0, 1)
    res = H.run_predictor(H.RunSpec("dev1", S), "stack_lr", seeds={"model_seed": 1}, score_fn=oracle_backdoor,
                          save=False, ledger=False)
    print("oracle", res.metrics["E_fixed_dev"]["auprc"])
    assert res.metrics["E_fixed_dev"]["auprc"] > 0.99


@pytest.mark.slow
@pytest.mark.parametrize("method", ["stack_lr", "stack_lr_wr2", "m1_hgb"])
def test_predictors_run_with_zero_sl_screened_rows_n10(dev_b, method):
    """n = 10 screened genes including 5 zero-SL rows (DESIGN v3 assertion); m1_hgb with prop_folds = 5."""
    sp = views.dev_split(1)
    cnt = views.query_sl_counts(dev_b, sp.pool)
    zero = sp.pool[cnt == 0][:5]
    other = sp.pool[cnt > 0][:5]
    assert len(zero) == 5
    S = tuple(np.sort(np.r_[zero, other]))
    res = H.run_predictor(H.RunSpec("dev1", S), method, {"drop": ["genome"]} if method == "m1_hgb" else None,
                          seeds={"model_seed": 1}, save=False, ledger=False)
    assert np.isfinite(res.scores["E_fixed"]).all()
