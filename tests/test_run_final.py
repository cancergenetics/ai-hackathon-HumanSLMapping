"""Final-evaluation machinery (DESIGN v3.1 §8, §10; hsl.final, scripts/run_final.py, scripts/make_frozen.py).

Never touches T: refusals are tested on throw-away git repositories and temporary output directories; the unblinded
loader is only exercised through its own refusal (no HSL_FINAL) or replaced by a mock. Crash-resume is tested on the
dev bundle with a stub predictor (no model fits)."""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate)

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest
import yaml

from hsl import config, eval as E, final as F, harness as H, views

ROOT = F.ROOT


def _script(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------------- plan
@pytest.fixture(scope="module")
def final_plan():
    cfg, P = F.mode_params("final")
    return cfg, P, F.build_plan(cfg, P)


def _units(plan):
    return [F.Unit(c, int(r), int(n), p, bool(e)) for j in plan["jobs"] for c, r, n, p, e in j["units"]]


def test_plan_matches_config(final_plan):
    cfg, P, plan = final_plan
    U = _units(plan)
    assert len(U) == len({u.uid for u in U}) == plan["n_units"]
    bud = config.budgets()
    q2 = [int(n) for n in bud["q2_checkpoints"]]
    R_s, R_r = cfg["R"]["strategy"], cfg["R"]["random"]
    for fid in cfg["finalists"]:                                   # confirmatory arms: R_s x Q2 checkpoints, P*
        got = sorted((u.traj, u.n) for u in U if u.cid == fid and u.pred == "pstar")
        assert got == sorted((r, n) for r in range(R_s) for n in q2)
    pool = P["pool_n"]
    grid = sorted({min(int(n), pool) for n in bud["random_grid"]})
    for n in grid:                                                 # Q1: every grid budget, R_r (pool max: R_by_n)
        R = cfg["random"]["R_by_n"].get(n, cfg["random"]["R_by_n"].get(2627, R_r)) if n == pool else R_r
        for pk in cfg["random"]["predictors"]:
            assert sorted(u.traj for u in U if u.cid == "R0" and u.n == n and u.pred == pk) == list(range(R))
    m1 = sorted({min(int(n), pool) for n in bud["m1_hgb_q1_grid"]})
    assert sorted({u.n for u in U if u.cid == "R0" and u.pred == "m1_hgb"}) == m1
    for fid in cfg["finalists"]:
        assert sorted({u.n for u in U if u.cid == fid and u.pred == "m1_hgb"}) == [160, 490, 660]   # P4 review
        assert len({u.traj for u in U if u.cid == fid and u.pred == "m1_hgb"}) == bud["m1_hgb_q2_trajectories"]
    # descriptive arms exist with the right roles; dev-only scenarios never planned
    roles = {c: d["role"] for c, d in plan["arms"].items()}
    assert roles["R1"] == "R1" and roles["H-static"] == "static" and roles["F-static_r16"] == "static"
    assert {c for c, r in roles.items() if r == "det"} == {f"{f}__det" for f in cfg["descriptive"]["deterministic"]}
    assert {u.traj for u in U if roles[u.cid] == "det"} == {cfg["deterministic_traj"]}
    ws = [c for c, r in roles.items() if r in ("ws", "ws_random")]
    assert all(u.n > bud["warm_start_n0"] for u in U if u.cid in ws)
    assert 440 in {u.n for u in U if u.cid == "R1"}                # R1-continuation reference of the warm start
    assert not any("feas" in c or "min" in c for c in roles)
    assert plan["arms"]["F_r16"]["knobs"] == {"r": 16} and plan["arms"]["F-static_r16"]["knobs"] == {"r": 16}
    assert plan["arms"]["C"]["knobs"] == {"normalize": False}
    assert set(cfg["descriptive"]["deterministic_caveat"]) == {"C", "F_r16"}
    assert "minimal_features" not in cfg["dev_only"] and "no_AP-MS" in cfg["dev_only"]


def test_jobs_are_single_trajectory_and_size_classed(final_plan):
    cfg, P, plan = final_plan
    for j in plan["jobs"]:
        assert {(c, int(r)) for c, r, *_ in j["units"]} == {(j["cid"], j["traj"])}
        assert {F.size_class(int(n)) for _, _, n, *_ in j["units"]} == {j["cls"]}
        assert plan["sel_max_n"][j["cid"]] >= j["n_max"]
    jobs = F.job_entries(plan, P)
    assert jobs[0]["name"] == f"final_{F.LABELS_JOB}"
    assert all(j["env"].get("HSL_FINAL") == "1" for j in jobs)
    _, Pd = F.mode_params("dryrun")
    jd = F.job_entries(F.build_plan(cfg, Pd), Pd)
    assert not any("HSL_FINAL" in j["env"] for j in jd)
    assert all(j["env"]["OMP_NUM_THREADS"] == "2" for j in jd)


def test_plan_hash_is_deterministic(final_plan):
    cfg, P, plan = final_plan
    assert F.build_plan(cfg, P)["sha256"] == plan["sha256"]


# ----------------------------------------------------------------------------- seeds
def test_seed_streams_disjoint_and_crn():
    cfg = F.load_cfg()
    s = config.seeds()
    keys = [cfg["seeds"][m][k] for m in ("final", "dryrun") for k in ("trajectories", "bootstrap")]
    vals = [s[k] for k in keys] + [s["dev_traj_base_seed"], s["bootstrap_seed"]]
    assert len(set(vals)) == len(vals)                              # final, dry-run and dev bases all distinct
    _, P = F.mode_params("final")
    t3 = F.traj_seed(P, 3)
    sel, ev = H.stream_seeds(t3)
    a = np.random.default_rng(sel).integers(0, 2**31, 8)
    b = np.random.default_rng(ev).integers(0, 2**31, 8)
    assert not np.array_equal(a, b)                                 # selection and evaluation streams separate
    assert F.traj_seed(P, 3).entropy == config.seed("final_traj_base_seed") and F.traj_seed(P, 3).spawn_key == (3,)
    assert H.model_seed_of(H.stream_seeds(F.traj_seed(P, 3))[1]) == H.model_seed_of(ev)   # CRN: arm-independent
    dev = config.trajectory_seeds("dev", 4)[3]
    assert H.model_seed_of(H.stream_seeds(dev)[1]) != H.model_seed_of(ev)


# ----------------------------------------------------------------------------- refusals on a throw-away repo
def _g(root, *a):
    r = subprocess.run(["git", "-C", str(root), "-c", "user.email=t@t", "-c", "user.name=t", *a],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


SPEC = {"design_s8": ["src/*.py", "config/protocol.yaml"], "final_plan": [F.FINALISTS], "extended": ["scripts/*.py"]}


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    r = tmp_path / "repo"
    for p, t in {"src/a.py": "x = 1\n", "config/protocol.yaml": "a: 1\n", F.FINALISTS: "status: test\n",
                 "scripts/run.py": "print(1)\n", ".gitignore": "/out/\n"}.items():
        (r / p).parent.mkdir(parents=True, exist_ok=True)
        (r / p).write_text(t)
    _g(r, "init", "-q")
    monkeypatch.setattr(F, "FROZEN_SPEC", SPEC)
    return r


def test_refusal_finalists_not_committed(repo):
    assert not F.check_finalists_committed(repo).ok
    _g(repo, "add", "-A")
    _g(repo, "commit", "-qm", "init")
    assert F.check_finalists_committed(repo).ok
    (repo / F.FINALISTS).write_text("status: edited\n")
    c = F.check_finalists_committed(repo)
    assert not c.ok and "uncommitted" in c.detail


def test_refusal_dirty_tree_but_ignored_outputs_allowed(repo):
    _g(repo, "add", "-A")
    _g(repo, "commit", "-qm", "init")
    assert F.check_tree_clean(repo).ok
    (repo / "out").mkdir()
    (repo / "out" / "big.npz").write_text("ignored")
    assert F.check_tree_clean(repo).ok                              # ignored outputs allowed
    (repo / "notes.txt").write_text("untracked")
    assert not F.check_tree_clean(repo).ok                          # untracked counts as dirty
    (repo / "notes.txt").unlink()
    (repo / "src" / "a.py").write_text("x = 2\n")
    assert not F.check_tree_clean(repo).ok                          # modified tracked file


def test_refusal_frozen_missing_mismatch_added(repo):
    mf = _script("make_frozen.py")
    _g(repo, "add", "-A")
    _g(repo, "commit", "-qm", "init")
    c = F.check_frozen(repo)
    assert not c.ok and "missing" in c.detail
    mf.make(repo, repo / F.FROZEN)
    assert not F.check_frozen(repo).ok                              # written but not committed
    _g(repo, "add", F.FROZEN)
    _g(repo, "commit", "-qm", "freeze")
    c = F.check_frozen(repo)
    assert c.ok, c.detail
    fz = yaml.safe_load((repo / F.FROZEN).read_text())
    assert set(fz["files"]) == set(SPEC) and "src/a.py" in fz["files"]["design_s8"]
    (repo / "src" / "a.py").write_text("x = 3\n")                   # change a frozen file (even if committed)
    _g(repo, "commit", "-qam", "edit")
    c = F.check_frozen(repo)
    assert not c.ok and "changed src/a.py" in c.detail
    _g(repo, "revert", "--no-edit", "HEAD")
    assert F.check_frozen(repo).ok
    (repo / "src" / "b.py").write_text("y = 1\n")                   # new file under a frozen glob
    _g(repo, "add", "-A")
    _g(repo, "commit", "-qm", "add")
    c = F.check_frozen(repo)
    assert not c.ok and "added src/b.py" in c.detail


def test_make_frozen_refuses_dirty_or_uncommitted_plan(repo):
    mf = _script("make_frozen.py")
    with pytest.raises(SystemExit):
        mf.make(repo, repo / F.FROZEN)                              # finalists.yaml not committed
    _g(repo, "add", "-A")
    _g(repo, "commit", "-qm", "init")
    (repo / "x.txt").write_text("dirty")
    with pytest.raises(SystemExit):
        mf.make(repo, repo / F.FROZEN)                              # dirty tree


def test_refusal_lock_and_exclusive_create(tmp_path):
    lock = tmp_path / "final" / ".lock"
    assert F.check_lock_absent(lock).ok
    F.write_lock(lock, {"git_head": "abc"})
    assert not F.check_lock_absent(lock).ok
    with pytest.raises(FileExistsError):
        F.write_lock(lock, {"git_head": "def"})                     # never overwritten
    assert F.read_lock(lock)["git_head"] == "abc"


def _redirect(monkeypatch, tmp_path, mode="final"):
    real = F.mode_params

    def mp(m, root=ROOT):
        cfg, P = real(m, root)
        out = tmp_path / m
        P.update(out=out, run=out / "run", lock=out / ".lock")
        return cfg, P
    monkeypatch.setattr(F, "mode_params", mp)
    return mp(mode)[1]


def test_start_refuses_in_real_repo_state_and_writes_no_lock(monkeypatch, tmp_path):
    rf = _script("run_final.py")
    P = _redirect(monkeypatch, tmp_path)
    fails = {"finalists_committed": False, "tree_clean": False, "frozen_hashes": False}
    for name in fails:                                              # each refusal alone blocks the start
        for nm, fn in (("finalists_committed", "check_finalists_committed"), ("tree_clean", "check_tree_clean"),
                       ("frozen_hashes", "check_frozen")):
            monkeypatch.setattr(F, fn, (lambda nm_=nm: (lambda root=ROOT: F.Check(nm_, nm_ != name, "test")))())
        with pytest.raises(F.Refusal, match=name):
            rf.start("final")
        assert not P["lock"].exists()
    for fn, nm in (("check_finalists_committed", "a"), ("check_tree_clean", "b"), ("check_frozen", "c")):
        monkeypatch.setattr(F, fn, (lambda nm_=nm: (lambda root=ROOT: F.Check(nm_, True, "")))())
    lk = rf.start("final")                                          # all pass -> lock + plan + jobs
    assert P["lock"].exists() and lk["git_head"] == views.git_state()["head"]
    plan = json.loads((P["run"] / "plan.json").read_text())
    assert plan["sha256"] == lk["plan_sha256"]
    jobs = json.loads((P["run"] / "jobs.json").read_text())
    assert jobs[0]["name"].endswith(F.LABELS_JOB) and all(j["env"]["HSL_FINAL"] == "1" for j in jobs)
    with pytest.raises(F.Refusal, match="lock_absent"):
        rf.start("final")                                           # second start refused by the lock


def test_dryrun_start_enforces_only_the_lock(monkeypatch, tmp_path):
    rf = _script("run_final.py")
    P = _redirect(monkeypatch, tmp_path, "dryrun")
    monkeypatch.setattr(F, "check_tree_clean", lambda root=ROOT: F.Check("tree_clean", False, "dirty (test)"))
    lk = rf.start("dryrun")
    assert P["lock"].exists() and any(c["name"] == "tree_clean" and not c["ok"] for c in lk["preflight"])
    with pytest.raises(F.Refusal, match="lock_absent"):
        rf.start("dryrun")


def test_guard_final_and_unblinded_loader(monkeypatch, tmp_path):
    rf = _script("run_final.py")
    P = _redirect(monkeypatch, tmp_path)
    monkeypatch.delenv("HSL_FINAL", raising=False)
    with pytest.raises(F.Refusal, match="HSL_FINAL"):
        F.guard_final(P)
    with pytest.raises(PermissionError):                            # the guard's own refusal (no HSL_FINAL)
        views.final_bundle()
    monkeypatch.setenv("HSL_FINAL", "1")
    with pytest.raises(F.Refusal, match="no lock"):
        F.guard_final(P)
    F.write_lock(P["lock"], {"git_head": "0" * 40, "plan_sha256": "x"})
    monkeypatch.setattr(views, "git_state", lambda: {"head": "1" * 40, "dirty": False, "diff_sha256": ""})
    with pytest.raises(F.Refusal, match="locked commit"):
        F.guard_final(P)
    monkeypatch.setattr(views, "git_state", lambda: {"head": "0" * 40, "dirty": True, "diff_sha256": ""})
    with pytest.raises(F.Refusal, match="dirty"):
        F.guard_final(P)
    calls = []
    monkeypatch.setattr(views, "final_bundle", lambda: calls.append(1) or "MOCK_BUNDLE")
    with pytest.raises(F.Refusal):
        rf.get_bundle(P)                                            # guard fails -> loader never called
    assert calls == []
    monkeypatch.setattr(views, "git_state", lambda: {"head": "0" * 40, "dirty": False, "diff_sha256": ""})
    assert rf.get_bundle(P) == "MOCK_BUNDLE" and calls == [1]
    with pytest.raises(F.Refusal, match="plan differs"):
        F.guard_final(P, {"sha256": "y"})
    monkeypatch.delenv("HSL_FINAL")
    _, Pd = F.mode_params("dryrun")
    assert rf.get_bundle(Pd) is views.dev_bundle()                  # dry run: dev bundle only


# ----------------------------------------------------------------------------- crash-resume (stub predictor)
def test_crash_resume_never_reselects(monkeypatch, tmp_path, dev_b):
    P = _redirect(monkeypatch, tmp_path, "dryrun")
    cfg, P = F.mode_params("dryrun")
    plan = F.build_plan(cfg, P)
    job = next(j for j in plan["jobs"] if j["cid"] == "R0" and j["group"] == "main" and j["cls"] == "small")
    calls = {"n": 0}

    def stub(P_, u, sel, bundle, crash_at=2):
        calls["n"] += 1
        if calls["n"] == crash_at:
            raise RuntimeError("simulated crash")
        d = P_["run"] / "predictions" / "stub"
        d.mkdir(parents=True, exist_ok=True)
        f = d / f"{u.uid.replace('|', '_')}.npz"
        np.savez(f, E_fixed=np.zeros(3, np.float32))
        return {"uid": u.uid, "cid": u.cid, "traj": u.traj, "n": u.n, "pred": u.pred,
                "selected_head": sel["selected"][: u.n][:3], "pred_file": str(f.relative_to(ROOT))
                if str(f).startswith(str(ROOT)) else str(f)}
    monkeypatch.setattr(F, "ROOT", Path("/"))                       # stub files live in tmp (absolute paths)
    counter = {}
    with pytest.raises(RuntimeError, match="simulated crash"):
        F.work(P, plan, job["name"], dev_b, runner=stub, counter=counter)
    sel_file = F.sel_path(P, "R0", job["traj"])
    assert sel_file.exists() and counter == {"R0": 1}              # selection saved before any evaluation
    mtime = sel_file.stat().st_mtime_ns
    rows1 = F.read_parts(P, [job["name"]])
    assert len(rows1) == 1
    res = F.work(P, plan, job["name"], dev_b, runner=lambda *a: stub(*a, crash_at=-1), counter=counter)
    assert counter == {"R0": 1} and sel_file.stat().st_mtime_ns == mtime     # loaded, never recomputed
    rows = F.read_parts(P, [job["name"]])
    assert res["ran"] == len(job["units"]) - 1 and len(rows) == len(job["units"])
    assert len({x["uid"] for x in rows}) == len(rows)               # no duplicates
    assert [x["selection_loaded"] for x in rows] == [False] + [True] * (len(rows) - 1)
    again = F.work(P, plan, job["name"], dev_b, runner=stub)
    assert again == {"ran": 0, "skipped": len(job["units"])}
    sel = json.loads(sel_file.read_text())                          # = harness RandomSelector permutation (seeded)
    rs = H.RandomSelector()
    rs.reset(views.get_split(P["split"]).pool, np.random.default_rng(H.stream_seeds(F.traj_seed(P, job["traj"]))[0]))
    assert sel["selected"] == rs._perm[: len(sel["selected"])]
    ev = [json.loads(x) for x in (P["run"] / "LEDGER_final.jsonl").read_text().splitlines()]
    assert [e["event"] for e in ev].count("job_start") == 3


# ----------------------------------------------------------------------------- labels / engine
def test_chromosome_rule_matches_pair_store():
    assert list(F.chromosome_of(["YAL001C", "YPR200C", "YLR312C-B", "Q0045", "YBR020W"])) == ["A", "P", "L", "", "B"]
    pdir = Path(hsl.paths.yeast_sls_root()) / "data" / "processed" / "features" / "pairs"
    genes = np.array(pdir.joinpath("genes.txt").read_text().split())
    M = np.load(pdir / "gen_same_chr.npy", mmap_mode="r")
    rng = np.random.default_rng(0)
    i, j = rng.integers(0, len(genes), 4000), rng.integers(0, len(genes), 4000)
    ci, cj = F.chromosome_of(genes[i]), F.chromosome_of(genes[j])
    m = np.asarray(M[i, j])                                         # -1 = coordinates missing in the store
    ok = (ci != "") & (cj != "") & (i != j) & (m >= 0)
    assert ok.mean() > 0.95
    assert np.array_equal((ci[ok] == cj[ok]), m[ok] == 1)


def test_boot_engine_matches_eval_and_two_way(tmp_path):
    rng = np.random.default_rng(1)
    nT, nA, ncell = 30, 50, 1200
    rows, cols = rng.integers(0, nT, ncell).astype(np.int32), rng.integers(0, nA, ncell).astype(np.int32)
    y = (rng.random(ncell) < 0.1).astype(np.int8)
    s = (rng.random(ncell) + y * 0.3).round(2).astype(np.float32)          # ties on purpose
    col_target = np.full(nA, -1)
    col_target[:5] = np.arange(5)
    run = tmp_path / "run"
    run.mkdir()
    np.savez(run / "labels.npz", y=y, rows=rows, cols=cols, targets=np.array([f"T{i}" for i in range(nT)]),
             AO=cols >= 5, TT=cols < 5, query_strata=np.arange(nT) % 3, on_array_rows=np.ones(nT, bool),
             lockbox_rows=np.arange(nT) > 2, same_chr=np.zeros(ncell, bool), col_target=col_target,
             chr_unresolved=np.zeros(ncell, bool))
    (run / "label_manifest.json").write_text(json.dumps({"eset": "toy"}))
    pdir = tmp_path / "pred"
    pdir.mkdir()
    np.savez(pdir / "x.npz", E_fixed=s)
    P = {"run": run, "B": 50, "B_two_way": 20, "seeds": {"bootstrap": "final_dryrun_bootstrap_seed"}}
    L = F.Labels(P)
    eng = F.BootEngine(P, L)
    old = F.ROOT
    F.ROOT = Path("/")
    try:
        v = eng.vector(str(pdir / "x.npz"))
        assert abs(v[0] - E.average_precision(y, s)) < 1e-12 and len(v) == 51
        b = 7
        assert abs(v[1 + b] - E.average_precision(y, s, eng.W[b][rows])) < 1e-12
        m = L.mask("lockbox")
        vl = eng.vector(str(pdir / "x.npz"), "lockbox")
        assert abs(vl[0] - E.average_precision(y[m], s[m])) < 1e-12
        tw = eng.two_way_vector(str(pdir / "x.npz"))
        Wq, V = eng.two_way_weights()
        assert np.array_equal(V[:, :5], Wq[:, :5])                 # target columns share the query weight
        assert abs(tw[0] - v[0]) < 1e-12
        assert abs(tw[1 + 3] - E.average_precision(y, s, Wq[3][rows] * V[3][cols])) < 1e-12
        x = {"traj": 0, "n": 10, "pred_file": str(pdir / "x.npz"), "auprc_inline": float(v[0])}
        assert F.check_inline([x], eng) == 1
        with pytest.raises(AssertionError):
            F.check_inline([{**x, "auprc_inline": float(v[0]) + 1e-6}], eng)
        A, est, tr = F.arm_matrix([x], [10, 20], eng)
        assert A.shape == (50, 1, 2) and np.isnan(est[0, 1]) and est[0, 0] == v[0]
    finally:
        F.ROOT = old


def test_final_code_never_calls_the_unblinded_loader_outside_run_final():
    pat = "final_bundle("
    assert pat not in (ROOT / "src" / "hsl" / "final.py").read_text()
    assert pat not in (ROOT / "scripts" / "score_final.py").read_text()
    assert pat not in (ROOT / "scripts" / "make_frozen.py").read_text()
    assert "views.final_bundle()" in (ROOT / "scripts" / "run_final.py").read_text()


def test_frozen_spec_covers_design_list():
    s8 = F.FROZEN_SPEC["design_s8"]
    for must in ("src/hsl/eval.py", "src/hsl/views.py", "src/hsl/guard.py", "src/hsl/harness.py",
                 "config/protocol.yaml", "config/split/*", "config/seeds.yaml", "config/budgets.yaml",
                 "config/predictor.yaml", "src/hsl/predictors/*.py", "third_party/ysl/PROVENANCE.md",
                 "config/feature_provenance.yaml", "DATA.md", "config/yeast_sls_commit.txt",
                 "results/data_checksums/checksums.json"):
        assert must in s8
    for grp, pats in F.FROZEN_SPEC.items():
        if grp == "final_plan":
            continue                                                  # finalists.yaml is committed at the freeze
        for p in pats:
            assert F.tracked_files(ROOT, p) or p in ("scripts/run_final.py", "scripts/score_final.py",
                                                      "src/hsl/final.py"), p
    assert os.path.exists(ROOT / ".gitignore")
    ign = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "results/final/.lock", "results/final/run/x"],
                         capture_output=True, text=True).stdout.split()
    assert ign == ["results/final/.lock", "results/final/run/x"]
