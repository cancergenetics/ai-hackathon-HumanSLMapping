"""Final evaluation machinery (DESIGN v3.1 §5, §6.2–§6.4, §7.1, §8, §10): plan, refusals, lock, selections,
evaluation units, labels and the replicate-AP engine. Driven by config/finalists.yaml; used by scripts/run_final.py
(the ONLY place that loads the unblinded bundle) and scripts/score_final.py (no bundle; reads the run directory).

Modes (finalists.yaml `modes`)
  final    split 'final' (targets T, pool D), out results/final, refusals enforced, HSL_FINAL=1 in every job.
  dryrun   pseudo-final on a dev split (dev1 stands in for T), reduced R / grids, out results/final_dryrun;
           refusals are evaluated and reported, only the lock is enforced. Same code path otherwise.

Phases (each resumable; nothing is ever re-selected)
  selection   one per (arm, trajectory): harness.Trajectory with NO predictors, run to the arm's maximum budget;
              the ordered selection (+ per-round timing and selector diagnostics) is written atomically to
              run/selections/<cid>_r<r>.json BEFORE any evaluation. An existing file is loaded, never recomputed.
  evaluation  one unit per (arm, trajectory, budget n, predictor): harness.run_predictor on the first n genes of
              the saved selection (model seed = the trajectory's evaluation stream), predictions SAVED (float32,
              run/predictions/) before the harness scores them; one JSON row per unit in run/parts/<job>.jsonl.
              A unit is done iff its row exists and its prediction file exists.
  scoring     scripts/score_final.py (separate step): every statistic recomputed from the saved predictions and
              run/labels.npz (written by the `_labels` job), with the replicate-AP engine below.
Working files live in run/ (gitignored), so the tree stays clean during the run (views.final_bundle requires it).
Seeds: trajectory r of every arm = SeedSequence(<trajectories key>).spawn(.)[r]; selection stream = child 0,
evaluation stream = child 1 (harness.stream_seeds); bootstrap = named children of SeedSequence(<bootstrap key>).
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import hashlib
import json
import os
import re
import socket
import subprocess
import time
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

import hsl.predictors  # noqa: F401  (registers the frozen P* method)
from hsl import config, eval as E, harness as H, p1, paths, views

ROOT = paths.PROJECT_ROOT
FINALISTS = "config/finalists.yaml"
FROZEN = "config/frozen.yaml"
MODES = ("final", "dryrun")
PY = ".venv/bin/python"
LABELS_JOB = "_labels"

# Freeze list (DESIGN §8) + the final plan + an EXTENDED set (code that produces the final numbers but that §8 does
# not name: selectors, features, upstream adapter, this module and the final scripts, the vendored predictor code).
FROZEN_SPEC = {
    "design_s8": ["src/hsl/eval.py", "src/hsl/views.py", "src/hsl/guard.py", "src/hsl/harness.py",
                  "config/protocol.yaml", "config/split/*", "config/seeds.yaml", "config/budgets.yaml",
                  "config/predictor.yaml", "src/hsl/predictors/*.py", "third_party/ysl/PROVENANCE.md",
                  "config/feature_provenance.yaml", "DATA.md", "config/yeast_sls_commit.txt",
                  "results/data_checksums/checksums.json"],
    "final_plan": [FINALISTS],
    "extended": ["config/strategies.yaml", "config/paths.yaml", "src/hsl/__init__.py", "src/hsl/config.py",
                 "src/hsl/paths.py", "src/hsl/features.py", "src/hsl/upstream.py", "src/hsl/p1.py",
                 "src/hsl/q2dev.py", "src/hsl/final.py", "src/hsl/strategies/*.py", "scripts/run_final.py",
                 "scripts/score_final.py", "scripts/launch.py", "third_party/ysl/**/*.py"],
}


class Refusal(RuntimeError):
    """A final-evaluation precondition is not met (DESIGN §10)."""


# ============================================================================ config / modes
def load_cfg(root: Path = ROOT) -> dict:
    return yaml.safe_load((Path(root) / FINALISTS).read_text())


def _budget(v):
    """A budgets.yaml key (dotted path allowed) or a literal list / int."""
    if isinstance(v, str):
        d = config.budgets()
        for k in v.split("."):
            d = d[k]
        return d
    return v


def predictor_table(cfg: dict) -> dict:
    pred = yaml.safe_load((ROOT / "config" / "predictor.yaml").read_text())
    out = {}
    for k, v in cfg["predictors"].items():
        if "from" in v:
            d = pred[v["from"]]
            out[k] = (str(d["method"]), dict(d.get("params") or {}))
        else:
            out[k] = (str(v["method"]), dict(v.get("params") or {}))
    return out


@dataclass(frozen=True)
class Arm:
    cid: str
    role: str                 # random | finalist | R1 | static | det | ws | ws_random | batch
    selector: str             # hsl.strategies registry name, or "R0"
    knobs: tuple = ()         # sorted (key, value) overrides of config/strategies.yaml
    deterministic: bool = False
    batch_size: int | None = None
    warm: int = 0
    base: str = ""            # finalist the arm derives from

    @property
    def knob_dict(self) -> dict:
        return dict(self.knobs)


@dataclass(frozen=True, order=True)
class Unit:
    cid: str
    traj: int
    n: int
    pred: str
    e_unobs: bool = False

    @property
    def uid(self) -> str:
        return f"{self.cid}|r{self.traj}|n{self.n}|{self.pred}"


@dataclass
class Job:
    name: str
    cid: str
    traj: int
    group: str
    cls: str
    units: list = field(default_factory=list)

    @property
    def n_max(self) -> int:
        return max(u.n for u in self.units) if self.units else 0


def resolve(cfg: dict, mode: str, root: Path = ROOT) -> dict:
    """Concrete parameters of a mode (budgets resolved, dry-run overrides applied, pool-max mapping)."""
    if mode not in MODES:
        raise ValueError(mode)
    m = cfg["modes"][mode]
    ov = m.get("overrides") or {}
    R = ov.get("R") or cfg["R"]
    split = views.get_split(m["split"])
    pool_n = int(len(split.pool))
    cap = lambda ns: sorted({min(int(n), pool_n) for n in ns})  # noqa: E731  (n above the pool -> pool max)
    sec = cfg["secondary"]
    ws = cfg["descriptive"]["warm_start"]
    pr = config.protocol()
    out = Path(root) / m["out"]
    P = {
        "mode": mode, "split": m["split"], "kind": split.kind, "out": out, "run": out / "run",
        "lock": out / ".lock", "enforce": bool(m["enforce_refusals"]), "threads": m.get("threads"),
        "pool_n": pool_n, "n_targets": int(len(split.targets)),
        "R_s": int(R["strategy"]), "R_r": int(R["random"]),
        "q2": cap(ov.get("q2_checkpoints") or _budget("q2_checkpoints")),
        "random_grid": cap(ov.get("random_grid") or _budget(cfg["random"]["grid"])),
        "random_R_by_n": {min(int(k), pool_n): int(v) for k, v in
                          (ov["random_R_by_n"] if "random_R_by_n" in ov else cfg["random"].get("R_by_n") or {}).items()},
        "random_unobs": cap(ov.get("q2_checkpoints") or _budget(cfg["random"]["e_unobs_at"])),
        "ws_n0": int(_budget(ws["n0"])),
        "ws_cps": cap(ov.get("warm_start_checkpoints") or _budget(ws["checkpoints"])),
        "batch_cps": cap(ov.get("batch_checkpoints") or cfg["descriptive"]["batch"]["checkpoints"]),
        "m1_random": {"cps": cap((ov.get("m1_random") or {}).get("checkpoints")
                                 or _budget(sec["m1_hgb"]["random"]["checkpoints"])),
                      "R": int((ov.get("m1_random") or {}).get("R") or _budget(sec["m1_hgb"]["random"]["R"]))},
        "m1_strat": {"cps": cap((ov.get("m1_strategies") or {}).get("checkpoints")
                                or _budget(sec["m1_hgb"]["strategies"]["checkpoints"])),
                     "R": int((ov.get("m1_strategies") or {}).get("R") or _budget(sec["m1_hgb"]["strategies"]["R"]))},
        "ref": {"cps": cap((ov.get("reference") or {}).get("checkpoints") or _budget(sec["reference"]["checkpoints"])),
                "R": int((ov.get("reference") or {}).get("R") or sec["reference"]["R"]),
                "arms": list(sec["reference"]["arms"])},
        "B": int(ov.get("B_query_bootstrap") or pr["bootstrap"]["B_query_bootstrap"]),
        "B_two_way": int(ov.get("B_two_way") or pr["bootstrap"]["B_two_way"]),
        "seeds": dict(cfg["seeds"][mode]),
        "preds": predictor_table(cfg),
        "c1": dict(config.budgets()["contrast_c1"]), "brief": dict(config.budgets()["brief_example"]),
        "e_unobs_save": bool(cfg["e_unobs"].get("save", True)),
    }
    P["ws_cps"] = [n for n in P["ws_cps"] if n > P["ws_n0"]]
    return P


def check_global_knobs(cfg: dict) -> None:
    """The finalists were tuned with these global selector settings; refuse if config drifted."""
    from hsl.strategies import base
    exp = cfg["expected_global_knobs"]
    got = {"offarray_calibration": base.knobs().get("offarray_calibration"),
           "randomised_greedy_top_k": config.protocol()["selector"]["randomised_greedy_top_k"],
           "cluster_cap_per_batch": config.protocol()["selector"]["cluster_cap_per_batch"]}
    bad = {k: (exp[k], got[k]) for k in exp if exp[k] != got[k]}
    if bad:
        raise Refusal(f"global selector settings differ from finalists.yaml expected_global_knobs: {bad}")


# ============================================================================ arms / plan
def _static_cid(fid: str, fin_sel: str, static_sel: str) -> str:
    return static_sel + fid[len(fin_sel):] if fid.startswith(fin_sel) else f"{fid}-static"


def arms(cfg: dict, P: dict) -> dict:
    out = {"R0": Arm("R0", "random", "R0"), "R1": Arm("R1", "R1", cfg["descriptive"]["R1"]["selector"])}
    fins = cfg["finalists"]
    for fid, d in fins.items():
        out[fid] = Arm(fid, "finalist", d["selector"], tuple(sorted((d.get("knobs") or {}).items())), base=fid)
    for fid, d in (cfg["descriptive"].get("static") or {}).items():
        if d:
            cid = _static_cid(fid, fins[fid]["selector"], d["selector"])
            out[cid] = Arm(cid, "static", d["selector"], tuple(sorted((d.get("knobs") or {}).items())), base=fid)
    for fid in cfg["descriptive"].get("deterministic") or []:
        a = out[fid]
        out[f"{fid}__det"] = Arm(f"{fid}__det", "det", a.selector, a.knobs, deterministic=True, base=fid)
    ws = cfg["descriptive"]["warm_start"]
    n0 = P["ws_n0"]
    for fid in ws.get("finalists") or []:
        a = out[fid]
        out[f"{fid}__ws{n0}"] = Arm(f"{fid}__ws{n0}", "ws", a.selector, a.knobs, warm=n0, base=fid)
    if ws.get("random_continuation") and ws.get("finalists"):
        out[f"R0__ws{n0}"] = Arm(f"R0__ws{n0}", "ws_random", "R0", warm=n0)
    bt = cfg["descriptive"]["batch"]
    for fid in bt.get("finalists") or []:
        a = out[fid]
        for b in bt["sizes"]:
            out[f"{fid}__b{b}"] = Arm(f"{fid}__b{b}", "batch", a.selector, a.knobs, batch_size=int(b), base=fid)
    return out


def units(cfg: dict, P: dict, A: dict) -> list[Unit]:
    U = set()
    R_s, R_r = P["R_s"], P["R_r"]
    for r in range(R_r):                                              # Q1 random curves (+ Q2 random arm)
        for n in P["random_grid"]:
            if r < P["random_R_by_n"].get(n, R_r):
                for pk in cfg["random"]["predictors"]:
                    U.add(Unit("R0", r, n, pk, pk == "pstar" and n in P["random_unobs"]))
    for r in range(P["m1_random"]["R"]):
        for n in P["m1_random"]["cps"]:
            U.add(Unit("R0", r, n, "m1_hgb"))
    fins = list(cfg["finalists"])
    for cid, a in A.items():
        if a.role == "random":
            continue
        trajs = [int(cfg["deterministic_traj"])] if a.role == "det" else range(R_s)
        cps = {"ws": P["ws_cps"], "ws_random": P["ws_cps"], "batch": P["batch_cps"]}.get(a.role, P["q2"])
        if a.role == "R1":
            cps = sorted(set(cps) | {min(int(n), P["pool_n"]) for n in cfg["descriptive"]["R1"].get("extra_checkpoints", [])})
        for r in trajs:
            for n in cps:
                U.add(Unit(cid, r, n, "pstar", True))
    for fid in fins:
        for r in range(P["m1_strat"]["R"]):
            for n in P["m1_strat"]["cps"]:
                U.add(Unit(fid, r, n, "m1_hgb"))
    ref_arms = [x for a in P["ref"]["arms"] for x in (fins if a == "finalists" else [a])]
    for cid in ref_arms:
        for r in range(P["ref"]["R"]):
            for n in P["ref"]["cps"]:
                U.add(Unit(cid, r, n, "ref"))
    return sorted(U)


def size_class(n: int) -> str:
    c = config.protocol()["compute"]
    if n >= int(c["alone_min_n"]):
        return "alone"
    return "small" if n <= int(c["small_job_max_n"]) else "large"


_GROUP = {"pstar": "main", "B1_array_rate": "main", "B2_product": "main", "stack_lr": "main", "m1_hgb": "m1",
          "ref": "ref"}
_ROLE_RANK = {"finalist": 0, "random": 1, "R1": 2, "static": 3, "det": 4, "ws": 5, "ws_random": 6, "batch": 7}
_CLS_RANK = {"small": 0, "large": 1, "alone": 2}
_GROUP_RANK = {"main": 0, "ref": 1, "m1": 2}


def build_plan(cfg: dict, P: dict) -> dict:
    A = arms(cfg, P)
    U = units(cfg, P, A)
    jobs: dict = {}
    for u in U:
        g, c = _GROUP.get(u.pred, "main"), size_class(u.n)
        nm = f"{u.cid}_r{u.traj}_{g}_{c}"
        jobs.setdefault(nm, Job(nm, u.cid, u.traj, g, c)).units.append(u)
    order = sorted(jobs.values(), key=lambda j: (_CLS_RANK[j.cls], _GROUP_RANK[j.group], _ROLE_RANK[A[j.cid].role],
                                                  j.cid, j.traj))
    sel_max = {}
    for u in U:
        sel_max[u.cid] = max(sel_max.get(u.cid, 0), u.n)
    if any(a.warm for a in A.values()):
        sel_max["R1"] = max(sel_max.get("R1", 0), P["ws_n0"])
    plan = {"mode": P["mode"], "split": P["split"], "pool_n": P["pool_n"], "R_s": P["R_s"], "R_r": P["R_r"],
            "arms": {cid: {"role": a.role, "selector": a.selector, "knobs": a.knob_dict,
                           "deterministic": a.deterministic, "batch_size": a.batch_size, "warm": a.warm,
                           "base": a.base} for cid, a in A.items()},
            "sel_max_n": sel_max, "predictors": {k: list(v) for k, v in P["preds"].items()},
            "jobs": [{"name": j.name, "cid": j.cid, "traj": j.traj, "group": j.group, "cls": j.cls,
                      "n_max": j.n_max, "units": [[u.cid, u.traj, u.n, u.pred, u.e_unobs] for u in j.units]}
                     for j in order]}
    plan["n_units"] = len(U)
    plan["sha256"] = hashlib.sha256(json.dumps({k: v for k, v in plan.items() if k != "sha256"}, sort_keys=True,
                                               default=str).encode()).hexdigest()
    return plan


def plan_arms(plan: dict) -> dict:
    return {cid: Arm(cid, d["role"], d["selector"], tuple(sorted(d["knobs"].items())), d["deterministic"],
                     d["batch_size"], int(d["warm"]), d["base"]) for cid, d in plan["arms"].items()}


def plan_job(plan: dict, name: str) -> Job:
    for j in plan["jobs"]:
        if j["name"] == name:
            return Job(j["name"], j["cid"], int(j["traj"]), j["group"], j["cls"],
                       [Unit(c, int(r), int(n), p, bool(e)) for c, r, n, p, e in j["units"]])
    raise KeyError(f"job {name!r} not in the plan")


def job_entries(plan: dict, P: dict, names=None) -> list[dict]:
    """launch.py job list: the labels job first, then the plan order. Final jobs carry HSL_FINAL=1."""
    env = {}
    if P["mode"] == "final":
        env["HSL_FINAL"] = "1"
    if P.get("threads"):                          # same variables as scripts/launch.py THREAD_VARS
        env.update({k: str(P["threads"]) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "HSL_THREADS")})
    mk = lambda nm, n: {"name": f"{P['mode']}_{nm}", "n": int(n), "env": dict(env),  # noqa: E731
                        "cmd": [PY, "scripts/run_final.py", "work", "--mode", P["mode"], "--job", nm]}
    out = [mk(LABELS_JOB, 0)] if names is None or LABELS_JOB in names else []
    out += [mk(j["name"], j["n_max"]) for j in plan["jobs"] if names is None or j["name"] in names]
    return out


# ============================================================================ refusals / frozen / lock
def _git(root, *args) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


def check_finalists_committed(root: Path = ROOT) -> Check:
    if not (Path(root) / FINALISTS).exists():
        return Check("finalists_committed", False, f"{FINALISTS} missing")
    if _git(root, "ls-files", "--error-unmatch", FINALISTS).returncode != 0:
        return Check("finalists_committed", False, f"{FINALISTS} is not tracked by git")
    st = _git(root, "status", "--porcelain", "--", FINALISTS).stdout.strip()
    if st:
        return Check("finalists_committed", False, f"{FINALISTS} has uncommitted changes: {st}")
    return Check("finalists_committed", True, _git(root, "log", "-1", "--format=%h", "--", FINALISTS).stdout.strip())


def check_tree_clean(root: Path = ROOT) -> Check:
    """`git status --porcelain` empty: no modified, staged or untracked (non-ignored) files. Ignored outputs allowed."""
    r = _git(root, "status", "--porcelain")
    if r.returncode != 0:
        return Check("tree_clean", False, f"git status failed: {r.stderr.strip()}")
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return Check("tree_clean", not lines, f"{len(lines)} dirty paths: " + "; ".join(lines[:8]) if lines else "clean")


def tracked_files(root: Path, pattern: str) -> list[str]:
    r = _git(root, "ls-files", "-z", "--", f":(glob){pattern}")
    return sorted(p for p in r.stdout.split("\0") if p)


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frozen_hashes(root: Path = ROOT, spec: dict | None = None) -> dict:
    spec = FROZEN_SPEC if spec is None else spec
    out = {}
    for grp, pats in spec.items():
        files = {}
        for pat in pats:
            got = tracked_files(root, pat)
            if not got:
                raise FileNotFoundError(f"frozen pattern {pat!r} matches no tracked file")
            for f in got:
                files[f] = sha256_file(Path(root) / f)
        out[grp] = dict(sorted(files.items()))
    return out


def check_frozen(root: Path = ROOT) -> Check:
    p = Path(root) / FROZEN
    if not p.exists():
        return Check("frozen_hashes", False, f"{FROZEN} missing (scripts/make_frozen.py writes it at the freeze)")
    if _git(root, "ls-files", "--error-unmatch", FROZEN).returncode != 0:
        return Check("frozen_hashes", False, f"{FROZEN} is not committed")
    fz = yaml.safe_load(p.read_text()) or {}
    try:
        now = frozen_hashes(root, fz.get("spec") or {})
    except FileNotFoundError as e:
        return Check("frozen_hashes", False, str(e))
    bad = []
    for grp, files in (fz.get("files") or {}).items():
        cur = now.get(grp, {})
        for f, h in files.items():
            if f not in cur:
                bad.append(f"missing {f}")
            elif cur[f] != h:
                bad.append(f"changed {f}")
        bad += [f"added {f}" for f in cur if f not in files]
    if not fz.get("files"):
        bad.append("no file hashes recorded")
    for grp in FROZEN_SPEC:
        if grp not in (fz.get("spec") or {}):
            bad.append(f"spec group {grp!r} absent")
    return Check("frozen_hashes", not bad, "; ".join(bad[:10]) if bad else f"{sum(map(len, now.values()))} files match")


def check_lock_absent(lock: Path) -> Check:
    return Check("lock_absent", not Path(lock).exists(), f"{lock} {'exists' if Path(lock).exists() else 'absent'}")


def preflight(root: Path, lock: Path) -> list[Check]:
    return [check_finalists_committed(root), check_tree_clean(root), check_frozen(root), check_lock_absent(lock)]


def write_lock(lock: Path, info: dict) -> dict:
    """Create the lock atomically (O_EXCL): a second start can never overwrite it."""
    lock = Path(lock)
    lock.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    with os.fdopen(fd, "w") as f:
        f.write(json.dumps(info, indent=1, sort_keys=True) + "\n")
    return info


def read_lock(lock: Path) -> dict | None:
    lock = Path(lock)
    return json.loads(lock.read_text()) if lock.exists() else None


def guard_final(P: dict, plan: dict | None = None) -> None:
    """Every final job: HSL_FINAL=1, lock present, HEAD == the locked commit, clean tree, same plan."""
    if P["mode"] != "final":
        return
    if os.environ.get("HSL_FINAL") != "1":
        raise Refusal("final jobs need HSL_FINAL=1 (set by the job file)")
    lk = read_lock(P["lock"])
    if lk is None:
        raise Refusal(f"no lock at {P['lock']}: run `scripts/run_final.py start --mode final` first")
    gs = views.git_state()
    if gs["dirty"]:
        raise Refusal("git tree is dirty")
    if gs["head"] != lk["git_head"]:
        raise Refusal(f"HEAD {gs['head'][:10]} != locked commit {lk['git_head'][:10]}")
    if plan is not None and plan["sha256"] != lk["plan_sha256"]:
        raise Refusal("plan differs from the locked plan")


def ledger_event(P: dict, event: str, **kw) -> None:
    rec = {"t": time.strftime("%Y-%m-%dT%H:%M:%S"), "event": event, "mode": P["mode"], "pid": os.getpid(),
           "host": socket.gethostname(), **kw}
    H.append_ledger(rec, P["run"] / "LEDGER_final.jsonl")


# ============================================================================ seeds
def traj_seed(P: dict, r: int) -> np.random.SeedSequence:
    return np.random.SeedSequence(config.seed(P["seeds"]["trajectories"])).spawn(int(r) + 1)[int(r)]


def boot_seed(P: dict, name: str) -> np.random.SeedSequence:
    return np.random.SeedSequence(config.seed(P["seeds"]["bootstrap"]), spawn_key=(zlib.crc32(name.encode()),))


# ============================================================================ selections
def make_selector(a: Arm):
    if a.selector == "R0":
        return H.RandomSelector()
    from hsl import strategies
    return strategies.make(a.selector, deterministic=a.deterministic, **a.knob_dict)


def sel_path(P: dict, cid: str, r: int) -> Path:
    return P["run"] / "selections" / f"{cid}_r{int(r)}.json"


def _atomic_write(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f".{p.name}.tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, p)


def ensure_selection(P: dict, plan: dict, cid: str, r: int, bundle, counter: dict | None = None) -> tuple[dict, bool]:
    """(selection record, loaded_from_disk). Existing selections are loaded and never recomputed."""
    p = sel_path(P, cid, r)
    need = int(plan["sel_max_n"][cid])
    if p.exists():
        rec = json.loads(p.read_text())
        if rec["cid"] != cid or int(rec["traj"]) != int(r) or len(rec["selected"]) < need:
            raise AssertionError(f"selection file {p} does not match the plan")
        return rec, True
    A = plan_arms(plan)
    a = A[cid]
    s0 = []
    if a.warm:
        r1, _ = ensure_selection(P, plan, "R1", r, bundle, counter)
        s0 = list(r1["selected"][:a.warm])
    sel = make_selector(a)
    t = H.Trajectory(sel, P["split"], traj_seed(P, r), int(r), arm=cid, batch_size=a.batch_size, checkpoints=(),
                     max_n=need, predictors=(), s0=tuple(s0), save=False, ledger=False, write_log=False,
                     bundle=bundle)
    t0 = time.perf_counter()
    log = t.run()
    wall = time.perf_counter() - t0
    if counter is not None:
        counter[cid] = counter.get(cid, 0) + 1
    rec = {"cid": cid, "traj": int(r), "role": a.role, "mode": P["mode"], "split": P["split"],
           "traj_id": log["traj_id"], "model_seed": int(log["model_seed"]), "batch_size": log["batch_size"],
           "s0_n": len(s0), "max_n": need, "selected": [str(g) for g in log["selected"]],
           "rounds": [{k: x[k] for k in ("round", "n_before", "view_s", "select_s")} for x in log["rounds"]],
           "selection_wall_s": wall, "selector_params": getattr(sel, "params", {}),
           "seeds": {"trajectories": P["seeds"]["trajectories"], "index": int(r)}, "git": log["git"],
           "diagnostics": getattr(sel, "diagnostics", [])}
    if len(rec["selected"]) < need:
        raise AssertionError(f"{cid} r{r}: selection stopped at {len(rec['selected'])} < {need}")
    _atomic_write(p, json.dumps(rec, default=H._json_default))
    return rec, False


# ============================================================================ evaluation units
def part_path(P: dict, job: str) -> Path:
    return P["run"] / "parts" / f"{job}.jsonl"


def read_parts(P: dict, jobs=None) -> list[dict]:
    d = P["run"] / "parts"
    files = sorted(d.glob("*.jsonl")) if jobs is None else [part_path(P, j) for j in jobs]
    return p1.read_rows(files)


def done_uids(P: dict, job: str) -> set:
    return {x["uid"] for x in p1.read_rows([part_path(P, job)])
            if x.get("pred_file") and (ROOT / x["pred_file"]).exists()}


_SL_TOTAL: dict = {}


def _pool_sl_total(split_name: str, bundle) -> int:
    """Strong SLs in the candidate pool's rows (denominator of 'fraction of SLs detected by the screened rows')."""
    k = (split_name, id(bundle))
    if k not in _SL_TOTAL:
        _SL_TOTAL[k] = int(views.query_sl_counts(bundle, views.get_split(split_name).cand_full).sum())
    return _SL_TOTAL[k]


def run_unit(P: dict, u: Unit, sel: dict, bundle) -> dict:
    method, params = P["preds"][u.pred]
    S = list(sel["selected"][:u.n])
    if len(S) != u.n:
        raise AssertionError(f"{u.uid}: selection shorter than n")
    ev = ("E_fixed", "E_unobs") if u.e_unobs else ("E_fixed",)
    sv = ev if (u.e_unobs and P["e_unobs_save"]) else ("E_fixed",)
    H.PRED_DIR = P["run"] / "predictions"
    with p1.PeakRSS() as pk:
        res = H.run_predictor(H.RunSpec(P["split"], tuple(S), ev, sv, label=f"{P['mode']}:{u.cid}"), method, params,
                              {"model_seed": int(sel["model_seed"]), "traj_index": int(u.traj)}, bundle=bundle,
                              parent=sel["traj_id"], save=True, ledger=False, tag=f"{P['mode']}:{u.cid}")
    es = res.esets["E_fixed"]
    m = res.metrics[es.name]
    preds = res.manifest.get("predictions", {})
    rt = res.manifest["runtime_s"]
    on = bundle.orf_to_col
    s0 = int(sel.get("s0_n", 0))
    used = [x for x in sel["rounds"] if x["n_before"] < u.n]
    k = int(views.query_sl_counts(bundle, np.asarray(S, dtype=str)).sum())
    row = {"uid": u.uid, "cid": u.cid, "traj": u.traj, "n": u.n, "pred": u.pred, "method": method, "params": params,
           "model_seed": int(sel["model_seed"]), "traj_id": sel["traj_id"], "run_id": res.run_id,
           "eset": es.name, "eset_digest": es.digest(), "n_cells": m["n_cells"], "n_pos": m["n_pos"],
           "auprc_inline": m["auprc"],
           "pred_file": preds.get("E_fixed", {}).get("file"),
           "pred_file_unobs": preds.get("E_unobs", {}).get("file"),
           "sl_detected": k, "sl_detected_frac": k / _pool_sl_total(P["split"], bundle),
           "offarray_share": float(np.mean([g not in on for g in S])),
           "offarray_share_new": float(np.mean([g not in on for g in S[s0:]])) if len(S) > s0 else None,
           "s0_n": s0, "batch_size": sel["batch_size"], "n_rounds": len(used),
           "select_s_cum": float(sum(x["select_s"] for x in used)),
           "view_s_cum": float(sum(x["view_s"] for x in used)),
           "fit_predict_s": rt.get("fit_predict_s"), "eval_total_s": rt.get("total_s"),
           "peak_rss_fit_gb": pk.peak, "n_threads": res.manifest["n_threads"], "git_head": res.manifest["git"]["head"]}
    if "E_unobs" in res.metrics:
        mu = res.metrics["E_unobs"]
        row["e_unobs"] = {k2: mu.get(k2) for k2 in ("n_cells", "n_pos", "prevalence", "auprc", "lift", "macro_auprc",
                                                     "auroc")}
    return row


def work(P: dict, plan: dict, job_name: str, bundle, runner=None, counter: dict | None = None) -> dict:
    """Run every not-yet-done unit of one job. Returns {'ran': k, 'skipped': m}."""
    guard_final(P, plan)
    job = plan_job(plan, job_name)
    done = done_uids(P, job.name)
    todo = [u for u in job.units if u.uid not in done]
    ledger_event(P, "job_start", job=job.name, n_units=len(job.units), n_todo=len(todo))
    if not todo:
        ledger_event(P, "job_end", job=job.name, ran=0, skipped=len(job.units))
        return {"ran": 0, "skipped": len(job.units)}
    sel, loaded = ensure_selection(P, plan, job.cid, job.traj, bundle, counter)
    runner = run_unit if runner is None else runner
    part = part_path(P, job.name)
    for u in todo:
        t0 = time.perf_counter()
        row = runner(P, u, sel, bundle)
        row.update(job=job.name, selection_loaded=bool(loaded), unit_wall_s=time.perf_counter() - t0,
                   t_done=time.strftime("%Y-%m-%dT%H:%M:%S"))
        p1.append_row(part, row)
        print(f"[final] {job.name}: {u.uid} done ({row['unit_wall_s']:.0f} s)", flush=True)   # no accuracy printed
    ledger_event(P, "job_end", job=job.name, ran=len(todo), skipped=len(job.units) - len(todo),
                 selection_loaded=bool(loaded))
    return {"ran": len(todo), "skipped": len(job.units) - len(todo)}


def progress(P: dict, plan: dict) -> dict:
    """Done / total units per job (rows whose prediction file exists)."""
    rows = read_parts(P)
    ok = {x["uid"] for x in rows if x.get("pred_file") and (ROOT / x["pred_file"]).exists()}
    out = {}
    for j in plan["jobs"]:
        uids = [Unit(c, int(r), int(n), p, bool(e)).uid for c, r, n, p, e in j["units"]]
        out[j["name"]] = (sum(u in ok for u in uids), len(uids))
    return out


# ============================================================================ labels (written by the `_labels` job)
_CHR = re.compile(r"^Y([A-P])[LR]\d{3}[WC]")


def chromosome_of(names) -> np.ndarray:
    """Chromosome letter A..P of a systematic yeast ORF name (current SGD name); '' if not a standard nuclear ORF."""
    out = []
    for g in np.asarray(names, dtype=str):
        m = _CHR.match(g)
        out.append(m.group(1) if m else "")
    return np.array(out, dtype="<U1")


def read_lockbox() -> np.ndarray:
    return np.array(sorted(paths.lockbox_file().read_text().split()), dtype=str)


def write_labels(P: dict, bundle) -> dict:
    """E_fixed cell labels + strata + bootstrap strata for the score step (label-only manifest, DESIGN §3)."""
    split = views.get_split(P["split"])
    es = H.get_eset(bundle, split, "E_fixed")
    tg = split.targets
    tcol = "TT" if split.kind == "final" else "VV"
    row_cur = views.current_of(bundle, tg)
    chr_r, chr_c = chromosome_of(row_cur), chromosome_of(bundle.array_orfs_current)
    same_chr = (chr_r[es.rows] != "") & (chr_r[es.rows] == chr_c[es.cols])
    col_target = np.full(bundle.n_arrays, -1, dtype=np.int64)
    for i, t in enumerate(tg):
        j = bundle.orf_to_col.get(t)
        if j is not None:
            col_target[j] = i
    lb = read_lockbox()
    arrays = dict(y=es.y, rows=es.rows, cols=es.cols, targets=tg, AO=np.asarray(es.strata["AO"], bool),
                  TT=np.asarray(es.strata[tcol], bool), query_strata=views.query_strata(tg, bundle),
                  on_array_rows=np.array([t in bundle.orf_to_col for t in tg]), lockbox_rows=np.isin(tg, lb),
                  same_chr=same_chr, col_target=col_target, chr_unresolved=(chr_r[es.rows] == "") | (chr_c[es.cols] == ""))
    P["run"].mkdir(parents=True, exist_ok=True)
    tmp = P["run"] / f".labels.tmp{os.getpid()}.npz"
    np.savez(tmp, **arrays)
    os.replace(tmp, P["run"] / "labels.npz")
    man = {"eset": es.name, "eset_digest": es.digest(), "split": P["split"], "split_digest": split.digest(),
           "counts": es.counts(), "identity_excluded_arrays": es.meta.get("identity_excluded_arrays"),
           "strata_name_T": tcol, "n_same_chr_cells": int(same_chr.sum()),
           "n_same_chr_pos": int(es.y[same_chr].sum()), "n_chr_unresolved_cells": int(arrays["chr_unresolved"].sum()),
           "n_lockbox_rows": int(arrays["lockbox_rows"].sum()), "n_non_lockbox_rows": int((~arrays["lockbox_rows"]).sum()),
           "n_on_array_rows": int(arrays["on_array_rows"].sum()),
           "query_strata_counts": {int(k): int(v) for k, v in zip(*np.unique(arrays["query_strata"], return_counts=True))},
           "written": time.strftime("%Y-%m-%dT%H:%M:%S")}
    _atomic_write(P["run"] / "label_manifest.json", json.dumps(man, indent=1))
    return man


# ============================================================================ replicate-AP engine (score step)
class Labels:
    def __init__(self, P: dict):
        with np.load(P["run"] / "labels.npz") as z:
            self.z = {k: z[k] for k in z.files}
        self.man = json.loads((P["run"] / "label_manifest.json").read_text())
        self.y = self.z["y"].astype(np.int8)
        self.rows = self.z["rows"]
        self.cols = self.z["cols"]
        self.n_targets = len(self.z["targets"])

    def mask(self, stratum: str):
        z = self.z
        if stratum == "all":
            return None
        if stratum == "AO":
            return z["AO"]
        if stratum in ("TT", "VV"):
            return z["TT"]
        if stratum == "onrows":
            return z["on_array_rows"][self.rows]
        if stratum == "lockbox":
            return z["lockbox_rows"][self.rows]
        if stratum == "nosamechr":
            return ~z["same_chr"]
        raise KeyError(stratum)


class BootEngine:
    """[1 + B] vectors (unit-weight AP, then B replicate APs under the shared stratified query weights) per
    (prediction file, stratum), cached in run/boot/. Two-way vectors use per-cell weights w_t · v_a."""

    def __init__(self, P: dict, L: Labels):
        self.P, self.L = P, L
        self.W = E.bootstrap_weights(L.z["query_strata"], P["B"], np.random.default_rng(boot_seed(P, "query_weights")))
        self.cache = P["run"] / "boot"
        self.chunk = int(config.protocol()["bootstrap"]["chunk"])
        self._tw = None
        self.n_vectors = 0

    def scores(self, pred_file: str, key: str = "E_fixed") -> np.ndarray:
        with np.load(ROOT / pred_file) as z:
            s = np.asarray(z[key], dtype=np.float32)
        if len(s) != len(self.L.y):
            raise AssertionError(f"{pred_file}: {len(s)} cells != labels {len(self.L.y)}")
        return s

    def vector(self, pred_file: str, stratum: str = "all") -> np.ndarray:
        cp = self.cache / f"{Path(pred_file).stem}__{stratum}__B{self.P['B']}.npy"
        if cp.exists():
            return np.load(cp)
        s = self.scores(pred_file)
        m = self.L.mask(stratum)
        y, r = (self.L.y, self.L.rows) if m is None else (self.L.y[m], self.L.rows[m])
        e = E.WeightedAP(y, s if m is None else s[m], r, self.L.n_targets)
        A = np.concatenate([e.ap(self.W[i:i + self.chunk]) for i in range(0, len(self.W), self.chunk)])
        v = np.r_[float(e.ap(np.ones(self.L.n_targets))[0]), A]
        cp.parent.mkdir(parents=True, exist_ok=True)
        np.save(cp, v)
        self.n_vectors += 1
        return v

    def two_way_weights(self):
        """(Wq [B2 x n_targets], V [B2 x n_arrays]): stratified query weights and array-column weights; an array
        column that is a target gene takes that gene's query weight (DESIGN §6.4 sensitivity)."""
        if self._tw is None:
            B2 = self.P["B_two_way"]
            Wq = E.bootstrap_weights(self.L.z["query_strata"], B2, np.random.default_rng(boot_seed(self.P, "two_way_q")))
            ct = self.L.z["col_target"]
            cols = np.unique(self.L.cols[ct[self.L.cols] < 0])
            rng = np.random.default_rng(boot_seed(self.P, "two_way_a"))
            V = np.zeros((B2, len(ct)), dtype=np.float64)
            V[:, cols] = rng.multinomial(len(cols), np.full(len(cols), 1.0 / len(cols)), size=B2)
            tc = np.flatnonzero(ct >= 0)
            V[:, tc] = Wq[:, ct[tc]]
            self._tw = (Wq, V)
        return self._tw

    def two_way_vector(self, pred_file: str) -> np.ndarray:
        cp = self.cache / f"{Path(pred_file).stem}__twoway__B{self.P['B_two_way']}.npy"
        if cp.exists():
            return np.load(cp)
        Wq, V = self.two_way_weights()
        s = self.scores(pred_file)
        sa = E.SortedAP(self.L.y, s)
        r, c = self.L.rows, self.L.cols
        vals = [sa.ap(Wq[b][r] * V[b][c]) for b in range(len(Wq))]
        v = np.r_[sa.ap(np.ones(len(s))), vals]
        cp.parent.mkdir(parents=True, exist_ok=True)
        np.save(cp, v)
        return v


def arm_matrix(rows: list[dict], grid, eng: BootEngine, stratum: str = "all", trajs=None, two_way: bool = False):
    """(A [B x R x J], est [R x J], trajs) for one (arm, predictor) row set (NaN where missing)."""
    tr = sorted({int(x["traj"]) for x in rows}) if trajs is None else list(trajs)
    B = (eng.P["B_two_way"] if two_way else eng.P["B"])
    A = np.full((B, len(tr), len(grid)), np.nan)
    est = np.full((len(tr), len(grid)), np.nan)
    gi = {int(n): j for j, n in enumerate(grid)}
    ti = {t: i for i, t in enumerate(tr)}
    for x in rows:
        if int(x["n"]) in gi and int(x["traj"]) in ti:
            v = eng.two_way_vector(x["pred_file"]) if two_way else eng.vector(x["pred_file"], stratum)
            A[:, ti[int(x["traj"])], gi[int(x["n"])]] = v[1:]
            est[ti[int(x["traj"])], gi[int(x["n"])]] = v[0]
    return A, est, tr


def check_inline(rows: list[dict], eng: BootEngine, tol: float = 1e-9) -> int:
    """The engine's unit-weight AP must equal the AP the harness computed at run time (same tie-grouped AP)."""
    k = 0
    for x in rows:
        v = eng.vector(x["pred_file"])
        if abs(v[0] - float(x["auprc_inline"])) > tol:
            raise AssertionError(f"engine AP {v[0]} != inline {x['auprc_inline']} ({x['pred_file']})")
        k += 1
    return k


def mode_params(mode: str, root: Path = ROOT) -> tuple[dict, dict]:
    cfg = load_cfg(root)
    return cfg, resolve(cfg, mode, root)


def load_plan(P: dict) -> dict:
    return json.loads((P["run"] / "plan.json").read_text())


__all__ = ["Refusal", "load_cfg", "resolve", "arms", "units", "build_plan", "preflight", "check_frozen",
           "frozen_hashes", "write_lock", "read_lock", "guard_final", "ensure_selection", "work", "write_labels",
           "Labels", "BootEngine", "arm_matrix", "check_inline", "progress", "mode_params", "load_plan"]
