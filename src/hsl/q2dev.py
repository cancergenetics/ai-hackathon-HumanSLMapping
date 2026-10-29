"""P3 (oracle ceilings) and P4 (strategy development) helpers on DEV splits only (DESIGN v3.1 §6, §8, §9).

Lives OUTSIDE hsl.strategies because it builds the evaluation-only oracles (hsl.oracles) and reads dev labels for
diagnostics; selectors themselves only ever receive a SelectorView (harness). Refuses any non-dev split.

Arm configurations are named by a config id: '<arm>[_<knob><value>...][__<scenario>...]', e.g.
    R0 | R1 | H | H-static | F_r8 | F-static_r8 | C | M_lam0.3 | oracle_hub | oracle_greedy
    scenario tokens: det (temperature 0) | b50 / b110 (batch size) | ws220 (warm start: S_0 = first 220 of R1
    trajectory r, then this arm) | feas (HAP1-feasible pool: pool minus the bottom array-SMF decile of on-array pool
    genes; off-array kept) | min (HAP1-minimal selector features: no AP-MS in the degree model or selector kernels).
Trajectory r of every arm uses dev trajectory seed r (SeedSequence(dev_traj_base_seed).spawn(.)[r]): the same model
seed as random permutation r (common random numbers), as DESIGN §5 requires.

Runs: run_config(cid, split, r, checkpoints, ...) -> harness.Trajectory with P* (config/predictor.yaml) at every
checkpoint; E_fixed_dev cell predictions saved (float32, results/predictions/p34/, gitignored); per-checkpoint rows
(accuracy + strata + discovery + runtime) appended to a JSON-lines part file. No per-run LEDGER lines
(ledger=False); one line per configuration is written by the analysis scripts.

Statistics: replicate APs under the split's shared stratified query weights (p1.query_weights(split, W_NAME)),
one cached [1 + B] vector per prediction file (results/p{3,4}/boot/, gitignored); arm matrices A[b, r, j].
"""
from __future__ import annotations

import hsl  # noqa: F401  (guard + label gate; must be first)

import json
import time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml

import hsl.predictors  # noqa: F401  (registers the frozen P* method)
from hsl import config, eval as E, harness as H, p1, paths, views

RESULTS = paths.output_root() / "results"
PRED_ROOT = RESULTS / "predictions" / "p34"
_FROZEN = yaml.safe_load((paths.PROJECT_ROOT / "config" / "predictor.yaml").read_text())
PSTAR = str(_FROZEN["primary"]["method"])
PSTAR_PARAMS = dict(_FROZEN["primary"]["params"] or {})
PREDICTORS = ((PSTAR, PSTAR_PARAMS),)
W_NAME = "q2_dev"                       # one shared weight matrix per dev split for every P3/P4 arm (incl. random)


def use_pred_dir() -> Path:
    H.PRED_DIR = PRED_ROOT
    return PRED_ROOT


def dev_split(name: str) -> views.Split:
    sp = views.get_split(name)
    if sp.kind != "dev" or name == "final":
        raise PermissionError("P3/P4 helpers run on dev splits only")
    return sp


# ============================================================================ config ids
@dataclass(frozen=True)
class ArmSpec:
    cid: str
    arm: str
    knobs: dict = field(default_factory=dict)
    deterministic: bool = False
    batch_size: int | None = None
    warm: int = 0
    feasible: bool = False
    minimal: bool = False


_KNOB_PREFIX = {"r": int, "lam": float}


def parse(cid: str) -> ArmSpec:
    toks = cid.split("__")
    head, scen = toks[0], toks[1:]
    parts = head.split("_")
    arm = parts[0]
    if arm == "oracle":
        arm, parts = "_".join(parts[:2]), parts[2:]
    else:
        parts = parts[1:]
    kn = {}
    for p in parts:
        for k, typ in _KNOB_PREFIX.items():
            if p.startswith(k) and p[len(k):].replace(".", "", 1).replace("e-", "", 1).isdigit():
                kn[k] = typ(p[len(k):])
                break
        else:
            raise ValueError(f"unknown knob token {p!r} in {cid!r}")
    spec = dict(cid=cid, arm=arm, knobs=kn)
    for s in scen:
        if s == "det":
            spec["deterministic"] = True
        elif s.startswith("b") and s[1:].isdigit():
            spec["batch_size"] = int(s[1:])
        elif s.startswith("ws") and s[2:].isdigit():
            spec["warm"] = int(s[2:])
        elif s == "feas":
            spec["feasible"] = True
        elif s == "min":
            spec["minimal"] = True
        else:
            raise ValueError(f"unknown scenario token {s!r} in {cid!r}")
    return ArmSpec(**spec)


def family(arm: str) -> str:
    a = arm.split("_")[0]
    return {"H": "H", "H-static": "H", "F": "F", "F-static": "F", "C": "C/M", "M": "C/M"}.get(a, a)


# ============================================================================ selectors
def minimal_knobs(k: dict) -> dict:
    """HAP1-minimal selector features (DESIGN §6.3): fitness + coexpression + protein embeddings; no AP-MS."""
    import copy
    k = copy.deepcopy(k)
    k["degree_model"]["features"]["apms_degree"] = None
    k["kernels"]["external"] = {kk: v for kk, v in k["kernels"]["external"].items() if kk not in ("apms", "apmse")}
    return k


def _patch_minimal():
    from hsl.strategies import base
    if getattr(base, "_hsl_orig_knobs", None) is not None:
        return
    orig_fn = base.knobs
    mk = minimal_knobs(orig_fn())
    base._hsl_orig_knobs = orig_fn
    base.knobs = lambda: mk                          # every call site uses base.knobs() (module attribute lookup)


def _unpatch_minimal():
    from hsl.strategies import base
    if getattr(base, "_hsl_orig_knobs", None) is not None:
        base.knobs = base._hsl_orig_knobs
        base._hsl_orig_knobs = None


class _PoolView:
    """SelectorView proxy whose candidate list is restricted to `allowed` (feasible-pool scenario). Every other
    attribute is the underlying SelectorView's."""

    def __init__(self, v, allowed: frozenset):
        object.__setattr__(self, "_v", v)
        keep = np.array([g in allowed for g in v.candidates], dtype=bool)
        object.__setattr__(self, "candidates", np.asarray(v.candidates)[keep])
        object.__setattr__(self, "candidate_col", np.asarray(v.candidate_col)[keep])

    def __getattr__(self, k):
        return getattr(object.__getattribute__(self, "_v"), k)


class RestrictPool:
    """Wraps a selector so it only ever sees / picks candidates in `allowed` (label-free pool restriction)."""

    def __init__(self, inner, allowed):
        self.inner, self.allowed = inner, frozenset(str(g) for g in allowed)
        self.needs_view = getattr(inner, "needs_view", True)
        self.static = getattr(inner, "static", False)
        self.name = getattr(inner, "name", type(inner).__name__)
        self.params = {**getattr(inner, "params", {}), "pool": "feasible", "n_allowed": len(self.allowed)}
        self.diagnostics = getattr(inner, "diagnostics", [])

    def reset(self, candidates, rng):
        self.inner.reset(np.array([g for g in candidates if g in self.allowed]), rng)
        self.diagnostics = getattr(self.inner, "diagnostics", [])

    def select(self, view, k, rng):
        out = self.inner.select(None if view is None else _PoolView(view, self.allowed), k, rng)
        self.diagnostics = getattr(self.inner, "diagnostics", [])
        return out


def feasible_pool(split: str) -> np.ndarray:
    """Pool minus the bottom decile of ARRAY SMF among on-array pool genes (off-array genes kept: no SMF)."""
    sp = dev_split(split)
    b = views.dev_bundle()
    col = np.array([b.orf_to_col.get(o, -1) for o in sp.pool])
    smf = np.full(len(col), np.nan)
    smf[col >= 0] = np.asarray(b.array_smf, dtype=float)[col[col >= 0]]
    cut = np.nanquantile(smf, 0.1)
    drop = np.isfinite(smf) & (smf <= cut)
    return sp.pool[~drop]


def make_selector(spec: ArmSpec, split: str, n_threads: int = 2):
    from hsl import oracles, strategies
    if spec.minimal:
        _patch_minimal()
    if spec.arm == "R0":
        sel = H.RandomSelector()
    elif spec.arm == "oracle_hub":
        sel = oracles.OracleHubFirst(split)
    elif spec.arm == "oracle_greedy":
        sel = oracles.GreedyOracle(split, n_threads=n_threads)
    else:
        sel = strategies.make(spec.arm, deterministic=spec.deterministic, **spec.knobs)
    if spec.feasible:
        sel = RestrictPool(sel, feasible_pool(split))
    return sel


def traj_seed(r: int) -> np.random.SeedSequence:
    return config.trajectory_seeds("dev", int(r) + 1)[int(r)]


def r1_prefix(split: str, r: int, n: int) -> list[str]:
    """First n genes of R1 trajectory r (selection only: identical to the R1 arm's picks; nothing is fitted)."""
    from hsl import strategies
    t = H.Trajectory(strategies.make("R1"), split, traj_seed(r), r, arm="R1", checkpoints=(), max_n=n,
                     predictors=(), save=False, ledger=False, write_log=False)
    return list(t.run()["selected"][:n])


# ============================================================================ runs
@lru_cache(maxsize=4)
def _eset(split: str):
    return H.get_eset(views.dev_bundle(), dev_split(split), "E_fixed")


@lru_cache(maxsize=1)
def _on_array() -> frozenset:
    return frozenset(views.dev_bundle().orf_to_col)


def checkpoint_rows(log: dict, cid: str, split: str, r: int, phase: str, wall_s: float, s0_n: int = 0) -> list[dict]:
    es = _eset(split)
    on = _on_array()
    rows = []
    sel = log["selected"]
    rounds = log.get("rounds", [])
    for c, d in sorted(log["checkpoints"].items(), key=lambda kv: int(kv[0])):
        c = int(c)
        dd = d[PSTAR]
        m = dd["metrics"]["E_fixed_dev"]
        rt = dd["runtime_s"]
        pred = PRED_ROOT / split / f"{dd['run_id']}.npz"
        s = p1.load_scores(str(pred.relative_to(paths.PROJECT_ROOT)))
        S = sel[:c]
        new = S[s0_n:] if s0_n else S
        used = [x for x in rounds if x["n_before"] < c]
        row = {"phase": phase, "config": cid, "arm": parse(cid).arm, "family": family(parse(cid).arm),
               "split": split, "traj": int(r), "n": c, "predictor": PSTAR, "traj_id": log["traj_id"],
               "model_seed": log["model_seed"], "run_id": dd["run_id"],
               "pred_file": str(pred.relative_to(paths.PROJECT_ROOT)),
               **{k: m.get(k) for k in ("n_cells", "n_pos", "prevalence", "auprc", "lift", "macro_auprc", "auroc",
                                        "p_at_10", "p_at_50", "auprc_AO", "prev_AO", "auprc_VV", "prev_VV")},
               **p1.extra_metrics(split, es, s), **p1.sl_detected(S),
               "offarray_share": float(np.mean([g not in on for g in S])) if len(S) else None,
               "offarray_share_new": float(np.mean([g not in on for g in new])) if len(new) else None,
               "fit_predict_s": rt.get("fit_predict_s"), "eval_total_s": rt.get("total_s"),
               "select_s_cum": float(sum(x["select_s"] for x in used)),
               "view_s_cum": float(sum(x["view_s"] for x in used)), "n_rounds": len(used),
               "traj_wall_s": wall_s, "batch_size": log.get("batch_size"), "s0_n": s0_n}
        rows.append(row)
    return rows


def run_config(cid: str, split: str, r: int, checkpoints, part: Path, *, max_n: int | None = None,
               phase: str = "p4", n_threads: int = 2, extra_log: Path | None = None,
               write_log: bool = True) -> list[dict]:
    from hsl import strategies
    dev_split(split)
    use_pred_dir()
    spec = parse(cid)
    sel = make_selector(spec, split, n_threads=n_threads)
    s0 = r1_prefix(split, r, spec.warm) if spec.warm else []
    cps = tuple(sorted(int(c) for c in checkpoints))
    traj = H.Trajectory(sel, split, traj_seed(r), int(r), arm=cid, batch_size=spec.batch_size, checkpoints=cps,
                        max_n=int(max_n or max(cps)), predictors=PREDICTORS, s0=tuple(s0), save=True, ledger=False,
                        write_log=write_log)
    t0 = time.perf_counter()
    try:
        log = strategies.run_trajectory(traj) if hasattr(sel, "diagnostics") else traj.run()
    finally:
        _unpatch_minimal()                          # the minimal-features knob patch is scoped to this trajectory
    wall = time.perf_counter() - t0
    if extra_log is not None and hasattr(sel, "log"):
        extra_log.parent.mkdir(parents=True, exist_ok=True)
        extra_log.write_text(json.dumps({"traj_id": log["traj_id"], "rounds": sel.log}, indent=0, default=float))
    rows = checkpoint_rows(log, cid, split, r, phase, wall, len(s0))
    for row in rows:
        p1.append_row(part, row)
    return rows


def done(parts, cid: str, split: str, r: int, n: int) -> bool:
    return any(x["config"] == cid and x["split"] == split and x["traj"] == r and x["n"] == n
               for x in p1.read_rows(parts))


# ============================================================================ random arm (P1 rows, frozen P*)
def random_rows() -> list[dict]:
    """Frozen-P* R0 rows produced by P1 (n = 30/160/330/660 reused from P2; the rest from P1 batch4)."""
    parts = sorted((p1.P1_DIR / "parts").glob("p1lb_*.jsonl")) + [p1.P1_DIR / "parts" / "p2_reuse_frozen_pstar.jsonl"]
    out = {}
    for x in p1.read_rows(parts):
        if x.get("predictor") != PSTAR or x.get("kind", "curve") != "curve" or not x.get("pred_file"):
            continue
        key = (x["split"], int(x["perm"]), int(x["n"]))
        out[key] = {"config": "R0", "arm": "R0", "family": "R0", "split": x["split"], "traj": int(x["perm"]),
                    "n": int(x["n"]), **{k: x.get(k) for k in ("auprc", "auprc_AO", "auprc_VV", "auprc_on_array_rows",
                                                           "auprc_off_array_rows", "macro_auprc", "sl_detected",
                                                           "fit_predict_s", "pred_file", "model_seed", "run_id",
                                                           "prevalence", "prec_at_recall10")}}
    return list(out.values())


# ============================================================================ replicate APs
_W: dict = {}


def weights(split: str) -> np.ndarray:
    if split not in _W:
        _W[split] = p1.query_weights(split, W_NAME)
    return _W[split]


def boot_vector(pred_file: str, split: str, cache_dir: Path) -> np.ndarray:
    """[1 + B]: unit-weight AP, then the B replicate APs under the split's shared query weights. Cached."""
    cp = Path(cache_dir) / split / f"{Path(pred_file).stem}.npy"
    if cp.exists():
        return np.load(cp)
    es = _eset(split)
    s = p1.load_scores(pred_file)
    W = weights(split)
    e = E.WeightedAP(es.y, s, es.rows, W.shape[1])
    ch = int(config.protocol()["bootstrap"]["chunk"])
    A = np.concatenate([e.ap(W[i:i + ch]) for i in range(0, len(W), ch)])
    est = float(e.ap(np.ones(W.shape[1]))[0])
    v = np.r_[est, A]
    cp.parent.mkdir(parents=True, exist_ok=True)
    np.save(cp, v)
    return v


def arm_matrix(rows: list[dict], split: str, grid, cache_dir: Path, trajs=None):
    """(A [B x R x J], est [R x J], traj ids) for one arm's rows on one split (NaN where a checkpoint is missing)."""
    rows = [x for x in rows if x["split"] == split]
    tr = sorted({int(x["traj"]) for x in rows}) if trajs is None else list(trajs)
    B = weights(split).shape[0]
    A = np.full((B, len(tr), len(grid)), np.nan)
    est = np.full((len(tr), len(grid)), np.nan)
    gi = {int(n): j for j, n in enumerate(grid)}
    ti = {t: i for i, t in enumerate(tr)}
    for x in rows:
        if int(x["n"]) in gi and int(x["traj"]) in ti:
            v = boot_vector(x["pred_file"], split, cache_dir)
            A[:, ti[int(x["traj"])], gi[int(x["n"])]] = v[1:]
            est[ti[int(x["traj"])], gi[int(x["n"])]] = v[0]
    return A, est, tr


def cols(A, est, j):
    """Replicates/estimates at grid index j restricted to trajectories that have it."""
    ok = np.isfinite(est[:, j])
    return A[:, ok, j], est[ok, j]


def check_est(rows: list[dict], split: str, cache_dir: Path, tol: float = 1e-9) -> int:
    """Unit-weight engine AP must equal the stored harness AUPRC for every row (same tie-grouped AP)."""
    k = 0
    for x in rows:
        if x["split"] != split or x.get("auprc") is None:
            continue
        v = boot_vector(x["pred_file"], split, cache_dir)
        if abs(v[0] - float(x["auprc"])) > tol:
            raise AssertionError(f"engine AP {v[0]} != stored {x['auprc']} for {x['pred_file']}")
        k += 1
    return k


# ============================================================================ DESIGN §6.4 statistics
def ci(x) -> tuple[float, float]:
    x = np.asarray(x, float)
    return E.percentile_ci(x[np.isfinite(x)]) if np.isfinite(x).any() else (float("nan"), float("nan"))


def gap_closure_joint(pairs, rng) -> dict:
    """G = Σ_k (s160 − r160) / Σ_k (r660 − r160) over splits k (one split: the DESIGN G), joint bootstrap
    (query weights w^b of each split + whole trajectories resampled within arm; splits independent).
    pairs: list of (A_s [B x R_s x 1], A_r [B x R_r x 2] (lo, hi), est_s [R_s x 1], est_r [R_r x 2])."""
    num_b = den_b = 0.0
    num = den = 0.0
    for As, Ar, es, er in pairs:
        Ms = E.joint_resample(As, rng)[:, 0]
        Mr = E.joint_resample(Ar, rng)
        num_b = num_b + (Ms - Mr[:, 0])
        den_b = den_b + (Mr[:, 1] - Mr[:, 0])
        num += float(np.nanmean(es[:, 0]) - np.nanmean(er[:, 0]))
        den += float(np.nanmean(er[:, 1]) - np.nanmean(er[:, 0]))
    with np.errstate(divide="ignore", invalid="ignore"):
        Gb = num_b / den_b
    lo_den = float(np.quantile(den_b, 0.025))
    return {"G": num / den if den else float("nan"), "G_ci": ci(Gb), "num": num, "den": den,
            "den_ci": ci(den_b), "defined": bool(lo_den > 0), "G_boot_median": float(np.nanmedian(Gb))}


def welch_at(As, es, Ar, er, margin: float = 0.0, scale_r: float = 1.0) -> dict:
    """E.welch on (strategy at budget a) vs (scale_r · random at budget b) using only trajectories present."""
    oks, okr = np.isfinite(es), np.isfinite(er)
    return E.welch(As[:, oks], es[oks], scale_r * Ar[:, okr], scale_r * er[okr], margin=margin)


def contrast_suite(As, es, Ar, er, grid, tag: str, n_s: int = 160) -> dict:
    """Per split: Δ(n) at every common checkpoint; c1-style s160 vs r490; brief s160 vs r660; matches contrast
    s160 − 0.95·r660 > 0; G (joint bootstrap) and G*; relative Δ(n)."""
    g = list(grid)
    out = {"delta": {}}
    for n in g:
        j = g.index(n)
        if np.isfinite(es[:, j]).sum() >= 1 and np.isfinite(er[:, j]).sum() >= 2:
            w = welch_at(As[:, :, j], es[:, j], Ar[:, :, j], er[:, j])
            w["rel"] = w["delta"] / float(np.nanmean(er[:, j]))
            w["s_mean"], w["r_mean"] = float(np.nanmean(es[:, j])), float(np.nanmean(er[:, j]))
            w["s_sd"] = float(np.nanstd(es[:, j], ddof=1)) if np.isfinite(es[:, j]).sum() > 1 else float("nan")
            w["r_sd"] = float(np.nanstd(er[:, j], ddof=1))
            out["delta"][n] = w
    js = g.index(n_s)
    m = float(config.protocol()["stats"]["match_margin_frac"])
    if 490 in g and np.isfinite(er[:, g.index(490)]).sum() >= 2:
        out["c1"] = welch_at(As[:, :, js], es[:, js], Ar[:, :, g.index(490)], er[:, g.index(490)])
    if 660 in g:
        j6 = g.index(660)
        out["brief"] = welch_at(As[:, :, js], es[:, js], Ar[:, :, j6], er[:, j6])
        out["matches"] = welch_at(As[:, :, js], es[:, js], Ar[:, :, j6], er[:, j6], scale_r=1.0 - m)
        r160, r660 = float(np.nanmean(er[:, js])), float(np.nanmean(er[:, j6]))
        out["G_star"] = float(E.g_star(r160, r660))
        rng = np.random.default_rng(p1.boot_seed(f"q2dev_G_{tag}"))
        out["G"] = gap_closure_joint([(As[:, :, [js]], Ar[:, :, [js, j6]], es[:, [js]], er[:, [js, j6]])], rng)
    return out


def n_req_joint(As, es, Ar, er, grid_r, j_s: int, rng) -> dict:
    """Random-equivalent budget of the strategy mean at one checkpoint: point (trajectory means) + joint bootstrap
    (frozen PAVA / log-n inversion of the random mean curve on the dense grid; censored replicates 0/+inf)."""
    target = float(np.nanmean(es[:, j_s]))
    curve = np.nanmean(er, axis=0)
    pt = E.n_req(target, grid_r, curve)
    Ms = E.joint_resample(As[:, :, [j_s]], rng)[:, 0]
    Mr = E.joint_resample(Ar, rng)
    reps = [E.n_req(Ms[b], grid_r, Mr[b]) for b in range(len(Ms))]
    v = np.array([x.n for x in reps])
    return {"n_req": pt.n, "censored": pt.censored, "ci": ci(v),
            "frac_left": float(np.mean([x.censored == "left" for x in reps])),
            "frac_right": float(np.mean([x.censored == "right" for x in reps]))}


def ledger(rec: dict) -> None:
    """One LEDGER line per configuration (DESIGN §8), with git state, predictor and bootstrap-weight name."""
    rec = dict(rec)
    rec.setdefault("manifest", {}).update({"git": H._git(), "weights": W_NAME, "predictor": PSTAR,
                                           "predictor_params": PSTAR_PARAMS})
    H.append_ledger(rec)
