#!/usr/bin/env python
"""Independent rescorer of the final evaluation (DESIGN v3.1 §11 "independent rescoring").

Written from the DESIGN text (§3, §6.4, §7.1, revision-history deviations) and the documented run-directory formats
only. It imports numpy / pandas / scipy / sklearn and NOTHING from the project (no hsl, no ysl, nothing from src/ or
third_party/), so a bug shared with scripts/score_final.py cannot pass silently.

Inputs (read-only), with RUN = <run-dir>/run:
  RUN/labels.npz          E_fixed cells: y, rows, cols, targets, AO, TT, query_strata (11 split strata),
                          on_array_rows, lockbox_rows, same_chr, sl_tertile, col_target
  RUN/label_manifest.json counts written with the labels (cross-checked here)
  RUN/parts/*.jsonl       one row per evaluation unit (cid, traj, n, pred, pred_file, auprc_inline, ...)
  RUN/predictions/...     float32 score vector per unit, key 'E_fixed', aligned with labels.npz cells
  RUN/predictions/_cells/<eset digest>.npz   cell index (rows, cols) written by the harness (alignment check)
  RUN/plan.json           arms / roles / pool size
  RUN/boot/*.npy          pipeline cache: element 0 = the pipeline's unit-weight AP per (prediction, stratum)
  <run-dir>/summary.json, confirmatory.csv, contrasts.csv   pipeline outputs to compare against

What is recomputed
  * every unit's AUPRC (tie-grouped average precision, sklearn semantics) on E_fixed and the strata AO, TT,
    on-array rows (plus offrows, SL tertiles, lockbox-only rows and no-same-chromosome cells, for the contrast table)
  * an own stratified query bootstrap (11 strata, B = 2000 by default, own seed) driving exact weighted AP replicates
  * the confirmatory family (Welch SE, nu, one-sided p, Holm at alpha = 0.025), the brief contrast s160 - r660,
    gap closure G with a joint-bootstrap CI, G*, the 'matches' linear contrast and its 3-condition rule,
    the headline sign rule, Q1 n_80 / n_90 (PAVA + log-n interpolation, joint bootstrap) and n_req(n)
  * a comparison table against the pipeline (exact for point values, approximate for CIs and p-values)

    .venv/bin/python analysis/rescore/rescore.py --run-dir results/final_dryrun --out analysis/rescore/out/dryrun
"""
from __future__ import annotations

import os

# <= 2 threads (compute constraint). Must precede the numpy import.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "2"

import argparse  # noqa: E402
import glob  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import zlib  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

DEFAULT_SEED = 20261030          # own bootstrap seed (rescorer; not a project seed, recorded in the output)
DEFAULT_B = 2000                 # DESIGN §6.4
ALPHA = 0.025                    # DESIGN §6.4 familywise, one-sided
MATCH_MARGIN = 0.05              # DESIGN §6.4: H0 s160 - 0.95 r660 <= 0
G_STAR_MIN = 0.5                 # DESIGN §6.4
N_FRACS = (0.80, 0.90)           # DESIGN §0 / §7.1
CI = 0.95
SETS = ("all", "AO", "TT", "onrows", "offrows", "tert1", "tert2", "tert3", "lockbox", "nosamechr")
EXACT_TOL = 1e-9
ROOT = Path(__file__).resolve().parents[2]


# ====================================================================================== core statistics
def ap_point(y, s) -> float:
    """Average precision with tied scores grouped at one threshold (sklearn average_precision_score semantics):
    AP = sum over distinct thresholds (descending) of (R_k - R_{k-1}) * P_k."""
    y = np.asarray(y).astype(np.float64)
    s = np.asarray(s, dtype=np.float64)
    o = np.argsort(-s, kind="mergesort")
    return ap_sorted(y[o], s[o])


def ap_sorted(yy, ss) -> float:
    """ap_point on cells already sorted by descending score (a subset of a sorted vector stays sorted)."""
    yy = np.asarray(yy, dtype=np.float64)
    if yy.size == 0 or yy.sum() == 0:
        return float("nan")
    last = np.r_[ss[1:] != ss[:-1], True]           # last element of each tie group
    tp = np.cumsum(yy)[last]
    nn = np.cumsum(np.ones_like(yy))[last]
    dtp = np.diff(np.r_[0.0, tp])
    return float(np.sum(dtp * (tp / nn)) / tp[-1])


class ReplicateAP:
    """Exact weighted AP for many query-weight vectors at once. Cell weight = weight of its query row.
    With w = cell weights, AP(w) = sum_g dTP_g(w) * TPcum_g(w) / Ncum_g(w) / TP(w) over tie groups g (only groups
    with a positive contribute). Ncum at positive groups = cumulative sum of per-row cell counts in the blocks between
    consecutive positive groups; both counts are integer, so float32 matmuls are exact (< 2^24)."""

    def __init__(self, y, s, rows, n_rows: int):
        y = np.asarray(y).astype(bool)
        s = np.asarray(s, dtype=np.float64)
        rows = np.asarray(rows, dtype=np.int64)
        if not np.all(np.isfinite(s)):
            raise ValueError("non-finite scores")
        o = np.argsort(-s, kind="mergesort")
        ss, yy, rr = s[o], y[o], rows[o]
        g = np.cumsum(np.r_[True, ss[1:] != ss[:-1]]) - 1        # tie-group id, descending score
        posg = np.unique(g[yy])
        K = len(posg)
        k = np.searchsorted(posg, g, side="left")                  # block = first positive group at or after g
        keep = k < K
        self.K, self.n_rows = K, n_rows
        self.D = np.bincount(rr[keep] * K + k[keep], minlength=n_rows * K).reshape(n_rows, K).astype(np.float32)
        self.Pg = np.bincount(rr[yy] * K + k[yy], minlength=n_rows * K).reshape(n_rows, K).astype(np.float32)
        if n_rows * K * 4 * 2 > 2.5e8:
            raise MemoryError(f"ReplicateAP matrices too large ({n_rows} x {K})")

    def ap(self, W, chunk: int = 250) -> np.ndarray:
        W = np.atleast_2d(np.asarray(W, dtype=np.float32))
        out = np.empty(len(W))
        for i in range(0, len(W), chunk):
            w = W[i:i + chunk]
            ncum = np.cumsum((w @ self.D).astype(np.float64), axis=1)
            dtp = (w @ self.Pg).astype(np.float64)
            tpc = np.cumsum(dtp, axis=1)
            tot = tpc[:, -1]
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.where(dtp > 0, dtp * tpc / np.where(ncum > 0, ncum, 1.0), 0.0).sum(axis=1) / tot
            out[i:i + chunk] = np.where(tot > 0, r, np.nan)
        return out


def stratified_query_weights(strata, B: int, rng) -> np.ndarray:
    """[B x n_queries] multiplicities: within each stratum, resample its queries with replacement (DESIGN §6.4)."""
    strata = np.asarray(strata)
    W = np.zeros((B, len(strata)), dtype=np.float64)
    for h in np.unique(strata):
        idx = np.flatnonzero(strata == h)
        W[:, idx] = rng.multinomial(len(idx), np.full(len(idx), 1.0 / len(idx)), size=B)
    return W


def welch(rep_a, pt_a, rep_b, pt_b, scale_b: float = 1.0) -> dict:
    """DESIGN §6.4. Estimand mean_a - scale_b * mean_b over trajectories (point APs).
    SE^2 = V_Q + s2_a/R_a + scale_b^2 s2_b/R_b; V_Q = variance over bootstrap replicates of the trajectory-mean
    contrast (same query weights for both arms); nu = SE^4 / [(s2_a/R_a)^2/(R_a-1) + (scale_b^2 s2_b/R_b)^2/(R_b-1)];
    one-sided p = 1 - T_nu(delta/SE); 95% CI delta +- t_{nu,0.975} SE. A term with R = 1 has no variance estimate:
    it is set to 0 and dropped from nu (flagged)."""
    pt_a, pt_b = np.asarray(pt_a, float), np.asarray(pt_b, float)
    Ra, Rb = len(pt_a), len(pt_b)
    delta = pt_a.mean() - scale_b * pt_b.mean()
    out = {"delta": float(delta), "mean_a": float(pt_a.mean()), "mean_b": float(pt_b.mean()), "R_a": Ra, "R_b": Rb,
           "sd_a": float(pt_a.std(ddof=1)) if Ra > 1 else float("nan"),
           "sd_b": float(pt_b.std(ddof=1)) if Rb > 1 else float("nan")}
    if rep_a is None or rep_b is None:
        return out
    c = np.asarray(rep_a).mean(axis=1) - scale_b * np.asarray(rep_b).mean(axis=1)
    VQ = float(np.var(c, ddof=1))
    ta = pt_a.var(ddof=1) / Ra if Ra > 1 else 0.0
    tb = scale_b ** 2 * pt_b.var(ddof=1) / Rb if Rb > 1 else 0.0
    se2 = VQ + ta + tb
    den = (ta ** 2 / (Ra - 1) if Ra > 1 else 0.0) + (tb ** 2 / (Rb - 1) if Rb > 1 else 0.0)
    nu = se2 ** 2 / den if den > 0 else float("inf")
    se = float(np.sqrt(se2))
    p = float(stats.t.sf(delta / se, nu)) if se > 0 else float("nan")
    q = float(stats.t.ppf(0.5 + CI / 2, nu))
    out.update(V_Q=VQ, se=se, df=float(nu), p=p, lo=float(delta - q * se), hi=float(delta + q * se),
               traj_var_missing=bool(Ra < 2 or Rb < 2))
    return out


def holm(p, alpha: float = ALPHA) -> tuple[np.ndarray, np.ndarray]:
    """Holm step-down: adjusted p_(i) = max_{j<=i} min(1, (m-j+1) p_(j)); reject iff adjusted p <= alpha."""
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p, kind="mergesort")
    adj_sorted = np.minimum(1.0, np.maximum.accumulate((m - np.arange(m)) * p[o]))
    adj = np.empty(m)
    adj[o] = adj_sorted
    return adj, adj <= alpha


def pava(v, w=None) -> np.ndarray:
    """Non-decreasing isotonic regression (pool adjacent violators), equal weights unless given."""
    v = np.asarray(v, float)
    w = np.ones_like(v) if w is None else np.asarray(w, float)
    blocks = []                                    # [mean, weight, count]
    for x, wt in zip(v, w):
        blocks.append([x, wt, 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            m2, w2, c2 = blocks.pop()
            m1, w1, c1 = blocks.pop()
            blocks.append([(m1 * w1 + m2 * w2) / (w1 + w2), w1 + w2, c1 + c2])
    return np.concatenate([np.full(c, m) for m, _, c in blocks])


def invert_curve(target: float, grid, fitted) -> tuple[float, str]:
    """Smallest n with f(n) >= target, f = the isotonic curve linearly interpolated in log n between grid points
    (never extrapolated). Below the first grid point -> 0 ('left'); never reached -> +inf ('right'). On a flat
    segment this is its left endpoint."""
    grid = np.asarray(grid, float)
    f = np.asarray(fitted, float)
    if not np.isfinite(target) or not np.all(np.isfinite(f)):
        return float("nan"), "nan"
    if f[0] >= target:
        return 0.0, "left"
    if f[-1] < target:
        return float("inf"), "right"
    k = int(np.argmax(f >= target))
    x0, x1 = np.log(grid[k - 1]), np.log(grid[k])
    x = x0 + (target - f[k - 1]) / (f[k] - f[k - 1]) * (x1 - x0)
    return float(np.exp(x)), ""


def pct_ci(vals, level: float = CI) -> tuple[float, float]:
    """Percentile interval; +inf / 0 codes of censored replicates are kept (one-sided bounds); NaN dropped."""
    v = np.asarray(vals, float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return float("nan"), float("nan")
    big = 1e300
    v = np.where(np.isposinf(v), big, v)
    lo, hi = np.quantile(v, [(1 - level) / 2, (1 + level) / 2])
    f = lambda x: float("inf") if x > 1e299 else float(x)  # noqa: E731
    return f(lo), f(hi)


# ====================================================================================== data
class Run:
    def __init__(self, run_dir: Path, root: Path):
        run_dir = Path(run_dir)
        if (run_dir / "run").is_dir():
            self.out_dir, self.run = run_dir, run_dir / "run"
        else:
            self.out_dir, self.run = run_dir.parent, run_dir
        self.root = root
        z = np.load(self.run / "labels.npz")
        self.z = {k: z[k] for k in z.files}
        self.y = self.z["y"].astype(np.int8)
        self.rows = self.z["rows"].astype(np.int64)
        self.cols = self.z["cols"].astype(np.int64)
        self.n_rows = len(self.z["targets"])
        self.man = json.loads((self.run / "label_manifest.json").read_text())
        self.plan = json.loads((self.run / "plan.json").read_text())
        self.units = self._read_parts()
        self.checks = []

    def _read_parts(self) -> pd.DataFrame:
        recs, bad = [], 0
        for f in sorted(glob.glob(str(self.run / "parts" / "*.jsonl"))):
            lines = Path(f).read_text().split("\n")
            for i, ln in enumerate(lines):
                if not ln.strip():
                    continue
                try:
                    recs.append(json.loads(ln))
                except json.JSONDecodeError:
                    bad += 1          # a truncated last line of a killed job (its unit is re-run)
        df = pd.DataFrame(recs)
        self.n_bad_lines = bad
        self.n_dup = int(df.duplicated("uid", keep="last").sum())
        df = df.drop_duplicates("uid", keep="last").reset_index(drop=True)
        df["pred_path"] = [str(self.root / p) if isinstance(p, str) else None for p in df.pred_file]
        df["has_pred"] = [p is not None and Path(p).exists() for p in df.pred_path]
        return df

    def mask(self, st: str):
        z, r = self.z, self.rows
        return {"all": None, "AO": z["AO"].astype(bool), "TT": z["TT"].astype(bool),
                "onrows": z["on_array_rows"].astype(bool)[r], "offrows": ~z["on_array_rows"].astype(bool)[r],
                "tert1": z["sl_tertile"][r] == 1, "tert2": z["sl_tertile"][r] == 2, "tert3": z["sl_tertile"][r] == 3,
                "lockbox": z["lockbox_rows"].astype(bool)[r], "nosamechr": ~z["same_chr"].astype(bool)}[st]

    def scores(self, path: str) -> np.ndarray:
        with np.load(path) as z:
            s = np.asarray(z["E_fixed"])
        if len(s) != len(self.y):
            raise ValueError(f"{path}: {len(s)} scores != {len(self.y)} cells")
        return s

    # ---------------------------------------------------------------- integrity checks
    def integrity(self) -> list[dict]:
        out = []

        def chk(name, ok, detail=""):
            out.append({"check": name, "ok": bool(ok), "detail": str(detail)})

        y, z, m = self.y, self.z, self.man
        c = m["counts"]
        chk("labels binary", set(np.unique(y)) <= {0, 1}, np.unique(y))
        chk("manifest n_cells", len(y) == c["n_cells"], f"{len(y)} vs {c['n_cells']}")
        chk("manifest n_pos", int(y.sum()) == c["n_pos"], f"{int(y.sum())} vs {c['n_pos']}")
        chk("manifest n_rows", self.n_rows == c["n_rows"], f"{self.n_rows} vs {c['n_rows']}")
        chk("manifest n_cols", len(np.unique(self.cols)) == c["n_cols"], f"{len(np.unique(self.cols))} vs {c['n_cols']}")
        AO, TT = z["AO"].astype(bool), z["TT"].astype(bool)
        tname = m.get("strata_name_T", "TT")
        chk("AO/TT partition E_fixed", not np.any(AO & TT) and np.all(AO | TT), f"overlap {int((AO & TT).sum())}")
        chk("manifest AO counts", int(AO.sum()) == c["AO"]["n_cells"] and int(y[AO].sum()) == c["AO"]["n_pos"])
        chk(f"manifest {tname} counts", int(TT.sum()) == c[tname]["n_cells"] and int(y[TT].sum()) == c[tname]["n_pos"])
        # TT columns are target genes on the array; AO columns are not target genes
        ct = z["col_target"]
        chk("TT columns are target genes", np.all(ct[self.cols[TT]] >= 0) and np.all(ct[self.cols[AO]] < 0))
        # no self pairs
        chk("no self pairs", not np.any(ct[self.cols] == self.rows), int(np.sum(ct[self.cols] == self.rows)))
        # unique cells
        key = self.rows * (int(self.cols.max()) + 1) + self.cols
        chk("cells unique", len(np.unique(key)) == len(key))
        pos_row = np.bincount(self.rows, weights=y, minlength=self.n_rows)
        chk("manifest rows with zero positives", int((pos_row == 0).sum()) == c["n_rows_zero_pos"],
            f"{int((pos_row == 0).sum())} vs {c['n_rows_zero_pos']}")
        # 11 split strata (DESIGN §2: zero stratum + 10 quantile bins)
        qs = z["query_strata"]
        chk("11 bootstrap strata (0..10)", set(np.unique(qs)) == set(range(11)), np.unique(qs))
        chk("stratum 0 rows have no E_fixed positives", np.all(pos_row[qs == 0] == 0))
        means = [pos_row[qs == h].mean() for h in range(1, 11) if np.any(qs == h)]
        rho = stats.spearmanr(np.arange(len(means)), means).statistic
        chk("strata ordered by SL count (Spearman of stratum means >= 0.9)", rho >= 0.9, f"rho = {rho:.3f}")
        tm = [pos_row[z["sl_tertile"] == t].mean() for t in (1, 2, 3)]
        chk("SL tertiles ordered by E_fixed positives", tm[0] < tm[1] < tm[2], np.round(tm, 2))
        chk("manifest on-array rows", int(z["on_array_rows"].sum()) == m["n_on_array_rows"])
        chk("manifest lockbox rows", int(z["lockbox_rows"].sum()) == m["n_lockbox_rows"])
        chk("manifest same-chr cells", int(z["same_chr"].sum()) == m["n_same_chr_cells"])
        # harness cell index == labels cell index
        cp = self.run / "predictions" / "_cells" / f"{m['eset_digest']}.npz"
        if cp.exists():
            with np.load(cp) as cz:
                ok = (np.array_equal(cz["rows"], self.rows) and np.array_equal(cz["cols"], self.cols)
                      and np.array_equal(cz["targets"], z["targets"]))
            chk("harness cell index == labels cells", ok, cp.name)
        else:
            chk("harness cell index == labels cells", False, f"missing {cp}")
        U = self.units
        chk("unit rows: eset digest", (U.eset_digest == m["eset_digest"]).all())
        chk("unit rows: n_cells / n_pos", (U.n_cells == c["n_cells"]).all() and (U.n_pos == c["n_pos"]).all())
        chk("all prediction files present", U.has_pred.all(), f"{int((~U.has_pred).sum())} missing")
        chk("parts: no malformed / duplicate rows", self.n_bad_lines == 0 and self.n_dup == 0,
            f"bad {self.n_bad_lines}, dup {self.n_dup}")
        chk("units done == plan n_units", len(U) == self.plan["n_units"], f"{len(U)} / {self.plan['n_units']}")
        self.checks = out
        return out


# ====================================================================================== engine
class Engine:
    """Replicate-AP vectors [unit-weight AP, B replicate APs] per prediction file under ONE weight matrix W shared by
    all arms (DESIGN §6.4). Own weights: SeedSequence(seed, spawn_key=(0,)). W may be given explicitly (diagnostic)."""

    def __init__(self, R: Run, B: int, seed: int, cache: Path, W=None, tag: str | None = None):
        self.R, self.B, self.seed = R, B, seed
        if W is None:
            W = stratified_query_weights(R.z["query_strata"], B,
                                         np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(0,))))
        self.W = np.asarray(W, dtype=np.float64)
        self.B = len(self.W)
        self.tag = tag or f"B{self.B}_s{seed}"
        self.cache = cache
        cache.mkdir(parents=True, exist_ok=True)
        self.n_new = 0

    def replicate(self, path: str, stratum: str = "all") -> np.ndarray:
        cp = self.cache / f"{Path(path).stem}__{stratum}__{self.tag}.npy"
        if cp.exists():
            return np.load(cp)
        s = self.R.scores(path)
        m = self.R.mask(stratum)
        y, rows = (self.R.y, self.R.rows) if m is None else (self.R.y[m], self.R.rows[m])
        e = ReplicateAP(y, s if m is None else s[m], rows, self.R.n_rows)
        v = np.r_[e.ap(np.ones(self.R.n_rows)), e.ap(self.W)]
        np.save(cp, v)
        self.n_new += 1
        return v


def point_table(R: Run) -> pd.DataFrame:
    """Unit-weight AP of every unit on every set (one sort per unit)."""
    masks = {st: R.mask(st) for st in SETS}
    recs = []
    for i, u in R.units[R.units.has_pred].iterrows():
        s = R.scores(u.pred_path).astype(np.float64)
        o = np.argsort(-s, kind="mergesort")
        ss, yy = s[o], R.y[o]
        rec = {"uid": u.uid, "cid": u.cid, "pred": u.pred, "n": int(u.n), "traj": int(u.traj),
               "pred_file": u.pred_file, "auprc_inline": u.auprc_inline,
               "n_nonfinite": int((~np.isfinite(s)).sum())}
        for st, mk in masks.items():
            if mk is None:
                rec[f"ap_{st}"] = ap_sorted(yy, ss)
            else:
                mo = mk[o]
                rec[f"ap_{st}"] = ap_sorted(yy[mo], ss[mo])
        recs.append(rec)
    return pd.DataFrame(recs)


# ====================================================================================== analyses
class Rescorer:
    def __init__(self, R: Run, E: Engine, PT: pd.DataFrame):
        self.R, self.E, self.PT = R, E, PT
        self.fins = [c for c, d in R.plan["arms"].items() if d["role"] == "finalist"]

    def pts(self, cid, pred, n, st="all") -> pd.DataFrame:
        d = self.PT[(self.PT.cid == cid) & (self.PT.pred == pred) & (self.PT.n == n)].sort_values("traj")
        return d

    def mat(self, cid, pred, n):
        """(replicates [B x R], points [R], trajs) for the full set."""
        d = self.pts(cid, pred, n)
        if d.empty:
            return None, None, []
        V = np.stack([self.E.replicate(self.R.root / f) for f in d.pred_file], axis=1)   # (1+B) x R
        return V[1:], V[0], list(d.traj)

    def w(self, a, pa, na, b, pb, nb, scale=1.0, boot=True, st="all"):
        if boot:
            Ra, ea, _ = self.mat(a, pa, na)
            Rb, eb, _ = self.mat(b, pb, nb)
            if ea is None or eb is None:
                return None
            return welch(Ra, ea, Rb, eb, scale)
        da, db = self.pts(a, pa, na), self.pts(b, pb, nb)
        if da.empty or db.empty:
            return None
        return welch(None, da[f"ap_{st}"].to_numpy(), None, db[f"ap_{st}"].to_numpy(), scale)

    # ---------------------------------------------------------------- family
    def family(self) -> pd.DataFrame:
        recs = []
        for f in self.fins:
            for h, (nb) in (("c1", 490), ("a", 160)):
                r = self.w(f, "pstar", 160, "R0", "pstar", nb)
                recs.append({"finalist": f, "hypothesis": h, "n_a": 160, "n_b": nb, **(r or {"p": np.nan})})
        df = pd.DataFrame(recs)
        ok = np.isfinite(df.p.to_numpy(float))
        df["p_holm"], df["reject"] = np.nan, False
        adj, rej = holm(df.p.to_numpy(float)[ok])
        df.loc[ok, "p_holm"], df.loc[ok, "reject"] = adj, rej
        df["m_hypotheses"] = int(ok.sum())
        return df

    # ---------------------------------------------------------------- random-arm resampling
    def random_block(self, pred, grid):
        """Replicates [B x R x J] (NaN where a trajectory lacks n), points [R x J], trajectory list."""
        mats = [self.mat("R0", pred, n) for n in grid]
        trajs = sorted({t for _, _, tr in mats for t in tr})
        ti = {t: i for i, t in enumerate(trajs)}
        A = np.full((self.E.B, len(trajs), len(grid)), np.nan)
        P = np.full((len(trajs), len(grid)), np.nan)
        for j, (Rp, ep, tr) in enumerate(mats):
            for k, t in enumerate(tr):
                A[:, ti[t], j] = Rp[:, k]
                P[ti[t], j] = ep[k]
        return A, P, trajs

    @staticmethod
    def resample_random(A, P, rng):
        """Joint bootstrap curve per replicate: query weights (already in A) x trajectories resampled with
        replacement within the arm. Budgets observed on different trajectory sets (R by n, D10) get their own draw
        per distinct set. Returns [B x J] trajectory-mean curves."""
        B, R, J = A.shape
        avail = [tuple(np.flatnonzero(np.isfinite(P[:, j]))) for j in range(J)]
        out = np.empty((B, J))
        for sset in set(avail):
            js = [j for j in range(J) if avail[j] == sset]
            idx = np.asarray(sset)
            draw = idx[rng.integers(0, len(idx), size=(B, len(idx)))]      # B x |set|
            for j in js:
                out[:, j] = np.take_along_axis(A[:, :, j], draw, axis=1).mean(axis=1)
        return out

    # ---------------------------------------------------------------- brief, G, G*, matches
    def brief_g(self, fam: pd.DataFrame) -> list[dict]:
        out = []
        Rr, er, trr = self.random_block("pstar", [160, 660])
        for i, f in enumerate(self.fins):
            br = self.w(f, "pstar", 160, "R0", "pstar", 660)
            mt = self.w(f, "pstar", 160, "R0", "pstar", 660, scale=1 - MATCH_MARGIN)
            Rs, es, _ = self.mat(f, "pstar", 160)
            r160, r660 = np.nanmean(er[:, 0]), np.nanmean(er[:, 1])
            s160 = float(np.mean(es))
            G = (s160 - r160) / (r660 - r160)
            Gstar = 1 - MATCH_MARGIN * r660 / (r660 - r160)
            rng = np.random.default_rng(np.random.SeedSequence(self.E.seed, spawn_key=(1, i)))
            cur = self.resample_random(Rr, er, rng)                         # B x 2 (r160, r660)
            sb = np.take_along_axis(Rs, rng.integers(0, Rs.shape[1], size=Rs.shape), axis=1).mean(axis=1)
            den = cur[:, 1] - cur[:, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                Gb = (sb - cur[:, 0]) / den
            # DESIGN: "G undefined if the 97.5% interval of r660 - r160 includes 0". Two readings, both reported:
            # one-sided 97.5% (lower 2.5% quantile > 0) and central two-sided 97.5% (1.25% quantile > 0).
            q0125, q025, q975, q9875 = np.quantile(den, [0.0125, 0.025, 0.975, 0.9875])
            glo, ghi = pct_ci(Gb[np.isfinite(Gb)])
            c1_rej = bool(fam[(fam.finalist == f) & (fam.hypothesis == "c1")].reject.any())
            mt_rej = bool(mt["p"] <= ALPHA)
            out.append({"finalist": f, "s160": s160, "r160": r160, "r660": r660,
                        "brief_delta": br["delta"], "brief_lo": br["lo"], "brief_hi": br["hi"], "brief_p": br["p"],
                        "matches_delta": mt["delta"], "matches_lo": mt["lo"], "matches_hi": mt["hi"],
                        "matches_p": mt["p"], "matches_se": mt["se"], "matches_df": mt["df"],
                        "matches_rejected": mt_rej, "c1_rejected_holm": c1_rej,
                        "G": G, "G_lo": glo, "G_hi": ghi, "G_defined": bool(q025 > 0),
                        "G_defined_two_sided975": bool(q0125 > 0 or q9875 < 0),
                        "gap_r660_r160": r660 - r160, "gap_ci95": [float(q025), float(q975)],
                        "gap_ci975_two_sided": [float(q0125), float(q9875)],
                        "G_star": Gstar, "G_star_discriminating": bool(Gstar >= G_STAR_MIN),
                        "matches_claim_allowed": bool(mt_rej and c1_rej and Gstar >= G_STAR_MIN)})
        return out

    # ---------------------------------------------------------------- headline rule
    def headline(self) -> dict:
        res = {}
        for f in self.fins:
            signs = {}
            for pk in ("pstar", "m1_hgb"):
                for st in ("all", "AO", "TT"):
                    r = self.w(f, pk, 160, "R0", pk, 160, boot=False, st=st)
                    signs[f"{pk}:{st}"] = None if r is None else float(np.sign(r["delta"]))
            vals = [v for v in signs.values() if v is not None]
            res[f] = {"signs": signs, "complete": len(vals) == 6,
                      "headline": len(vals) == 6 and all(v > 0 for v in vals)}
        return res

    # ---------------------------------------------------------------- Q1 n_80 / n_90
    def q1(self, boot_preds=("pstar",)) -> tuple[list, list]:
        nx, curves = [], []
        pool_n = int(self.R.plan["pool_n"])
        preds = sorted(set(self.PT[self.PT.cid == "R0"].pred) - {"ref"})    # ref: paired sensitivity, not a Q1 curve
        for pk in preds:
            grid = sorted(set(self.PT[(self.PT.cid == "R0") & (self.PT.pred == pk)].n))
            P = np.full((0, len(grid)), np.nan)
            d = self.PT[(self.PT.cid == "R0") & (self.PT.pred == pk)]
            piv = d.pivot(index="traj", columns="n", values="ap_all").reindex(columns=grid)
            P = piv.to_numpy()
            curve = np.nanmean(P, axis=0)
            fit = pava(curve)
            rec = {"predictor": pk, "grid": "|".join(map(str, grid)), "anchor_n": grid[-1],
                   "anchor_is_pool_max": grid[-1] == pool_n, "auprc_anchor_measured": float(curve[-1]),
                   "auprc_anchor_pava": float(fit[-1]), "R_at_anchor": int(np.isfinite(P[:, -1]).sum())}
            Mb = None
            if pk in boot_preds:
                A, P2, _ = self.random_block(pk, grid)
                Mb = self.resample_random(A, P2, np.random.default_rng(np.random.SeedSequence(self.E.seed,
                                                                                               spawn_key=(2, 0))))
            for fr in N_FRACS:
                key = f"n{int(round(fr * 100))}"
                for anc, lab in ((curve[-1], ""), (fit[-1], "_pava_anchor")):
                    v, cen = invert_curve(fr * anc, grid, fit)
                    rec[f"{key}{lab}"], rec[f"{key}{lab}_censored"] = v, cen
                if Mb is not None:
                    vals = np.array([invert_curve(fr * Mb[b, -1], grid, pava(Mb[b]))[0] for b in range(len(Mb))])
                    rec[f"{key}_lo"], rec[f"{key}_hi"] = pct_ci(vals)
                    rec[f"{key}_frac_censored"] = float(np.mean((vals == 0) | np.isinf(vals)))
            nx.append(rec)
            for j, n in enumerate(grid):
                dd = d[d.n == n]
                c = {"predictor": pk, "n": n, "R": len(dd), "auprc_mean": float(dd.ap_all.mean()),
                     "auprc_sd": float(dd.ap_all.std(ddof=1)) if len(dd) > 1 else float("nan"),
                     **{f"auprc_{st}": float(dd[f"ap_{st}"].mean()) for st in ("AO", "TT", "onrows")}}
                if Mb is not None:
                    c["ci_lo"], c["ci_hi"] = pct_ci(Mb[:, j])
                curves.append(c)
        return nx, curves

    # ---------------------------------------------------------------- n_req(n)
    def n_req(self, ci_arms) -> list[dict]:
        d = self.PT[(self.PT.cid == "R0") & (self.PT.pred == "pstar")]
        grid = sorted(set(d.n))
        piv = d.pivot(index="traj", columns="n", values="ap_all").reindex(columns=grid).to_numpy()
        fit = pava(np.nanmean(piv, axis=0))
        roles = {c: v["role"] for c, v in self.R.plan["arms"].items()}
        arms = self.fins + [c for c, r in roles.items() if r in ("R1", "static", "det")]
        out = []
        A = Pr = None
        for ia, cid in enumerate(arms):
            ns = sorted(set(self.PT[(self.PT.cid == cid) & (self.PT.pred == "pstar")].n))
            for jn, n in enumerate(ns):
                s = self.pts(cid, "pstar", n)
                v, cen = invert_curve(float(s.ap_all.mean()), grid, fit)
                rec = {"arm": cid, "n": n, "n_req": v, "censored": cen}
                if cid in ci_arms:
                    if A is None:
                        A, Pr, _ = self.random_block("pstar", grid)
                    rng = np.random.default_rng(np.random.SeedSequence(self.E.seed, spawn_key=(3, ia, jn)))
                    Rs, es, _ = self.mat(cid, "pstar", n)
                    cur = self.resample_random(A, Pr, rng)
                    sb = np.take_along_axis(Rs, rng.integers(0, Rs.shape[1], size=Rs.shape), axis=1).mean(axis=1)
                    vals = np.array([invert_curve(sb[b], grid, pava(cur[b]))[0] for b in range(len(sb))])
                    rec["lo"], rec["hi"] = pct_ci(vals)
                    rec["frac_left"] = float(np.mean(vals == 0))
                    rec["frac_right"] = float(np.mean(np.isinf(vals)))
                out.append(rec)
        return out


# ====================================================================================== comparison
class Cmp:
    def __init__(self):
        self.rows = []

    def exact(self, section, item, qty, mine, pipe, tol=EXACT_TOL, rel=False):
        mine, pipe = _num(mine), _num(pipe)
        if pipe is None or mine is None:
            st, d = "MISSING", np.nan
        elif np.isnan(mine) and np.isnan(pipe):
            st, d = "OK", 0.0
        elif np.isinf(mine) or np.isinf(pipe):
            st, d = ("OK" if mine == pipe else "MISMATCH"), (0.0 if mine == pipe else np.inf)
        else:
            d = abs(mine - pipe)
            t = tol * max(1.0, abs(pipe)) if rel else tol
            st = "OK" if d <= t else "MISMATCH"
        self.rows.append(dict(section=section, item=item, quantity=qty, kind="exact", mine=mine, pipeline=pipe,
                              abs_diff=d, tolerance=tol, status=st))

    def variant(self, section, item, qty, mine, pipe):
        """An alternative reading of an ambiguous DESIGN sentence: a difference is shown, not counted as a bug."""
        a, b = _num(mine), _num(pipe)
        if a is not None and b is not None:
            same = (a == b) or (np.isfinite(a) and np.isfinite(b) and abs(a - b) <= EXACT_TOL * max(1.0, abs(b)))
        else:
            same = _canon(mine) == _canon(pipe)
        self.rows.append(dict(section=section, item=item, quantity=qty, kind="variant", mine=mine, pipeline=pipe,
                              abs_diff=np.nan, tolerance=np.nan, status="OK" if same else "INTERPRETATION_DIFF"))

    def decision(self, section, item, qty, mine, pipe):
        st = "MISSING" if pipe is None else ("OK" if _canon(mine) == _canon(pipe) else "DECISION_MISMATCH")
        self.rows.append(dict(section=section, item=item, quantity=qty, kind="decision", mine=str(mine),
                              pipeline=str(pipe), abs_diff=np.nan, tolerance=np.nan, status=st))

    def approx(self, section, item, qty, mine, pipe, tol, how="abs"):
        """how: 'abs' (|diff| <= tol), 'rel' (|diff| <= tol*|pipe|), 'z' (one-sided p on the normal-quantile scale)."""
        mine, pipe = _num(mine), _num(pipe)
        if pipe is None or mine is None:
            st, d = "MISSING", np.nan
        elif np.isinf(mine) or np.isinf(pipe) or np.isnan(mine) or np.isnan(pipe):
            same = (mine == pipe) or (np.isnan(mine) and np.isnan(pipe))
            st, d = ("APPROX_OK" if same else "APPROX_OUTSIDE"), (0.0 if same else np.inf)
        else:
            if how == "z":
                c = lambda p: stats.norm.isf(min(max(p, 1e-15), 1 - 1e-15))  # noqa: E731
                d = abs(c(mine) - c(pipe))
            else:
                d = abs(mine - pipe)
            t = tol * abs(pipe) if how == "rel" else tol
            st = "APPROX_OK" if d <= t else "APPROX_OUTSIDE"
        self.rows.append(dict(section=section, item=item, quantity=qty, kind=f"approx_{how}", mine=mine,
                              pipeline=pipe, abs_diff=d, tolerance=tol, status=st))

    def ci(self, section, item, qty, mine_lo, mine_hi, pipe_lo, pipe_hi, frac=0.25):
        """CI endpoints agree to within frac x the pipeline's CI width (Monte-Carlo difference only)."""
        pl, ph = _num(pipe_lo), _num(pipe_hi)
        wdt = (ph - pl) if (pl is not None and ph is not None and np.isfinite(ph - pl)) else None
        if wdt:
            tol = frac * wdt
        else:                                   # one-sided (censored) interval: 25% of the finite endpoint
            fin = [abs(v) for v in (pl, ph) if v is not None and np.isfinite(v)]
            tol = frac * max(fin) if fin else 0.0
        self.approx(section, item, qty + "_lo", mine_lo, pipe_lo, tol)
        self.approx(section, item, qty + "_hi", mine_hi, pipe_hi, tol)

    def frame(self):
        return pd.DataFrame(self.rows)


def _num(x):
    if x is None:
        return None
    try:
        if isinstance(x, str) and x.strip() == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def _canon(x):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, str) and x in ("True", "False"):
        return x == "True"
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    return x


def compare(R: Run, S: Rescorer, PT, fam, bg, head, nx, curves, nreq) -> pd.DataFrame:
    C = Cmp()
    od = R.out_dir
    # ---- 1. per-unit AUPRC vs the harness value recorded at run time
    for _, u in PT.iterrows():
        C.exact("unit_ap", u.uid, "ap_all_vs_auprc_inline", u.ap_all, u.auprc_inline)
    # ---- 2. per-unit strata AP vs the pipeline's cached unit-weight AP (element 0 of run/boot vectors)
    boot = R.run / "boot"
    if boot.is_dir():
        idx = {}
        for p in boot.glob("*.npy"):
            stem, st, b = p.stem.split("__")
            if st != "twoway":
                idx.setdefault((stem, st), p)
        for _, u in PT.iterrows():
            stem = Path(u.pred_file).stem
            for st in SETS:
                key = (stem, "VV" if st == "TT" and (stem, "VV") in idx else st)
                if key in idx:
                    C.exact("unit_ap_strata", u.uid, f"ap_{st}", u[f"ap_{st}"], float(np.load(idx[key])[0]))
    # ---- 3. confirmatory family
    cf = od / "confirmatory.csv"
    if cf.exists():
        P = pd.read_csv(cf)
        for _, r in fam.iterrows():
            q = P[(P.finalist == r.finalist) & (P.hypothesis == r.hypothesis)]
            if q.empty:
                C.exact("family", f"{r.finalist}:{r.hypothesis}", "row", 0, None)
                continue
            q = q.iloc[0]
            it = f"{r.finalist}:{r.hypothesis}"
            C.exact("family", it, "s_a", r.mean_a, q.s_a)
            C.exact("family", it, "s_b", r.mean_b, q.s_b)
            C.exact("family", it, "delta", r.delta, q.delta)
            C.exact("family", it, "R_a", r.R_a, q.R_a)
            C.exact("family", it, "R_b", r.R_b, q.R_b)
            C.approx("family", it, "V_Q", r.V_Q, q.V_Q, 0.25, "rel")
            C.approx("family", it, "se", r.se, q.se, 0.15, "rel")
            C.approx("family", it, "df", r.df, q.df, 0.5, "rel")
            C.approx("family", it, "p", r.p, q.p, 0.25, "z")
            C.ci("family", it, "ci", r.lo, r.hi, q.lo, q.hi)
            C.decision("family", it, "reject_holm", bool(r.reject), bool(q.reject))
            C.exact("family", it, "m_hypotheses", r.m_hypotheses, q.m_hypotheses)
        ok = np.isfinite(P.p.to_numpy(float))
        adj, rej = holm(P.p.to_numpy(float)[ok])
        for (_, q), a_, r_ in zip(P[ok].iterrows(), adj, rej):
            it = f"{q.finalist}:{q.hypothesis}"
            C.exact("family_holm_on_pipeline_p", it, "p_holm", a_, q.p_holm)
            C.decision("family_holm_on_pipeline_p", it, "reject", bool(r_), bool(q.reject))
    # ---- 4. contrasts.csv: point values of every contrast on the query/two-way bootstrap rows
    cc = od / "contrasts.csv"
    if cc.exists():
        P = pd.read_csv(cc)
        for _, q in P.iterrows():
            st = "TT" if q.set == "VV" else q.set
            if st not in SETS:
                continue
            r = S.w(q.a, q.pred_a, int(q.n_a), q.b, q.pred_b, int(q.n_b), scale=float(q.scale_b), boot=False, st=st)
            it = f"{q.a}/{q.pred_a}@{q.n_a} - {q.scale_b:g}*{q.b}/{q.pred_b}@{q.n_b} [{q.set},{q.bootstrap}]"
            if r is None:
                C.exact("contrasts", it, "row", None, 0)
                continue
            C.exact("contrasts", it, "s_a", r["mean_a"], q.s_a)
            C.exact("contrasts", it, "s_b", r["mean_b"], q.s_b)
            C.exact("contrasts", it, "delta", r["delta"], q.delta)
            C.exact("contrasts", it, "sd_a", r["sd_a"], q.sd_a)
            C.exact("contrasts", it, "sd_b", r["sd_b"], q.sd_b)
    # ---- 5. summary.json: brief, G, G*, matches rule, headline, n_req, Q1
    sj = od / "summary.json"
    if sj.exists():
        J = json.loads(sj.read_text())
        pb = {x["finalist"]: x for x in J.get("brief_G", [])}
        for x in bg:
            f, q = x["finalist"], pb.get(x["finalist"], {})
            for k in ("brief_delta", "matches_delta", "G", "G_star", "gap_r660_r160"):
                C.exact("brief_G", f, k, x[k], q.get(k))
            C.approx("brief_G", f, "brief_p", x["brief_p"], q.get("brief_p"), 0.25, "z")
            C.approx("brief_G", f, "matches_p", x["matches_p"], q.get("matches_p"), 0.25, "z")
            C.ci("brief_G", f, "brief", x["brief_lo"], x["brief_hi"], q.get("brief_lo"), q.get("brief_hi"))
            C.ci("brief_G", f, "matches", x["matches_lo"], x["matches_hi"], q.get("matches_lo"), q.get("matches_hi"))
            C.ci("brief_G", f, "G", x["G_lo"], x["G_hi"], q.get("G_lo"), q.get("G_hi"))
            gc = q.get("gap_ci") or [None, None]
            C.ci("brief_G", f, "gap_ci95", x["gap_ci95"][0], x["gap_ci95"][1], gc[0], gc[1])
            for k in ("matches_rejected", "c1_rejected_holm", "G_defined", "G_star_discriminating",
                      "matches_claim_allowed"):
                C.decision("brief_G", f, k, x[k], q.get(k))
            C.variant("brief_G", f, "G_defined (two-sided 97.5% reading)", x["G_defined_two_sided975"],
                      q.get("G_defined"))
        ph = (J.get("strata_sensitivities") or {}).get("headline", {})
        for f, h in head.items():
            q = ph.get(f, {})
            for k, v in h["signs"].items():
                C.decision("headline", f, f"sign {k}", v, (q.get("signs") or {}).get(k))
            C.decision("headline", f, "headline", h["headline"], q.get("headline"))
        pn = {(x["arm"], int(x["n"])): x for x in J.get("n_req", [])}
        for x in nreq:
            q = pn.get((x["arm"], int(x["n"])), {})
            it = f"{x['arm']}@{x['n']}"
            C.exact("n_req", it, "n_req", x["n_req"], q.get("n_req"), rel=True)
            C.decision("n_req", it, "censored", x["censored"], q.get("censored", ""))
            if "lo" in x:
                C.ci("n_req", it, "ci", x["lo"], x["hi"], q.get("lo"), q.get("hi"))
        px = {x["predictor"]: x for x in J.get("q1_nx", [])}
        for x in nx:
            q = px.get(x["predictor"], {})
            it = x["predictor"]
            C.exact("q1_nx", it, "anchor_n", x["anchor_n"], q.get("anchor_n"))
            C.exact("q1_nx", it, "auprc_anchor_pava", x["auprc_anchor_pava"], q.get("auprc_anchor_pava"))
            for key in ("n80", "n90"):
                # like-for-like: target = frac x PAVA value at the anchor (the pipeline's choice); DESIGN §0 says
                # "measured AUPRC at the full dev pool" -> the measured-anchor reading is reported as a variant
                C.exact("q1_nx", it, key, x[f"{key}_pava_anchor"], q.get(key), rel=True)
                C.variant("q1_nx", it, f"{key} (measured-anchor reading)", x[key], q.get(key))
                C.decision("q1_nx", it, f"{key}_censored", x[f"{key}_censored"], q.get(f"{key}_censored", ""))
                if f"{key}_lo" in x:
                    C.ci("q1_nx", it, key, x[f"{key}_lo"], x[f"{key}_hi"], q.get(f"{key}_lo"), q.get(f"{key}_hi"))
        pc = {(x["predictor"], int(x["n"])): x for x in J.get("q1_curves", [])}
        for x in curves:
            q = pc.get((x["predictor"], int(x["n"])), {})
            it = f"{x['predictor']}@{x['n']}"
            C.exact("q1_curves", it, "R", x["R"], q.get("R"))
            for k in ("auprc_mean", "auprc_sd"):
                C.exact("q1_curves", it, k, x[k], q.get(k))
            if x["predictor"] == "pstar":
                for k in ("auprc_AO", "auprc_TT", "auprc_onrows"):
                    C.exact("q1_curves", it, k, x[k], q.get(k))
            if "ci_lo" in x:
                C.ci("q1_curves", it, "ci", x["ci_lo"], x["ci_hi"], q.get("ci_lo"), q.get("ci_hi"))
    return C.frame()


# ====================================================================================== diagnostic
def pipeline_weights_check(R: Run, PT: pd.DataFrame, cache: Path, n_sample: int = 40, seed: int = 0):
    """DIAGNOSTIC (not part of the independent estimate): rebuild the pipeline's query-weight matrix from its
    documented seed scheme (config/seeds.yaml: bootstrap = SeedSequence(<mode bootstrap seed>, spawn_key=(crc32(name),)),
    name 'query_weights', B from summary.json) with THIS module's stratified_query_weights, then
    (1) compare our replicate APs to the pipeline's cached replicate vectors (run/boot) element by element and
    (2) recompute the family under those weights: V_Q, SE, nu, p, CI and Holm should then match the pipeline exactly.
    Separates Monte-Carlo differences (own seed) from definitional ones."""
    sj = R.out_dir / "summary.json"
    if not sj.exists() or not (R.run / "boot").is_dir():
        return None, None, "no summary.json or run/boot"
    B = int(json.loads(sj.read_text())["B"])
    key = {"final": "final_bootstrap_seed", "dryrun": "final_dryrun_bootstrap_seed"}.get(R.plan.get("mode"))
    m = re.search(rf"^{key}:\s*(\d+)", (R.root / "config" / "seeds.yaml").read_text(), flags=re.M) if key else None
    if not m:
        return None, None, f"bootstrap seed for mode {R.plan.get('mode')} not found"
    ss = np.random.SeedSequence(int(m.group(1)), spawn_key=(zlib.crc32(b"query_weights"),))
    W = stratified_query_weights(R.z["query_strata"], B, np.random.default_rng(ss))
    E2 = Engine(R, B, int(m.group(1)), cache, W=W, tag=f"pipeW_B{B}")
    S2 = Rescorer(R, E2, PT)
    fam2 = S2.family()
    # replicate vectors: every vector the family uses + a random sample of the other cached ones
    by_stem = {Path(f).stem: f for f in PT.pred_file}
    need = set()
    for f in S2.fins:
        for n in (160,):
            need |= {(Path(x).stem, "all") for x in S2.pts(f, "pstar", n).pred_file}
    for n in (160, 490, 660):
        need |= {(Path(x).stem, "all") for x in S2.pts("R0", "pstar", n).pred_file}
    cached = []
    for p in sorted((R.run / "boot").glob(f"*__B{B}.npy")):
        stem, st, _ = p.stem.split("__")
        if stem in by_stem and st != "twoway":
            cached.append((stem, st, p))
    rng = np.random.default_rng(seed)
    others = [c for c in cached if (c[0], c[1]) not in need]
    pick = [others[i] for i in rng.choice(len(others), size=min(n_sample, len(others)), replace=False)] if others else []
    recs = []
    for stem, st, p in [c for c in cached if (c[0], c[1]) in need] + pick:
        mine = E2.replicate(R.root / by_stem[stem], "TT" if st == "VV" else st)
        theirs = np.load(p)
        recs.append({"vector": p.name, "B": len(theirs) - 1, "max_abs_diff": float(np.max(np.abs(mine - theirs)))})
    return fam2, pd.DataFrame(recs), f"seed {key} = {m.group(1)}, spawn_key crc32('query_weights'), B = {B}"


# ====================================================================================== main
def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if np.isfinite(f) else ("inf" if f > 0 else ("-inf" if f < 0 else None))
    return o


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2 ** 30 if sys.platform == "darwin" else r / 2 ** 20


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-dir", required=True, help="results/final or results/final_dryrun (or its run/ subdir)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=DEFAULT_B)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--root", default=str(ROOT), help="project root that unit pred_file paths are relative to")
    ap.add_argument("--cache", default=None, help="replicate-vector cache (default <out>/cache)")
    ap.add_argument("--nreq-ci-arms", default="finalists", help="'finalists', 'none' or comma list")
    ap.add_argument("--no-pipeline-weights-check", action="store_true",
                    help="skip the diagnostic that rebuilds the pipeline's weights and reproduces its replicates")
    a = ap.parse_args(argv)
    t0 = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if any(m.split(".")[0] in ("hsl", "ysl") for m in sys.modules):
        raise RuntimeError("project package imported: the rescorer must stay independent")

    R = Run(Path(a.run_dir), Path(a.root))
    checks = R.integrity()
    print(f"[rescore] {len(R.units)} units; integrity: {sum(c['ok'] for c in checks)}/{len(checks)} checks pass",
          flush=True)
    PT = point_table(R)
    PT.to_csv(out / "unit_ap.csv", index=False)
    print(f"[rescore] point APs: {len(PT)} units x {len(SETS)} sets ({time.time() - t0:.0f} s)", flush=True)
    # sklearn spot-check of the point AP on a sample of units
    from sklearn.metrics import average_precision_score
    rng = np.random.default_rng(0)
    sk = []
    for i in rng.choice(len(PT), size=min(10, len(PT)), replace=False):
        s = R.scores(R.root / PT.pred_file.iloc[i])
        sk.append(abs(average_precision_score(R.y, s) - PT.ap_all.iloc[i]))
    checks.append({"check": "point AP == sklearn average_precision_score (10 random units)",
                   "ok": bool(max(sk) <= 1e-12), "detail": f"max |diff| = {max(sk):.2e}"})

    E = Engine(R, a.B, a.seed, Path(a.cache) if a.cache else out / "cache")
    S = Rescorer(R, E, PT)
    fam = S.family()
    print(f"[rescore] family done ({time.time() - t0:.0f} s, {E.n_new} replicate vectors)", flush=True)
    bg = S.brief_g(fam)
    head = S.headline()
    nx, curves = S.q1()
    ci_arms = {"finalists": S.fins, "none": []}.get(a.nreq_ci_arms, a.nreq_ci_arms.split(","))
    nreq = S.n_req(set(ci_arms))
    print(f"[rescore] analyses done ({time.time() - t0:.0f} s, {E.n_new} replicate vectors)", flush=True)

    fam.to_csv(out / "confirmatory_rescored.csv", index=False)
    pd.DataFrame(bg).to_csv(out / "brief_G_rescored.csv", index=False)
    pd.DataFrame(nx).to_csv(out / "q1_nx_rescored.csv", index=False)
    pd.DataFrame(curves).to_csv(out / "q1_curves_rescored.csv", index=False)
    pd.DataFrame(nreq).to_csv(out / "n_req_rescored.csv", index=False)
    cmp = compare(R, S, PT, fam, bg, head, nx, curves, nreq)
    pw_note = "skipped"
    if not a.no_pipeline_weights_check:
        fam2, vecs, pw_note = pipeline_weights_check(R, PT, E.cache)
        if fam2 is not None:
            fam2.to_csv(out / "confirmatory_pipeline_weights.csv", index=False)
            vecs.to_csv(out / "pipeline_weights_vectors.csv", index=False)
            C2 = Cmp()
            P = pd.read_csv(R.out_dir / "confirmatory.csv")
            for _, r in fam2.iterrows():
                q = P[(P.finalist == r.finalist) & (P.hypothesis == r.hypothesis)].iloc[0]
                it = f"{r.finalist}:{r.hypothesis}"
                for k_m, k_p in (("V_Q", "V_Q"), ("se", "se"), ("df", "df"), ("p", "p"), ("lo", "lo"), ("hi", "hi"),
                                 ("p_holm", "p_holm")):
                    C2.exact("family_pipeline_weights", it, k_p, r[k_m], q[k_p], tol=1e-9, rel=True)
                C2.decision("family_pipeline_weights", it, "reject", bool(r.reject), bool(q.reject))
            for _, v in vecs.iterrows():
                C2.exact("replicates_pipeline_weights", v.vector, "max_abs_diff_over_B", v.max_abs_diff, 0.0, tol=1e-12)
            cmp = pd.concat([cmp, C2.frame()], ignore_index=True)
        print(f"[rescore] pipeline-weights diagnostic: {pw_note} ({time.time() - t0:.0f} s)", flush=True)
    cmp.to_csv(out / "comparison.csv", index=False)
    counts = cmp.groupby(["section", "status"]).size().unstack(fill_value=0)
    flagged = cmp[~cmp.status.isin(["OK", "APPROX_OK"])]
    summ = {"run_dir": str(R.out_dir), "mode": R.plan.get("mode"), "split": R.plan.get("split"),
            "B": a.B, "seed": a.seed, "seed_note": "own SeedSequence; children: query weights, then spawn_key (1,i) "
            "G finalist i, (2,0) Q1, (3,arm,n) n_req", "n_units": len(R.units), "n_units_scored": len(PT),
            "integrity_checks": checks, "family": fam.to_dict("records"), "brief_G": bg, "headline": head,
            "q1_nx": nx, "n_req": nreq,
            "comparison_counts": {s: counts.loc[s].to_dict() for s in counts.index},
            "n_flagged": int(len(flagged)), "flagged": flagged.to_dict("records"),
            "pipeline_weights_diagnostic": pw_note,
            "replicate_vectors_computed": E.n_new, "wall_s": time.time() - t0, "peak_rss_gb": peak_rss_gb(),
            "imports": sorted({m.split(".")[0] for m in sys.modules if m.split(".")[0] in ("hsl", "ysl")})}
    (out / "rescore_summary.json").write_text(json.dumps(_jsonable(summ), indent=1))
    print(counts.to_string())
    bad_checks = [c for c in checks if not c["ok"]]
    for c in bad_checks:
        print(f"[rescore] INTEGRITY FAIL: {c['check']}: {c['detail']}")
    if len(flagged):
        print(f"[rescore] {len(flagged)} flagged comparison rows:")
        print(flagged.to_string(max_colwidth=60))
    print(f"[rescore] family (Holm, alpha={ALPHA}):")
    print(fam[["finalist", "hypothesis", "delta", "se", "df", "p", "p_holm", "reject"]].to_string())
    print(f"[rescore] wall {time.time() - t0:.0f} s, peak RSS {peak_rss_gb():.2f} GB -> {out}")
    return 0 if not bad_checks and not len(flagged[flagged.status.isin(["MISMATCH", "DECISION_MISMATCH"])]) else 1


if __name__ == "__main__":
    sys.exit(main())
