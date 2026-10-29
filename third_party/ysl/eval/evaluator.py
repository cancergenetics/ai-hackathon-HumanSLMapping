"""Evaluator: the ONLY code that reads held-out labels (DESIGN.md §8–§10; r1 agent-failure §2.1). PI-owned, frozen.

    ev = Evaluator(bundle, split)          # bundle = ysl.data.dense.load_bundle(dataset[, unblind=True for final])
    m  = ev.evaluate(S)                    # S: float [len(split.test) x n_arrays] gene-level FULL-GRID scores for the
                                           #    sorted held-out ORFs; broadcast to every strain of an ORF; scored only on
                                           #    observed cells. Rejects wrong shape, NaN/inf anywhere, constant scores.
    df = ev.keyed_predictions(S)           # [query_strain_id, array_strain_id, score] over observed test cells
    ev.evaluate_keyed(df)                  # same metrics from keyed predictions (exact index-set check)
    ev.paired(S_a, S_b)                    # paired query-ORF cluster bootstrap of ΔAUPRC
    ev.oracle_diagnostic(S, array_rate)    # DEV ONLY: true held-out query positive rate as propensity

Evaluation unit: every observed (query strain, array strain) cell of every held-out query ORF. Bootstrap clusters =
query ORFs (all strains of an ORF together); arrays are not resampled. Labels never leave this module: there is no
public accessor returning them, and nothing in ysl.models may import this module.

Strata (DESIGN §8; boolean masks over cells; each reported with prevalence and lift):
  recip_available / recip_unavailable : x has an array strain, a's ORF is a screened query, and that training cell
                                        (a's ORF as query, x's array strain) is observed
  x_on_array / x_off_array            : x's ORF has an array strain
  x_on_array_no_recip                 : on array AND no reciprocal (the dominant stratum)
  same_chromosome / diff_chromosome   : SGD chromosome of x and a (current names); pairs with unknown coordinates are
                                        in neither (count reported)
  adjacent_screened_query / no_adjacent_screened_query : a screened query ORF lies within ADJ_BP of x (edge-to-edge,
                                        same chromosome; SGD coordinates)
  paralog_pair                        : (x, a) is an Ensembl Compara paralog pair (any paralog type; features table)
  lockbox / non_lockbox               : x is a lockbox ORF (final splits)
  cellmap only: core_strain / cellmap_only_strain (query strain present in S1 or not), supp_strain / non_supp_strain
"""
from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd

from ysl.eval import metrics as M

log = logging.getLogger("ysl.eval")
ADJ_BP = 3000
ORACLE_SMOOTH = 0.5  # Jeffreys-style smoothing of the oracle query rate: (k + 0.5) / (n + 1)


class PredictionError(ValueError):
    pass


def _logit(p):
    p = np.asarray(p, dtype=np.float64)
    return np.log(p) - np.log1p(-p)


class Evaluator:
    def __init__(self, bundle, split, n_boot: int = 1000, boot_seed: int = 0, strata: bool = True):
        if split.dataset and split.dataset != bundle.dataset:
            raise ValueError(f"split dataset {split.dataset} != bundle dataset {bundle.dataset}")
        lb_test = np.isin(split.test, split.lockbox)
        if split.kind == "dev" and lb_test.any():
            raise PermissionError("dev split contains lockbox ORFs in its test set")
        if lb_test.any() and not bundle.unblinded:
            raise PermissionError("lockbox ORFs are evaluated only with an unblinded bundle (final evaluation)")
        self.split = split
        self.dataset = bundle.dataset
        self.n_boot, self.boot_seed = int(n_boot), int(boot_seed)
        self.n_arrays = bundle.n_arrays
        self.array_strains = bundle.array_strains
        # explicit alignment contract: prediction rows are split.test in this exact (strictly increasing) order
        if len(split.test) > 1 and not (split.test[1:] > split.test[:-1]).all():
            raise AssertionError("split.test must be strictly increasing (prediction row order)")
        R = bundle.strain_rows(split.test)                 # sorted strain rows of the held-out ORFs
        tpos = np.searchsorted(split.test, bundle.strain_orf[R])
        if not np.array_equal(split.test[tpos], bundle.strain_orf[R]):
            raise AssertionError("strain -> test ORF alignment failed")
        if len(np.unique(tpos)) != len(split.test):
            raise AssertionError("some held-out ORFs have no strain rows")
        Yr = bundle.Y[R]
        obs = Yr >= 0
        ri, ci = np.nonzero(obs)
        self.strain_ids = bundle.strain_ids[R]
        self._cell_strain = ri.astype(np.int32)
        self._cell_t = tpos[ri].astype(np.int32)
        self._cell_a = ci.astype(np.int32)
        self.__y = Yr[obs] == 1
        self.n_cells = int(len(self.__y))
        # global ORF code (position in bundle.all_query_orfs) per cell, for joint multi-split bootstraps
        self.all_query_orfs = bundle.all_query_orfs
        gpos = np.searchsorted(bundle.all_query_orfs, split.test)
        if not np.array_equal(bundle.all_query_orfs[np.minimum(gpos, len(bundle.all_query_orfs) - 1)], split.test):
            raise AssertionError("held-out ORFs missing from bundle.all_query_orfs")
        self._cell_orf_global = gpos[self._cell_t].astype(np.int32)
        self._cell_lockbox = np.isin(split.test, split.lockbox)[self._cell_t]
        self.skipped_strata: dict = {}
        self._strata = self._build_strata(bundle, split, R) if strata else {}

    # ------------------------------------------------------------------ strata
    def _build_strata(self, bundle, split, R) -> dict:
        st = {}
        orf_to_col = bundle.orf_to_col
        col_x = np.array([orf_to_col.get(o, -1) for o in split.test], dtype=np.int64)
        on = col_x[self._cell_t] >= 0
        # reciprocal availability from TRAINING rows only
        tr_rows = bundle.gene_rows(split.train)
        arr_tr = np.full(bundle.n_arrays, -1, dtype=np.int64)
        for i, o in enumerate(split.train):
            j = orf_to_col.get(o, -1)
            if j >= 0:
                arr_tr[j] = i
        obs_train = bundle.Y[tr_rows] >= 0
        ta = arr_tr[self._cell_a]
        cx = col_x[self._cell_t]
        rec = (ta >= 0) & (cx >= 0)
        rec[rec] = obs_train[ta[rec], cx[rec]]
        st["recip_available"], st["recip_unavailable"] = rec, ~rec
        st["x_on_array"], st["x_off_array"] = on, ~on
        st["x_on_array_no_recip"] = on & ~rec
        # genome coordinates (current SGD names)
        try:
            from ysl.data.sgd import orf_coordinates
            co = orf_coordinates()
            cur = bundle.current_name
            chrom_x = pd.Series([cur[o] for o in split.test]).map(co["chromosome"]).to_numpy(dtype=float)
            chrom_a = pd.Series(bundle.array_orfs_current).map(co["chromosome"]).to_numpy(dtype=float)
            cxv, cav = chrom_x[self._cell_t], chrom_a[self._cell_a]
            known = ~np.isnan(cxv) & ~np.isnan(cav)
            st["same_chromosome"] = known & (cxv == cav)
            st["diff_chromosome"] = known & (cxv != cav)
            self.skipped_strata["n_cells_unknown_chromosome"] = int((~known).sum())
            adj = self._adjacent(split, co, cur)
            st["adjacent_screened_query"] = adj[self._cell_t]
            st["no_adjacent_screened_query"] = ~adj[self._cell_t]
        except Exception as e:  # noqa: BLE001 — logged and recorded, never silent
            msg = f"coordinate strata skipped: {type(e).__name__}: {e}"
            warnings.warn(msg)
            log.warning(msg)
            self.skipped_strata["coordinates"] = msg
        try:
            st["paralog_pair"] = self._paralog_cells(bundle, split)
        except Exception as e:  # noqa: BLE001
            msg = f"paralog stratum skipped: {type(e).__name__}: {e}"
            warnings.warn(msg)
            log.warning(msg)
            self.skipped_strata["paralog_pair"] = msg
        lb = np.isin(split.test, split.lockbox)[self._cell_t]
        if lb.any():
            st["lockbox"], st["non_lockbox"] = lb, ~lb
        if bundle.dataset == "cellmap":
            from ysl.data import costanzo as C
            s1_strains = set(pd.read_parquet(C.require(C.queries_path("s1")))["query_strain"].astype(str))
            core = np.array([s in s1_strains for s in self.strain_ids])[self._cell_strain]
            st["core_strain"], st["cellmap_only_strain"] = core, ~core
            supp = bundle.strain_is_supp[R][self._cell_strain]
            st["supp_strain"], st["non_supp_strain"] = supp, ~supp
        return st

    @staticmethod
    def _adjacent(split, co, cur) -> np.ndarray:
        """Per held-out ORF: does a screened query ORF lie within ADJ_BP (edge-to-edge, same chromosome)?"""
        tr = co.reindex([cur[o] for o in split.train]).dropna(subset=["chromosome"])
        out = np.zeros(len(split.test), dtype=bool)
        by_chr = {c: g[["lo", "hi"]].to_numpy() for c, g in tr.groupby("chromosome")}
        te = co.reindex([cur[o] for o in split.test])
        for i, (c, lo, hi) in enumerate(zip(te["chromosome"], te["lo"], te["hi"])):
            if pd.isna(c) or c not in by_chr:
                continue
            iv = by_chr[c]
            gap = np.maximum(0, np.maximum(iv[:, 0] - hi, lo - iv[:, 1]))
            out[i] = bool((gap < ADJ_BP).any())
        return out

    def _paralog_cells(self, bundle, split) -> np.ndarray:
        from ysl.data.features import load_paralogs
        para = load_paralogs()
        partners: dict = {}
        for a, b in zip(para["orf_a"].astype(str), para["orf_b"].astype(str)):
            partners.setdefault(a, set()).add(b)
            partners.setdefault(b, set()).add(a)
        cur = bundle.current_name
        acur = bundle.array_orfs_current
        col_of_cur: dict = {}
        for j, c in enumerate(acur):
            col_of_cur.setdefault(c, []).append(j)
        rows, cols = [], []
        for i, o in enumerate(split.test):
            for p in partners.get(cur[o], ()):
                for j in col_of_cur.get(p, ()):
                    rows.append(i)
                    cols.append(j)
        mat = np.zeros((len(split.test), bundle.n_arrays), dtype=bool)
        if rows:
            mat[np.array(rows), np.array(cols)] = True
        return mat[self._cell_t, self._cell_a]

    # ------------------------------------------------------------------ predictions
    def _scores(self, S) -> np.ndarray:
        S = np.asarray(S)
        if S.shape != (len(self.split.test), self.n_arrays):
            raise PredictionError(f"predictions must be the full grid {(len(self.split.test), self.n_arrays)}; "
                                  f"got {S.shape}")
        if not np.issubdtype(S.dtype, np.floating):
            raise PredictionError(f"predictions must be float, got {S.dtype}")
        if not np.isfinite(S).all():
            raise PredictionError("predictions contain NaN/inf (the full grid must be finite; no silent fallbacks)")
        s = S[self._cell_t, self._cell_a].astype(np.float64)
        if s.min() == s.max():
            raise PredictionError("predictions are constant over the scored cells")
        return s

    def flat_scores(self, S) -> np.ndarray:
        """Validated scores of the full-grid matrix S at the scored cells (cell order of this evaluator)."""
        return self._scores(S)

    def keyed_predictions(self, S) -> pd.DataFrame:
        s = self._scores(S)
        # saved at the method's own precision (float32 stays float32) so rescoring sees exactly the evaluated scores
        keep32 = np.asarray(S).dtype == np.float32
        return pd.DataFrame({
            "query_strain_id": pd.Categorical.from_codes(self._cell_strain, categories=self.strain_ids),
            "array_strain_id": pd.Categorical.from_codes(self._cell_a, categories=self.array_strains),
            "score": s.astype(np.float32) if keep32 else s})

    def evaluate_keyed(self, df: pd.DataFrame, **kw) -> dict:
        """Metrics from keyed predictions; the key set must equal the observed held-out cells exactly."""
        return self._metrics(self.scores_from_keyed(df), **kw)

    def scores_from_keyed(self, df: pd.DataFrame) -> np.ndarray:
        """Flat scores (in this evaluator's cell order) from keyed predictions, with an exact index-set check."""
        if len(df) != self.n_cells:
            raise PredictionError(f"{len(df)} keyed predictions for {self.n_cells} scored cells")
        si = pd.Index(self.strain_ids).get_indexer(df["query_strain_id"].astype(str))
        ai = pd.Index(self.array_strains).get_indexer(df["array_strain_id"].astype(str))
        if (si < 0).any() or (ai < 0).any():
            raise PredictionError("keyed predictions contain unknown strain IDs")
        key = si.astype(np.int64) * self.n_arrays + ai
        if len(np.unique(key)) != len(key):
            raise PredictionError("keyed predictions contain duplicate (query strain, array strain) keys")
        want = self._cell_strain.astype(np.int64) * self.n_arrays + self._cell_a
        pos = pd.Index(key).get_indexer(want)
        if (pos < 0).any():
            raise PredictionError("keyed prediction index != observed held-out cells (missing/extra/duplicate)")
        s = df["score"].to_numpy(dtype=np.float64)[pos]
        if not np.isfinite(s).all() or s.min() == s.max():
            raise PredictionError("keyed predictions are non-finite or constant")
        return s

    # ------------------------------------------------------------------ metrics
    def evaluate(self, S, n_boot: int | None = None, macro: bool = True) -> dict:
        return self._metrics(self._scores(S), n_boot=n_boot, macro=macro)

    def _metrics(self, s, n_boot: int | None = None, macro: bool = True) -> dict:
        y = self.__y
        n_boot = self.n_boot if n_boot is None else int(n_boot)
        out = {"auprc": M.pooled_auprc(y, s), "auroc": M.pooled_auroc(y, s), "prevalence": float(y.mean()),
               "n": self.n_cells, "n_pos": int(y.sum()), "n_test_orfs": int(len(self.split.test)),
               "n_test_strains": int(len(self.strain_ids)), "precision_at_k": M.precision_at_k(y, s)}
        out["lift"] = out["auprc"] / out["prevalence"]
        if macro:
            out.update(M.macro_auprc(y, s, self._cell_strain))
        if self._strata:
            out["strata"] = M.stratified(y, s, self._strata)
        if self.skipped_strata:
            out["strata_notes"] = dict(self.skipped_strata)
        if n_boot > 0:
            out["auprc_bootstrap"] = M.cluster_bootstrap(y, s, self._cell_t, n_boot=n_boot, seed=self.boot_seed)
        return out

    def paired(self, S_a, S_b, n_boot: int | None = None) -> dict:
        n_boot = self.n_boot if n_boot is None else int(n_boot)
        return M.paired_bootstrap(self.__y, self._scores(S_a), self._scores(S_b), self._cell_t, n_boot=n_boot,
                                  seed=self.boot_seed)

    def paired_flat(self, s_a, s_b, n_boot: int | None = None) -> dict:
        """Paired bootstrap from flat scores in cell order (e.g. from scores_from_keyed)."""
        n_boot = self.n_boot if n_boot is None else int(n_boot)
        return M.paired_bootstrap(self.__y, s_a, s_b, self._cell_t, n_boot=n_boot, seed=self.boot_seed)

    def stratum_sizes(self) -> dict:
        return {k: int(v.sum()) for k, v in self._strata.items()}

    # ------------------------------------------------------------------ oracle (dev only)
    def oracle_diagnostic(self, S, array_rate, betas=(1.0,)) -> dict:
        """DEV ONLY. True held-out query positive rate r_x (smoothed, all strains of the ORF) used as propensity:
          oracle_degree_product : array_rate(a) * r_x                       (B2-style product)
          method                : S
          method_plus_oracle    : S + beta * logit(r_x)                     (logit-sum; S should be on a logit scale)
          method_recentred      : S - rowmean_x(S) + logit(r_x)             (method's per-query level replaced)
        Roadmap rule (DESIGN §8): gain = AUPRC(method_plus_oracle, beta=1) - AUPRC(method) >= 0.01 -> prioritise
        propensity modelling."""
        if self.split.kind != "dev":
            raise PermissionError("oracle diagnostics are dev-only")
        y = self.__y
        s = self._scores(S)
        k = np.bincount(self._cell_t, weights=y.astype(float), minlength=len(self.split.test))
        n = np.bincount(self._cell_t, minlength=len(self.split.test))
        r = (k + ORACLE_SMOOTH) / (n + 2 * ORACLE_SMOOTH)
        lr = _logit(r)[self._cell_t]
        ar = np.asarray(array_rate, dtype=np.float64)
        if ar.shape != (self.n_arrays,):
            raise ValueError("array_rate must have one value per array")
        S = np.asarray(S, dtype=np.float64)
        rowmean = S.mean(1)[self._cell_t]
        res = {"method": M.pooled_auprc(y, s),
               "oracle_degree_product": M.pooled_auprc(y, ar[self._cell_a] * r[self._cell_t]),
               "method_recentred": M.pooled_auprc(y, s - rowmean + lr)}
        for b in betas:
            res[f"method_plus_oracle_beta{b:g}"] = M.pooled_auprc(y, s + b * lr)
        res["gain_logit_sum"] = res["method_plus_oracle_beta1"] - res["method"] if 1.0 in betas else None
        res["roadmap_propensity_first"] = bool(res["gain_logit_sum"] is not None and res["gain_logit_sum"] >= 0.01)
        res["note"] = ("dev-only oracle using held-out labels; not a method result. method_plus_oracle assumes S is "
                       "on a logit scale (LR decision function).")
        return res


# ----------------------------------------------------------------------------- multi-split (joint) bootstraps
def _same_universe(evs) -> int:
    u = evs[0].all_query_orfs
    if any(not np.array_equal(e.all_query_orfs, u) for e in evs):
        raise ValueError("evaluators come from different datasets")
    return len(u)


def _subset(e, subset):
    if subset is None:
        return np.ones(e.n_cells, dtype=bool)
    if subset == "lockbox":
        if not e._cell_lockbox.any():
            raise ValueError("no lockbox test cells in this split")
        return e._cell_lockbox
    raise ValueError(f"unknown subset {subset!r}")


def joint_bootstrap(evaluators: list, flat_scores: list, n_boot: int = 1000, seed: int = 0,
                    subset: str | None = None) -> dict:
    """Mean pooled AUPRC over splits with ONE ORF-weight vector per replicate applied in every split (DESIGN §4/§9).
    flat_scores[i] are cell-order scores for evaluators[i] (Evaluator.flat_scores / scores_from_keyed).
    subset="lockbox": only cells of lockbox query ORFs (the confirmatory endpoint); the resampled pool is then the
    lockbox ORFs that are held out in at least one split."""
    n = _same_universe(evaluators)
    parts = []
    for e, s in zip(evaluators, flat_scores):
        m = _subset(e, subset)
        parts.append((e._Evaluator__y[m], np.asarray(s, dtype=np.float64)[m], e._cell_orf_global[m]))
    return M.joint_bootstrap(parts, n, n_boot=n_boot, seed=seed) | {"subset": subset or "all"}


def joint_paired_bootstrap(evaluators: list, scores_a: list, scores_b: list, n_boot: int = 1000,
                           seed: int = 0, subset: str | None = None) -> dict:
    """Joint multi-split paired bootstrap of Δ mean AUPRC (method a - method b); subset as in joint_bootstrap."""
    n = _same_universe(evaluators)
    parts = []
    for e, a, b in zip(evaluators, scores_a, scores_b):
        m = _subset(e, subset)
        parts.append((e._Evaluator__y[m], np.asarray(a, dtype=np.float64)[m], np.asarray(b, dtype=np.float64)[m],
                      e._cell_orf_global[m]))
    return M.joint_paired_bootstrap(parts, n, n_boot=n_boot, seed=seed) | {"subset": subset or "all"}


# ----------------------------------------------------------------------------- confirmatory endpoint (DESIGN §4, §9)
def confirmatory_analysis(evaluators: list, scores: dict, candidates: list, comparator: str, n_boot: int = 1000,
                          seed: int = 0, alpha: float = 0.05) -> dict:
    """Pre-registered final analysis over the final splits (one Evaluator per split, unblinded bundle).

    scores = {method: [flat cell-order scores per split]} for every method run in the final batch.
    PRIMARY (confirmatory): pooled AUPRC on LOCKBOX-ORF test rows, mean over splits, with one joint ORF bootstrap
      (lockbox ORFs resampled once per replicate, same weights in every split); each of <= 3 pre-declared candidates
      vs the pre-declared comparator by a joint paired bootstrap of Δ mean AUPRC; two-sided bootstrap p-values,
      Holm-corrected across the candidates.
    SECONDARY: the same on ALL test rows (the brief's metric; ~80% of these queries were used for dev selection).
    Prevalence and lift are reported next to every AUPRC."""
    candidates = list(candidates)
    if not 1 <= len(candidates) <= 3:
        raise ValueError("1-3 pre-declared candidates")
    if comparator in candidates or comparator not in scores or not set(candidates) <= set(scores):
        raise ValueError("comparator/candidates must be distinct methods present in `scores`")
    if any(e.split.kind != "final" for e in evaluators):
        raise PermissionError("the confirmatory analysis runs on final splits only")
    if any(len(v) != len(evaluators) for v in scores.values()):
        raise ValueError("every method needs scores for every split")
    out = {"candidates": candidates, "comparator": comparator, "n_splits": len(evaluators), "n_boot": n_boot,
           "alpha": alpha, "seeds": [e.split.seed for e in evaluators]}
    for subset, key in (("lockbox", "primary_lockbox_rows"), (None, "secondary_all_test_rows")):
        masks = [_subset(e, subset) for e in evaluators]
        prev = [float(e._Evaluator__y[m].mean()) for e, m in zip(evaluators, masks)]
        res = {"prevalence_mean": float(np.mean(prev)), "prevalence_per_split": prev,
               "n_cells_per_split": [int(m.sum()) for m in masks],
               "n_orfs_per_split": [int(len(np.unique(e._cell_orf_global[m]))) for e, m in zip(evaluators, masks)],
               "methods": {}, "vs_comparator": {}}
        for name, ss in scores.items():
            per = [M.pooled_auprc(e._Evaluator__y[m], np.asarray(s, dtype=np.float64)[m])
                   for e, s, m in zip(evaluators, ss, masks)]
            jb = joint_bootstrap(evaluators, ss, n_boot=n_boot, seed=seed, subset=subset)
            res["methods"][name] = {"auprc_mean": float(np.mean(per)), "auprc_sd_over_splits": float(np.std(per, ddof=1))
                                    if len(per) > 1 else 0.0, "auprc_per_split": per,
                                    "lift": float(np.mean(per)) / res["prevalence_mean"], "joint_bootstrap": jb}
        for c in candidates:
            res["vs_comparator"][c] = joint_paired_bootstrap(evaluators, scores[c], scores[comparator], n_boot=n_boot,
                                                             seed=seed, subset=subset)
        adj = M.holm({c: res["vs_comparator"][c]["p_two_sided"] for c in candidates})
        for c in candidates:
            res["vs_comparator"][c]["p_holm"] = adj[c]
            res["vs_comparator"][c]["significant_holm"] = bool(adj[c] < alpha)
        out[key] = res
    return out
