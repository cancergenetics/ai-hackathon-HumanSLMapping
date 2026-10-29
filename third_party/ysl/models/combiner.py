"""Integrated combiners (DESIGN.md §6 rung 3): `stack_lr_all` (cross-fitted LR) and `m1_hgb` (M1, gradient boosting)
over ONE feature set:

  core    stack_lr inputs (baselines.pair_features): log array rate, log query propensity (x's array-side rate),
          on-array flag, reciprocal (-eps, label, -signed log10 P), 3 kNN scores, array SMF
  wr      within-row completion scores = stack_lr_wr2's set (foldin2.wr2_features: multi-kernel fold-in 'mk',
          tied low-rank, refined kNN, 2-hop pp/pg, external-only fold-in for off-array targets 'ext_off'), each
          divided by its SD over the view's targets (as in stack_lr_wr / stack_lr_wr2). params["wr"] may name any
          completion.WR_FEATURES / foldin2.WR2_FEATURES entry (the stack_lr_wr set was used for the stage-2/3 tuning)
  prop    predicted query propensity (propensity.propensity_predict, nested cross-fitted sub-model; LR also gets
          prop x off-array)
  coex    coex_z, coex_mr                      apms    apms_edge, apms_shared, apms_jaccard
  para    para_any, para_within, para_identity, esm_cos
  genome  gen_same_chr, gen_dist_kb            rest    phylo_ident_corr, dom_ipr_jaccard, abund_sum, abund_absdiff
          (external pair features, src/ysl/features/pairs.py; static, identical at train and test time)
`params["drop"]` removes groups (ablation); `core` is always kept.

One feature function (`view_features` + `design`) for training and test. Training cells = every observed cell of
every screened query; the features of fold k's rows come from view.without_queries(fold k) (K = ctx.n_folds = 10).
Additive statistics (array rate, query propensity): `additive="loo"` overrides them with exact leave-one-query-out
values (LR only; stack_lr_wr convention); `additive="kfold"` keeps the fold view's K-fold values — REQUIRED for
trees (decision log 2026-10-26: exact-LOO rates leak the training label to a tree). m1_hgb refuses "loo".

Test-time views (`test_view`):
  full     test features from the full view (T screened queries)
  matched  test features from the K matched views T \\ fold k (same targets; propensity.matched_target_view), so the
           count-based features (rates, kNN, reciprocal availability, completion scores) see exactly as many screened
           queries as the training features did; predictions (log-odds) averaged over k.
Model pooling (`ensemble`, m1_hgb): pooled = one model on all folds' rows; fold = matched-view ensemble (DESIGN §6):
model k fitted on fold k's rows only and scored on matched view k's test features; log-odds averaged (requires
test_view="matched").
The propensity feature at test time always comes from propensity_predict(full view), whose own test features are
already matched-view averages (results/m4_propensity).

Missing values:
  m1_hgb       NaN passed to HistGradientBoosting (native missing-value routing); no imputation.
  stack_lr_all the core/wr columns use stack_lr_wr's design (unavailable -> 0 next to the on-array / reciprocal /
               SMF availability flags). External pair columns (after a fixed monotone transform: log mutual rank,
               log1p of AP-MS edge/shared counts and of distance) are z-scored with the TRAINING mean/SD over available
               entries and only then set to 0 where missing (= imputing the training mean), next to one availability
               indicator per mask group (coex, apms_bait, apms_obs, compara, esm, coords, uniprot, compara_profile,
               abundance; gen_dist_kb's same-chromosome mask is carried by gen_same_chr). This is the one imputation.
HGB: no class weights; no negative subsampling (optional `subsample` = uniform fraction of ALL training rows, recorded);
early stopping off by default, or (`es_val_frac` > 0) on a held-out fraction of training QUERIES passed as X_val
(never validation_fraction). Optional monotonic constraints on the GI-derived scores (`monotone`).
"""
from __future__ import annotations

import warnings

import numpy as np

from ysl.features import pairs as PF
from ysl.models import baselines as B
from ysl.models import completion as CP
from ysl.models import foldin2 as FD
from ysl.models import propensity as PR
from ysl.models.registry import register
from ysl.models.view import MethodContext, ScreenView

KNN = B._knn_needed(B.LR_FEATURES["stack_lr"])
STACK_BASE = B.LR_FEATURES["stack_lr"]
WR = tuple(FD.WR2_STACK_NAMES)       # = stack_lr_wr2's within-row features (accepted reference)
WR1 = tuple(CP.WR_STACK_NAMES)      # stack_lr_wr's set

# external pair features: name -> (group, LR transform)
PAIRS = {
    "coex_z": ("coex", None), "coex_mr": ("coex", "log"),
    "apms_edge": ("apms", "log1p"), "apms_shared": ("apms", "log1p"), "apms_jaccard": ("apms", None),
    "para_any": ("para", None), "para_within": ("para", None), "para_identity": ("para", None),
    "esm_cos": ("para", None),
    "gen_same_chr": ("genome", None), "gen_dist_kb": ("genome", "log1p"),
    "phylo_ident_corr": ("rest", None), "dom_ipr_jaccard": ("rest", None), "abund_sum": ("rest", None),
    "abund_absdiff": ("rest", None),
}
GROUPS = ("core", "wr", "prop", "coex", "apms", "para", "genome", "rest")
NO_FLAG_MASKS = {"coords_same_chr"}

# HGB columns (raw, NaN = unavailable) of the core group
HGB_CORE = ("log_array_rate", "array_smf", "log_query_prop", "on_array", "recip_negeps", "recip_label",
            "recip_negsp", "log_knn_qs_col", "log_knn_qs_mixed", "log_knn_as_col")
# +1: larger = more interaction-like (used only when params["monotone"]); SMF: lower fitness -> more GIs
MONOTONE = {"log_array_rate": 1, "log_query_prop": 1, "recip_negeps": 1, "recip_label": 1, "recip_negsp": 1,
            "log_knn_qs_col": 1, "log_knn_qs_mixed": 1, "log_knn_as_col": 1, "array_smf": -1, "foldin": 1,
            "mk": 1, "ext_off": 1, "knn_lab": 1, "prop_logit": 1}

COMBINER_DEFAULTS = B._LR_DEFAULTS | B.KNN_DEFAULTS | PR.PROP_DEFAULTS | {
    "prop_fset": "all", "wr": list(WR), "drop": [], "test_view": "matched"}
# HGB recipe selected on dev seeds 0-2 only (scripts/43_tune.py stages 2-4, results/m4b_combiner/tune_hgb.csv), then
# frozen: a-priori default (lr 0.05, 31 leaves, min leaf 1000, l2 1) peaked at <= 50 iterations; lr 0.02 + 50% feature
# subsampling per split (max_features) was best on every seed; 150 iterations (100-150 equal on the wr2 feature set).
HGB_DEFAULTS = COMBINER_DEFAULTS | {
    "additive": "kfold", "ensemble": "pooled", "learning_rate": 0.02, "max_iter": 150, "max_leaf_nodes": 31,
    "min_samples_leaf": 1000, "l2_regularization": 1.0, "max_bins": 255, "max_features": 0.5, "monotone": False,
    "subsample": 1.0, "es_val_frac": 0.0, "n_iter_no_change": 20}


# ----------------------------------------------------------------------------- one feature function
def groups_of(params: dict) -> tuple[str, ...]:
    drop = set(params.get("drop") or [])
    bad = drop - set(GROUPS)
    if bad or "core" in drop:
        raise ValueError(f"cannot drop {sorted(bad | ({'core'} & drop))}")
    return tuple(g for g in GROUPS if g not in drop)


def pair_names(groups) -> tuple[str, ...]:
    return tuple(n for n, (g, _) in PAIRS.items() if g in groups)


def view_features(v: ScreenView, ctx: MethodContext, groups, prop: np.ndarray | None = None) -> dict:
    """Per-target / per-cell features of v's targets (reads only v and external feature files)."""
    f = B.pair_features(v, ctx.params, KNN)
    f["n_t"], f["n_a"] = v.n_targets, v.n_arrays
    if "wr" in groups:
        seed = int(ctx.rng("tied").integers(2**31))
        names = tuple(ctx.get("wr", WR))
        fw, _ = CP._standardise_wr(FD.wr2_features(v, ctx.params, names, seed), None)
        f["wr"] = fw
        f["wr_names"] = names
    if "prop" in groups:
        if prop is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                prop = PR.propensity_predict(v, ctx)
        f["prop"] = np.asarray(prop, dtype=np.float64)
    f["tidx"] = PF.index_of(v.target_orfs_current)
    f["aidx"] = PF.index_of(v.array_orfs_current)
    return f


def _pair_block(f: dict, name: str, rows, transform: bool) -> np.ndarray:
    M = PF.block(name, f["tidx"][rows], f["aidx"])
    if transform:
        tr = PAIRS[name][1]
        if tr == "log":
            M = np.log(M)
        elif tr == "log1p":
            M = np.log1p(np.maximum(M, 0))
    return M


def _cell(x: np.ndarray, n_a: int) -> np.ndarray:
    """per-target vector -> [n_rows x n_a]"""
    return np.broadcast_to(np.asarray(x, dtype=np.float32)[:, None], (len(x), n_a))


def hgb_columns(groups, wr_names=WR) -> list[str]:
    cols = list(HGB_CORE)
    if "wr" in groups:
        cols += list(wr_names)
    if "prop" in groups:
        cols.append("prop_logit")
    return cols + list(pair_names(groups))


def hgb_design(f: dict, groups, rows=slice(None)) -> np.ndarray:
    """HGB inputs [n_rows*n_a x d], NaN = unavailable (no imputation)."""
    n_a = f["n_a"]
    sub = lambda k: np.asarray(f[k])[rows]  # noqa: E731
    on = sub("cols") >= 0
    n_r = len(on)
    out = []
    with np.errstate(divide="ignore", invalid="ignore"):
        out.append(np.log(sub("array_rate")))
        out.append(sub("array_smf"))
        out.append(np.log(sub("query_prop")))                          # NaN off-array
        out.append(_cell(on.astype(np.float32), n_a))
        out.append(-np.clip(sub("recip_eps"), *B.EPS_CLIP))            # NaN where unavailable
        out.append(sub("recip_label"))
        out.append(-np.clip(sub("recip_sp"), -B.SP_CLIP, B.SP_CLIP))
        for k in KNN:
            out.append(np.log(sub(k)))                                 # NaN off-array
    if "wr" in groups:
        for nm in f["wr_names"]:
            out.append(np.asarray(f["wr"][nm])[rows])
    if "prop" in groups:
        out.append(_cell(np.asarray(f["prop"])[rows], n_a))
    for nm in pair_names(groups):
        out.append(_pair_block(f, nm, rows, transform=False))
    X = np.empty((n_r * n_a, len(out)), dtype=np.float32)
    for j, c in enumerate(out):
        X[:, j] = np.broadcast_to(np.asarray(c, dtype=np.float32), (n_r, n_a)).ravel()
    return X


def lr_design_raw(f: dict, groups, rows=slice(None)) -> np.ndarray:
    """LR inputs before the pair-column standardisation: stack_lr_wr columns (+ prop, prop x off-array) followed by
    the transformed external pair columns (NaN where unavailable)."""
    n_a = f["n_a"]
    parts = [B.design_matrix(f, STACK_BASE, rows)]
    if "wr" in groups:
        parts.append(CP._wr_design(f["wr"], f["wr_names"], rows))
    if "prop" in groups:
        parts.append(PR._prop_design(f["prop"], f["cols"], n_a, rows))
    pn = pair_names(groups)
    if pn:
        n_r = len(np.asarray(f["cols"])[rows])
        P = np.empty((n_r * n_a, len(pn)), dtype=np.float32)
        for j, nm in enumerate(pn):
            P[:, j] = _pair_block(f, nm, rows, transform=True).ravel()
        parts.append(P)
    return np.hstack(parts).astype(np.float32, copy=False)


def lr_columns(groups, wr_names=WR) -> list[str]:
    cols = list(STACK_BASE) + (list(wr_names) if "wr" in groups else []) + \
        (list(PR.PROP_COLS) if "prop" in groups else [])
    return cols + list(pair_names(groups))


class PairStandardiser:
    """Pair columns (the last n_pair columns): z-score with the TRAINING mean/SD over available entries, THEN missing
    -> 0, plus one availability indicator per mask group. Other columns pass through unchanged."""

    def __init__(self, pair_cols: tuple[str, ...]):
        self.pair_cols = tuple(pair_cols)
        masks = []
        for nm in self.pair_cols:
            m = PF.FEATURES[nm].mask
            if m not in NO_FLAG_MASKS and m not in [x[0] for x in masks]:
                masks.append((m, nm))
        self.flag_from = masks          # (mask name, representative column)

    def fit(self, X: np.ndarray):
        p = len(self.pair_cols)
        P = X[:, X.shape[1] - p:].astype(np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.mean_ = np.nan_to_num(np.nanmean(P, 0)) if p else np.zeros(0)
            sd = np.nanstd(P, 0) if p else np.zeros(0)
        self.sd_ = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        p = len(self.pair_cols)
        if p == 0:
            return X
        n0 = X.shape[1] - p
        P = X[:, n0:]
        Z = (P - self.mean_[None, :].astype(np.float32)) / self.sd_[None, :].astype(np.float32)
        miss = np.isnan(Z)
        Z[miss] = 0.0
        F = np.stack([~np.isnan(P[:, self.pair_cols.index(nm)]) for _, nm in self.flag_from], axis=1) \
            if self.flag_from else np.zeros((len(X), 0), bool)
        return np.hstack([X[:, :n0], Z, F.astype(np.float32)]).astype(np.float32, copy=False)

    def flag_names(self) -> list[str]:
        return [f"avail_{m}" for m, _ in self.flag_from]


# ----------------------------------------------------------------------------- training set (cross-fitted)
def _fold_rows_features(v: ScreenView, ctx: MethodContext, k: int, folds: np.ndarray, groups, loo: dict | None):
    rows = np.flatnonzero(folds == k)
    fv = v.without_queries(v.train_orfs[rows])
    f = view_features(fv, ctx, groups)
    if loo is not None:
        f = B._loo_override(f, loo, rows)
    return rows, f


def crossfit_training_set(v: ScreenView, ctx: MethodContext, groups, kind: str):
    """(X, y, fold id, training-row id) over all observed training cells; X from `design` of the fold views."""
    folds = ctx.folds(v)
    additive = ctx.get("additive", "loo")
    if kind == "hgb" and additive != "kfold":
        raise ValueError("trees need K-fold additive statistics (exact LOO leaks the label)")
    loo = B.loo_additive(v, ctx.get("rate_prior", B.DEFAULT_RATE_PRIOR)) if additive == "loo" else None
    Xs, ys, fs, qs = [], [], [], []
    for k in range(ctx.n_folds):
        rows, f = _fold_rows_features(v, ctx, k, folds, groups, loo)
        X = hgb_design(f, groups) if kind == "hgb" else lr_design_raw(f, groups)
        y = v.Y[rows].ravel()
        keep = ~np.isnan(y)
        Xs.append(X[keep])
        ys.append(y[keep])
        fs.append(np.full(int(keep.sum()), k, dtype=np.int16))
        qs.append(np.repeat(rows, v.n_arrays)[keep].astype(np.int32))
        del X, f
    return np.concatenate(Xs), np.concatenate(ys).astype(np.int8), np.concatenate(fs), np.concatenate(qs)


def test_feature_views(v: ScreenView, ctx: MethodContext):
    """[(k or None, view)]: the full view, or the K matched views T \\ fold k with v's targets."""
    mode = ctx.get("test_view", "matched")
    if mode == "full":
        return [(None, v)]
    if mode != "matched":
        raise ValueError(mode)
    folds = ctx.folds(v)
    return [(k, PR.matched_target_view(v, v.train_orfs[folds == k])) for k in range(ctx.n_folds)]


def iter_test_designs(v: ScreenView, ctx: MethodContext, groups, kind: str, prop_full: np.ndarray | None):
    """Yield (k, row slice, raw design) chunks of the test grid for every test feature view."""
    for k, tv in test_feature_views(v, ctx):
        f = view_features(tv, ctx, groups, prop=prop_full)
        for i in range(0, v.n_targets, ctx.chunk_targets):
            sl = slice(i, min(i + ctx.chunk_targets, v.n_targets))
            yield k, sl, (hgb_design(f, groups, sl) if kind == "hgb" else lr_design_raw(f, groups, sl))
        del f


def _full_prop(v: ScreenView, ctx: MethodContext, groups):
    if "prop" not in groups:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return PR.propensity_predict(v, ctx)


# ----------------------------------------------------------------------------- drift diagnostic
def drift_check(Xtr: np.ndarray, qtr: np.ndarray, Xte: np.ndarray, qte: np.ndarray, names, seed: int,
                n_iter: int = 100) -> dict:
    """Adversarial AUROC (train rows vs test rows; HGB, 2-fold CV with folds split by QUERY so query identity cannot
    be memorised) + standardised mean differences per feature (NaN-aware) + missing-rate differences."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    X = np.vstack([Xtr, Xte]).astype(np.float32)
    lab = np.r_[np.zeros(len(Xtr)), np.ones(len(Xte))].astype(np.int8)
    grp = np.r_[qtr, qte + 10_000_000]
    ug = np.unique(grp)
    half = set(rng.permutation(ug)[: len(ug) // 2].tolist())
    fold = np.array([g in half for g in grp])
    pred = np.empty(len(X))
    for f in (True, False):
        tr, te = fold != f, fold == f
        m = HistGradientBoostingClassifier(max_iter=n_iter, learning_rate=0.1, max_leaf_nodes=31,
                                           min_samples_leaf=200, early_stopping=False, random_state=0)
        m.fit(X[tr], lab[tr])
        pred[te] = m.decision_function(X[te])
    auc = float(roc_auc_score(lab, pred))
    smd, miss = {}, {}
    for j, nm in enumerate(names):
        a, b = Xtr[:, j].astype(np.float64), Xte[:, j].astype(np.float64)
        a_, b_ = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a_) > 1 and len(b_) > 1:
            s = np.sqrt((a_.var() + b_.var()) / 2)
            smd[nm] = float((b_.mean() - a_.mean()) / s) if s > 0 else 0.0
        miss[nm] = float(np.isnan(b).mean() - np.isnan(a).mean())
    worst = sorted(smd, key=lambda k: -abs(smd[k]))[:8]
    return {"adversarial_auroc": round(auc, 4), "max_abs_smd": round(max(abs(x) for x in smd.values()), 4),
            "worst_smd": {k: round(smd[k], 3) for k in worst},
            "max_abs_missing_rate_diff": round(max(abs(x) for x in miss.values()), 4),
            "smd": {k: round(x, 4) for k, x in smd.items()}, "n_train_rows": int(len(Xtr)),
            "n_test_rows": int(len(Xte))}


def _sample_rows(n: int, m: int, rng: np.random.Generator) -> np.ndarray:
    return np.sort(rng.choice(n, size=min(n, m), replace=False))


# ----------------------------------------------------------------------------- stack_lr_all
# stack_lr_all test view: 'full' selected on dev seeds 0-2 (full vs matched +0.0005 [+0.0003, +0.0006]; matched-view
# averaging helps the HGB, +0.0013 [+0.0010, +0.0016], but not the LR) -> results/m4b_combiner/SUMMARY.md.
@register("stack_lr_all", defaults=COMBINER_DEFAULTS | {"test_view": "full"}, logit_scale=True,
          description="stack_lr_wr2 + propensity + 15 external pair features (z-scored, 0 + flag); K=10 cross-fitted "
                      "LR; full test view")
def stack_lr_all(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    groups = groups_of(ctx.params)
    X, y, _, qtr = crossfit_training_set(v, ctx, groups, "lr")
    ps = PairStandardiser(pair_names(groups)).fit(X)
    X = ps.transform(X)
    sc, clf = B.fit_lr(X, y, ctx)
    names = lr_columns(groups, tuple(ctx.get("wr", WR))) + ps.flag_names()
    ctx.info["lr_coef"] = dict(zip(names, np.round(clf.coef_[0], 5).tolist()))
    ctx.info["n_training_cells"] = int(len(y))
    rng = ctx.rng("drift")
    do_drift = bool(ctx.get("drift", True))
    if do_drift:
        itr = _sample_rows(len(X), 50_000, rng)
        Xtr_s, qtr_s = X[itr], qtr[itr]
    del X
    S, te_s, qte_s, n_views = _predict_views(v, ctx, groups, "lr", lambda k, D: clf.decision_function(
        sc.transform(ps.transform(D))), ps.transform if do_drift else None, rng)
    if do_drift:
        ctx.info["drift_full_grid_gi"] = _gi_drift(Xtr_s, qtr_s, te_s, qte_s, names, groups, rng)
    ctx.info["test_views"] = n_views
    return S


def _gi_drift(Xtr, qtr, Xte, qte, names, groups, rng) -> dict:
    """In-method drift check on the GI-derived (+ propensity) columns only, over the FULL test grid. External pair
    columns are static functions of the gene pair (no train/test feature drift by construction); over the full grid
    their distribution differs only because linkage-masked cells (close same-chromosome pairs) are never observed in
    training. The observed-cell drift check over all columns is scripts/43_drift.py."""
    pn = set(pair_names(groups))
    keep = [j for j, n in enumerate(names) if n not in pn and not n.startswith("avail_")]
    return drift_check(Xtr[:, keep], qtr, Xte[:, keep], qte, [names[j] for j in keep], int(rng.integers(2**31)))


def _predict_views(v, ctx, groups, kind, score_fn, drift_tf, rng, n_drift: int = 50_000):
    """Average the per-view log-odds over the test feature views; also return a uniform drift sample of test-grid
    rows (~n_drift in total, spread over the views)."""
    n_t, n_a = v.n_targets, v.n_arrays
    S = np.zeros((n_t, n_a), dtype=np.float64)
    seen = {}
    te, qte = [], []
    prop_full = _full_prop(v, ctx, groups)
    n_exp = 1 if ctx.get("test_view", "matched") == "full" else ctx.n_folds
    for k, sl, D in iter_test_designs(v, ctx, groups, kind, prop_full):
        S[sl] += np.asarray(score_fn(k, D), dtype=np.float64).reshape(-1, n_a)
        seen[k] = True
        if drift_tf is not None:
            m = max(1, int(round(n_drift * (sl.stop - sl.start) / n_t / n_exp)))
            idx = _sample_rows(len(D), m, rng)
            te.append(drift_tf(D[idx]))
            qte.append((sl.start + idx // n_a).astype(np.int32))
    n_views = len(seen)
    S /= n_views
    te_s = np.vstack(te) if te else None
    qte_s = np.concatenate(qte) if qte else None
    return S.astype(np.float32), te_s, qte_s, n_views


# ----------------------------------------------------------------------------- m1_hgb
def _hgb(ctx: MethodContext, names, n_iter=None):
    from sklearn.ensemble import HistGradientBoostingClassifier
    mono = [MONOTONE.get(n, 0) for n in names] if ctx.get("monotone", False) else None
    es = float(ctx.get("es_val_frac", 0.0)) > 0
    return HistGradientBoostingClassifier(
        learning_rate=float(ctx.get("learning_rate")), max_iter=int(n_iter or ctx.get("max_iter")),
        max_leaf_nodes=int(ctx.get("max_leaf_nodes")), min_samples_leaf=int(ctx.get("min_samples_leaf")),
        l2_regularization=float(ctx.get("l2_regularization")), max_bins=int(ctx.get("max_bins")),
        max_features=float(ctx.get("max_features", 1.0)), monotonic_cst=mono, early_stopping=es,
        n_iter_no_change=int(ctx.get("n_iter_no_change", 20)), scoring="loss", class_weight=None,
        random_state=int(ctx.rng("hgb").integers(2**31 - 1)))


class FittedHGB:
    """A fitted HGB plus the training columns that were entirely unavailable (all NaN; e.g. external features on toy
    genes): those are set to 0 (constant -> never split on) in training and in every scored design."""

    def __init__(self, clf, dead: np.ndarray):
        self.clf, self.dead = clf, dead
        self.n_iter_ = None

    def _x(self, D: np.ndarray) -> np.ndarray:
        if not self.dead.any():
            return D
        D = np.array(D, copy=True)
        D[:, self.dead] = 0.0
        return D

    def decision_function(self, D: np.ndarray) -> np.ndarray:
        return self.clf.decision_function(self._x(D))


def fit_hgb(X: np.ndarray, y: np.ndarray, q: np.ndarray, ctx: MethodContext, names, tag: str = "") -> FittedHGB:
    """Fit one HGB; optional uniform row subsample (over ALL training rows) and query-grouped early-stopping split
    (X_val = the rows of a random fraction of training queries)."""
    frac = float(ctx.get("subsample", 1.0))
    keep = ctx.rng("subsample" + tag).random(len(y)) < frac if frac < 1.0 else None
    dead = np.isnan(X).all(0)
    model = FittedHGB(_hgb(ctx, names), dead)
    if keep is not None:
        X, y, q = X[keep], y[keep], q[keep]
    X = model._x(X)
    es = float(ctx.get("es_val_frac", 0.0))
    if es > 0:
        uq = np.unique(q)
        vq = ctx.rng("es_val" + tag).random(len(uq)) < es
        is_val = np.isin(q, uq[vq])
        model.clf.fit(X[~is_val], y[~is_val], X_val=X[is_val], y_val=y[is_val])
    else:
        model.clf.fit(X, y)
    model.n_iter_ = model.clf.n_iter_
    return model


@register("m1_hgb", defaults=dict(HGB_DEFAULTS), logit_scale=True,
          description="M1: HistGradientBoosting on core + within-row + propensity + external pair features (NaN "
                      "native); K-fold additive stats; pooled model, matched test views")
def m1_hgb(v: ScreenView, ctx: MethodContext) -> np.ndarray:
    groups = groups_of(ctx.params)
    names = hgb_columns(groups, tuple(ctx.get("wr", WR)))
    X, y, fold, qtr = crossfit_training_set(v, ctx, groups, "hgb")
    ens = ctx.get("ensemble", "pooled")
    if ens == "pooled":
        models = {None: fit_hgb(X, y, qtr, ctx, names)}
        score = lambda k, D: models[None].decision_function(D)  # noqa: E731
    elif ens == "fold":
        if ctx.get("test_view", "matched") != "matched":
            raise ValueError("the fold (matched-view) ensemble needs test_view='matched'")
        models = {k: fit_hgb(X[fold == k], y[fold == k], qtr[fold == k], ctx, names, tag=str(k))
                  for k in range(ctx.n_folds)}
        score = lambda k, D: models[k].decision_function(D)  # noqa: E731
    else:
        raise ValueError(ens)
    ctx.info["n_training_cells"] = int(len(y))
    ctx.info["subsample"] = float(ctx.get("subsample", 1.0))
    ctx.info["n_iter"] = {str(k): int(m.n_iter_) for k, m in models.items()}
    ctx.info["features"] = names
    rng = ctx.rng("drift")
    do_drift = bool(ctx.get("drift", True))
    if do_drift:
        itr = _sample_rows(len(X), 50_000, rng)
        Xtr_s, qtr_s = X[itr], qtr[itr]
    del X
    S, te_s, qte_s, n_views = _predict_views(v, ctx, groups, "hgb", score, (lambda D: D) if do_drift else None, rng)
    if do_drift:
        ctx.info["drift_full_grid_gi"] = _gi_drift(Xtr_s, qtr_s, te_s, qte_s, names, groups, rng)
    ctx.info["test_views"] = n_views
    return S
