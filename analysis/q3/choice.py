"""Q3 analysis 3 — choice of query lines (DESIGN v3.1 §7.4): what do the adaptive arms pick early, on
human-available features, and how much of their gain does a static, interpretable rank rule recover?

    .venv/bin/python -m analysis.q3.choice profile --source dev|final|dryrun     (no model fits)
    .venv/bin/python -m analysis.q3.choice evaluate --source dev                 (P* fits, 1 job, 2 threads)
    .venv/bin/python -m analysis.q3.choice report --source dev

profile   For each arm, the first `n_early` (160) picks of every stored trajectory vs its candidate pool:
          standardised mean difference SMD = (mean over picks − pool mean) / pool SD, per feature, averaged over
          trajectories; 95% CI by resampling trajectories within split. Random (R0) is the null control. The
          true query SL degree (dev labels, evaluation only) is shown for context, never used by a rule.
rule      Target: per split, the fraction of finalist trajectories (H, C, F_r16) that picked each pool gene in
          their first n_early. Features: the curated human-available set below (standardised within pool; missing
          -> pool median, so off-array genes are route-neutral, DESIGN §1). Models: forward-selected logistic
          score with 1, 2, 3 features (binomial deviance, fractional response), and an L1 logistic on all
          features. Leave-split-out: the rule applied to split k is fitted on the other splits' pools minus split
          k's D_val genes (their pick labels never enter).
evaluate  Static selection of n in {160, 330, 660} on each dev split from the leave-split-out rule: randomised
          top-k picks (protocol selector.randomised_greedy_top_k) with seeds 0..S-1, then P* (analysis.q3.fit).
          Recovered gain = (rule − random) / (strategy − random) at the same n (random = P1 R0 mean, strategy =
          P4 mean of H, and of the three finalists).
Final mode: profile only (rule fitting on final selections is descriptive; no T information is used: features
and selections only; final evaluation of the rule would need T labels and is out of scope).
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from analysis.q3 import common as C

# human-available features (DESIGN §1 analogues): name -> (gene_features column | 'smf', label, human analogue)
FEATURES = {
    "smf": ("smf", "single-mutant fitness (array)", "HAP1 single-KO fitness"),
    "abundance": ("log10_ppm_mueller2020", "protein abundance (PaxDb)", "PaxDb human"),
    "expr_mean": ("expr_mean_rnaseq", "mRNA expression level", "HAP1/GTEx expression"),
    "expr_sd": ("expr_sd_rnaseq", "mRNA expression variability", "ARCHS4 variability"),
    "coexp_hub": ("coexp_top20_mean", "coexpression hub-ness (top-20 mean)", "COXPRESdb/ARCHS4"),
    "apms_degree": ("log1p_ppi_degree", "AP-MS degree (log1p)", "BioPlex/OpenCell degree"),
    "apms_in_network": ("ppi_in_network", "in AP-MS network", "BioPlex coverage"),
    "conservation": ("phylo_breadth", "phylogenetic breadth", "ortholog breadth"),
    "identity": ("phylo_mean_identity", "mean ortholog identity", "ortholog identity"),
    "paralogs": ("n_paralogs_wsp", "paralog count", "Ensembl Compara"),
    "paralog_identity": ("max_paralog_identity", "max paralog identity", "Ensembl Compara"),
    "length": ("log_protein_length", "protein length (log)", "same"),
    "plddt": ("plddt_mean", "AlphaFold mean pLDDT", "AlphaFold human"),
    "disorder": ("fraction_disordered", "fraction disordered", "AlphaFold human"),
    "domains": ("n_interpro", "InterPro domain count", "InterPro human"),
    "dnds": ("dnds_yn00_mean", "dN/dS (population)", "gnomAD constraint"),
    "lof_freq": ("lof_freq", "LoF variant frequency", "gnomAD LOEUF"),
}
ARMS = ("R0", "H", "C", "F_r16", "H-static", "R1")
N_EARLY = 160


def feature_table(genes, smf_map: dict) -> pd.DataFrame:
    from hsl import features as HF
    g = np.asarray(sorted(set(map(str, genes))), dtype=str)
    X = HF.gene_features(g)
    X.index = g
    out = pd.DataFrame(index=g)
    for k, (col, _, _) in FEATURES.items():
        out[k] = [smf_map.get(x, np.nan) for x in g] if col == "smf" else X[col].to_numpy(float)
    return out


def smf_map() -> dict:
    from analysis.q3 import cache
    d = cache.load()
    return {o: float(v) for o, v in zip(d["array_orfs"], d["array_smf"]) if np.isfinite(v)}


def degree_map() -> dict:
    """True strong-SL count of each dev query row (evaluation-only context)."""
    from analysis.q3 import cache
    d = cache.load()
    keep = ~np.isin(d["array_orfs"], d["test"])
    k = (d["Y"][:, keep] == 1).sum(1)
    return dict(zip(d["query_orfs"], k.astype(float)))


def pools(source: str) -> dict:
    from analysis.q3 import cache
    if source == "dev":
        return {sp: cache.eset(sp)["pool"] for sp in ("dev1", "dev2", "dev3")}
    from hsl import views
    if source == "dryrun":
        return {"dev1": cache.eset("dev1")["pool"]}
    return {"final": np.asarray(views.final_split().pool, dtype=str)}


def standardise(F: pd.DataFrame, pool) -> pd.DataFrame:
    P = F.loc[pool]
    mu, sd = P.mean(), P.std(ddof=0).replace(0, 1)
    Z = (F - mu) / sd
    return Z.fillna(0.0)          # pool median ~ 0 after standardising; route-neutral imputation


# ============================================================================ profile
def profile(source: str, out=None) -> pd.DataFrame:
    t0 = time.perf_counter()
    O = C.out_dir(source, out)
    sels = C.selections(source, arms=ARMS)
    P = pools(source)
    allg = np.unique(np.concatenate(list(P.values())))
    F = feature_table(allg, smf_map())
    deg = degree_map() if source != "final" else {}
    if deg:
        F["true_sl_degree"] = [deg.get(g, np.nan) for g in F.index]
    recs = []
    for (arm, sp, r), s in sels.items():
        if sp not in P:
            continue
        pool = P[sp]
        picks = [g for g in s["selected"][:N_EARLY]]
        Fp = F.loc[pool]
        mu, sd = Fp.mean(), Fp.std(ddof=1)
        m = F.loc[picks].mean()
        rec = {"arm": arm, "split": sp, "traj": r}
        for k in F.columns:
            rec[k] = float((m[k] - mu[k]) / sd[k]) if sd[k] > 0 else np.nan
        rec["offarray_share"] = float(np.mean([g not in smf_map_cached() for g in picks]))
        recs.append(rec)
    T = pd.DataFrame(recs)
    T.to_csv(O / "choice_profile_traj.csv", index=False)
    rng = np.random.default_rng(C.hap1_constants()["bootstrap"]["seed"])
    rows = []
    feats = [c for c in T.columns if c not in ("arm", "split", "traj")]
    for arm, g in T.groupby("arm"):
        bysplit = [x for _, x in g.groupby("split")]
        for f in feats:
            pt = float(np.mean([x[f].mean() for x in bysplit]))
            reps = [np.mean([x[f].to_numpy()[rng.integers(0, len(x), len(x))].mean() for x in bysplit])
                    for _ in range(1000)]
            lo, hi = np.nanpercentile(reps, [2.5, 97.5])
            rows.append({"arm": arm, "feature": f, "smd": pt, "lo": lo, "hi": hi, "R": len(g),
                         "label": FEATURES.get(f, (None, f, None))[1], "human_analogue": FEATURES.get(f, (0, 0, ""))[2]})
    prof = pd.DataFrame(rows)
    prof.to_csv(O / "choice_profile.csv", index=False)
    fig_profile(prof, O / "fig_choice_profile.png", source)
    C.write_json(O / "choice_profile_meta.json", {"source": source, "n_early": N_EARLY,
                                                   "runtime_s": time.perf_counter() - t0})
    return prof


_SMF = {}


def smf_map_cached():
    if not _SMF:
        _SMF.update(smf_map())
    return _SMF


# ============================================================================ rule
def pick_freq(sels, pool, sp, arms=C.FINALISTS, n_early=N_EARLY) -> pd.Series:
    cnt = pd.Series(0.0, index=pool)
    R = 0
    for (arm, s2, r), s in sels.items():
        if s2 != sp or arm not in arms:
            continue
        R += 1
        cnt.loc[[g for g in s["selected"][:n_early] if g in cnt.index]] += 1
    return cnt / max(R, 1), R


def fit_logit(Z: np.ndarray, y: np.ndarray, w: np.ndarray, C_reg: float = 1e4, l1: bool = False):
    """Fractional-response logistic via duplicated rows (y and 1−y with weights w·y, w·(1−y))."""
    from sklearn.linear_model import LogisticRegression
    X2 = np.vstack([Z, Z])
    y2 = np.r_[np.ones(len(y)), np.zeros(len(y))]
    w2 = np.r_[w * y, w * (1 - y)]
    m = LogisticRegression(C=C_reg, l1_ratio=1.0 if l1 else 0.0, solver="liblinear" if l1 else "lbfgs",
                           max_iter=2000)
    m.fit(X2, y2, sample_weight=w2 + 1e-12)
    return m


def deviance(m, Z, y, w):
    p = np.clip(m.predict_proba(Z)[:, 1], 1e-9, 1 - 1e-9)
    return float(-2 * np.sum(w * (y * np.log(p) + (1 - y) * np.log(1 - p))))


def forward_select(Z: pd.DataFrame, y, w, k_max=3):
    chosen, path = [], []
    for _ in range(k_max):
        best = None
        for f in Z.columns:
            if f in chosen:
                continue
            cols = chosen + [f]
            m = fit_logit(Z[cols].to_numpy(), y, w)
            d = deviance(m, Z[cols].to_numpy(), y, w)
            if best is None or d < best[1]:
                best = (f, d, m)
        chosen.append(best[0])
        path.append({"k": len(chosen), "features": list(chosen), "deviance": best[1],
                     "coef": dict(zip(chosen, best[2].coef_[0].round(4).tolist()))})
    return path


def train_rule(source: str, eval_split: str | None, F: pd.DataFrame, sels, P, k_max=3):
    """Rule fitted on all splits != eval_split, rows = their pools minus eval_split's D_val (dev)."""
    from analysis.q3 import cache
    excl = set(cache.eset(eval_split)["targets"]) if eval_split else set()
    Zs, ys, ws = [], [], []
    for sp, pool in P.items():
        if sp == eval_split:
            continue
        pool = np.array([g for g in pool if g not in excl])
        y, R = pick_freq(sels, pool, sp)
        Z = standardise(F[list(FEATURES)], pool).loc[pool]
        Zs.append(Z)
        ys.append(y.loc[pool].to_numpy())
        ws.append(np.full(len(pool), float(R)))
    Z = pd.concat(Zs)
    y, w = np.concatenate(ys), np.concatenate(ws)
    path = forward_select(Z, y, w, k_max)
    l1 = fit_logit(Z.to_numpy(), y, w, C_reg=0.05, l1=True)
    return path, {"features": list(Z.columns), "coef": dict(zip(Z.columns, l1.coef_[0].round(4).tolist()))}


def rule_scores(F, pool, coef: dict, seed: int = 0) -> np.ndarray:
    """Linear rule score. Route-stratified for genes without SMF (off-array in yeast; every HAP1 gene has a
    fitness value): they are ranked among themselves by the rule minus its SMF term (ties broken at random) and
    their percentile is mapped onto the on-array score distribution (as the strategies' D5 calibration), so a
    missing SMF is never a score term in either direction (DESIGN §1)."""
    Z = standardise(F[list(FEATURES)], pool).loc[pool]
    if not coef:
        return np.zeros(len(pool))
    w = np.array(list(coef.values()))
    s = Z[list(coef)].to_numpy() @ w
    off = ~np.isfinite(F.loc[pool, "smf"].to_numpy(float))
    if "smf" in coef and off.any() and (~off).any():
        rest = {k: v for k, v in coef.items() if k != "smf"}
        s_off = (Z.loc[off, list(rest)].to_numpy() @ np.array(list(rest.values()))) if rest else np.zeros(off.sum())
        rng = np.random.default_rng(seed)
        u = (pd.Series(s_off + 1e-9 * rng.standard_normal(off.sum())).rank(method="first").to_numpy() - 0.5) / off.sum()
        s[off] = np.quantile(s[~off], u)
    return s


def rules(source: str, out=None) -> dict:
    """Leave-split-out rules + their held-out ranking AUROC against the finalists' pick frequency."""
    from sklearn.metrics import roc_auc_score
    from scipy.stats import spearmanr
    O = C.out_dir(source, out)
    sels = C.selections(source, arms=ARMS)
    P = pools(source)
    allg = np.unique(np.concatenate(list(P.values())))
    F = feature_table(allg, smf_map())
    res, rows = {}, []
    deg = degree_map() if source != "final" else {}
    splits = list(P) if source == "dev" else [None]
    for ev in splits:
        path, l1 = train_rule(source, ev, F, sels, P)
        res[str(ev)] = {"forward": path, "l1": l1}
        if ev is None:
            continue
        pool = P[ev]
        y, _ = pick_freq(sels, pool, ev)
        for k, p in enumerate(path, 1):
            sc = pd.Series(rule_scores(F, pool, p["coef"]), index=pool)
            d = np.array([deg.get(g, np.nan) for g in pool])
            rows.append({"eval_split": ev, "rule": f"forward_k{k}", "features": "+".join(p["features"]),
                         "coef": json.dumps(p["coef"]),
                         "auroc_vs_picked_any": roc_auc_score(y.to_numpy() > 0, sc.to_numpy()),
                         "spearman_vs_pickfreq": spearmanr(sc, y).statistic,
                         "spearman_vs_true_degree": spearmanr(sc, d, nan_policy="omit").statistic})
        sc = pd.Series(rule_scores(F, pool, {k: v for k, v in l1["coef"].items() if v != 0}), index=pool)
        rows.append({"eval_split": ev, "rule": "l1_all", "features": "+".join(k for k, v in l1["coef"].items() if v),
                     "coef": json.dumps(l1["coef"]),
                     "auroc_vs_picked_any": roc_auc_score(y.to_numpy() > 0, sc.to_numpy()),
                     "spearman_vs_pickfreq": spearmanr(sc, y).statistic,
                     "spearman_vs_true_degree": spearmanr(sc, [deg.get(g, np.nan) for g in pool],
                                                          nan_policy="omit").statistic})
    if source == "dev":     # rule on all three splits (the one to report / apply to HAP1)
        path, l1 = train_rule(source, None, F, sels, P)
        res["all"] = {"forward": path, "l1": l1}
    pd.DataFrame(rows).to_csv(O / "choice_rules_heldout.csv", index=False)
    C.write_json(O / "choice_rules.json", res)
    return res


# ============================================================================ evaluate (P* fits)
def static_pick(scores: pd.Series, n: int, seed: int, top_k: int) -> list:
    """Randomised top-k static selection (each pick uniform among the top-k remaining by score); seed -1 =
    deterministic top-n."""
    order = scores.sort_values(ascending=False, kind="mergesort").index.tolist()
    if seed < 0:
        return order[:n]
    rng = np.random.default_rng(np.random.SeedSequence([int(C.hap1_constants()["bootstrap"]["seed"]), int(seed)]))
    rem, out = list(order), []
    while len(out) < n:
        j = int(rng.integers(0, min(top_k, len(rem))))
        out.append(rem.pop(j))
    return out


def evaluate(source: str, out=None, budgets=(160, 330, 660), seeds=(0, 1, 2), rule="forward_k2") -> None:
    if source != "dev":
        raise PermissionError("rule evaluation needs dev labels: dev only")
    from analysis.q3 import fit
    O = C.out_dir(source, out)
    part = O / "choice_eval_rows.jsonl"
    done = set()
    if part.exists():
        done = {(x["split"], x["rule"], x["seed"], x["n"]) for x in map(json.loads, part.read_text().splitlines())}
    res = json.loads((O / "choice_rules.json").read_text())
    sels = C.selections(source, arms=ARMS)
    P = pools(source)
    F = feature_table(np.unique(np.concatenate(list(P.values()))), smf_map())
    top_k = int(C.protocol()["selector"]["randomised_greedy_top_k"])
    b = C.dev_bundle()
    rand = C.rows("dev")
    for sp in P:
        k = int(rule.split("_k")[1]) if "_k" in rule else None
        coef = res[sp]["forward"][k - 1]["coef"] if k else {a: v for a, v in res[sp]["l1"]["coef"].items() if v}
        for seed in seeds:
            sc = pd.Series(rule_scores(F, P[sp], coef, seed=seed), index=P[sp])
            sel = static_pick(sc, 660, seed, top_k)        # nested prefixes: same order for every budget
            ms = int(C.rows("dev").query("arm == 'R0' and split == @sp and traj == @seed").model_seed.iloc[0])
            for n in budgets:
                if (sp, rule, seed, n) in done:
                    continue
                r, _ = fit.fit_score(b, sp, sel[:n], ms)
                r.update({"rule": rule, "seed": seed, "coef": coef,
                          "offarray_share": float(np.mean([g not in smf_map_cached() for g in sel[:n]]))})
                with open(part, "a") as fh:
                    fh.write(json.dumps(r, default=float) + "\n")
                print(f"{sp} {rule} seed{seed} n{n}: AUPRC {r['auprc']:.4f} ({r['fit_s']:.0f} s)", flush=True)
    del rand


def report(source: str, out=None) -> pd.DataFrame:
    O = C.out_dir(source, out)
    ev = pd.DataFrame([json.loads(x) for x in (O / "choice_eval_rows.jsonl").read_text().splitlines() if x.strip()])
    rows = C.rows("dev")
    out_rows = []
    for (rule, n), g in ev.groupby(["rule", "n"]):
        rec = {"rule": rule, "n": n}
        for col in ("auprc", "auprc_VV"):
            per = []
            for sp, gs in g.groupby("split"):
                rr = rows[(rows.arm == "R0") & (rows.split == sp) & (rows.n == n)][col].mean()
                hh = rows[(rows.arm == "H") & (rows.split == sp) & (rows.n == n)][col].mean()
                ff = rows[(rows.arm.isin(C.FINALISTS)) & (rows.split == sp) & (rows.n == n)][col].mean()
                hs = rows[(rows.arm == "H-static") & (rows.split == sp) & (rows.n == n)][col].mean()
                per.append({"rule": gs[col].mean(), "random": rr, "H": hh, "fin": ff, "Hstatic": hs})
            D = pd.DataFrame(per)
            rec.update({f"{col}_rule": D.rule.mean(), f"{col}_random": D.random.mean(), f"{col}_H": D.H.mean(),
                        f"{col}_finalists": D.fin.mean(), f"{col}_H_static": D.Hstatic.mean(),
                        f"{col}_rel_gain_rule": (D.rule / D.random - 1).mean(),
                        f"{col}_rel_gain_H": (D.H / D.random - 1).mean(),
                        f"{col}_recovered_vs_H": ((D.rule - D.random) / (D.H - D.random)).mean(),
                        f"{col}_recovered_vs_finalists": ((D.rule - D.random) / (D.fin - D.random)).mean(),
                        f"{col}_recovered_H_static_vs_H": ((D.Hstatic - D.random) / (D.H - D.random)).mean()})
        rec["R_seeds"] = int(g.seed.nunique())
        rec["fit_s_mean"] = float(g.fit_s.mean())
        out_rows.append(rec)
    R = pd.DataFrame(out_rows)
    R.to_csv(O / "choice_rule_gain.csv", index=False)
    return R


# ============================================================================ figure
def fig_profile(prof, path, source):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
                         "axes.spines.right": False})
    style = {"R0": "#6f6e69", "H": "#2a78d6", "C": "#eb6834", "F_r16": "#1baf7a"}
    feats = [f for f in list(FEATURES) + ["true_sl_degree"] if f in set(prof.feature)]
    order = (prof[prof.arm == "H"].set_index("feature").smd.reindex(feats).abs().sort_values().index.tolist())
    f, ax = plt.subplots(figsize=(7, 0.32 * len(order) + 1.2), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.axvline(0, color="#9a9992", lw=0.8)
    off = {"R0": -0.27, "H": -0.09, "C": 0.09, "F_r16": 0.27}
    for a, col in style.items():
        d = prof[prof.arm == a].set_index("feature").reindex(order)
        y = np.arange(len(order)) + off[a]
        ax.errorbar(d.smd, y, xerr=[d.smd - d.lo, d.hi - d.smd], fmt="o", ms=4, color=col, elinewidth=1,
                    label={"R0": "random (null)"}.get(a, a))
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([FEATURES.get(o, (0, "true SL degree (eval only)"))[1] for o in order])
    ax.set_xlabel(f"standardised mean difference, first {N_EARLY} picks vs candidate pool")
    ax.grid(True, axis="x", color="#e8e7e2", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title(f"What the adaptive arms pick early ({source})", fontsize=9)
    f.tight_layout()
    f.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(f)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("profile", "rules", "evaluate", "report"))
    ap.add_argument("--source", choices=C.SOURCES, default="dev")
    ap.add_argument("--out", default=None)
    ap.add_argument("--rule", default="forward_k2")
    ap.add_argument("--budgets", default="160,330,660")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--rss-cap-gb", type=float, default=3.0)
    a = ap.parse_args(argv)
    C.set_threads(2)
    if a.cmd == "profile":
        p = profile(a.source, a.out)
        print(p[p.arm.isin(["H", "R0"])].round(3).to_string())
    elif a.cmd == "rules":
        r = rules(a.source, a.out)
        print(json.dumps(r, indent=1, default=str)[:3000])
    elif a.cmd == "evaluate":
        from analysis.q3.lbo import rss_watchdog
        rss_watchdog(a.rss_cap_gb - 0.05)
        evaluate(a.source, a.out, rule=a.rule, budgets=tuple(int(x) for x in a.budgets.split(",")),
                 seeds=tuple(int(x) for x in a.seeds.split(",")))
    else:
        print(report(a.source, a.out).round(4).to_string())


if __name__ == "__main__":
    main()
