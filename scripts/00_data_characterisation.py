#!/usr/bin/env python
"""Aggregate characterisation of the Costanzo 2016 S1 NxN data (read in place from YEAST_SLs via the vendored ysl
bundle loader). Only aggregate statistics; allowed before the test split (no per-gene selection is made here).

Outputs: results/data_characterisation/summary.json and summary.md

Gene level = upstream's definition: one canonical query strain per query ORF (ysl.data.costanzo.canonical_query_strains;
S1 has exactly one strain per ORF), array columns = array ORFs (one DMA strain each). Label = stringent negative
(P < 0.05 and eps < -0.12; ysl.eval.labels.stringent_negative), computed by upstream when the bundle was built.
"""
from __future__ import annotations

import hsl  # noqa: F401  (YEAST_SLs write guard + paths; must be first)

import json
import time
from pathlib import Path

import numpy as np
from scipy import stats

from hsl import paths
from ysl.data.dense import load_bundle

OUT = paths.output_root() / "results" / "data_characterisation"


def qdist(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    qs = [0, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]
    return {"n": int(len(x)), "mean": float(x.mean()), "sd": float(x.std(ddof=1)),
            "quantiles": {f"q{int(q * 100):02d}": float(np.quantile(x, q)) for q in qs},
            "n_zero": int((x == 0).sum()), "frac_zero": float((x == 0).mean())}


def top_share(counts: np.ndarray, fracs=(0.01, 0.05, 0.10, 0.20)) -> dict:
    c = np.sort(np.asarray(counts, dtype=float))[::-1]
    tot = c.sum()
    out = {}
    for f in fracs:
        k = max(1, int(round(f * len(c))))
        out[f"top{int(f * 100)}pct_n{k}"] = float(c[:k].sum() / tot)
    # Gini coefficient
    s = np.sort(np.asarray(counts, dtype=float))
    n = len(s)
    out["gini"] = float((2 * np.arange(1, n + 1) - n - 1).dot(s) / (n * s.sum()))
    return out


def corr(a, b) -> dict:
    return {"n": int(len(a)), "pearson": float(stats.pearsonr(a, b)[0]), "spearman": float(stats.spearmanr(a, b)[0])}


def main() -> None:
    t0 = time.time()
    ds = paths.costanzo_dataset()
    b = load_bundle(ds, unblind=True)   # aggregate statistics over ALL queries (incl. upstream lockbox)
    rows = b.gene_rows(b.all_query_orfs)
    Y = b.Y[rows]                       # int8: 1 / 0 / -1 (unobserved)
    E = b.E[rows]
    obs = Y >= 0
    pos = Y == 1
    qorfs = b.all_query_orfs
    aorfs = b.array_orfs
    qset, aset = set(qorfs), set(aorfs)
    both = sorted(qset & aset)

    n_obs, n_pos = int(obs.sum()), int(pos.sum())
    q_sl = pos.sum(1)
    q_obs = obs.sum(1)
    a_sl = pos.sum(0)
    a_obs = obs.sum(0)

    # query-side vs array-side SL count for genes that are both
    o2c = b.orf_to_col
    qi = {o: i for i, o in enumerate(qorfs)}
    bi = np.array([qi[g] for g in both])
    bc = np.array([o2c[g] for g in both])
    side = {"counts": corr(q_sl[bi], a_sl[bc]),
            "rates": corr(q_sl[bi] / q_obs[bi], a_sl[bc] / np.maximum(a_obs[bc], 1))}

    # reciprocal (AB/BA) pairs among genes that are both query and array
    Yb = Y[np.ix_(bi, bc)]              # [g as query, h as array], g, h in `both` (same order)
    Eb = E[np.ix_(bi, bc)]
    ob = Yb >= 0
    iu = np.triu_indices(len(both), k=1)
    ab_obs, ba_obs = ob[iu], ob.T[iu]
    ab_y, ba_y = Yb[iu], Yb.T[iu]
    m_any = ab_obs | ba_obs
    m_both = ab_obs & ba_obs
    yab, yba = ab_y[m_both] == 1, ba_y[m_both] == 1
    eab, eba = Eb[iu][m_both], Eb.T[iu][m_both]
    n_union_pos = int((yab | yba).sum())
    recip = {
        "n_genes_query_and_array": len(both),
        "n_unordered_pairs_measured_any_direction": int(m_any.sum()),
        "n_unordered_pairs_measured_both_directions": int(m_both.sum()),
        "frac_of_measured_pairs_in_both_directions": float(m_both.sum() / m_any.sum()),
        "frac_of_all_observed_cells_that_are_part_of_a_reciprocal_pair": float(2 * m_both.sum() / n_obs),
        "prevalence_AB": float(yab.mean()), "prevalence_BA": float(yba.mean()),
        "n_pos_both_directions": int((yab & yba).sum()),
        "n_pos_exactly_one_direction": int((yab ^ yba).sum()),
        "P_BA_pos_given_AB_pos": float((yab & yba).sum() / yab.sum()),
        "P_AB_pos_given_BA_pos": float((yab & yba).sum() / yba.sum()),
        "jaccard_pos": float((yab & yba).sum() / n_union_pos),
        "cohen_kappa": float(_kappa(yab, yba)),
        "overall_label_agreement": float((yab == yba).mean()),
        "epsilon_pearson": float(stats.pearsonr(eab, eba)[0]),
        "epsilon_spearman": float(stats.spearmanr(eab, eba)[0]),
        "self_pairs_observed(query ORF == array ORF)": int(np.diag(ob).sum()),
    }

    # strain / allele structure
    stypes = np.unique([s.split("_")[1].rstrip("0123456789") for s in b.strain_ids], return_counts=True)
    atypes = np.unique([s.split("_")[1].rstrip("0123456789") for s in b.array_strains], return_counts=True)
    strains = {
        "n_query_strains": int(len(b.strain_ids)), "n_query_orfs": int(len(qorfs)),
        "query_strains_per_orf_max": int(np.unique(b.strain_orf, return_counts=True)[1].max()),
        "query_strain_types": dict(zip(stypes[0].tolist(), stypes[1].tolist())),
        "query_strains_with_suppressor_allele": int(b.strain_is_supp.sum()),
        "screen_temperatures_C": sorted(set(b.strain_temperature.tolist())),
        "n_array_strains": int(len(b.array_strains)), "n_array_orfs": int(len(set(aorfs))),
        "array_strain_types": dict(zip(atypes[0].tolist(), atypes[1].tolist())),
        "array_strain_type_note": "dma = deletion mutant array (non-essential gene deletions); no TS/DAmP alleles in NxN",
        "array_smf": {"n_missing": int(np.isnan(b.array_smf).sum()), **qdist(b.array_smf[~np.isnan(b.array_smf)])},
    }

    summary = {
        "dataset": ds, "source": "Costanzo et al. 2016 Science Data File S1, SGA_NxN.txt (via upstream bundle)",
        "bundle_meta": {k: b.meta[k] for k in ("raw_nxn_sha256", "processed_parquet_sha256", "bundle_sha256", "shape",
                                               "n_observed", "n_positive")},
        "label": "stringent negative: P < 0.05 and epsilon < -0.12",
        "genes": {
            "n_query_orfs": int(len(qset)), "n_array_orfs": int(len(aset)),
            "n_query_and_array": len(both), "n_array_only": int(len(aset - qset)),
            "n_query_only_off_array": int(len(qset - aset)),
            "frac_queries_on_array": float(len(both) / len(qset)),
            "n_union": int(len(qset | aset)),
        },
        "cells": {"n_grid": int(Y.size), "n_observed": n_obs, "frac_observed": float(n_obs / Y.size),
                  "n_positive": n_pos, "prevalence": float(n_pos / n_obs),
                  "observed_per_query": qdist(q_obs), "observed_per_array": qdist(a_obs)},
        "per_query_sl_counts": {**qdist(q_sl), "hub_share_of_all_SLs": top_share(q_sl),
                                "per_query_prevalence": qdist(q_sl / q_obs)},
        "per_array_sl_counts": {**qdist(a_sl), "hub_share_of_all_SLs": top_share(a_sl)},
        "query_vs_array_side_sl_count_same_gene": side,
        "reciprocal_AB_BA": recip,
        "strains_and_alleles": strains,
        "runtime_s": round(time.time() - t0, 1),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    (OUT / "summary.md").write_text(render_md(summary))
    print(json.dumps({k: summary[k] for k in ("genes", "cells")}, indent=1)[:2000])


def _kappa(a: np.ndarray, b: np.ndarray) -> float:
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe)


def render_md(s: dict) -> str:
    g, c, q, a, r, st = (s["genes"], s["cells"], s["per_query_sl_counts"], s["per_array_sl_counts"],
                         s["reciprocal_AB_BA"], s["strains_and_alleles"])
    qq, aq = q["quantiles"], a["quantiles"]
    side = s["query_vs_array_side_sl_count_same_gene"]
    return f"""# Data characterisation — Costanzo 2016 S1 NxN (aggregate statistics)

Generated by `scripts/00_data_characterisation.py` (numbers in `summary.json`). Label: {s['label']}.
Gene level as upstream: one canonical strain per query ORF (S1: exactly one `_sn` strain per ORF, 30 °C, no
suppressor alleles), arrays = {st['n_array_orfs']} DMA (non-essential deletion) strains, one per ORF.

| quantity | value |
|---|---|
| query ORFs | {g['n_query_orfs']:,} |
| array ORFs | {g['n_array_orfs']:,} |
| query ∩ array | {g['n_query_and_array']:,} ({g['frac_queries_on_array']:.1%} of queries) |
| array-only | {g['n_array_only']:,} |
| query-only (off-array) | {g['n_query_only_off_array']:,} |
| observed cells | {c['n_observed']:,} of {c['n_grid']:,} ({c['frac_observed']:.1%}) |
| strong SL (positives) | {c['n_positive']:,} |
| prevalence | {c['prevalence']:.4f} |

**Per-query SL counts** (query rows): median {qq['q50']:.0f}, mean {q['mean']:.1f}, IQR {qq['q25']:.0f}–{qq['q75']:.0f},
90th pct {qq['q90']:.0f}, 99th pct {qq['q99']:.0f}, max {qq['q100']:.0f}; zero-SL queries {q['n_zero']} ({q['frac_zero']:.1%}).
Hub concentration: top 1% / 5% / 10% / 20% of queries hold
{_v(q, 1):.1%} / {_v(q, 5):.1%} / {_v(q, 10):.1%} / {_v(q, 20):.1%} of all SLs (Gini {q['hub_share_of_all_SLs']['gini']:.2f}).

**Per-array SL counts**: median {aq['q50']:.0f}, mean {a['mean']:.1f}, 90th pct {aq['q90']:.0f}, max {aq['q100']:.0f};
zero-SL arrays {a['n_zero']} ({a['frac_zero']:.1%}); top 10% of arrays hold {_v(a, 10):.1%} of SLs (Gini {a['hub_share_of_all_SLs']['gini']:.2f}).

**Same gene, query side vs array side** ({side['counts']['n']} genes on both sides): SL counts Pearson
{side['counts']['pearson']:.2f}, Spearman {side['counts']['spearman']:.2f}; SL rates Spearman {side['rates']['spearman']:.2f}.

**Reciprocal pairs (AB/BA)** among the {r['n_genes_query_and_array']:,} genes that are both query and array:
{r['n_unordered_pairs_measured_both_directions']:,} of {r['n_unordered_pairs_measured_any_direction']:,} measured unordered
pairs ({r['frac_of_measured_pairs_in_both_directions']:.1%}) are measured in both directions
({r['frac_of_all_observed_cells_that_are_part_of_a_reciprocal_pair']:.1%} of all observed cells). Among them:
P(BA SL | AB SL) = {r['P_BA_pos_given_AB_pos']:.3f}, P(AB SL | BA SL) = {r['P_AB_pos_given_BA_pos']:.3f},
Jaccard of SL calls {r['jaccard_pos']:.3f}, Cohen's κ {r['cohen_kappa']:.3f}, ε Pearson {r['epsilon_pearson']:.3f};
{r['n_pos_both_directions']:,} pairs SL in both directions, {r['n_pos_exactly_one_direction']:,} in exactly one.
Self pairs (query ORF = array ORF) observed: {r['self_pairs_observed(query ORF == array ORF)']}.

**Essential genes / other allele types**: none in S1 NxN (queries all `sn`, arrays all `dma`). The upstream raw zips
also contain SGA_ExN(_NxE), SGA_ExE and SGA_DAmP tables, but upstream only extracted NxN + the strain/SMF table;
those tables are not on disk unzipped and are not used.
"""


def _v(d: dict, pct: int) -> float:
    return next(v for k, v in d["hub_share_of_all_SLs"].items() if k.startswith(f"top{pct}pct"))


if __name__ == "__main__":
    main()
