"""P* fit on an arbitrary screened set of a dev split, scored on that split's E_fixed_dev cells — the same calls
hsl.harness.run_predictor makes (upstream.construct_view + upstream.run_method with the frozen P* params,
n_folds = min(10, max(n, 2)), split_seed 0), WITHOUT the harness side effects (no write to the tracked
results/LEDGER.jsonl, no prediction files under results/predictions). Dev splits only.

`check` reproduces a stored P1 random-arm AUPRC exactly (same screened set and model seed).
    .venv/bin/python -m analysis.q3.fit check
"""
from __future__ import annotations

import sys
import time

import numpy as np
import yaml

from analysis.q3 import common as C

C.set_threads(2)


def pstar_spec():
    import hsl.predictors  # noqa: F401  (registers P*)
    pz = yaml.safe_load((C.ROOT / "config" / "predictor.yaml").read_text())["primary"]
    return pz["method"], dict(pz["params"] or {})


def fit_score(b, split: str, screened, model_seed: int, n_threads: int = 2) -> dict:
    from hsl import eval as E, upstream as U, views
    from analysis.q3 import cache
    sp = views.get_split(split)
    if sp.kind != "dev":
        raise PermissionError("analysis.q3.fit runs on dev splits only")
    S = np.sort(np.asarray(list(screened), dtype=str))
    if not np.isin(S, sp.pool).all():
        raise ValueError("screened genes must be in the split's pool")
    es = cache.eset(split)
    if not np.array_equal(es["targets"], np.sort(sp.targets)):
        raise AssertionError("cache targets differ from the split targets")
    method, params = pstar_spec()
    pr = C.protocol()["predictor"]
    t0 = time.perf_counter()
    view = U.construct_view(b, S, sp.targets)
    Sm, _ = U.run_method(view, method, model_seed=int(model_seed), split_seed=0, n_threads=n_threads,
                         n_folds=int(min(int(pr["n_folds"]), max(len(S), 2))),
                         chunk_targets=int(pr["chunk_targets"]), params=params)
    s = np.asarray(Sm, dtype=np.float32)[es["rows"], es["cols"]]
    y = es["y"]
    out = {"split": split, "n": int(len(S)), "model_seed": int(model_seed), "fit_s": time.perf_counter() - t0,
           "auprc": E.average_precision(y, s),
           "auprc_VV": E.average_precision(y[es["VV"]], s[es["VV"]]),
           "auprc_AO": E.average_precision(y[es["AO"]], s[es["AO"]])}
    return out, s


def check() -> None:
    import pandas as pd
    b = C.dev_bundle()
    sels = C.selections("dev", arms=[C.RANDOM])
    d = pd.read_csv(C.ROOT / "results" / "p1" / "curves.csv")
    for sp, r, n in (("dev1", 0, 160),):
        sel = sels[(C.RANDOM, sp, r)]
        x = d[(d.predictor == C.PSTAR) & (d.split == sp) & (d.perm == r) & (d.n == n)].iloc[0]
        out, _ = fit_score(b, sp, sel["selected"][:n], sel["model_seed"])
        print(f"{sp} r{r} n{n}: ours {out['auprc']:.10f} stored {x.auprc:.10f} diff {out['auprc'] - x.auprc:.2e} "
              f"fit {out['fit_s']:.0f} s")
        if abs(out["auprc"] - x.auprc) > 1e-6:
            sys.exit("MISMATCH")


if __name__ == "__main__":
    if sys.argv[1:] == ["check"]:
        check()
