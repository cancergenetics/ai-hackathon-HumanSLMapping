"""SGD chromosomal feature coordinates (data/raw/sgd/SGD_features.tab; scripts/01b_download_sgd_features.py).

Only coordinates are used (feature name, chromosome, start, stop). Description/annotation columns are ignored.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ysl._paths import UPSTREAM_ROOT as ROOT  # humanSLmapping PATCH (was Path(__file__).resolve().parents[N])
SGD_FEATURES = ROOT / "data" / "raw" / "sgd" / "SGD_features.tab"
COLS = ["sgdid", "feature_type", "qualifier", "feature_name", "gene_name", "alias", "parent", "sgdid2",
        "chromosome", "start", "stop", "strand", "genetic_pos", "coord_version", "seq_version", "description"]


def orf_coordinates(path: Path | str = SGD_FEATURES) -> pd.DataFrame:
    """One row per named chromosomal feature with coordinates: orf, chromosome (int; 17 = mito), start, end,
    mid (bp). start/end are min/max of the two coordinates (strand-independent)."""
    df = pd.read_csv(path, sep="\t", header=None, names=COLS, dtype=str, keep_default_na=False)
    df = df[(df["feature_name"] != "") & (df["chromosome"] != "") & (df["start"] != "")]
    df = df[df["parent"].str.startswith("chromosome")]  # top-level features only (not CDS/intron sub-features)
    out = pd.DataFrame({
        "orf": df["feature_name"].values,
        "feature_type": df["feature_type"].values,
        "chromosome": pd.to_numeric(df["chromosome"]).astype(int).values,
        "start": pd.to_numeric(df["start"]).astype(int).values,
        "stop": pd.to_numeric(df["stop"]).astype(int).values,
    })
    out["lo"] = out[["start", "stop"]].min(axis=1)
    out["hi"] = out[["start", "stop"]].max(axis=1)
    out["mid"] = (out["lo"] + out["hi"]) / 2
    return out.drop_duplicates("orf").set_index("orf")[["feature_type", "chromosome", "lo", "hi", "mid"]]
