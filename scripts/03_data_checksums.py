#!/usr/bin/env python
"""SHA-256 of every upstream data file we may read in place (YEAST_SLs/data/raw/** and data/processed/** except the
upstream prediction files), cross-checked against the checksums upstream recorded (YEAST_SLs/DATA.md, every
MANIFEST.json under data/raw, the bundle sidecars, the pair/co-fitness store manifests).

Output: results/data_checksums/checksums.json  (path relative to YEAST_SLs, bytes, sha256, recorded_upstream: bool)
Read-only on YEAST_SLs (write guard active).
"""
from __future__ import annotations

import hsl  # noqa: F401  (YEAST_SLs write guard + paths; must be first)

import hashlib
import json
import re
import time

from hsl import paths

OUT = paths.output_root() / "results" / "data_checksums"
SKIP_DIRS = ("data/processed/predictions",)
SHA_RE = re.compile(r"\b[0-9a-f]{64}\b")


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 24):
            h.update(b)
    return h.hexdigest()


def recorded_hashes(root) -> set[str]:
    texts = [(root / "DATA.md").read_text(encoding="utf-8")]
    for pat in ("data/raw/*/MANIFEST.json", "data/processed/*.json", "data/processed/features/*/manifest.json",
                "data/processed/features/*.json", "config/frozen.yaml"):
        texts += [p.read_text(encoding="utf-8", errors="replace") for p in sorted(root.glob(pat))]
    return set(SHA_RE.findall("\n".join(texts)))


def main() -> None:
    t0 = time.time()
    root = paths.yeast_sls_root()
    rec = recorded_hashes(root)
    files = []
    for top in ("data/raw", "data/processed"):
        for p in sorted((root / top).rglob("*")):
            rel = str(p.relative_to(root))
            if not p.is_file() or p.name == ".DS_Store" or any(rel.startswith(s) for s in SKIP_DIRS):
                continue
            h = sha256(p)
            files.append({"path": rel, "bytes": p.stat().st_size, "sha256": h, "recorded_upstream": h in rec})
    OUT.mkdir(parents=True, exist_ok=True)
    out = {"yeast_sls_root": "YEAST_SLs (config/paths.yaml)", "commit": (paths.PROJECT_ROOT / "config" /
                                                                         "yeast_sls_commit.txt").read_text().strip(),
           "computed": time.strftime("%Y-%m-%dT%H:%M:%S"), "n_files": len(files),
           "n_matching_an_upstream_record": sum(f["recorded_upstream"] for f in files),
           "skipped": list(SKIP_DIRS) + ["*.DS_Store"], "runtime_s": round(time.time() - t0, 1), "files": files}
    (OUT / "checksums.json").write_text(json.dumps(out, indent=1))
    print(f"{len(files)} files, {out['n_matching_an_upstream_record']} match an upstream record, {out['runtime_s']} s")


if __name__ == "__main__":
    main()
