"""P4 review (independent): real-data selection replay of a stored dev trajectory (H, dev1, trajectory 0) under bundle
poisonings (DESIGN §6.1). Selection only, no predictor fit (~5 s per case). DEV ONLY. Skips if the stored trajectory
is absent. Full matrix (H/C/F_r16 x dev1/dev2 x 4 poisonings + 2 positive controls, n = 160):
scripts/review/p4_poison_replay.py -> results/review_p4/poison_replay.jsonl."""
from __future__ import annotations

import hsl  # noqa: F401

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
N = 60


def _mod():
    spec = importlib.util.spec_from_file_location("p4_poison_replay", ROOT / "scripts" / "review" / "p4_poison_replay.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _replay(poison, cfg="H", split="dev1", r=0):
    from hsl import harness as H, q2dev
    part = q2dev.RESULTS / "p4" / "parts" / f"{cfg}_{split}_r{r}.jsonl"
    if not part.exists():
        pytest.skip("stored P4 trajectory not present")
    tid = json.loads(part.read_text().splitlines()[0])["traj_id"]
    stored = json.loads((H.TRAJ_DIR / split / f"{tid}.json").read_text())["selected"]
    m = _mod()
    b, _ = m.poisoned_bundle(poison, q2dev.dev_split(split), stored, N)
    spec = q2dev.parse(cfg)
    t = H.Trajectory(q2dev.make_selector(spec, split), split, q2dev.traj_seed(r), r, arm=cfg, checkpoints=(), max_n=N,
                     predictors=(), save=False, ledger=False, write_log=False, bundle=b)
    return t.run()["selected"][:N], stored[:N]


@pytest.mark.parametrize("poison", ["none", "dval_rows", "dval_cols", "unscreened_rows", "dval_smf"])
def test_h_selection_invariant_to_masked_channels(poison):
    got, ref = _replay(poison)
    assert got == ref


def test_h_selection_changes_when_visible_columns_poisoned():
    got, ref = _replay("pc_vis_cols")
    assert got != ref
