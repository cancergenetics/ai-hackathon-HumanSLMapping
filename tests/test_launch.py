"""scripts/launch.py: concurrency policy, timeouts, logs, thread env."""
from __future__ import annotations

import hsl  # noqa: F401

import importlib.util
import json
import sys

from hsl import config, paths


def _mod():
    spec = importlib.util.spec_from_file_location("launch", paths.PROJECT_ROOT / "scripts" / "launch.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_policy():
    L = _mod()
    c = config.protocol()["compute"]
    small = [{"cls": "small", "rss": 0.5}] * 3
    assert L.job_class({"n": 660}, c) == "small" and L.job_class({"n": 1310}, c) == "large"
    assert L.job_class({"n": 2627}, c) == "alone" and L.job_class({}, c) == "large"
    assert L.can_start("small", small, c, None)                       # 4th small job
    assert not L.can_start("small", small + [{"cls": "small", "rss": 0.5}], c, None)
    assert not L.can_start("large", [{"cls": "large"}, {"cls": "small"}], c, None)
    assert L.can_start("large", [{"cls": "small"}], c, None)
    assert not L.can_start("alone", [{"cls": "small"}], c, None)
    assert not L.can_start("small", [{"cls": "small", "rss": 9.0}], c, None)   # total RSS cap
    assert L.threads_for("small", c) == 2 and L.threads_for("large", c) == 3


def test_run_jobs_timeout_and_logs(tmp_path, monkeypatch):
    monkeypatch.setenv("HSL_LOGS_DIR", str(tmp_path))
    L = _mod()
    jobs = [{"name": "ok", "cmd": [sys.executable, "-c", "import os;print(os.environ['OMP_NUM_THREADS'])"], "n": 30},
            {"name": "slow", "cmd": [sys.executable, "-c", "import time;time.sleep(30)"], "timeout_s": 1},
            {"name": "bad", "cmd": [sys.executable, "-c", "raise SystemExit(3)"]}]
    jf = tmp_path / "jobs.json"
    jf.write_text(json.dumps(jobs))
    rc = L.main([str(jf)])
    assert rc == 2
    assert (tmp_path / "ok.log").read_text().strip() == "2"
    recs = {r["name"]: r for f in tmp_path.glob("launch_*.jsonl") for r in map(json.loads, f.read_text().splitlines())}
    assert recs["ok"]["status"] == "ok" and recs["slow"]["status"] == "timeout" and recs["bad"]["status"] == "exit_3"


def test_memory_gate_opt_in():
    L = _mod()
    c = config.protocol()["compute"]
    j = {"n": 660, "rss_gb": 3.0}
    assert L.gate_ok(j, [], c, 2.0, 6, 2, free_fn=lambda: 0.0)                    # nothing running: progress
    run = [{"n": 160, "rss_hint": 3.0, "rss": 1.0}]                                # 2 GB not yet claimed
    assert L.gate_ok(j, run, c, 2.0, 6, 2, free_fn=lambda: 7.0)                    # 7 - 2 >= 3 + 2
    assert not L.gate_ok(j, run, c, 2.0, 6, 2, free_fn=lambda: 6.9)
    assert not L.gate_ok(j, run * 6, c, 2.0, 6, 2, free_fn=lambda: 100.0)          # job cap
    big = {"n": 2627, "rss_gb": 9.0}
    two_big = [{"n": 1970, "rss_hint": 7.0, "rss": 7.0}] * 2
    assert not L.gate_ok(big, two_big, c, 2.0, 6, 2, free_fn=lambda: 100.0)        # <= 2 jobs at n >= 1970
    assert L.gate_ok(big, two_big[:1], c, 2.0, 6, 2, free_fn=lambda: 11.0)
    assert not L.gate_ok(j, run, c, 2.0, 6, 2, free_fn=lambda: None)               # unknown memory: wait
    assert L.job_hint({}, c) == float(c["rss_cap_gb_per_job"])
    free = L.unused_gb()
    assert free is None or free >= 0


def test_gate_runs_jobs_and_job_rss_cap(tmp_path, monkeypatch):
    monkeypatch.setenv("HSL_LOGS_DIR", str(tmp_path))
    L = _mod()
    jobs = [{"name": f"j{i}", "cmd": [sys.executable, "-c", "pass"], "n": 30, "rss_gb": 0.1, "rss_cap_gb": 7.5}
            for i in range(3)]
    jf = tmp_path / "jobs.json"
    jf.write_text(json.dumps(jobs))
    assert L.main([str(jf), "--min-free-gb", "0.1"]) == 0
    recs = [r for f in tmp_path.glob("launch_*.jsonl") for r in map(json.loads, f.read_text().splitlines())]
    assert len(recs) == 3 and all(r["status"] == "ok" and r["rss_cap_gb"] == 7.5 for r in recs)
