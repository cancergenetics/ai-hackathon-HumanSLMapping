"""The YEAST_SLs write guard (hsl.guard) must refuse writes BEFORE the OS call.

Safety of the test itself: for the real YEAST_SLs root we never call a real writing function. We only fire the audit
events with `sys.audit(...)` (no OS call is made by sys.audit) and assert the hook raises. Real write calls are only
exercised against a throw-away protected root inside pytest's tmp_path.
"""
from __future__ import annotations

import os
import subprocess
import sys

import pytest

import hsl  # noqa: F401  (installs the guard)
from hsl import guard, paths

YEAST = paths.yeast_sls_root()
NEVER = YEAST / "results" / "__hsl_guard_test_must_never_exist__.txt"


def test_guard_installed_for_yeast_root():
    assert os.path.realpath(YEAST) in guard.protected_roots()
    assert YEAST.name == "YEAST_SLs"


@pytest.mark.parametrize("mode", ["w", "wb", "a", "x", "r+", "w+b"])
def test_open_write_event_on_yeast_is_refused(mode):
    with pytest.raises(guard.ReadOnlyViolation):
        sys.audit("open", str(NEVER), mode, os.O_WRONLY | os.O_CREAT)
    assert not NEVER.exists()


def test_os_open_create_flags_on_yeast_refused():
    with pytest.raises(guard.ReadOnlyViolation):
        sys.audit("open", str(NEVER), None, os.O_CREAT | os.O_EXCL | os.O_WRONLY)   # what os.open / importlib emit


@pytest.mark.parametrize("event,args", [
    ("os.remove", (str(NEVER), -1)),
    ("os.rmdir", (str(YEAST / "results"), -1)),
    ("os.rename", (str(NEVER), "/tmp/x", -1, -1)),
    ("os.rename", ("/tmp/x", str(NEVER), -1, -1)),
    ("os.mkdir", (str(YEAST / "newdir"), 0o777, -1)),
    ("os.symlink", ("/tmp/x", str(NEVER), -1)),
    ("os.truncate", (str(NEVER), 0)),
    ("os.chmod", (str(NEVER), 0o644, -1)),
    ("os.utime", (str(NEVER), None, None, -1)),
    ("shutil.rmtree", (str(YEAST / "results"), None)),
    ("shutil.copyfile", ("/tmp/x", str(NEVER))),
    ("shutil.move", (str(NEVER), "/tmp/x")),
    ("subprocess.Popen", ("git", ["git", "status"], str(YEAST), None)),
    ("subprocess.Popen", ("rm", ["rm", str(NEVER)], "/tmp", None)),
])
def test_mutating_events_on_yeast_refused(event, args):
    with pytest.raises(guard.ReadOnlyViolation):
        sys.audit(event, *args)
    assert not NEVER.exists()


def test_reads_of_yeast_still_allowed():
    sys.audit("open", str(YEAST / "README.md"), "r", os.O_RDONLY)   # no exception
    with open(YEAST / "README.md", encoding="utf-8") as f:            # real read-only open works
        assert f.read(10)
    sys.audit("subprocess.Popen", "ls", ["ls", str(YEAST)], "/tmp", None)  # read-only tool mentioning the root


def test_writes_elsewhere_allowed(tmp_path):
    p = tmp_path / "ok.txt"
    p.write_text("x")
    assert p.read_text() == "x"


def test_real_write_calls_blocked_before_os_call(tmp_path):
    """Real Python write APIs against a throw-away protected root: refused and nothing is created."""
    fake = tmp_path / "fake_yeast"
    fake.mkdir()
    (fake / "keep.txt").write_text("keep")
    guard.install([fake], guard_pyarrow=True)
    target = fake / "new.txt"
    with pytest.raises(guard.ReadOnlyViolation):
        open(target, "w")
    with pytest.raises(guard.ReadOnlyViolation):
        target.write_text("x")
    with pytest.raises(guard.ReadOnlyViolation):
        target.touch()
    with pytest.raises(guard.ReadOnlyViolation):
        os.open(target, os.O_CREAT | os.O_WRONLY)
    with pytest.raises(guard.ReadOnlyViolation):
        (fake / "sub").mkdir()
    with pytest.raises(guard.ReadOnlyViolation):
        os.remove(fake / "keep.txt")
    with pytest.raises(guard.ReadOnlyViolation):
        os.replace(fake / "keep.txt", tmp_path / "moved.txt")
    with pytest.raises(guard.ReadOnlyViolation):
        open(fake / "keep.txt", "a")
    import numpy as np
    with pytest.raises(guard.ReadOnlyViolation):
        np.save(fake / "a.npy", np.zeros(3))
    import pyarrow as pa
    import pyarrow.parquet as pq
    with pytest.raises(guard.ReadOnlyViolation):
        pq.write_table(pa.table({"a": [1]}), str(fake / "t.parquet"))
    import pandas as pd
    with pytest.raises(guard.ReadOnlyViolation):
        pd.DataFrame({"a": [1]}).to_parquet(fake / "t2.parquet")
    with pytest.raises(guard.ReadOnlyViolation):
        subprocess.run(["true"], cwd=fake)
    assert sorted(os.listdir(fake)) == ["keep.txt"]
    assert (fake / "keep.txt").read_text() == "keep"   # reading is fine


def test_vendored_ysl_resolves_upstream_root_and_writable_output_elsewhere():
    from ysl import _paths, runner
    assert _paths.UPSTREAM_ROOT == YEAST.resolve()
    assert not str(_paths.OUTPUT_ROOT).startswith(str(YEAST.resolve()))
    assert not str(runner.RESULTS.resolve()).startswith(str(YEAST.resolve()))
    assert not str(runner.PRED_DIR.resolve()).startswith(str(YEAST.resolve()))
    assert guard.is_protected(runner.ROOT / "data" / "processed" / "x")
