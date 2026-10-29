"""Write guard: make the upstream YEAST_SLs checkout read-only for this Python process.

`install(roots)` adds a `sys.addaudithook` hook (cannot be removed for the life of the process) that raises
`ReadOnlyViolation` (a PermissionError) BEFORE the operating-system call when Python code tries to modify anything
under a protected root:

  * open / io.open / os.open / pathlib writes / np.save / torch.save / tarfile / zipfile ... ("open" audit event)
    in a write, append, create, truncate or update mode (mode contains w, a, x or +, or flags contain O_WRONLY,
    O_RDWR, O_CREAT, O_TRUNC or O_APPEND). Read-only opens are allowed.
  * os.remove/unlink, os.rmdir, os.rename/replace (source or destination), os.mkdir/makedirs, os.symlink, os.link,
    os.truncate, os.chmod, os.chown, os.utime, os.chflags, os.lchflags, os.mkfifo, os.mknod, os.setxattr,
    os.removexattr, shutil.rmtree/copyfile/copytree/move/make_archive/unpack_archive (destination),
    tempfile.mkstemp/mkdtemp, sqlite3.connect.
  * subprocess.Popen / os.system / os.exec* / os.posix_spawn whose working directory (cwd) is inside a protected
    root, or whose command line mentions a protected root together with an obviously mutating tool
    (git, uv, pip, rm, mv, cp, touch, tee, sed -i, ...). Read-only `find`/`ls`/`cat` are allowed.
  * pyarrow writers (which open files in C++ and raise no Python audit event) are wrapped: pyarrow.parquet
    write_table / ParquetWriter / write_to_dataset, pyarrow.feather.write_feather, pyarrow.csv.write_csv,
    pyarrow.dataset.write_dataset, pyarrow.output_stream / OSFile(mode != 'r') / create_memory_map.

Python's bytecode cache writes (importlib) go through os.open and are therefore also refused; importlib swallows the
OSError, so importing from a protected tree still works without writing __pycache__.

Limits (documented, not silently assumed): C extensions that write via their own C file APIs (other than the pyarrow
writers above) and child processes started with a cwd outside the root are not audited. Every entry point imports
this module (via `import hsl`), and the final read-only check (`find -newer .yeast_readonly_stamp`, upstream
`git status`) is the backstop.
"""
from __future__ import annotations

import os
import sys
import threading

__all__ = ["ReadOnlyViolation", "install", "protected_roots", "is_protected", "check_open"]


class ReadOnlyViolation(PermissionError):
    """Raised when code tries to modify a protected (read-only) tree."""


_STATE = {"roots": [], "hook": False, "pyarrow": False}
_LOCK = threading.Lock()

_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND
_PATH_EVENTS = {  # event -> indices of path arguments that must not be protected
    "os.remove": (0,), "os.rmdir": (0,), "os.rename": (0, 1), "os.mkdir": (0,), "os.symlink": (1,),
    "os.link": (1,), "os.truncate": (0,), "os.chmod": (0,), "os.chown": (0,), "os.utime": (0,),
    "os.chflags": (0,), "os.lchflags": (0,), "os.mkfifo": (0,), "os.mknod": (0,), "os.setxattr": (0,),
    "os.removexattr": (0,), "shutil.rmtree": (0,), "shutil.copyfile": (1,), "shutil.copymode": (1,),
    "shutil.copystat": (1,), "shutil.copytree": (1,), "shutil.move": (0, 1), "shutil.chown": (0,),
    "shutil.make_archive": (0, 1), "shutil.unpack_archive": (1,), "tempfile.mkstemp": (0,),
    "tempfile.mkdtemp": (0,), "sqlite3.connect": (0,),
}
_MUTATING_TOOLS = ("git", "uv", "pip", "pip3", "rm", "rmdir", "mv", "cp", "touch", "tee", "mkdir", "ln", "chmod",
                   "chown", "truncate", "rsync", "dd", "install", "pytest", "python", "python3", "make")


def protected_roots() -> list[str]:
    return list(_STATE["roots"])


def _norm(path) -> str | None:
    if isinstance(path, int) or path is None:
        return None
    try:
        p = os.fsdecode(os.fspath(path))
    except TypeError:
        return None
    if p == "" or p.startswith(("file:", "http:", "https:", "s3:")):
        return None
    return os.path.realpath(os.path.abspath(os.path.expanduser(p)))


def is_protected(path) -> bool:
    p = _norm(path)
    if p is None:
        return False
    return any(p == r or p.startswith(r + os.sep) for r in _STATE["roots"])


def _deny(what: str, path) -> None:
    raise ReadOnlyViolation(f"YEAST_SLs write guard: {what} on protected path {os.fsdecode(os.fspath(path))!r} "
                            f"(protected roots: {_STATE['roots']})")


def _is_write_mode(mode, flags) -> bool:
    if isinstance(mode, str) and any(c in mode for c in "wax+"):
        return True
    if isinstance(flags, int) and flags & _WRITE_FLAGS:
        return True
    return False


def check_open(path, mode=None, flags=None) -> None:
    """The check run for an 'open' audit event (exposed for tests)."""
    if _is_write_mode(mode, flags) and is_protected(path):
        _deny(f"open(mode={mode!r}, flags={flags!r})", path)


def _cmd_tokens(args) -> list[str]:
    if isinstance(args, (str, bytes, os.PathLike)):
        return os.fsdecode(args).split()
    try:
        return [os.fsdecode(a) if isinstance(a, (str, bytes, os.PathLike)) else str(a) for a in args]
    except TypeError:
        return []


def _check_command(tokens: list[str], cwd) -> None:
    if cwd is not None and is_protected(cwd):
        _deny("subprocess with working directory", cwd)
    if not tokens:
        return
    text = " ".join(tokens)
    if not any(r in text for r in _STATE["roots"]) and "YEAST_SLs" not in text:
        return
    words = [os.path.basename(t) for t in " ".join(tokens).replace(";", " ").replace("&&", " ").split()]
    if any(w in _MUTATING_TOOLS for w in words) or "-i" in words or ">" in text:
        _deny(f"subprocess {text[:200]!r} mentioning", next((r for r in _STATE["roots"]), "YEAST_SLs"))


def _hook(event: str, args: tuple) -> None:
    if not _STATE["roots"]:
        return
    if event == "open":
        if args:
            check_open(args[0], args[1] if len(args) > 1 else None, args[2] if len(args) > 2 else None)
        return
    idx = _PATH_EVENTS.get(event)
    if idx is not None:
        for i in idx:
            if i < len(args) and is_protected(args[i]):
                _deny(event, args[i])
        return
    if event == "subprocess.Popen":           # (executable, args, cwd, env)
        _check_command(_cmd_tokens(args[1]) if len(args) > 1 else [], args[2] if len(args) > 2 else None)
    elif event == "os.system":                # (command,)
        _check_command(_cmd_tokens(args[0]), None)
    elif event in ("os.exec", "os.posix_spawn", "os.spawn"):   # (path, args, env) / (mode, path, args, env)
        toks = [a for a in args if isinstance(a, (list, tuple))]
        _check_command(_cmd_tokens(toks[0]) if toks else [], None)


# ----------------------------------------------------------------------------- pyarrow writers (C++ file I/O)
def _wrap_first_path(fn, what: str, path_pos: int = 1, path_kw=("where", "path", "root_path", "base_dir",
                                                                  "output_file", "dest", "sink")):
    def wrapper(*a, **k):
        cands = [a[path_pos]] if len(a) > path_pos else []
        cands += [k[x] for x in path_kw if x in k]
        for c in cands:
            if isinstance(c, (str, bytes, os.PathLike)) and is_protected(c):
                _deny(what, c)
        return fn(*a, **k)
    wrapper.__wrapped__ = fn
    wrapper.__name__ = getattr(fn, "__name__", what)
    return wrapper


def _wrap_osfile(cls):
    class GuardedOSFile(cls):  # type: ignore[misc, valid-type]
        # the Cython base opens the file in __cinit__ (tp_new), so the check must run in __new__, not __init__
        def __new__(klass, path, mode="r", *a, **k):
            if mode != "r" and is_protected(path):
                _deny(f"pyarrow.OSFile(mode={mode!r})", path)
            return super().__new__(klass, path, mode, *a, **k)
    GuardedOSFile.__name__ = "OSFile"
    return GuardedOSFile


def _guard_pyarrow() -> None:
    if _STATE["pyarrow"]:
        return
    try:
        import pyarrow as pa
        import pyarrow.csv as pcsv
        import pyarrow.dataset as pds
        import pyarrow.feather as pf
        import pyarrow.parquet as pq
    except ImportError:
        return
    pq.write_table = _wrap_first_path(pq.write_table, "pyarrow.parquet.write_table")
    pq.write_to_dataset = _wrap_first_path(pq.write_to_dataset, "pyarrow.parquet.write_to_dataset")
    pq.ParquetWriter = _wrap_parquet_writer(pq.ParquetWriter)
    pf.write_feather = _wrap_first_path(pf.write_feather, "pyarrow.feather.write_feather")
    pcsv.write_csv = _wrap_first_path(pcsv.write_csv, "pyarrow.csv.write_csv")
    pds.write_dataset = _wrap_first_path(pds.write_dataset, "pyarrow.dataset.write_dataset")
    pa.output_stream = _wrap_first_path(pa.output_stream, "pyarrow.output_stream", path_pos=0)
    pa.create_memory_map = _wrap_first_path(pa.create_memory_map, "pyarrow.create_memory_map", path_pos=0)
    try:
        pa.OSFile = _wrap_osfile(pa.OSFile)
    except TypeError:  # final class: cannot subclass; leave as is (documented limit)
        pass
    _STATE["pyarrow"] = True


def _wrap_parquet_writer(cls):
    class GuardedParquetWriter(cls):  # type: ignore[misc, valid-type]
        def __init__(self, where, *a, **k):
            if isinstance(where, (str, bytes, os.PathLike)) and is_protected(where):
                _deny("pyarrow.parquet.ParquetWriter", where)
            super().__init__(where, *a, **k)
    GuardedParquetWriter.__name__ = "ParquetWriter"
    return GuardedParquetWriter


# ----------------------------------------------------------------------------- install
def install(roots=None, guard_pyarrow: bool = True) -> list[str]:
    """Protect `roots` (default: the YEAST_SLs root from config/paths.yaml). Idempotent; roots accumulate."""
    if roots is None:
        from hsl.paths import yeast_sls_root
        roots = [yeast_sls_root()]
    with _LOCK:
        for r in roots:
            rp = os.path.realpath(os.path.abspath(os.fspath(r)))
            if rp == os.sep:
                raise ValueError("refusing to protect '/'")
            if rp not in _STATE["roots"]:
                _STATE["roots"].append(rp)
        if not _STATE["hook"]:
            sys.addaudithook(_hook)
            _STATE["hook"] = True
    if guard_pyarrow:
        _guard_pyarrow()
    return protected_roots()
