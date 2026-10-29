#!/usr/bin/env bash
# YEAST_SLs read-only check (DESIGN v3 §11). Exit 0 = clean; non-zero = any finding.
#  1. stamp check: no file or directory under YEAST_SLs is newer than .yeast_readonly_stamp, except the `.git`
#     directory entry itself (its mtime moves when git reads the repo);
#  2. `GIT_OPTIONAL_LOCKS=0 git -C YEAST_SLs status --short` is empty (no modified/untracked files).
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
STAMP="$HERE/.yeast_readonly_stamp"
YEAST="$(cd "$HERE" && cd "$(.venv/bin/python -c 'import yaml;print(yaml.safe_load(open("config/paths.yaml"))["yeast_sls_root"])')" && pwd)"
fail=0
[ -f "$STAMP" ] || { echo "MISSING stamp $STAMP"; exit 2; }
echo "YEAST_SLs: $YEAST"
echo "stamp:     $STAMP ($(stat -f '%Sm' "$STAMP" 2>/dev/null || stat -c '%y' "$STAMP"))"
newer="$(find "$YEAST" -newer "$STAMP" ! -path "$YEAST/.git" 2>/dev/null)"
if [ -n "$newer" ]; then
  echo "FINDING: entries newer than the stamp:"; echo "$newer" | head -50; fail=1
else
  echo "stamp check: OK (no entry newer than the stamp other than the .git directory itself)"
fi
st="$(GIT_OPTIONAL_LOCKS=0 git -C "$YEAST" status --short 2>&1)"
if [ -n "$st" ]; then
  echo "FINDING: git status not clean:"; echo "$st" | head -50; fail=1
else
  echo "git status --short: clean (HEAD $(GIT_OPTIONAL_LOCKS=0 git -C "$YEAST" rev-parse --short HEAD))"
fi
exit $fail
