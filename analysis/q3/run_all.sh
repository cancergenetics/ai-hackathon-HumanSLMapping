#!/bin/bash
# Re-run every Q3 analysis on one source:   bash analysis/q3/run_all.sh dev|final|dryrun [OUT_DIR]
#  - analyses without model fits run directly (minutes);
#  - LBO-CV refits (stopping estimator, CV thresholds) are a separate, long queue (one P* job at a time):
#      SOURCE=final SPLITS=final TRAJS="0 1" ARMS="R0 H" bash analysis/q3/run_lbo.sh     (~30 min per trajectory)
#    then re-run this script so predcontrib/stopping pick the LBO rows up.
#  - final: run only after scripts/score_final.py has finished (reads results/final/*, incl. run/labels.npz = T labels).
set -e
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2 HSL_THREADS=2
SRC=${1:-dev}; OUTARG=${2:+--out $2}
PY=.venv/bin/python
$PY -m analysis.q3.cache >/dev/null 2>&1 || true   # dev label cache (no T rows); rebuilt if missing
$PY -m analysis.q3.hap1 --source $SRC $OUTARG
$PY -m analysis.q3.frontier --source $SRC $OUTARG
$PY -m analysis.q3.choice profile --source $SRC $OUTARG
if [ "$SRC" = "dev" ]; then
  $PY -m analysis.q3.choice rules --source dev $OUTARG
  [ -f "${2:-results/q3_dev}/choice_eval_rows.jsonl" ] && $PY -m analysis.q3.choice report --source dev $OUTARG
fi
$PY -m analysis.q3.predcontrib --source $SRC $OUTARG
if [ -f "${2:-results/q3_$SRC}/lbo/lbo_rows.jsonl" ]; then
  if [ "$SRC" = "dev" ]; then $PY -m analysis.q3.stopping evaluate --source dev $OUTARG
  else $PY -m analysis.q3.stopping evaluate --source $SRC $OUTARG; $PY -m analysis.q3.stopping apply --source $SRC $OUTARG; fi
fi
$PY -m analysis.q3.report --source $SRC $OUTARG
