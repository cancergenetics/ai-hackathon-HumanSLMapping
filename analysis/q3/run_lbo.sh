#!/bin/bash
# Sequential LBO-CV queue (one P* fit at a time, 2 BLAS threads), any source. Resumable.
#   SOURCE=final SPLITS=final TRAJS="0 1" ARMS="R0 H" bash analysis/q3/run_lbo.sh
# Defaults reproduce the dev run (analysis/q3/run_lbo_dev.sh). Final mode reads results/final/selections.json.gz and
# fits on the DEV bundle (screened genes are pool genes; no T row is ever loaded; T columns are excluded from cells).
set -e
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2 HSL_THREADS=2
SOURCE=${SOURCE:-dev}; SPLITS=${SPLITS:-"dev1 dev2 dev3"}; TRAJS=${TRAJS:-"0 1"}; ARMS=${ARMS:-"R0 H"}
CPS=${CPS:-30,60,100,160,220,330,490,660}
for r in $TRAJS; do for s in $SPLITS; do for arm in $ARMS; do
  for c in ${CPS//,/ }; do  # one process per checkpoint (peak RSS < 3 GB)
    .venv/bin/python -m analysis.q3.lbo run --source $SOURCE --arm $arm --split $s --traj $r --checkpoints $c ${OUT:+--out $OUT} ${RSS:+--rss-cap-gb $RSS}
  done
done; done; done
echo "LBO queue done ($SOURCE)"
