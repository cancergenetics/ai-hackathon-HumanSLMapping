#!/bin/bash
# Sequential LBO-CV queue on dev (one P* fit at a time, 2 BLAS threads). Resumable: finished (arm, split, traj, n)
# rows in results/q3_dev/lbo/lbo_rows.jsonl are skipped. Usage: bash analysis/q3/run_lbo_dev.sh
set -e
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2 HSL_THREADS=2
CPS=${CPS:-30,60,100,160,220,330,490,660}
for r in 0 1; do
  for s in dev1 dev2 dev3; do
    for arm in R0 H; do
      .venv/bin/python -m analysis.q3.lbo run --source dev --arm $arm --split $s --traj $r --checkpoints $CPS
    done
  done
done
echo "LBO queue done"
