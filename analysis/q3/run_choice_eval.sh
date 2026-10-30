#!/bin/bash
# Static rank-rule evaluation on dev (P* refits; one process per (rule, budget); RSS watchdog 3 GB). Resumable.
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2 HSL_THREADS=2
for rule in forward_k1 forward_k3; do for n in 160 330 660; do
  .venv/bin/python -m analysis.q3.choice evaluate --source dev --rule $rule --budgets $n || echo "FAILED $rule n=$n"
done; done
echo "choice eval done"
