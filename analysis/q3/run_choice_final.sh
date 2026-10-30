#!/bin/bash
# ONE-TIME evaluation on T of the rank rules fixed in bde1622 (post-final descriptive; not part of the frozen plan).
# One process per (rule, budget); RSS watchdog 4 GB (PI final budget). Resumable; never refit.
cd "$(dirname "$0")/../.."
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2 HSL_THREADS=2
for rule in forward_k3 forward_k1; do for n in 160 330 660; do
  .venv/bin/python -m analysis.q3.choice evaluate --source final --rule $rule --budgets $n --rss-cap-gb 4.0 || echo "FAILED $rule n=$n"
done; done
echo "choice final eval done"
