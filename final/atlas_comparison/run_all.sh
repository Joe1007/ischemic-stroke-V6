#!/bin/bash
# Launch per-atlas 4-task tables (36 methods) and V6 core (24 methods) runs.
cd "$(dirname "$0")"
PY=/home/joe1007/miniconda3/envs/AD/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MPLCONFIGDIR=/tmp/mpl-stroke-v6 PYTHONWARNINGS=ignore
for a in BA JHU JUELICH_25; do
  $PY run_tables_atlas.py --atlas $a --jobs 2 --workers 5 > tables_$a.log 2>&1 &
  for t in mRS NIHSS; do $PY run_core_atlas.py --atlas $a --task $t --jobs 2 > core_${a}_$t.log 2>&1 & done
done
wait
echo ALL_DONE
