#!/bin/sh
# waits for the "all" group, then: B3 profiles, clean timings (one job at a time), classical timing, quadrature check
cd ~/SPINN_Allen_Khan/code
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=~/spectral_pinn_test/venv/bin/python
until grep -q "^finished all" ../data_r1/all.log; do sleep 30; done
$PY r1_profiles.py > ../data_r1/profiles.log 2>&1
$PY r1_run.py timing 1 > ../data_r1/timing.log 2>&1
$PY r1_classical2.py timing > ../data_r1/classical2_timing.log 2>&1
$PY r1_check_quad.py > ../data_r1/check_quad.log 2>&1
echo FINALIZE_DONE >> ../data_r1/timing.log
