#!/bin/bash
set -e

n_rest=1
n_jump=100
n_steps=1000000

L=14400

NR0=1
NR1=29
NRD=4

for n_realiz in $(seq ${NR0} ${NRD} ${NR1}); do
  sbatch -J phi4_${n_realiz}_${L} meluxina-run_1.sh ${L} ${n_realiz} ${n_steps} ${n_jump} ${n_rest}
done
