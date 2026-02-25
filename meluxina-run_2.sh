#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --qos=default
#SBATCH --account=p200715
###p200243
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH -J PHI4_MMC
#SBATCH --mail-user=loris.dicairano@uni.lu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err


export OMP_NUM_THREADS=1


ulimit -s unlimited
set -e

module purge
module load env/release/2021.5
module load foss/2021a

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir: ${SLURM_SUBMIT_DIR}"
cd "${SLURM_SUBMIT_DIR}"

L=$1
n_realiz=$2
n_steps=$3
n_jump=$4
n_rest=$5




EXE="${SLURM_SUBMIT_DIR}/phi4_mmc"
INPUT="${SLURM_SUBMIT_DIR}/input.inp"

NSAMP_MIN=1000
NSAMP_MAX=2550
NSAMP_STEP=50

NSAMP_PAR=32     # 32 nsamp in parallel
NREAL_PAR=4      # 4 realiz per nsamp  -> 32*4 = 128

active=0
for nsamp in $(seq "${NSAMP_MIN}" "${NSAMP_STEP}" "${NSAMP_MAX}"); do
  for k in $(seq 0 $((NREAL_PAR-1))); do
    nr=$((n_realiz + k))
    srun --exclusive -n 1 -c 1 "${EXE}" -i "${INPUT}" "${nsamp}" "${L}" "${n_steps}" "${n_jump}" "${nr}" "${n_rest}" &
  done

  active=$((active + 1))
  if [ "${active}" -ge "${NSAMP_PAR}" ]; then
    wait
    active=0
  fi
done
wait