#!/usr/bin/env bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00
#SBATCH --job-name compare_ANI
#SBATCH -o logs/slurm/ANI_%j.out
#SBATCH -e logs/slurm/ANI_%j.err

set -euo pipefail

mkdir -p logs/slurm logs/nextflow
# Singularity bind sources (work/ANI/pip_cache, work/ANI/python_packages) are
# created at config-parse time by the ani profile (see _bind_sources in
# conf/profile_ANI.config) — no launcher-script mkdir needed here.

source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow
module load singularity

cd "$(dirname "$(readlink -f "$0")")/.."

nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile ani \
    --pipeline compare_ani \
    -params-file nextflow/param_files/params_ani.yaml \
    -resume \
    "$@"
