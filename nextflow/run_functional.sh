#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_functional
#SBATCH --output=logs/nextflow/functional_launch.%j.log

# Launch the Nextflow functional annotation pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch pipeline/nextflow/run_functional.sh
#
# To run only specific tools:
#   sbatch pipeline/nextflow/run_functional.sh --run_pfam true --run_cazy false
#
# To limit to first N samples for testing:
#   sbatch pipeline/nextflow/run_functional.sh --n_test 5

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run pipeline/nextflow/genome_functional.nf \
    -c pipeline/nextflow/nextflow.config \
    -resume \
    "$@"
