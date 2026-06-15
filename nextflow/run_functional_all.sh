#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00 -p batch
#SBATCH --job-name=nxf_functional
#SBATCH --output=logs/functional_launch.%j.log

# Launch the Nextflow functional annotation pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_functional.sh
#
# To run only specific tools:
#   sbatch nextflow/run_functional.sh --run_pfam true --run_cazy false
#
# To limit to first N samples for testing:
#   sbatch nextflow/run_functional.sh --n_test 5
#
# To skip symlink setup (input/ already populated from a prior run):
#   sbatch nextflow/run_functional.sh --run_setup false

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/BFD.nf \
    -c nextflow/nextflow.config \
    -profile BFD \
    --run_setup true \
    --run_pfam true \
    --run_cazy true \
    --run_merops true \
    --run_signalp true \
    --run_tmhmm true \
    --run_targetp true \
    --run_idp true \
    --run_wolfpsort true \
    --run_predgpi true \
    --run_busco_genome true \
    --run_busco_pep true \
    -resume \
    "$@"
