#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_seqstats
#SBATCH --output=logs/seqstats_launch.%j.log

# Launch the sequence statistics pipeline (AA freq, codon freq, gene stats, chrom info).
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_seqstats.sh
#
# To run only specific analyses:
#   sbatch nextflow/run_seqstats.sh --run_aa_freq false --run_codon_freq false
#
# To limit to first N samples for testing:
#   sbatch nextflow/run_seqstats.sh --n_test 5

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/genome_seqstats.nf \
    -c nextflow/nextflow.config \
    -profile BFD \
    -resume \
    "$@"
