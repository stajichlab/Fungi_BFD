#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_earlgrey
#SBATCH --output=logs/earlgrey_launch.%j.log

# Launch the EarlGrey repeat-masking pipeline for large genomes.
# Builds a curated TE library once per species (best representative > cutoff_mb)
# and applies it with RepeatMasker to every conspecific strain, writing
# input_clean_genomes/<asmid>.masked.fasta (consumed automatically by funannotate).
#
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_earlgrey.sh
#
# Common overrides:
#   sbatch nextflow/run_earlgrey.sh --n_test 1            # one species (testing)
#   sbatch nextflow/run_earlgrey.sh --cutoff_mb 150       # lower size cutoff
#   sbatch nextflow/run_earlgrey.sh --repeat_taxon eukarya

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/earlgrey_mask.nf \
    -c nextflow/nextflow.config \
    -profile earlgrey \
    -resume \
    "$@"
