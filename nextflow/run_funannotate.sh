#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_funannotate
#SBATCH --output=logs/funannotate_launch.%j.log

# Launch the funannotate genome prediction/annotation pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_funannotate.sh
#
# The pipeline is dispatched by --pipeline (not -entry) and all profile-specific
# params are in conf/profile_funannotate.config.
#
# Common overrides:
#   sbatch nextflow/run_funannotate.sh --run_annotate true
#   sbatch nextflow/run_funannotate.sh --n_test 2 --only_clean true
#   sbatch nextflow/run_funannotate.sh --stop_after_sra_fetch true
#   sbatch nextflow/run_funannotate.sh --run_repeatmasker false
#   sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Ascomycota
#   sbatch nextflow/run_funannotate.sh --asmid GCA_000001405.15

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile funannotate \
    --pipeline funannotate \
    -resume \
    "$@"
