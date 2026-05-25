#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_funannotate
#SBATCH --output=logs/funannotate_launch.%j.log

# Launch the funannotate genome prediction/annotation pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_funannotate.sh
#
# Common overrides:
#   sbatch nextflow/run_funannotate.sh --run_annotate true
#   sbatch nextflow/run_funannotate.sh --n_test 2 --only_clean true
#   sbatch nextflow/run_funannotate.sh --stop_after_sra_fetch true
#   sbatch nextflow/run_funannotate.sh --run_repeatmasker false

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/funannotate.nf \
    -c nextflow/nextflow.config \
    -profile funannotate \
    -resume \
    "$@"
