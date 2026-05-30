#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_comparative
#SBATCH --output=logs/comparative_launch.%j.log

# Launch the comparative genomics clustering pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_comparative.sh --project MyProject --group group.csv
#   sbatch nextflow/run_comparative.sh --project MyProject --taxon CLASS:Dothideomycetes
#   sbatch nextflow/run_comparative.sh --project MyProject \
#       --taxon PHYLUM:Ascomycota,ORDER:Hypocreales --ignore exclude.txt
#
# Common overrides:
#   --run_orthofinder true          enable OrthoFinder sub-workflow
#   --run_mmseqs2 false             skip MMseqs2 clustering
#   --run_mcl false                 skip DIAMOND/MCL clustering
#   --mmseqs_min_id 0.50            raise MMseqs2 identity threshold
#   --mcl_inflation 2.0             tighter MCL clusters
#   --orthofinder_msa true          use MSA-based species tree in OrthoFinder
#   --n_test 10                     restrict to first N rows of samples.csv

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/comparative_genomics.nf \
    -c nextflow/nextflow.config \
    -profile comparative \
    -resume \
    "$@"
