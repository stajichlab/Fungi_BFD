#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_phyling
#SBATCH --output=logs/phyling_launch.%j.log

# Launch the PHYling phylogenomics pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv / phylo.csv lives):
#   sbatch nextflow/run_phyling.sh
#
# Common overrides (append after the script name):
#   --taxon PHYLUM:Ascomycota          restrict to Ascomycota
#   --taxon CLASS:Sordariomycetes      restrict to a class
#   --seq_type cds                     use CDS instead of protein input
#   --markerset fungi_odb12            single markerset (default)
#   --markerset "fungi_odb12,ascomycota_odb12"   multiple markersets (quoted)
#   --tree_method ft                   FastTree (fast); iqtree (default); raxml
#   --top_n 50                         top-N markers by treeness/RCV (default 50)
#   --n_test 10                        limit to first 10 taxa (pilot run)
#   --samples phylo.csv                explicit sample sheet (default: phylo.csv or samples.csv)
#
# Examples:
#   sbatch nextflow/run_phyling.sh --taxon PHYLUM:Ascomycota
#   sbatch nextflow/run_phyling.sh --taxon CLASS:Dothideomycetes --seq_type cds
#   sbatch nextflow/run_phyling.sh --markerset "fungi_odb12,ascomycota_odb12" --tree_method iqtree
#   sbatch nextflow/run_phyling.sh --n_test 20 --tree_method ft   # fast pilot

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/phyling.nf \
    -c nextflow/nextflow.config \
    -profile phyling \
    -resume \
    "$@"
