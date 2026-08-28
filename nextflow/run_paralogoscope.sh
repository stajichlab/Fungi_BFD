#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_paralogoscope
#SBATCH --output=logs/paralogoscope_launch.%j.log

# Launch the paralogoscope (wgd duplication dating) pipeline.
# Submit from the PROJECT ROOT directory (where samples.csv + input/ live):
#   sbatch nextflow/run_paralogoscope.sh --taxon CLASS:Sordariomycetes
#   sbatch nextflow/run_paralogoscope.sh --group group.csv
#
# Common overrides (append after the script name):
#   --taxon PHYLUM:Ascomycota            restrict to a phylum
#   --group group.csv                    explicit LOCUSTAG,GROUP selection
#   --ignore exclude.txt                 exclude these LOCUSTAGs
#   --run_wgd_syn true                   enable wgd syn (i-ADHoRe synteny; heavy)
#   --wgd_sif /path/to/wgd-2.0.38.sif    use a local SIF instead of ghcr docker URI
#   --n_test 10                          restrict to first 10 sample-sheet rows

set -euo pipefail

module load nextflow

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile paralogoscope \
    --pipeline paralogoscope \
    -resume \
    "$@"
