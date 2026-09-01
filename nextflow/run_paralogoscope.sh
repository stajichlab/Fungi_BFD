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
#   --wgd_sif /path/to/other.sif          override image (default = cached SIF,
#                                         hyphaltip_wgd2_complete-2.0.38)
#   --n_test 10                          restrict to first 10 sample-sheet rows

set -euo pipefail

module load nextflow

# Container cache (env-resolved params.singularity_cache): backwards-compat
# default shared across this repo's launchers.
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export NXF_APPTAINER_CACHEDIR=${NXF_APPTAINER_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR}}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-$NXF_SINGULARITY_CACHEDIR}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$NXF_APPTAINER_CACHEDIR}

mkdir -p logs/nextflow

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile paralogoscope \
    --pipeline paralogoscope \
    -resume \
    "$@"
