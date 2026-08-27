#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00 -p batch
#SBATCH --job-name=nxf_functional
#SBATCH --output=logs/functional_launch.%j.log

# Launch the Nextflow functional annotation pipeline (BFD).
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_functional.sh
#
# Tool selection lives in nextflow/params/params_functional.yaml (booleans
# must go in a -params-file, not CLI --flag overrides, on this repo's
# installed nextflow/nf-schema version -- confirmed 2026-08-26, see that
# file's header comment). Edit the yaml to change which tools run.
#
# To limit to first N samples for testing (safe on the CLI, not a boolean):
#   sbatch nextflow/run_functional.sh --n_test 5

set -euo pipefail

module load nextflow

# Container cache (env-resolved params.singularity_cache): backwards-compat
# fallback to the canonical shared dir, mirroring run_funannotate.sh/run_unified.sh.
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

mkdir -p logs/nextflow

PARAMS_FILE=nextflow/params/params_functional.yaml
if [ ! -f "$PARAMS_FILE" ]; then
    echo "ERROR: $PARAMS_FILE not found in $(pwd) -- run from the project root (where samples.csv lives)" >&2
    exit 1
fi

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile BFD \
    -params-file "$PARAMS_FILE" \
    -resume \
    "$@"
