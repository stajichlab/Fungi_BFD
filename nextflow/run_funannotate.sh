#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00 -p batch
#SBATCH --job-name=nxf_funannotate
#SBATCH --output=logs/funannotate_launch.%j.log

# Launch the funannotate genome prediction/annotation pipeline, toggling
# between the representatives-only pass and the full (all-strains) pass.
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#
#   sbatch nextflow/run_funannotate.sh representatives
#   sbatch nextflow/run_funannotate.sh all
#
# representatives -> params_predict_representatives.yaml
#   Predicts ONLY the ANI-selected representative strain per species and
#   backfills genome_annotation/_reuse_assignments/ with shared ab-initio
#   params. Run this first.
# all -> params_predict_all.yaml
#   Predicts every strain; non-representative siblings are gated on their
#   species' shared params being available (pre-existing, or backfilled by
#   an earlier representatives pass in this same run). Run this second.
#
# Any extra arguments are passed straight through to `nextflow run` as
# param overrides, e.g.:
#
#   sbatch nextflow/run_funannotate.sh all --run_annotate true
#   sbatch nextflow/run_funannotate.sh representatives --n_test 2 --only_clean true
#   sbatch nextflow/run_funannotate.sh all --stop_after_sra_fetch true
#   sbatch nextflow/run_funannotate.sh all --run_repeatmasker false
#   sbatch nextflow/run_funannotate.sh all --taxon PHYLUM:Ascomycota
#   sbatch nextflow/run_funannotate.sh representatives --asmid GCA_000001405.15
#
# This supersedes the older single-mode run_funannotate_represenatives.sh /
# run_funannotate_all.sh (kept as thin wrappers around this script).

set -euo pipefail

MODE="${1:-}"
case "$MODE" in
    representatives)
        PARAMS_FILE=params_predict_representatives.yaml
        ;;
    all)
        PARAMS_FILE=params_predict_all.yaml
        ;;
    *)
        echo "Usage: sbatch nextflow/run_funannotate.sh <representatives|all> [extra nextflow args...]" >&2
        echo "  representatives -> params_predict_representatives.yaml (run first)" >&2
        echo "  all             -> params_predict_all.yaml (run second)" >&2
        exit 1
        ;;
esac
shift

if [ ! -f "$PARAMS_FILE" ]; then
    echo "ERROR: $PARAMS_FILE not found in $(pwd) -- run from the project root (where samples.csv lives)" >&2
    exit 1
fi

module load nextflow
module load singularity
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
# Modules invoke `singularity exec docker://...` directly, bypassing Nextflow's
# own cacheDir handling -- SINGULARITY_CACHEDIR/APPTAINER_CACHEDIR are what the
# singularity/apptainer binaries themselves read for pull/convert caching.
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

mkdir -p logs/nextflow
# Singularity bind sources (work/ANI/pip_cache, work/ANI/python_packages) are
# created at config-parse time by the funannotate profile (see _bind_sources in
# conf/profile_funannotate.config) — no launcher-script mkdir needed here.

echo "[INFO] funannotate mode=$MODE params-file=$PARAMS_FILE" >&2

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile funannotate \
    --pipeline funannotate \
    -params-file "$PARAMS_FILE" \
    -resume \
    "$@"
