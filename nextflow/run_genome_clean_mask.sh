#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00 -p batch
#SBATCH --job-name=nxf_genome_clean_mask
#SBATCH --output=logs/genome_clean_mask_launch.%j.log

# Run the funannotate pipeline through genome cleaning + tantan soft-masking,
# then stop -- no RNA-seq acquisition, training, prediction, or annotation.
# Params in nextflow/param_files/params_genome_clean_mask.yaml
# (stop_after_genome_prep: true, run_repeatmasker: true).
#
# Differs from run_clean_genome.sh's --only_clean: that stops BEFORE masking
# (unmasked genomes only). This stops AFTER GENOME_CLEAN(_BATCH) +
# MASKREPEAT_TANTAN_RUN both complete -- cleaned + masked genomes, nothing
# downstream queued.
#
# Cleaned assemblies land in input_clean_genomes/<asmid>.fa.gz; masked in
# input_clean_genomes/<asmid>.masked.fasta.
#
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_genome_clean_mask.sh
#
# Common overrides (non-boolean; safe to pass on the CLI -- see the params
# file's header comment for why boolean overrides must go in a -params-file
# instead of `--flag true` on the CLI with this repo's installed
# nextflow/nf-schema version):
#   sbatch nextflow/run_genome_clean_mask.sh --taxon 'PHYLUM:Basidiomycota'
#   sbatch nextflow/run_genome_clean_mask.sh --asmid GCA_000000000.1
#   sbatch nextflow/run_genome_clean_mask.sh --n_test 2
#   sbatch nextflow/run_genome_clean_mask.sh --clean_batch_size 500

set -euo pipefail

module load nextflow

# Container cache (env-resolved params.singularity_cache): backwards-compat
# fallback to the canonical shared dir, mirroring run_funannotate.sh/run_unified.sh.
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

mkdir -p logs/nextflow

# Literal project-root-relative path (not $(dirname "$0")/BASH_SOURCE -- those
# resolve incorrectly once this runs inside a SLURM job's environment; see
# ~/.claude/CLAUDE.md). Matches run_funannotate.sh's PARAMS_FILE convention:
# submit from the project root, where samples.csv lives.
PARAMS_FILE=nextflow/param_files/params_genome_clean_mask.yaml
if [ ! -f "$PARAMS_FILE" ]; then
    echo "ERROR: $PARAMS_FILE not found in $(pwd) -- run from the project root (where samples.csv lives)" >&2
    exit 1
fi

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile funannotate \
    -params-file "$PARAMS_FILE" \
    -resume \
    "$@"
