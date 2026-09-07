#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_clean_genome
#SBATCH --output=logs/clean_genome_launch.%j.log

# Run the funannotate pipeline ONLY through the genome-cleaning step, then stop.
# only_clean (nextflow/param_files/params_only_clean.yaml) gates every downstream
# process (masking, SRA fetch, train, predict, annotate) so just
# SETUP_TAXONDB + cleaning run. Cleaned assemblies land in
# input_clean_genomes/<asmid>.fa. For clean + tantan masking, then stop, see
# run_genome_clean_mask.sh instead.
#
# By default genomes are cleaned in BATCHES (GENOME_CLEAN_BATCH): each SLURM job
# stages the FCS-GX DB into /dev/shm once (~30 min) and then cleans up to
# clean_batch_size genomes (default 1000) sequentially. Already-cleaned .fa files
# are skipped, so a killed/retried batch resumes cheaply. Set --clean_batch_size 0
# to revert to one SLURM job per genome (GENOME_CLEAN).
#
# Submit from the PROJECT ROOT directory (where samples.csv lives):
#   sbatch nextflow/run_clean_genome.sh
#
# Common overrides (safe on the CLI -- non-boolean; booleans must go in
# nextflow/param_files/params_only_clean.yaml instead, confirmed 2026-08-26, see
# that file's header comment):
#   sbatch nextflow/run_clean_genome.sh --taxon 'PHYLUM:Basidiomycota'
#   sbatch nextflow/run_clean_genome.sh --asmid GCA_000000000.1
#   sbatch nextflow/run_clean_genome.sh --n_test 2
#   sbatch nextflow/run_clean_genome.sh --clean_batch_size 500

set -euo pipefail

module load nextflow
module load apptainer
# Container cache (env-resolved params.singularity_cache): backwards-compat
# fallback to the canonical shared dir, mirroring run_funannotate.sh/run_unified.sh.
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

mkdir -p logs/nextflow

PARAMS_FILE=nextflow/param_files/params_only_clean.yaml
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
