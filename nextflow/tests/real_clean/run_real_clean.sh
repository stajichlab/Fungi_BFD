#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 48:00:00
#SBATCH --job-name=nxf_real_clean
#SBATCH --output=logs/nextflow/real_clean_launch.%j.log

# Real-dataset smoke test for the AAFTF-SIF genome-cleaning step.
#
# Runs the funannotate pipeline ONLY through genome cleaning (--only_clean) on 3
# small real assemblies symlinked into tests/real_clean/source/ (Ascomycota:
# S. cerevisiae, S. pombe; Basidiomycota: C. neoformans). The AAFTF + taxonkit
# calls run inside the AAFTF v0.7.0-beta.2 SIF via in-script `singularity exec`.
#
# Launching from the nextflow/ dir keeps launchDir here, so outputs land in
# nextflow/input_clean_genomes/ and params.taxondb resolves to nextflow/lib/taxdump
# (reused/cached by SETUP_TAXONDB's storeDir), exactly like production.
#
# Submit with:
#   sbatch nextflow/tests/real_clean/run_real_clean.sh
#
# Rough expected timeline: SETUP_TAXONDB (real taxdump, ~70 MB download) -> the
# highmem GENOME_CLEAN_BATCH stages the 465 GB FCS-GX DB into /dev/shm (~30 min)
# then purges the 3 genomes (~15 min). Total ~45-60 min wall clock.
#
# If the test passes (real *.fa.gz in input_clean_genomes/ plus per-genome
# .purge.fasta under input_clean_genomes/clean/), the production AAFTF module load
# is replaced permanently by the SIF form. See .living/decisions.md (2026-08-10).

set -euo pipefail

module load nextflow

# Resolve the nextflow/ dir (repository tooling), wherever this launcher was submitted.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$HERE" || { echo "[ERR] cannot cd to $HERE" >&2; exit 1; }

# The sbatch job may lack the stajichlab supplementary group, so it cannot create
# new top-level dirs here even though interactive shells can. Guard the mkdir (the
# nextflow/logs dir must pre-exist -- see --output above) and proceed regardless.
mkdir -p "$HERE/logs/nextflow" 2>/dev/null || true

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run main.nf \
    -c nextflow.config \
    -c conf/test_clean_real.config \
    -profile funannotate \
    --pipeline funannotate \
    -resume
