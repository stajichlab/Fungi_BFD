#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 2 --mem 4gb --time 1:00:00
#SBATCH --job-name=nxf_functional_test
#SBATCH --output=logs/nextflow/functional_test.%j.log

# Stub-run test: exercises the full DAG without invoking any real bioinformatics
# tools. Each process runs its stub: block (creates minimal placeholder outputs).
#
# Submit from the project root:
#   sbatch nextflow/run_test.sh
#
# Or run interactively:
#   bash nextflow/run_test.sh

set -euo pipefail

module load nextflow

NXFDIR="nextflow"
mkdir -p logs/nextflow

# Container cache (env-resolved params.singularity_cache): backwards-compat
# fallback to the canonical shared dir, mirroring run_funannotate.sh/run_unified.sh.
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-/bigdata/stajichlab/shared/singularity_cache}

echo "=== Step 1: Syntax + channel wiring check (preview) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile test \
    --pipeline BFD \
    -preview 2>&1 | tee logs/nextflow/functional_preview.log

echo ""
echo "=== Step 2: Full stub-run (all 9 tool subworkflows) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile test \
    --pipeline BFD \
    -stub-run 2>&1 | tee logs/nextflow/functional_stubrun.log

echo ""
echo "=== Step 3: Validate stub outputs ==="
python3 ${NXFDIR}/tests/validate_outputs.py \
    --tables   ${NXFDIR}/tests/output/tables \
    --outdir   ${NXFDIR}/tests/output/function

echo ""
echo "=== Step 4: compare_ani syntax + wiring check (preview) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile ani,test,test_ani \
    --pipeline compare_ani --compare SPECIES \
    -preview 2>&1 | tee logs/nextflow/ani_preview.log

echo ""
echo "=== Step 5: compare_ani stub-run (incl. PICK_REPRESENTATIVE_STRAIN) ==="
# --run_ani_reuse true exercises CONCAT_ANI_TSVS + PICK_REPRESENTATIVE_STRAIN,
# reading conf/test.config's BUSCO_genome stub output (genome_stats_outdir) --
# the same genome_stats -> ani dependency fixed in run_funannotate_neurospora.sh.
# See conf/test_ani.config for why output paths must be overridden here: an
# earlier version of this stub-run, without those overrides, silently
# overwrote a real project's abinitio_reuse_assignments.csv via the default
# ${launchDir}/genome_annotation publishDir target.
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile ani,test,test_ani \
    --pipeline compare_ani --compare SPECIES --run_ani_reuse true \
    -stub-run 2>&1 | tee logs/nextflow/ani_stubrun.log

echo ""
echo "=== Step 6: funannotate syntax + wiring check (preview) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile funannotate,test,test_funannotate \
    --pipeline funannotate \
    -preview 2>&1 | tee logs/nextflow/funannotate_preview.log

echo ""
echo "=== Step 7: funannotate stub-run (smoke test) ==="
# KNOWN GAPS (see conf/test_funannotate.config header for detail, tracked in
# issue #3): run_sra_fetch/run_ani_reuse are forced off here, so this does NOT
# exercise FUNANNOTATE_RNASEQ's abinitioReuseMap-driven representative pick,
# only its "no assignment data" fallback path. Separately, GENOME_CLEAN_BATCH's
# stub output currently isn't picked up by FUNANNOTATE_GENOME_PREP's
# downstream channel check in a fresh launchDir (looks like the same class of
# channel-construction-time-vs-task-completion race documented elsewhere in
# .living/learnings.md) -- so this stub-run only confirms the pipeline compiles
# and starts cleanly without touching real project paths, not a full DAG walk.
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile funannotate,test,test_funannotate \
    --pipeline funannotate \
    -stub-run 2>&1 | tee logs/nextflow/funannotate_stubrun.log

echo ""
echo "=== Step 8: --asmid filtering (compare_ani) ==="
# Locks in the compare_ANI.nf/BFD.nf --asmid fix (issue #3): test_samples.csv
# has two single-strain species (GCF_TEST001, GCF_TEST002); restricting to one
# ASMID must leave exactly that one row in predict_input_for_ani.tsv, proving
# ANI/BUSCO representative-picking no longer runs against the full unfiltered
# sample set when --asmid is set.
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile ani,test,test_ani \
    --pipeline compare_ani --compare SPECIES --run_ani_reuse true \
    --asmid GCF_TEST001 \
    -stub-run 2>&1 | tee logs/nextflow/asmid_stubrun.log

ASMID_TSV="${NXFDIR}/tests/output/ani_genome_annotation/_reuse_assignments/predict_input_for_ani.tsv"
n_rows=$(tail -n +2 "${ASMID_TSV}" | grep -c . || true)
if [ "${n_rows}" -ne 1 ]; then
    echo "FAIL: expected 1 row in ${ASMID_TSV} with --asmid GCF_TEST001, got ${n_rows}" >&2
    exit 1
fi
if ! grep -q '^GCF_TEST001' "${ASMID_TSV}"; then
    echo "FAIL: ${ASMID_TSV} does not contain the requested ASMID GCF_TEST001" >&2
    exit 1
fi
if grep -q 'GCF_TEST002' "${ASMID_TSV}"; then
    echo "FAIL: ${ASMID_TSV} contains GCF_TEST002 -- --asmid filtering did not restrict the sample set" >&2
    exit 1
fi
echo "OK: --asmid GCF_TEST001 correctly restricted to 1 row"

echo ""
echo "All tests passed."
