#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 2 --mem 4gb --time 1:00:00
#SBATCH --job-name=nxf_functional_test
#SBATCH --output=logs/nextflow/functional_test.%j.log

# Stub-run test: exercises the full DAG without invoking any real bioinformatics
# tools. Each process runs its stub: block (creates minimal placeholder outputs).
#
# Submit from the project root:
#   sbatch pipeline/nextflow/run_test.sh
#
# Or run interactively:
#   bash pipeline/nextflow/run_test.sh

set -euo pipefail

module load nextflow

NXFDIR="pipeline/nextflow"
mkdir -p logs/nextflow

echo "=== Step 1: Syntax + channel wiring check (preview) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/BFD.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile test \
    -preview 2>&1 | tee logs/nextflow/functional_preview.log

echo ""
echo "=== Step 2: Full stub-run (all 9 tool subworkflows) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/BFD.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile test \
    -stub-run 2>&1 | tee logs/nextflow/functional_stubrun.log

echo ""
echo "=== Step 3: Validate stub outputs ==="
python3 ${NXFDIR}/tests/validate_outputs.py \
    --tables   ${NXFDIR}/tests/output/tables \
    --outdir   ${NXFDIR}/tests/output/function

echo ""
echo "All tests passed."
