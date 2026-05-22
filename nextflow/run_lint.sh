#!/usr/bin/bash
# Lint the genome_functional.nf workflow.
# Run from the project root: bash pipeline/nextflow/run_lint.sh

set -euo pipefail

NXFDIR="pipeline/nextflow"

module load nextflow 2>/dev/null || true

echo "=== Nextflow syntax check (-preview, local executor, 0 samples) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/genome_functional.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile test \
    -preview \
    2>&1

echo ""
echo "=== nf-core lint (if available) ==="
if module load nf-core 2>/dev/null; then
    nf-core lint ${NXFDIR}/genome_functional.nf \
        --fail-ignored --fail-warned \
        2>&1 || true   # nf-core lint warns about non-nf-core structure; treat as advisory
else
    echo "nf-core not available, skipping."
fi

echo ""
echo "=== Python bin/ script syntax check ==="
python3 -m py_compile ${NXFDIR}/bin/merge_cazy.py     && echo "  OK  merge_cazy.py"
python3 -m py_compile ${NXFDIR}/bin/merge_merops.py   && echo "  OK  merge_merops.py"
python3 -m py_compile ${NXFDIR}/bin/merge_signalp.py  && echo "  OK  merge_signalp.py"
python3 -m py_compile ${NXFDIR}/bin/merge_tmhmm.py    && echo "  OK  merge_tmhmm.py"
python3 -m py_compile ${NXFDIR}/bin/merge_targetp.py  && echo "  OK  merge_targetp.py"
python3 -m py_compile ${NXFDIR}/bin/merge_wolfpsort.py && echo "  OK  merge_wolfpsort.py"
python3 -m py_compile ${NXFDIR}/bin/merge_predgpi.py  && echo "  OK  merge_predgpi.py"
python3 -m py_compile ${NXFDIR}/tests/validate_outputs.py && echo "  OK  validate_outputs.py"

echo ""
echo "Lint complete."
