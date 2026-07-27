#!/usr/bin/bash
# Lint the migrated nextflow pipelines via main.nf dispatch.
# Run from the project root: bash nextflow/run_lint.sh

set -euo pipefail

NXFDIR="nextflow"

module load nextflow 2>/dev/null || true

echo "=== Nextflow syntax check: main.nf --pipeline BFD (-preview, test profile) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile BFD,test \
    --pipeline BFD \
    -preview \
    2>&1

echo ""
echo "=== Nextflow syntax check: main.nf --pipeline funannotate (-preview, funannotate,test) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile funannotate,test \
    --pipeline funannotate \
    -preview \
    2>&1

echo ""
echo "=== Nextflow syntax check: main.nf --pipeline compare_ani (-preview, ani,test) ==="
NXF_OPTS="-Xms256m -Xmx2g" \
nextflow run ${NXFDIR}/main.nf \
    -c ${NXFDIR}/nextflow.config \
    -profile ani,test \
    --pipeline compare_ani \
    -params-file ${NXFDIR}/params_ani.yaml \
    -preview \
    2>&1

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
