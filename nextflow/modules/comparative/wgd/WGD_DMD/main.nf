// wgd dmd — identify and date WGD events by running an all-vs-all DIAMOND
// BLASTp on the translated CDS + MCL clustering to delineate paralog families.
// Containerized: runs inside the wgd SIF (container_wgd2_complete) via
// `apptainer exec`, load apptainer+squashfuse inline (beforeScript module
// loads do not reliably survive into the task script on this cluster).
//
// Exports the per-species paralog families TSV for the downstream WGD_KSD /
// WGD_SYN steps.
//
// Output layout follows the BFD convention: storeDir fan-out into a SHA-1
// hash sub-folder keyed on the stable LOCUSTAG, so no per-genome directory
// accumulates more than a few hundred files. Files keep wgd's natural names
// (they embed the source cds id), which also means a re-annotated sample
// produces a new filename and can never be served a stale storeDir hit.
include { hashBucketForType } from '../../../common/utils.nf'

process WGD_DMD {
    tag        "${meta.id}"
    label      'comparative_wgd'

    storeDir   { "${params.outdir}/wgd_dmd/${hashBucketForType('wgd_dmd', meta.locustag)}" }

    input:
    tuple val(meta), path(cds)

    output:
    tuple val(meta), path(cds), path("*.tsv"), emit: families

    script:
    """
    module load apptainer squashfuse fuse
    export TMPDIR=\${SCRATCH:-/tmp}
    PWD_REAL=\$(realpath \${PWD})
    SING_BINDS="--bind \${PWD_REAL}:\${PWD_REAL},${projectDir}:${projectDir},/bigdata:/bigdata,/rhome:/rhome,\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} \\
        --env PYTHONNOUSERSITE=1 \\
        ${params.wgd_sif} wgd dmd \\
        --outdir wgd_dmd \\
        --eval ${params.wgd_mcl_evalue} \\
        --inflation ${params.wgd_mcl_inflation} \\
        ${cds}

    [ -f wgd_dmd/*.tsv ] || { echo "ERROR: wgd dmd produced no families file" >&2; exit 1; }
    mv wgd_dmd/*.tsv .
    """

    stub:
    """
    mkdir -p wgd_dmd
    printf 'family\\tgene_a\\tgene_b\\n' > wgd_dmd/${cds}.tsv
    mv wgd_dmd/*.tsv .
    """
}
