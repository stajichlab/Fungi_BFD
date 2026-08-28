// wgd ksd — Ka/Ks (dS) estimation for paralog families from wgd dmd, giving
// the Ks distribution used to date whole-genome duplication events.
// Also emits the Ks histogram plot (pdf/svg) alongside the TSV.
// Containerized: same apptainer + squashfuse inline-load pattern as WGD_DMD.
//
// Requires a writable TMPDIR: wgd shells out to mafft, whose wrapper mkdirs
// in $TMPDIR; on this cluster the host $TMPDIR (/scratch) is node-local and
// NOT apptainer-mounted by default, so we bind ${SCRATCH}/$TMPDIR in.
//
// Output layout follows the BFD convention: storeDir fan-out into a SHA-1
// hash sub-folder keyed on the stable LOCUSTAG (see WGD_DMD notes). Files are
// moved out of the wgd_ksd/ subdir to the workdir root so storeDir places
// them directly in the bucket.
include { hashBucketForType } from '../../../common/utils.nf'

process WGD_KSD {
    tag        "${meta.id}"
    label      'comparative_wgd'

    storeDir   { "${params.outdir}/wgd_ksd/${hashBucketForType('wgd_ksd', meta.locustag)}" }

    input:
    tuple val(meta), path(families), path(cds)

    output:
    path("*.ks.tsv"),  emit: ks
    path("*.ksd.pdf"), emit: ks_plot

    script:
    """
    module load apptainer squashfuse fuse
    export TMPDIR=\${SCRATCH:-/tmp}
    PWD_REAL=\$(realpath \${PWD})
    SING_BINDS="--bind \${PWD_REAL}:\${PWD_REAL},${projectDir}:${projectDir},/bigdata:/bigdata,/rhome:/rhome,\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} \\
        --env PYTHONNOUSERSITE=1 \\
        ${params.wgd_sif} wgd ksd \\
        --outdir wgd_ksd \\
        --nthreads ${task.cpus} \\
        ${families} ${cds}

    [ -f wgd_ksd/*.ks.tsv ] || { echo "ERROR: wgd ksd produced no ks.tsv output" >&2; exit 1; }
    mv wgd_ksd/*.ks.tsv wgd_ksd/*.ksd.pdf .
    """

    stub:
    """
    mkdir -p wgd_ksd
    printf 'pair\\tks\\n' > wgd_ksd/${families}.ks.tsv
    touch wgd_ksd/${families}.ksd.pdf
    mv wgd_ksd/* .
    """
}
