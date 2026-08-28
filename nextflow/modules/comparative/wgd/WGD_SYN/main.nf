// wgd syn — synteny-based inference of collinear blocks / anchor points
// (i-ADHoRe) for paralog families. Heavier step, gated behind --run_wgd_syn.
// Takes the dmd families TSV plus the species GFF3 (gene coordinates).
// Containerized: same apptainer + squashfuse inline-load pattern as WGD_DMD.
//
// Output layout follows the BFD convention: storeDir fan-out into a SHA-1
// hash sub-folder keyed on the stable LOCUSTAG (see WGD_DMD notes). Files are
// moved out of the wgd_syn/ subdir to the workdir root so storeDir places
// them directly in the bucket.
include { hashBucketForType } from '../../../common/utils.nf'

process WGD_SYN {
    tag        "${meta.id}"
    label      'comparative_wgd'

    storeDir   { "${params.outdir}/wgd_syn/${hashBucketForType('wgd_syn', meta.locustag)}" }

    input:
    tuple val(meta), path(families), path(gff)

    output:
    path("anchors.csv"), optional: true, emit: anchors
    path("*.dot.pdf"), optional: true, emit: dotplot

    script:
    """
    module load apptainer squashfuse fuse
    export TMPDIR=\${SCRATCH:-/tmp}
    PWD_REAL=\$(realpath \${PWD})
    SING_BINDS="--bind \${PWD_REAL}:\${PWD_REAL},${projectDir}:${projectDir},/bigdata:/bigdata,/rhome:/rhome,\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} \\
        --env TMPDIR=\$TMPDIR \\
        --env PYTHONNOUSERSITE=1 \\
        --env LD_LIBRARY_PATH=/opt/conda/lib \\
        --env OMPI_ALLOW_RUN_AS_ROOT=1,OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \\
        ${params.wgd_sif} bash -c \\
        'env -u SLURM_JOB_ID -u SLURM_JOBID -u SLURM_JOB_NODELIST -u SLURM_JOB_NUM_NODES -u SLURM_JOB_PARTITION -u SLURM_JOB_USER -u SLURM_JOB_UID -u SLURM_NODELIST -u SLURM_NODEID -u SLURM_NTASKS -u SLURM_PROCID -u SLURM_STEP_NUM_TASKS -u SLURM_TASKS_PER_NODE -u SLURM_TASK_PID -u SLURM_CPUS_ON_NODE -u SLURM_CPU_BIND -u SLURM_LOCALID -u SLURM_SUBMIT_DIR -u SLURM_MEM_PER_NODE -u SLURM_MEM_PER_CPU -u SLURM_CLUSTER_NAME -u SLURM_CONF -u SLURM_WORKING_DIR -u PMI_FD -u PMI_RANK -u PMI_SIZE -u PMI_JOBID -u PMI_ID -u PMI_PORT -u PMI_DEBUG bash -c "wgd syn --outdir wgd_syn -f mRNA -a ID \\${families} \\${gff}"'

    [ -f wgd_syn/anchors.csv ] || [ -s wgd_syn/iadhore-out/anchorpoints.txt ] \
        || { echo "ERROR: wgd syn produced no i-ADHoRe output" >&2; exit 1; }
    mv wgd_syn/anchors.csv wgd_syn/*.dot.pdf . 2>/dev/null || true
    """

    stub:
    """
    mkdir -p wgd_syn
    touch wgd_syn/anchors.csv
    mv wgd_syn/anchors.csv .
    """
}
