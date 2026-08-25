include { hashBucketForType } from '../../common/utils.nf'

process RUN_PFAM {
    tag        "${meta.locustag}"
    label      'pfam'
    storeDir   { "${params.outdir}/pfam_hmmscan/${hashBucketForType('pfam_hmmscan', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.domtblout.gz"),    emit: domtbl
        path("${meta.locustag}.tblout.gz"),  emit: tblout

    script:
    // --cpu-bind=none is required for multi-task srun on this cluster -- the
    // default cpu-bind mode conflicts with the job step's cgroup/allocation
    // mask and fails every task outright ("CPU binding outside of job step
    // allocation"), independent of --mpi=<plugin> choice. See T-015
    // (.living/learnings.md 2026-08-01) for the full diagnosis.
    def mpi_launch = params.pfam_tasks > 1 ? "srun --cpu-bind=none -N ${params.pfam_nodes} -n ${params.pfam_tasks}" : ""
    def mpi_flag   = params.pfam_tasks > 1 ? "--mpi" : ""
    // --cpu is REJECTED by hmmsearch when combined with --mpi ("Option --cpu
    // is incompatible with option(s) --mpi") -- confirmed directly, not
    // assumed. Only pass --cpu in the non-MPI (threaded) branch.
    def cpu_flag   = params.pfam_tasks > 1 ? "" : "--cpu ${task.cpus}"
    """
    module load db-pfam
    if [ ! -z "${mpi_flag}" ]; then
        module load hmmer/3.4-mpi
    else
        # Non-MPI (production default: pfam_tasks=1): hmmsearch comes from the
        # HMMER 3.4 SIF, not `module load hmmer/3.4` (mirrors the AAFTF SIF
        # pattern). The MPI branch keeps the host module: the SIF has no MPI
        # runtime and --mpi needs the cluster's OpenMPI/PMI wiring.
        #
        # db-pfam must be loaded (above) before this bind-mount is built --
        # PFAM_DB has to already be set or SING_BINDS captures an empty path
        # and the container never gets /srv/projects/db/pfam/... mounted,
        # failing every hmmsearch with "File existence/permissions problem"
        # even though the path is valid and readable on the host.
        module load apptainer
        PFAM_SIF=${params.pfam_sif}
        export TMPDIR=\${SCRATCH:-/tmp}
        SING_BINDS="--bind \${PWD}:\${PWD},\${PFAM_DB}:\${PFAM_DB},\$TMPDIR:\$TMPDIR"
        SING="apptainer exec \${SING_BINDS} \${PFAM_SIF}"
    fi
    ${mpi_launch} \${SING} hmmsearch ${mpi_flag} --cut_ga --noali ${cpu_flag} \\
        --domtbl    ${meta.locustag}.domtblout \\
        --tblout    ${meta.locustag}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${meta.locustag}.domtblout ${meta.locustag}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${meta.locustag}.domtblout.gz
    printf '' | gzip     > ${meta.locustag}.tblout.gz
    """
}
