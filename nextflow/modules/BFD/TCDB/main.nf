include { hashBucketForType } from '../../common/utils.nf'

// Transporter classification: blastp each proteome against TCDB
// (Transporter Classification Database, tcdb.org), same containerized-blastp
// pattern as RUN_MEROPS (blastp_sif, no host `module load ncbi-blast`).
process RUN_TCDB {
    tag        "${meta.locustag}"
    label      'tcdb'
    storeDir   { "${params.outdir}/tcdb/${hashBucketForType('tcdb', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)
        path(tcdb_fasta)
        path(tcdb_blastdb)

    output:
        path("${meta.locustag}.blasttab.gz"), emit: blasttab

    script:
    """
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    # ${proteins} and the tcdb_fasta/tcdb_blastdb inputs (SETUP_TCDB_DB's
    # storeDir-cached output, params.tcdb_dbdir) may be symlinks Nextflow
    # staged from outside \$PWD; bind both real parent dirs too, or apptainer
    # won't see the symlink targets inside the container.
    PROT_REAL_DIR=\$(dirname \$(readlink -f ${proteins}))
    SING_BINDS="--bind \${PWD}:\${PWD},${params.tcdb_dbdir}:${params.tcdb_dbdir},\${PROT_REAL_DIR}:\${PROT_REAL_DIR},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.blastp_sif} blastp -query ${proteins} \\
        -db tcdb \\
        -out ${meta.locustag}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 5 \\
        -evalue ${params.tcdb_evalue} \\
        -outfmt 6
    pigz ${meta.locustag}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${meta.locustag}.blasttab.gz
    """
}
