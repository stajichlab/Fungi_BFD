include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

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
    def proteins_plain = plainName(proteins)
    """
    ${gzInflate(proteins)}
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    # ${proteins_plain} and the tcdb_fasta/tcdb_blastdb inputs (SETUP_TCDB_DB's
    # storeDir-cached output, params.tcdb_dbdir) may be symlinks Nextflow
    # staged from outside \$PWD; bind both real parent dirs too, or apptainer
    # won't see the symlink targets inside the container.
    # Bind the real directory of the staged BLAST db files, not params.tcdb_dbdir:
    # if tcdb_dbdir is itself a symlink, its target is not visible in the
    # container and blastp fails with "Database memory map file error".
    PROT_REAL_DIR=\$(dirname \$(readlink -f ${proteins_plain}))
    DB_REAL_DIR=\$(dirname \$(readlink -f tcdb.pin))
    SING_BINDS="--bind \${PWD}:\${PWD},\${DB_REAL_DIR}:\${DB_REAL_DIR},\${PROT_REAL_DIR}:\${PROT_REAL_DIR},\$TMPDIR:\$TMPDIR"
    # The task shell is 'bash -l' without -e, so check blastp explicitly. Without
    # this, a failed blastp leaves an empty table that pigz compresses and
    # storeDir caches as a finished genome.
    apptainer exec \${SING_BINDS} ${params.blastp_sif} blastp -query ${proteins_plain} \\
        -db tcdb \\
        -out ${meta.locustag}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 5 \\
        -evalue ${params.tcdb_evalue} \\
        -outfmt 6 \\
        || { echo "[ERROR] blastp failed for ${meta.locustag} (exit \$?)" >&2; rm -f ${meta.locustag}.blasttab; exit 1; }
    [ -f ${meta.locustag}.blasttab ] || { echo "[ERROR] blastp wrote no output for ${meta.locustag}" >&2; exit 1; }
    pigz ${meta.locustag}.blasttab || { echo "[ERROR] pigz failed for ${meta.locustag}" >&2; exit 1; }
    """

    stub:
    """
    printf '' | gzip > ${meta.locustag}.blasttab.gz
    """
}
