include { hashBucketForType } from '../../common/utils.nf'

process RUN_MEROPS {
    tag        "${meta.locustag}"
    label      'merops'
    storeDir   { "${params.outdir}/merops/${hashBucketForType('merops', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.blasttab.gz"), emit: blasttab

    script:
    // Containerized blastp replaces `module load ncbi-blast` so the version is
    // pinned (blast 2.16.0 SIF, same image RUN_SWISSPROT uses for its
    // alternative engine). The MEROPS DB path still comes from the db-merops
    // module env; the DB directory is bound so the container can read it.
    """
    module load db-merops/124
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    # ${proteins} may be a symlink Nextflow staged from outside \$PWD (e.g. a
    # shared genome_annotation dir elsewhere under /bigdata); bind its real
    # parent dir too, or apptainer won't see the symlink target inside the
    # container.
    PROT_REAL_DIR=\$(dirname \$(readlink -f ${proteins}))
    SING_BINDS="--bind \${PWD}:\${PWD},\${MEROPS_DB}:\${MEROPS_DB},\${PROT_REAL_DIR}:\${PROT_REAL_DIR},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.blastp_sif}"
    \${SING} blastp -query ${proteins} \\
        -db \$MEROPS_DB/merops_scan.lib \\
        -out ${meta.locustag}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -outfmt 6 \\
        -use_sw_tback
    pigz ${meta.locustag}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${meta.locustag}.blasttab.gz
    """
}
