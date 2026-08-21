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
    module load singularity
    singularity exec ${params.blastp_sif} blastp -query ${proteins} \\
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
