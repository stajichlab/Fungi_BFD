include { hashBucketForType } from '../../common/utils.nf'

process RUN_SWISSPROT {
    tag        "${meta.locustag}"
    label      'swissprot'
    storeDir   { "${params.outdir}/swissprot/${hashBucketForType('swissprot', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.blasttab.gz"), emit: blasttab

    // SwissProt homology search. Two interchangeable engines, both run from
    // Singularity containers (never the host module, so versions are pinned):
    //   diamond (default)  — diamond blastp --sensitive (quay.io/biocontainers
    //                        diamond 2.2.5, newer than the host `diamond`
    //                        2.1.24 module)
    //   blastp             — NCBI blastp (blast 2.16.0 container, same as
    //                        RUN_MEROPS)
    // Both emit the same 18-column blasttab, so MERGE_SWISSPROT is engine-
    // agnostic. E-value/counts are deliberately permissive (raw evidence is
    // kept); the 80-80 functional-transfer flag is derived in the merge, not
    // pre-filtered here.
    script:
    """
    module load singularity
    SWISSPROT_DB=${params.swissprot_dbdir}
    SING_BINDS="--bind \${SWISSPROT_DB}:\${SWISSPROT_DB}"
    if [ "${params.swissprot_search}" = "blastp" ]; then
        SING="singularity exec \${SING_BINDS} ${params.blastp_sif}"
        \${SING} blastp -query ${proteins} \\
            -db \${SWISSPROT_DB}/uniprot_sprot \\
            -out ${meta.locustag}.blasttab \\
            -num_threads ${task.cpus} \\
            -seg yes -soft_masking true \\
            -max_target_seqs ${params.swissprot_max_targets} \\
            -evalue ${params.swissprot_evalue} \\
            -outfmt "6 qseqid sseqid pident positive nident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen stitle"
    else
        SING="singularity exec \${SING_BINDS} ${params.diamond_sif}"
        \${SING} diamond blastp --sensitive \\
            -d \${SWISSPROT_DB}/uniprot_sprot.dmnd \\
            -q ${proteins} \\
            -o ${meta.locustag}.blasttab \\
            --threads ${task.cpus} \\
            -k ${params.swissprot_max_targets} \\
            -e ${params.swissprot_evalue} \\
            -f 6 qseqid sseqid pident positive nident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen stitle
    fi
    pigz ${meta.locustag}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${meta.locustag}.blasttab.gz
    """
}
