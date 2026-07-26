process RUN_MEROPS {
    tag        "${meta.locustag}"
    label      'merops'
    storeDir   "${params.outdir}/merops"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.blasttab.gz"), emit: blasttab

    script:
    """
    module load ncbi-blast
    module load db-merops
    blastp -query ${proteins} \\
        -db \$MEROPS_DB/merops_scan.lib \\
        -out ${meta.id}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -outfmt 6 \\
        -use_sw_tback
    pigz ${meta.id}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${meta.id}.blasttab.gz
    """
}
