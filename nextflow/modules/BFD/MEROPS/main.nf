process RUN_MEROPS {
    tag        "${locustag}"
    label      'merops'
    storeDir   "${params.outdir}/merops"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.blasttab.gz"), emit: blasttab

    script:
    """
    module load ncbi-blast
    module load db-merops
    blastp -query ${proteins} \\
        -db \$MEROPS_DB/merops_scan.lib \\
        -out ${basename}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -outfmt 6 \\
        -use_sw_tback
    pigz ${basename}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${basename}.blasttab.gz
    """
}
