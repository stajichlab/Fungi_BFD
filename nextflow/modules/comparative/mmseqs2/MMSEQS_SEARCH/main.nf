process MMSEQS_SEARCH {
    label    'comparative_mmseqs'
    tag      "search"
    storeDir "${params.outdir}/${params.project}/mmseqs2/search"

    input:
    path db_files

    output:
    path "search_result*", emit: result_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mkdir -p tmp
    mmseqs search combined_db combined_db search_result tmp \\
        --min-seq-id ${params.mmseqs_min_id} \\
        -c ${params.mmseqs_cov} --cov-mode 0 \\
        -s ${params.mmseqs_sensitivity} \\
        --threads ${task.cpus}
    rm -rf tmp
    """

    stub:
    """
    touch search_result search_result.index search_result.dbtype
    """
}
