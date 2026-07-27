process MMSEQS_CLUST {
    label    'comparative_mmseqs'
    tag      "clust"

    input:
    path db_files
    path result_files

    output:
    path "cluster_result*", emit: cluster_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mmseqs clust combined_db search_result cluster_result \\
        --cluster-mode ${params.mmseqs_cluster_mode} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch cluster_result cluster_result.index cluster_result.dbtype
    """
}
