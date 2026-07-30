process MMSEQS_CREATETSV {
    label      'comparative_mmseqs'
    tag        "createtsv"
    publishDir "${params.outdir}/${params.project}/mmseqs2", mode: 'copy'

    input:
    path db_files
    path cluster_files

    output:
    path "${params.project}.mmseqs2_clusters.tsv", emit: tsv

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mmseqs createtsv combined_db combined_db cluster_result \\
        ${params.project}.mmseqs2_clusters.tsv
    """

    stub:
    """
    touch ${params.project}.mmseqs2_clusters.tsv
    """
}
