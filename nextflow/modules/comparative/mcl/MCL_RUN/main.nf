process MCL_RUN {
    label      'comparative_mcl'
    tag        "mcl"
    publishDir "${params.outdir}/${params.project}/mcl", mode: 'copy'

    input:
    path abc

    output:
    path "${params.project}.mcl_clusters.tsv", emit: tsv

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mcl
    mcl ${abc} --abc \\
        -I ${params.mcl_inflation} \\
        -te ${task.cpus} \\
        -o ${params.project}.mcl_clusters.tsv
    """

    stub:
    """
    touch ${params.project}.mcl_clusters.tsv
    """
}
