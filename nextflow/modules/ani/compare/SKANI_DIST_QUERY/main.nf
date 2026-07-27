process SKANI_DIST_QUERY {
    tag   "${group_name} q=${query_genomes.size()} r=${ref_genomes.size()}"
    label 'skani'

    cpus   { (ref_genomes.size() + query_genomes.size()) > 2000 ? 32 : 8 }
    memory { (ref_genomes.size() + query_genomes.size()) > 2000 ? '64 GB' : '16 GB' }

    publishDir { "${params.outdir}/skani_query/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name), path(query_genomes), path(ref_genomes)

    output:
        tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    def q_list = query_genomes.collect { it.toString() }.join('\n')
    def r_list = ref_genomes.collect   { it.toString() }.join('\n')
    """
    # skani 0.3.x: sketch query and reference genomes into separate directories.
    # skani dist -q accepts a sketch directory; -r accepts a list of fasta/sketch targets.
    printf '%s\\n' "${q_list}" > query_list.txt
    printf '%s\\n' "${r_list}" > ref_list.txt
    skani sketch -t ${task.cpus} -l query_list.txt -o query_sketch
    skani sketch -t ${task.cpus} -l ref_list.txt   -o ref_sketch

    skani dist -q query_sketch -rl ref_list.txt \\
        --min-af ${params.skani_min_af} -n ${params.query_top_n} -t ${task.cpus} \\
        -o skani_raw.tsv

    # skani dist cols: name1  name2  ANI  AF_ref  AF_query
    awk -F'\\t' 'NR>1 { print \$1"\\t"\$2"\\t"\$3 }' skani_raw.tsv > ${group_name}.query.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.query.ani.tsv
    """
}
