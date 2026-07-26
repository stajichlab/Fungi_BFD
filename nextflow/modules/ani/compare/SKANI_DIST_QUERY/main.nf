// Asymmetric query-vs-reference ANI: O(queries x references), not O(refs^2).
process SKANI_DIST_QUERY {
    tag   "${group_name} q=${query_sketches.size()} r=${ref_sketches.size()}"
    label 'skani'

    cpus   { ref_sketches.size() > 2000 ? 32 : ref_sketches.size() > 500 ? 16 : 8 }
    memory { ref_sketches.size() > 2000 ? '64 GB' : ref_sketches.size() > 500 ? '32 GB' : '16 GB' }

    storeDir "${params.outdir}/skani_query/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(query_sketches, stageAs: 'query/*'), path(ref_sketches, stageAs: 'ref/*')

    output:
        tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    """
    ls query/*.sketch > query_list.txt
    ls ref/*.sketch   > ref_list.txt
    skani dist --ql query_list.txt --rl ref_list.txt \\
        --min-af ${params.skani_min_af} -n ${params.query_top_n} -t ${task.cpus} \\
        -o skani_raw.tsv

    # skani cols: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query [names]
    awk -F'\\t' 'NR>1 {
        nr=split(\$1,b,"/"); r=b[nr]; sub(/\\.sketch\$/,"",r);
        nq=split(\$2,a,"/"); q=a[nq]; sub(/\\.sketch\$/,"",q);
        print q"\\t"r"\\t"\$3
    }' skani_raw.tsv > ${group_name}.query.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.query.ani.tsv
    """
}
