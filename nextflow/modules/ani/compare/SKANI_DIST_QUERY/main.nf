// Asymmetric query-vs-reference ANI: O(queries x references), not O(refs^2).
// skani 0.3.x uses consolidated sketch databases (sketches.db) instead of
// individual .sketch files. Query and reference chunk DBs are built at sketch
// time, then merged per (group, role) before dist comparison.
process SKANI_DIST_QUERY {
    tag   "${group_name} q=${query_dbs.size()} r=${ref_dbs.size()}"
    label 'skani'

    cpus   { (ref_dbs.size() + query_dbs.size()) > 2000 ? 32 : 8 }
    memory { (ref_dbs.size() + query_dbs.size()) > 2000 ? '64 GB' : '16 GB' }

    publishDir { "${params.outdir}/skani_query/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name), path(query_dbs), path(ref_dbs)

    output:
        tuple val(group_name), path("${group_name}.query.ani.tsv")

    script:
    """
    # Merge query chunks into one database
    if [ "${query_dbs.size()}" -eq 1 ]; then
        cp ${query_dbs[0]} query_db.sketches.db
        [ -f "${query_dbs[0]}.index"   ] && cp "${query_dbs[0]}.index"   query_db.sketches.db.index
        [ -f "${query_dbs[0]}.markers" ] && cp "${query_dbs[0]}.markers" query_db.sketches.db.markers
    else
        QD=\$(for d in ${query_dbs.join(' ')}; do echo -n "-d \$d "; done)
        skani sketch --merge \$QD -o query_db
    fi

    # Merge reference chunks into one database
    if [ "${ref_dbs.size()}" -eq 1 ]; then
        cp ${ref_dbs[0]} ref_db.sketches.db
        [ -f "${ref_dbs[0]}.index"   ] && cp "${ref_dbs[0]}.index"   ref_db.sketches.db.index
        [ -f "${ref_dbs[0]}.markers" ] && cp "${ref_dbs[0]}.markers" ref_db.sketches.db.markers
    else
        RD=\$(for d in ${ref_dbs.join(' ')}; do echo -n "-d \$d "; done)
        skani sketch --merge \$RD -o ref_db
    fi

    skani dist -d ref_db -d query_db \\
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
