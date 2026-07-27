// skani 0.3.x: merge per-chunk databases into one group-level database,
// then run skani dist (all-vs-all within the group).
process SKANI_COMPARE {
    tag   "${group_name} n=${sketch_dbs.size()} chunks"
    label 'skani'

    cpus   { sketch_dbs.size() > 500 ? 32 : sketch_dbs.size() > 200 ? 16 : 8 }
    memory { sketch_dbs.size() > 500 ? '64 GB' : sketch_dbs.size() > 200 ? '32 GB' : '16 GB' }

    publishDir { "${params.outdir}/${params.ani_method}/${params.compare}/${group_name}/batches" }, mode: 'copy'

    input:
        tuple val(group_name), path(sketch_dbs)

    output:
        tuple val(group_name), path("${group_name}.full.ani.tsv")

    script:
    """
    # skani 0.3.x consolidates all sketches into one database via --merge.
    # If only one chunk exists, no merge needed (cp the single db through).
    if [ "${sketch_dbs.size()}" -eq 1 ]; then
        cp ${sketch_dbs[0]} group_db.sketches.db
        [ -f "${sketch_dbs[0]}.index"   ] && cp "${sketch_dbs[0]}.index"   group_db.sketches.db.index
        [ -f "${sketch_dbs[0]}.markers" ] && cp "${sketch_dbs[0]}.markers" group_db.sketches.db.markers
    else
        # --merge takes multiple -d flags; collect all chunk .sketches.db files.
        DFLAGS=\$(for d in ${sketch_dbs.join(' ')}; do echo -n "-d \$d "; done)
        skani sketch --merge \$DFLAGS -o group_db
    fi

    skani dist -d group_db --sparse \\
        --min-af ${params.skani_min_af} -t ${task.cpus} \\
        -o skani_raw.tsv

    # skani dist cols: name1  name2  ANI  AF_ref  AF_query
    awk -F'\\t' 'NR>1 && \$1 != \$2 {
        print \$1"\\t"\$2"\\t"\$3
    }' skani_raw.tsv > ${group_name}.full.ani.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > ${group_name}.full.ani.tsv
    """
}
