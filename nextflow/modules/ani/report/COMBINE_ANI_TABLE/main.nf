// Merge every group's pair table + names lookup into one labeled, queryable
// table (CSV + SQLite) so pairwise identities can be looked up across the
// whole run without re-parsing per-group files.
process COMBINE_ANI_TABLE {
    label 'report'

    publishDir "${params.outdir}/${params.ani_method}/${params.compare}", mode: 'copy'

    input:
        path ani_tsvs,   stageAs: 'ani_in/*'
        path names_tsvs, stageAs: 'names_in/*'

    output:
        path("all_pairs.csv")
        path("ani.db")

    script:
    """
    python3 ${projectDir}/bin/combine_ani_table.py \\
        --ani-dir   ani_in \\
        --names-dir names_in \\
        --compare-level "${params.compare}" \\
        --csv-output all_pairs.csv \\
        --db-output  ani.db
    """

    stub:
    """
    touch all_pairs.csv ani.db
    """
}
