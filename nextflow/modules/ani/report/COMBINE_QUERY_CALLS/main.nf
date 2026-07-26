// Merge every group's calls into one master orphan-classification table.
process COMBINE_QUERY_CALLS {
    label 'report'

    publishDir "${params.outdir}/skani_query/${params.compare}", mode: 'copy'

    input:
        path calls_tsvs, stageAs: 'calls_in/*'

    output:
        path("orphan_calls.csv")

    script:
    """
    python3 ${projectDir}/bin/combine_query_calls.py --calls-dir calls_in --output orphan_calls.csv
    """

    stub:
    """
    touch orphan_calls.csv
    """
}
