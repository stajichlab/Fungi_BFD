// Concatenate multiple per-group ANI TSVs into one merged TSV.
// Each input TSV: query<TAB>ref<TAB>ANI[<TAB>AF_ref<TAB>AF_query]
// Output: same 3-column format (query, ref, ANI) with header row.
process CONCAT_ANI_TSVS {
    tag   "CONCAT_ANI_TSVS"
    label 'report'
    publishDir "${params.outdir}/${params.ani_method}/${params.compare}", mode: 'copy', overwrite: true

    input:
        path manifest  // file listing .ani.tsv paths, one per line

    output:
        path("all_pairs_merged.tsv"), emit: out

    script:
    """
    # Write header, then all data rows from every file.
    printf 'query\\tref\\tANI\\n' > all_pairs_merged.tsv
    xargs -a "${manifest}" -I{} sh -c '[ -s "{}" ] && cat "{}"' >> all_pairs_merged.tsv
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > all_pairs_merged.tsv
    """
}
