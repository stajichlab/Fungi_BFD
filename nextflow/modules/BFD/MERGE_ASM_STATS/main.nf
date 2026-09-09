include { tablesDir } from '../../common/utils.nf'

process MERGE_ASM_STATS {
    label      'merge'

    input:
    path manifest
    path samples

    output:
    path "asm_stats.parquet", emit: parquet

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/summarize_asm_stats.py \\
        --manifest ${manifest} \\
        --samples  ${samples} \\
        -o         asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         asm_stats \\
        --local-parquet asm_stats.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_ASM_STATS \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    printf 'ASMID\\tSPECIES\\tSTRAIN\\tcontig_count\\ttotal_length_bp\\tmin_contig_bp\\tmax_contig_bp\\tmedian_contig_bp\\tmean_contig_bp\\tL50\\tN50_bp\\tL90\\tN90_bp\\tgc_pct\\tn_gap_count\\ttotal_n_bases\\tmasked_bases\\tmasked_pct\\tt2t_scaffolds\\ttelomere_fwd\\ttelomere_rev\\nSTUB_ASM\\tFoo bar\\tCBS 1\\t10\\t1000000\\t500\\t200000\\t50000\\t100000\\t3\\t150000\\t8\\t50000\\t48.5\\t0\\t0\\t10000\\t1.0\\t2\\t2\\t2\\n' | gzip > asm_stats.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('asm_stats.tsv.gz', delim='\\t', sample_size=-1)) TO 'asm_stats.parquet' (FORMAT PARQUET);"
    rm -f asm_stats.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         asm_stats \\
        --local-parquet asm_stats.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_ASM_STATS \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
