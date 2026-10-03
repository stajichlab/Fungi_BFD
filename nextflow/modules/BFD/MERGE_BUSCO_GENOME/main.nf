include { tablesDir } from '../../common/utils.nf'

process MERGE_BUSCO_GENOME {
    label      'merge'

    input:
    path manifest

    output:
    path "busco_genome.parquet", emit: parquet

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/summarize_busco_stats.py \\
        --manifest ${manifest} \\
        -o         busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         busco_genome \\
        --local-parquet busco_genome.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_BUSCO_GENOME \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    printf 'ASMID\\tcomplete_pct\\tsingle_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\tn_markers\\tlineage\\nSTUB_ASM\\t99.0\\t98.0\\t1.0\\t0.5\\t0.5\\t758\\tfungi_odb10\\n' | gzip > busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         busco_genome \\
        --local-parquet busco_genome.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_BUSCO_GENOME \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
