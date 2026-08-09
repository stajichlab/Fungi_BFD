include { tablesDir } from '../../common/utils.nf'

process MERGE_BUSCO_GENOME {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "busco_genome.parquet", emit: parquet

    script:
    """
    python3 ${projectDir}/bin/summarize_busco_stats.py \\
        --manifest ${manifest} \\
        -o         busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    """

    stub:
    """
    printf 'ASMID\\tcomplete_pct\\tsingle_pct\\tduplicated_pct\\tfragmented_pct\\tmissing_pct\\tn_markers\\tlineage\\n' | gzip > busco_genome.tsv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('busco_genome.tsv.gz', delim='\\t', sample_size=-1)) TO 'busco_genome.parquet' (FORMAT PARQUET);"
    rm -f busco_genome.tsv.gz
    """
}
