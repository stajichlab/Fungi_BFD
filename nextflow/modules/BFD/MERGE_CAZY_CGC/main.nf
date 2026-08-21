include { tablesDir } from '../../common/utils.nf'

process MERGE_CAZY_CGC {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(cgc_tsvs)

    output:
        path("cazy_cgc.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_cazy_cgc.py -o cazy_cgc.csv ${cgc_tsvs}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy_cgc.csv', sample_size=-1)) TO 'cazy_cgc.parquet' (FORMAT PARQUET);"
    rm -f cazy_cgc.csv
    """

    stub:
    """
    printf 'species_prefix,cgc_id,contig,gene_ids\\n' > cazy_cgc.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy_cgc.csv', sample_size=-1)) TO 'cazy_cgc.parquet' (FORMAT PARQUET);"
    rm -f cazy_cgc.csv
    """
}
