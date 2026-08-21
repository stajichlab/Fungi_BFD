include { tablesDir } from '../../common/utils.nf'

process MERGE_TRNASCAN {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("trnascan.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_trnascan.py -o trnascan.csv ${gff3s}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('trnascan.csv', sample_size=-1)) TO 'trnascan.parquet' (FORMAT PARQUET);"
    rm -f trnascan.csv
    """

    stub:
    """
    printf 'species_prefix,contig,start,end,strand,product,note\\n' > trnascan.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('trnascan.csv', sample_size=-1)) TO 'trnascan.parquet' (FORMAT PARQUET);"
    rm -f trnascan.csv
    """
}
