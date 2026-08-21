include { tablesDir } from '../../common/utils.nf'

process MERGE_PHOBIUS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(short_txts)

    output:
        path("phobius.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_phobius.py -o phobius.csv ${short_txts}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('phobius.csv', sample_size=-1)) TO 'phobius.parquet' (FORMAT PARQUET);"
    rm -f phobius.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,tm_count,sp_predicted,prediction\\n' > phobius.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('phobius.csv', sample_size=-1)) TO 'phobius.parquet' (FORMAT PARQUET);"
    rm -f phobius.csv
    """
}
