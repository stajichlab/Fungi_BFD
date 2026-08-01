include { tablesDir } from '../../common/utils.nf'

process MERGE_PREDGPI {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("predgpi.parquet"), emit: parquet

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_predgpi.py -o predgpi.csv ${gff3s}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('predgpi.csv', sample_size=-1)) TO 'predgpi.parquet' (FORMAT PARQUET);"
    rm -f predgpi.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,source,feature,start,end,score,strand,phase,attributes\\n' > predgpi.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('predgpi.csv', sample_size=-1)) TO 'predgpi.parquet' (FORMAT PARQUET);"
    rm -f predgpi.csv
    """
}
