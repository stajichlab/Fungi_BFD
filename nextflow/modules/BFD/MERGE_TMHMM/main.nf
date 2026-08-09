include { tablesDir } from '../../common/utils.nf'

process MERGE_TMHMM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(tsvs)

    output:
        path("tmhmm.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_tmhmm.py -o tmhmm.csv ${tsvs}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tmhmm.csv', sample_size=-1)) TO 'tmhmm.parquet' (FORMAT PARQUET);"
    rm -f tmhmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,len,ExpAA,First60,PredHel,Topology\\n' > tmhmm.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tmhmm.csv', sample_size=-1)) TO 'tmhmm.parquet' (FORMAT PARQUET);"
    rm -f tmhmm.csv
    """
}
