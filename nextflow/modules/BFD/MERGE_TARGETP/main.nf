include { tablesDir } from '../../common/utils.nf'

process MERGE_TARGETP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(summaries)

    output:
        path("targetP.parquet"), emit: parquet

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_targetp.py -o targetP.csv ${summaries}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('targetP.csv', sample_size=-1)) TO 'targetP.parquet' (FORMAT PARQUET);"
    rm -f targetP.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,prediction,probability,cleavage_position_start,cleavage_position_end,cleavage_probability,motif\\n' > targetP.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('targetP.csv', sample_size=-1)) TO 'targetP.parquet' (FORMAT PARQUET);"
    rm -f targetP.csv
    """
}
