include { tablesDir } from '../../common/utils.nf'

process MERGE_DEEPTMHMM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("deeptmhmm.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_deeptmhmm.py -o deeptmhmm.csv ${gff3s}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('deeptmhmm.csv', sample_size=-1)) TO 'deeptmhmm.parquet' (FORMAT PARQUET);"
    rm -f deeptmhmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,feature,start,end\\n' > deeptmhmm.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('deeptmhmm.csv', sample_size=-1)) TO 'deeptmhmm.parquet' (FORMAT PARQUET);"
    rm -f deeptmhmm.csv
    """
}
