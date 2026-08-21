include { tablesDir } from '../../common/utils.nf'

process MERGE_TCDB {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("tcdb.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_tcdb.py -o tcdb.csv ${blasttabs}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tcdb.csv', sample_size=-1)) TO 'tcdb.parquet' (FORMAT PARQUET);"
    rm -f tcdb.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,tc_number,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' > tcdb.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tcdb.csv', sample_size=-1)) TO 'tcdb.parquet' (FORMAT PARQUET);"
    rm -f tcdb.csv
    """
}
