include { tablesDir } from '../../common/utils.nf'

process MERGE_MEROPS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("merops.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_merops.py -o merops.csv ${blasttabs}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('merops.csv', sample_size=-1)) TO 'merops.parquet' (FORMAT PARQUET);"
    rm -f merops.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,merops_id,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' > merops.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('merops.csv', sample_size=-1)) TO 'merops.parquet' (FORMAT PARQUET);"
    rm -f merops.csv
    """
}
