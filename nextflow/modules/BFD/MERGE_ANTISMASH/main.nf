include { tablesDir } from '../../common/utils.nf'

process MERGE_ANTISMASH {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(jsons)

    output:
        path("antismash.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_antismash.py -o antismash.csv ${jsons}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('antismash.csv', sample_size=-1)) TO 'antismash.parquet' (FORMAT PARQUET);"
    rm -f antismash.csv
    """

    stub:
    """
    printf 'species_prefix,record_id,cluster_num,product,contig_edge,start,end\\n' > antismash.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('antismash.csv', sample_size=-1)) TO 'antismash.parquet' (FORMAT PARQUET);"
    rm -f antismash.csv
    """
}
