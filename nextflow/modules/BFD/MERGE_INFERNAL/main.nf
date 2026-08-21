include { tablesDir } from '../../common/utils.nf'

process MERGE_INFERNAL {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(tblouts)

    output:
        path("infernal_rfam.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_infernal.py -o infernal_rfam.csv ${tblouts}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('infernal_rfam.csv', sample_size=-1)) TO 'infernal_rfam.parquet' (FORMAT PARQUET);"
    rm -f infernal_rfam.csv
    """

    stub:
    """
    printf 'species_prefix,target_name,target_acc,query_name,contig,mdl_from,mdl_to,seq_from,seq_to,strand,evalue,score\\n' > infernal_rfam.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('infernal_rfam.csv', sample_size=-1)) TO 'infernal_rfam.parquet' (FORMAT PARQUET);"
    rm -f infernal_rfam.csv
    """
}
