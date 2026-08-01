include { tablesDir } from '../../common/utils.nf'

process MERGE_IDP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path 'idp/*'
        path 'sum/*'

    output:
        path("idp.parquet"),         emit: idp
        path("idp_summary.parquet"), emit: summary

    script:
    """
    first=1
    for f in idp/*.idp.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > idp.csv.gz

    first=1
    for f in sum/*.idp_summary.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > idp_summary.csv.gz

    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp.csv.gz', sample_size=-1)) TO 'idp.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp_summary.csv.gz', sample_size=-1)) TO 'idp_summary.parquet' (FORMAT PARQUET);"
    rm -f idp.csv.gz idp_summary.csv.gz
    """

    stub:
    """
    printf 'species_prefix,protein_id,IDP_start,IDP_end,IDP_length,mean_score\\n' | gzip > idp.csv.gz
    printf 'species_prefix,protein_id,IDP_residues,IDP_fraction,length\\n'        | gzip > idp_summary.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp.csv.gz', sample_size=-1)) TO 'idp.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp_summary.csv.gz', sample_size=-1)) TO 'idp_summary.parquet' (FORMAT PARQUET);"
    rm -f idp.csv.gz idp_summary.csv.gz
    """
}
