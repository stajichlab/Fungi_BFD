include { tablesDir } from '../../common/utils.nf'

// Transcription-factor inventory, derived from RUN_PFAM's domtblout hits
// (same input files as MERGE_PFAM) filtered against a curated list of
// fungal TF-associated Pfam accessions -- see bin/merge_tf_domains.py for why
// (FTFD's host does not resolve, so there's no dedicated TF DB to search).
process MERGE_TF {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(domtbls)

    output:
        path("tf_inventory.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_tf_domains.py \\
        --tf-domains ${projectDir}/assets/fungal_tf_pfam_domains.csv \\
        -o tf_inventory.csv \\
        ${domtbls}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tf_inventory.csv', sample_size=-1)) TO 'tf_inventory.parquet' (FORMAT PARQUET);"
    rm -f tf_inventory.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,pfam_acc,pfam_name,tf_family,domain_i_evalue,domain_score\\n' > tf_inventory.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tf_inventory.csv', sample_size=-1)) TO 'tf_inventory.parquet' (FORMAT PARQUET);"
    rm -f tf_inventory.csv
    """
}
