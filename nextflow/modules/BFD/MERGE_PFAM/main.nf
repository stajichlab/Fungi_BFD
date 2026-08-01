include { tablesDir } from '../../common/utils.nf'

process MERGE_PFAM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(domtbls)

    output:
        path("pfam.parquet"), emit: parquet

    script:
    """
    python3 ${params.scripts}/pfamtbl_to_long.py \\
        --outfile pfam.csv \\
        ${domtbls}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('pfam.csv', sample_size=-1)) TO 'pfam.parquet' (FORMAT PARQUET);"
    rm -f pfam.csv
    """

    stub:
    """
    printf 'protein_id,pfam_id,pfam_acc,pfam_len,full_seq_e_value,full_seq_score,full_seq_bias,domain_num,domain_num_of,domain_c_evalue,domain_i_evalue,domain_score,domain_bias,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to\\n' > pfam.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('pfam.csv', sample_size=-1)) TO 'pfam.parquet' (FORMAT PARQUET);"
    rm -f pfam.csv
    """
}
