include { tablesDir } from '../../common/utils.nf'

process MERGE_SWISSPROT {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("swissprot.parquet"), emit: parquet

    script:
    """
    python3 ${projectDir}/bin/merge_swissprot.py -o swissprot.csv ${blasttabs}
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot.csv', sample_size=-1)) TO 'swissprot.parquet' (FORMAT PARQUET);"
    rm -f swissprot.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,swissprot_acc,swissprot_entry,swissprot_name,pident,positive,nident,aln_length,q_len,s_len,qcovhsp,query_cov,hit_cov,q_start,q_end,s_start,s_end,gapopen,mismatch,evalue,bitscore,func_transfer_80_80\\n' > swissprot.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('swissprot.csv', sample_size=-1)) TO 'swissprot.parquet' (FORMAT PARQUET);"
    rm -f swissprot.csv
    """
}
