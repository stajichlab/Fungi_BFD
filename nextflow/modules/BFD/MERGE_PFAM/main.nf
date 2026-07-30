include { tablesDir } from '../../common/utils.nf'

process MERGE_PFAM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(domtbls)

    output:
        path("pfam.csv.gz"), emit: csv

    script:
    """
    python3 ${params.scripts}/pfamtbl_to_long.py \\
        --outfile pfam.csv \\
        ${domtbls}
    pigz pfam.csv
    """

    stub:
    """
    printf 'protein_id,pfam_id,pfam_acc,pfam_len,full_seq_e_value,full_seq_score,full_seq_bias,domain_num,domain_num_of,domain_c_evalue,domain_i_evalue,domain_score,domain_bias,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to\\n' | gzip > pfam.csv.gz
    """
}
