include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_PFAM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(domtbls)

    output:
        path("pfam.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/pfam.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${domtbls} | LC_ALL=C sort > pfam.files
    ${dm} --sql ${projectDir}/sql/merge/pfam.sql --files pfam.files --out pfam.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('pfam.parquet')}
    """

    stub:
    """
    printf 'protein_id,pfam_id,pfam_acc,pfam_len,full_seq_e_value,full_seq_score,full_seq_bias,domain_num,domain_num_of,domain_c_evalue,domain_i_evalue,domain_score,domain_bias,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to\\n' > pfam.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('pfam.csv', sample_size=-1)) TO 'pfam.parquet' (FORMAT PARQUET);"
    rm -f pfam.csv
    """
}
