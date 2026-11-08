include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_CAZY {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(overviews)
        path(cazymes)

    output:
        path("cazy.overview.parquet"),    emit: overview_parquet
        path("cazy.cazymes_hmm.parquet"), emit: hmm_parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/cazy_overview/cazy_cazymes_hmm.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${overviews} | LC_ALL=C sort > cazy_overview.files
    ${dm} --sql ${projectDir}/sql/merge/cazy_overview.sql --files cazy_overview.files --out cazy.overview.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('cazy.overview.parquet')}
    printf '%s\\n' ${cazymes} | LC_ALL=C sort > cazy_cazymes_hmm.files
    ${dm} --sql ${projectDir}/sql/merge/cazy_cazymes_hmm.sql --files cazy_cazymes_hmm.files --out cazy.cazymes_hmm.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('cazy.cazymes_hmm.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,EC,cazyme_fam,sub_fam,diamond_fam,substrate,toolcount\\n' > cazy.overview.csv
    printf 'species_prefix,HMM_id,profile_length,protein_id,protein_length,evalue,q_start,q_end,s_start,s_end,coverage\\n' > cazy.cazymes_hmm.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.overview.csv', sample_size=-1)) TO 'cazy.overview.parquet' (FORMAT PARQUET);"
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('cazy.cazymes_hmm.csv', sample_size=-1)) TO 'cazy.cazymes_hmm.parquet' (FORMAT PARQUET);"
    rm -f cazy.overview.csv cazy.cazymes_hmm.csv
    """
}
