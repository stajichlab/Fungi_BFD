include { tablesDir; coverageCheck; codeMd5 } from '../../common/utils.nf'

process MERGE_MEROPS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("merops.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/merops.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${blasttabs} | LC_ALL=C sort > merops.files
    ${dm} --sql ${projectDir}/sql/merge/merops.sql --files merops.files --out merops.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('merops.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,merops_id,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' > merops.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('merops.csv', sample_size=-1)) TO 'merops.parquet' (FORMAT PARQUET);"
    rm -f merops.csv
    """
}
