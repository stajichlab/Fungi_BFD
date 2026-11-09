include { tablesDir; coverageCheck; codeMd5 } from '../../common/utils.nf'

process MERGE_WOLFPSORT {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(results)

    output:
        path("wolfpsort.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/wolfpsort.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${results} | LC_ALL=C sort > wolfpsort.files
    ${dm} --sql ${projectDir}/sql/merge/wolfpsort.sql --files wolfpsort.files --out wolfpsort.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('wolfpsort.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,localization,score\\n' > wolfpsort.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('wolfpsort.csv', sample_size=-1)) TO 'wolfpsort.parquet' (FORMAT PARQUET);"
    rm -f wolfpsort.csv
    """
}
