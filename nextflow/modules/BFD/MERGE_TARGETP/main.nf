include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_TARGETP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(summaries)

    output:
        path("targetP.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/targetp.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${summaries} | LC_ALL=C sort > targetp.files
    ${dm} --sql ${projectDir}/sql/merge/targetp.sql --files targetp.files --out targetP.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('targetP.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,prediction,probability,cleavage_position_start,cleavage_position_end,cleavage_probability,motif\\n' > targetP.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('targetP.csv', sample_size=-1)) TO 'targetP.parquet' (FORMAT PARQUET);"
    rm -f targetP.csv
    """
}
