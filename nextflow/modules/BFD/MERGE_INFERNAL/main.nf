include { tablesDir; coverageCheck; codeMd5 } from '../../common/utils.nf'

process MERGE_INFERNAL {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(tblouts)

    output:
        path("infernal_rfam.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/infernal.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${tblouts} | LC_ALL=C sort > infernal.files
    ${dm} --sql ${projectDir}/sql/merge/infernal.sql --files infernal.files --out infernal_rfam.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('infernal_rfam.parquet')}
    """

    stub:
    """
    printf 'species_prefix,target_name,target_acc,query_name,contig,mdl_from,mdl_to,seq_from,seq_to,strand,evalue,score\\n' > infernal_rfam.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('infernal_rfam.csv', sample_size=-1)) TO 'infernal_rfam.parquet' (FORMAT PARQUET);"
    rm -f infernal_rfam.csv
    """
}
