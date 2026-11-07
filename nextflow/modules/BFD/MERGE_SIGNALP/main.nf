include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_SIGNALP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("signalp.signal_peptide.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/signalp.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${gff3s} | LC_ALL=C sort > signalp.files
    ${dm} --sql ${projectDir}/sql/merge/signalp.sql --files signalp.files --out signalp.signal_peptide.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('signalp.signal_peptide.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,peptide_start,peptide_end,probability\\n' > signalp.signal_peptide.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('signalp.signal_peptide.csv', sample_size=-1)) TO 'signalp.signal_peptide.parquet' (FORMAT PARQUET);"
    rm -f signalp.signal_peptide.csv
    """
}
