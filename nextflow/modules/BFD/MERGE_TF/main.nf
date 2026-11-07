include { tablesDir; coverageCheck } from '../../common/utils.nf'

// Transcription-factor inventory, derived from RUN_PFAM's domtblout hits
// (same input files as MERGE_PFAM) filtered against a curated list of
// fungal TF-associated Pfam accessions -- see bin/merge_tf_domains.py for why
// (FTFD's host does not resolve, so there's no dedicated TF DB to search).
process MERGE_TF {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(domtbls)

    output:
        path("tf_inventory.parquet"), emit: parquet

    script:
    // DuckDB parses the per-genome files directly (bin/duckdb_merge.py +
    // sql/merge/tf.sql); same rows, columns and types as the former Python
    // merge, checked on real BFD files (sql/merge/README_agent*.md, 2026-11-06).
    // File lists are written with the printf builtin (no ARG_MAX limit) and
    // sorted like the Python scripts sorted their inputs.
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${domtbls} | LC_ALL=C sort > tf.files
    ${dm} --sql ${projectDir}/sql/merge/tf.sql --files tf.files --out tf_inventory.parquet --temp-dir \${SCRATCH:-\$PWD} --param tf_domains=${projectDir}/assets/fungal_tf_pfam_domains.csv
    ${coverageCheck('tf_inventory.parquet')}
    """

    stub:
    """
    printf 'species_prefix,protein_id,pfam_acc,pfam_name,tf_family,domain_i_evalue,domain_score\\n' > tf_inventory.csv
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('tf_inventory.csv', sample_size=-1)) TO 'tf_inventory.parquet' (FORMAT PARQUET);"
    rm -f tf_inventory.csv
    """
}
