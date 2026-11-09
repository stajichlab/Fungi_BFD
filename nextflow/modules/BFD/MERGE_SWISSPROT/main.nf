include { tablesDir; codeMd5 } from '../../common/utils.nf'

// One task per genome bucket (sha1(LOCUSTAG)[:2], the storeDir bucket of the
// per-genome DIAMOND tables). DuckDB reads the 18-column tables directly and
// writes the narrow hit table of analysis/swissprot_design/REVIEW.md
// (sql/merge/swissprot_hits.sql): sp_id instead of accession/entry/name,
// rank, scaled identity/coverage, fixed coverages, -log10 e-value, ZSTD.
// On 1,018 genomes: 532 s vs 1,507 s for the Python path, 1.72 GB vs 6.36 GB,
// same rows and values (2026-11-06). Published to tables/swissprot_hits/.
process MERGE_SWISSPROT {
    tag        "${bucket}"
    label      'merge'
    publishDir path: { "${tablesDir()}/swissprot_hits" }, mode: 'copy'

    input:
        tuple val(bucket), path(blasttabs), path(subject)

    output:
        path("swissprot_hits.${bucket}.parquet"), emit: hits

    script:
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${blasttabs} | LC_ALL=C sort > swissprot.files
    python3 ${projectDir}/bin/duckdb_merge.py --sql ${projectDir}/sql/merge/swissprot_hits.sql \\
        --files swissprot.files --out swissprot_hits.${bucket}.parquet \\
        --threads ${task.cpus} --memory-limit ${(task.memory.toGiga() * 0.75) as int}GB \\
        --temp-dir \${SCRATCH:-\$PWD} --param subject="'${subject}'"
    """

    stub:
    """
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT 'X_1' AS protein_id, 1::INTEGER AS sp_id, 1::UTINYINT AS rk) TO 'swissprot_hits.${bucket}.parquet' (FORMAT PARQUET);"
    """
}
