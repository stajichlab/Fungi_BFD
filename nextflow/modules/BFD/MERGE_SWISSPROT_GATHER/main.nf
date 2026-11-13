include { tablesDir; coverageCheck; codeMd5 } from '../../common/utils.nf'

// After all MERGE_SWISSPROT buckets: swissprot_best.parquet (rank-1 hit per
// protein with subject strings) and swissprot_subject_summary.parquet (per
// subject counts). Reads the staged bucket files, not the published copies.
process MERGE_SWISSPROT_GATHER {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(hits)
        path(subject)

    output:
        path("swissprot_best.parquet"),            emit: best
        path("swissprot_subject_summary.parquet"), emit: summary

    script:
    def dm = "python3 ${projectDir}/bin/duckdb_merge.py --threads ${task.cpus} " +
             "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    # code md5 (task cache key): ${codeMd5(['bin/duckdb_merge.py', 'sql/merge'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${hits} | LC_ALL=C sort > hits.files
    ${dm} --sql ${projectDir}/sql/merge/swissprot_best.sql --files hits.files \\
        --out swissprot_best.parquet --temp-dir \${SCRATCH:-\$PWD} --param subject="'${subject}'"
    ${coverageCheck('swissprot_best.parquet')}
    ${dm} --sql ${projectDir}/sql/merge/swissprot_subject_summary.sql --files hits.files \\
        --out swissprot_subject_summary.parquet --temp-dir \${SCRATCH:-\$PWD}
    """

    stub:
    """
    touch swissprot_best.parquet swissprot_subject_summary.parquet
    """
}
