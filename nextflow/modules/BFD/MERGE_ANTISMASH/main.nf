include { tablesDir; codeMd5 } from '../../common/utils.nf'

// antiSMASH-lite tables (schema v2, docs/BFD_ANTISMASH_DESIGN.md): per-genome
// <LOCUSTAG>.antismash.{regions,cds,known,known_pairs}.tsv.gz -> four Parquet tables
// bgc_cluster, bgc_cluster_gene, bgc_mibig_hit, bgc_mibig_pair (join key: locustag, record, region).
process MERGE_ANTISMASH {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(tables)

    output:
        path("bgc_cluster.parquet"),      emit: parquet
        path("bgc_cluster_gene.parquet"), emit: gene
        path("bgc_mibig_hit.parquet"),    emit: hit
        path("bgc_mibig_pair.parquet"),   emit: pair

    script:
    // Coverage check (as coverageCheck() in utils.nf) on bgc_cluster only; its key column is `locustag`.
    def cov = (!params.merge_all.toBoolean() || params.allow_shrink.toBoolean()) ? '' :
        "python3 ${projectDir}/bin/check_table_coverage.py --new bgc_cluster.parquet --old ${tablesDir()}/bgc_cluster.parquet --key locustag --min-frac ${params.coverage_min_frac} || exit 1"
    """
    # code md5 (task cache key): ${codeMd5(['bin/merge_antismash.py', 'sql/merge/antismash.sql'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' ${tables} | LC_ALL=C sort > antismash.files
    python3 ${projectDir}/bin/merge_antismash.py --threads ${task.cpus} -o . @antismash.files
    ${cov}
    """

    stub:
    """
    module load duckdb 2>/dev/null || true
    for t in bgc_cluster bgc_cluster_gene bgc_mibig_hit bgc_mibig_pair; do
        duckdb -c "COPY (SELECT NULL::VARCHAR AS locustag WHERE false) TO '\$t.parquet' (FORMAT PARQUET);"
    done
    """
}
