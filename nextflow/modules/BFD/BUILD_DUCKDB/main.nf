include { dbDir } from '../../common/utils.nf'

// Rebuilds db/BFD.duckdb from tables/*.parquet (bin/build_BFD_duckDB.sh).
// `sync` is a gating manifest, not a real input -- the script reads
// params.tables directly off disk (it globs the whole tables/ dir, not a
// per-table file list), so `sync` only exists to make this task wait until
// every MERGE_* process invoked this run has published its parquet, same
// gating idiom as gatedGlobIn/toManifest elsewhere in BFD_MERGE.nf.
process BUILD_DUCKDB {
    label      'duckdb'
    publishDir path: { dbDir() }, mode: 'copy'

    input:
    path(sync)

    output:
    path "BFD.duckdb", emit: db

    script:
    """
    TABLES=${params.tables} DBDIR=. bash ${projectDir}/bin/build_BFD_duckDB.sh
    """

    stub:
    """
    touch BFD.duckdb
    """
}
