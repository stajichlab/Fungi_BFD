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
    // The MERGE_* outputs in `sync` are published to params.tables by publishDir,
    // which copies in the background after each task ends; this task can start
    // before a large copy finishes (2026-11-06: BUILD_DUCKDB read half-copied
    // gene_info/gene_proteins.parquet, "No magic bytes found"). Wait until each
    // published copy has the staged file's size and a Parquet footer ('PAR1').
    """
    set -euo pipefail
    deadline=\$(( \$(date +%s) + 7200 ))
    for f in ${sync}; do
        case "\$f" in *.parquet) ;; *) continue ;; esac
        b=\$(basename "\$f")
        want=\$(stat -L -c %s "\$f")
        # published at the top of tables/ or one folder down (swissprot_hits/)
        until dst=\$( { ls -d "${params.tables}/\$b" "${params.tables}"/*/"\$b" 2>/dev/null || true; } | head -1) \\
              && [ -n "\$dst" ] && [ "\$(stat -c %s "\$dst")" = "\$want" ] \\
              && [ "\$(tail -c 4 "\$dst")" = "PAR1" ]; do
            [ \$(date +%s) -lt \$deadline ] || { echo "ERROR: \$b not published under ${params.tables} after 2 h" >&2; exit 1; }
            sleep 30
        done
    done
    TABLES=${params.tables} DBDIR=. \\
    BIN_DIR=${projectDir}/bin \\
    SCHEMA=${projectDir}/../sql/table_schema.json \\
    bash ${projectDir}/bin/build_BFD_duckDB.sh
    """

    stub:
    """
    touch BFD.duckdb
    """
}
