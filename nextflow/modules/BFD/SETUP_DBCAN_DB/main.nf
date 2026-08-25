// Download the full run_dbcan reference DB set (CAZyDB, dbCAN-PUL, dbCAN_sub,
// tf-1/tf-2, stp, tcdb-fa, etc.) needed by CGC/PUL clustering -- distinct from
// (and larger than) the lightweight dbcanlight HMM set RUN_CAZY already uses,
// since dbcanlight has no `cgc` mode (confirmed via `dbcanlight search --help`
// 2026-08-20: only cazyme/sub/diamond). storeDir caches this at
// params.dbcan_dbdir so it downloads at most once.
//
// `run_dbcan database --db_dir . --aws_s3` verified against the real dbCAN
// 5.2.9 image (params.dbcan_cgc_sif) via `run_dbcan database --help`
// 2026-08-20, and matches nf-core/funcscan's RUNDBCAN_DATABASE module
// (same image version) exactly.
process SETUP_DBCAN_DB {
    label 'setup'

    storeDir { params.dbcan_dbdir }

    output:
        path("*"), emit: db

    script:
    """
    set -euo pipefail
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind \${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.dbcan_cgc_sif} run_dbcan database --db_dir . --aws_s3
    """

    stub:
    """
    touch CAZyDB.stub.fa
    """
}
