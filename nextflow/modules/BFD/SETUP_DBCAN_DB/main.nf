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
    # Link the CAZy files from a local release (params.dbcan_cazy_src) under
    # the names run_dbcan expects, then download only what is still missing
    # (--no-overwrite). Link targets are absolute, so they survive the
    # storeDir move.
    SRC="${params.dbcan_cazy_src ?: ''}"
    if [ -n "\$SRC" ]; then
        CAZY_DMND=\$(ls "\$SRC"/CAZyDB.*.dmnd | sort | tail -1)
        HMM_TXT=\$(ls "\$SRC"/dbCAN-HMMdb-V*.txt | sort -V | tail -1)
        MAP_TSV=\$(ls "\$SRC"/fam-substrate-mapping-*.tsv | sort | tail -1)
        ln -s "\$(readlink -f "\$CAZY_DMND")" CAZy.dmnd
        for e in "" .h3f .h3i .h3m .h3p; do
            ln -s "\$(readlink -f "\$HMM_TXT\$e")" "dbCAN.hmm\$e"
            ln -s "\$(readlink -f "\$SRC/dbCAN_sub.hmm\$e")" "dbCAN_sub.hmm\$e"
        done
        ln -s "\$(readlink -f "\$MAP_TSV")" fam-substrate-mapping.tsv
    fi
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    SING_BINDS="--bind \${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.dbcan_cgc_sif} run_dbcan database --db_dir . --aws_s3 --no-overwrite
    """

    stub:
    """
    touch CAZyDB.stub.fa
    """
}
