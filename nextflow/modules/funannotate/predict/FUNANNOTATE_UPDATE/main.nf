process FUNANNOTATE_UPDATE {
    tag "$out"

    cpus   16
    memory '96 GB'
    time   '48h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(r1), path(r2)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    script:
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no reads (empty marker file from SRA_FETCH) ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate update"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity
    # ── Containerized funannotate ──────────────────────────────────────────────
    # Same swap as FUNANNOTATE_TRAIN/main.nf: `funannotate update` runs via
    # `\$SING` (singularity exec \${params.funannotate_sif}), kept WITHOUT a
    # `container =` directive on this process -- the mysqldb sidecar below is
    # started via `singularity instance start` from this same host-executed
    # script, and both .sifs stay siblings launched from the host shell rather
    # than nesting one Singularity container inside another.
    #
    # Same confirmations as FUNANNOTATE_TRAIN/main.nf: Singularity shares the
    # host network namespace by default, so \$SING reaches the mysqldb instance
    # via MYHOSTNAME:PORT with no change needed; `singularity exec` passes
    # through all host-shell env vars (PASACONF etc.) unless --cleanenv is
    # added, so don't add --cleanenv without converting those to explicit
    # `--env` flags first; and the `which` failure some tools show under
    # `bash -c` inside this image is a HOST env leak (Rocky's GNU-which-
    # flavored bash function via BASH_FUNC_which%%), not an image bug -- it
    # doesn't affect PASA/funannotate's own Perl system()/backtick calls
    # (those spawn dash, which never imports BASH_FUNC_*%%). Unset defensively
    # below anyway.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    # PASA's hook loader (Pasa_conf.pm::_get_hook) resolves
    # HOOK_EXISTING_GENE_ANNOTATION_LOADER (GFF3::GFF3_annot_retriever, used by
    # funannotate update's annotation-comparison step) by searching @INC for
    # GFF3/GFF3_annot_retriever.pm. The container's default PERL5LIB is empty
    # and SAMPLE_HOOKS isn't on the default @INC, so this crashes with
    # "Error, couldn't resolve path for GFF3::GFF3_annot_retriever" on every
    # real update run. singularity exec passes host env through by default, so
    # exporting PERL5LIB here (before the singularity exec call) reaches the
    # container. Confirmed
    # via direct `singularity exec ... perl -e '...@INC...'` test, 2026-08-15.
    export PERL5LIB="/venv/opt/pasa/src/SAMPLE_HOOKS:/venv/opt/pasa/src/PerlLib"
    pasa_db_arg="--pasa_db sqlite"
    SING_BINDS="--bind ${params.training_target}:${params.training_target},${params.target}:${params.target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},\$TMPDIR:\$TMPDIR"
    SING="singularity exec \${SING_BINDS} ${params.funannotate_sif}"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.training_target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # Bind the WHOLE MYSQL_SCRATCH dir (see FUNANNOTATE_TRAIN/main.nf for
        # why -- MYSQL_SCRATCH/mysql_db, a prior value, was never created by
        # anything above and was a dead bind source; \$PASACONF lives under
        # MYSQL_SCRATCH/conf, which this now correctly covers). Appended to
        # SING_BINDS rather than a separate SINGULARITY_BINDPATH env var, to
        # match this project's SING_BINDS/SING convention.
        SING_BINDS="\${SING_BINDS},\$MYSQL_SCRATCH:\$MYSQL_SCRATCH"
        SING="singularity exec \${SING_BINDS} ${params.funannotate_sif}"
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Link training data into work dir so funannotate update finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "${out}/training"
    fi

    # A prior crashed update run can leave an empty update_results/ dir (mkdir
    # happens before funannotate writes anything into it). funannotate annotate's
    # directory-detection is `if isdir(update_results): use it` -- unconditional,
    # no content check -- so a stale empty dir permanently poisons annotate even
    # though predict_results/ is intact. Clear it before retrying so a failed
    # attempt doesn't leave annotate stuck. Found via real run, 2026-08-15.
    if [ -d "${params.target}/${out}/update_results" ] && [ ! -f "${params.target}/${out}/update_results/${out}.gbk" ]; then
        rm -rf "${params.target}/${out}/update_results"
    fi

    # r1/r2 are pre-normalized reads from SRA_FETCH (fastp-trimmed + bbnorm-normalized).
    # funannotate update will still run its internal alignment step against these.
    echo "[INFO] Running funannotate update for ${out}"
    \$SING funannotate update -i ${params.target}/${out} \\
        --left ${r1} --right ${r2} \\
        --cpus ${task.cpus} \\
        \$pasa_db_arg
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    EXPECTED="${params.target}/${out}/update_results/${out}.gbk"
    if [ ! -f "\$EXPECTED" ]; then
        echo "ERROR: funannotate update did not produce expected GBK: \$EXPECTED" >&2
        exit 1
    fi
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_UPDATE stub for ${out} (r1=${r1}, r2=${r2})"
    mkdir -p ${params.target}/${out}/update_results
    touch ${params.target}/${out}/update_results/${out}.tbl
    touch ${params.target}/${out}/update_results/${out}.gbk
    touch ${params.target}/${out}/update_results/${out}.gff3
    """
}
