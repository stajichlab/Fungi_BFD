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
    module load apptainer
    # ── Containerized funannotate ──────────────────────────────────────────────
    # Same swap as FUNANNOTATE_TRAIN/main.nf: `funannotate update` runs via
    # `\$SING` (apptainer exec \${params.funannotate_sif}), kept WITHOUT a
    # `container =` directive on this process -- the mysqldb sidecar below is
    # started via `singularity instance start` from this same host-executed
    # script, and both .sifs stay siblings launched from the host shell rather
    # than nesting one Singularity container inside another.
    #
    # Same confirmations as FUNANNOTATE_TRAIN/main.nf: Singularity shares the
    # host network namespace by default, so \$SING reaches the mysqldb instance
    # via MYHOSTNAME:PORT with no change needed; `apptainer exec` passes
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

    # APPTAINERENV_ prefix is REQUIRED here, not a plain export: the
    # funannotate-live image bakes its own ENV defaults for these vars
    # (AUGUSTUS_CONFIG_PATH=/venv/config, FUNANNOTATE_DB=/opt/databases/...),
    # sourced from the image's /.singularity.d/env/ scripts AFTER the host
    # shell env is inherited, so a plain `export` is silently overwritten
    # inside the container (confirmed empirically 2026-08-24 against
    # FUNANNOTATE_PREDICT's identical setup).
    export APPTAINERENV_AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export APPTAINERENV_FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    # PASA's hook loader (Pasa_conf.pm::_get_hook) resolves
    # HOOK_EXISTING_GENE_ANNOTATION_LOADER (GFF3::GFF3_annot_retriever, used by
    # funannotate update's annotation-comparison step) by searching @INC for
    # GFF3/GFF3_annot_retriever.pm. The container's default PERL5LIB is empty
    # and SAMPLE_HOOKS isn't on the default @INC, so this crashes with
    # "Error, couldn't resolve path for GFF3::GFF3_annot_retriever" on every
    # real update run. apptainer exec passes host env through by default, so
    # exporting PERL5LIB here (before the apptainer exec call) reaches the
    # container. Confirmed
    # via direct `apptainer exec ... perl -e '...@INC...'` test, 2026-08-15.
    export PERL5LIB="/venv/opt/pasa/src/SAMPLE_HOOKS:/venv/opt/pasa/src/PerlLib"
    pasa_db_arg="--pasa_db sqlite"
    # \$PWD (the task workdir) holds the Nextflow-staged r1/r2 read files passed
    # directly into `funannotate update`, plus the \${out}/training symlink
    # created below -- neither is under target/training_target, so \$PWD needs
    # its own explicit bind. Same missing-bind bug class as GENEMARK_RUN/
    # FUNANNOTATE_PREDICT (confirmed 2026-08-24).
    SING_BINDS="--bind \$PWD:\$PWD,${params.training_target}:${params.training_target},${params.target}:${params.target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        # Lives under \$TMPDIR (== node-local \$SCRATCH under SLURM) rather
        # than under training_target on shared storage -- see
        # FUNANNOTATE_TRAIN/main.nf (confirmed 2026-09-01) for why: this
        # datadir is pure per-task sidecar infra that nothing downstream ever
        # reads back, so a persistent, species-keyed copy on /bigdata just
        # let a crashed prior attempt strand an orphaned InnoDB tablespace
        # file that broke PASA's "-r" drop-then-recreate ("Directory not
        # empty" / "database exists"). Each SLURM job gets its own fresh
        # node-local scratch dir, so that bug class is now impossible. The
        # \$TMPDIR:\$TMPDIR bind above already covers this path.
        MYSQL_SCRATCH=\$TMPDIR/mysql_db_${out}
        rm -rf \$MYSQL_SCRATCH
        mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
        rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
            { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
        cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
            { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # No separate SING_BINDS entry needed for \$MYSQL_SCRATCH (\$PASACONF
        # lives under it) -- it's now a subdirectory of \$TMPDIR, already
        # covered by the \$TMPDIR:\$TMPDIR bind above.
        # >/dev/null too, not just stderr: apptainer/singularity prints its
        # "Usage: ... instance stop ..." help text to STDOUT (not stderr) when
        # the named instance is already gone -- happens every run here since
        # the EXIT trap always re-fires this after the explicit stop_mysqldb
        # call near the end of the script already stopped it (see
        # FUNANNOTATE_TRAIN/main.nf, same pattern).
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} >/dev/null 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        # No trailing "/usr/bin/mysqld_safe" here: `instance start` only runs
        # trailing args THROUGH the image's %startscript (see --help), and
        # Docker's ENTRYPOINT/CMD populate %runscript, not %startscript --
        # mariadb_sif is a custom build (nextflow/docs/mariadb-10.3.9.def)
        # whose %startscript itself execs mysqld_safe, so it needs no args.
        # Same fix as FUNANNOTATE_TRAIN/main.nf, confirmed 2026-08-25: a plain
        # docker://mariadb:10.3.9 (or an `apptainer pull`-cached sif of it,
        # both lacking a real %startscript) starts an empty instance that
        # never launches mysqld at all -- PASA then fails to connect.
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid}
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
