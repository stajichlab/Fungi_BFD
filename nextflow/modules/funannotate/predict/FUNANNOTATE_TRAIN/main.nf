// For non-representative strains: funannotate train --trinity <shared_fasta> runs only
// PASA (skips Trimmomatic, normalization, HISAT2, and Trinity-GG assembly).
// Falls back to a full train when no shared Trinity is available (e.g. species with
// a single strain or when run_sra_fetch is false).
process FUNANNOTATE_TRAIN {
    tag "$out"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa), val(pasa_tier)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), emit: predict_input
    // Audit row emitted only when a composite/composite_fallback-tier run gracefully
    // degrades (see the TRAIN_STATUS handling near the end) -- collected by
    // FUNANNOTATE_RNASEQ.nf into one reviewable TSV, same convention as
    // repr_assignments.tsv/rnaseq_blacklist_candidates.csv elsewhere in this project.
    path("${out}.composite_train_failed.tsv"), optional: true, emit: composite_failed

    script:
    def pasa_db_arg = "--pasa_db sqlite"
    // Real, symlink-resolved locations of rnaseq_reads/ and rnaseq_data/ (both
    // are themselves top-level symlinks into ../rnaseq_reads, ../rnaseq_data
    // from launchDir) -- these are NOT in SING_BINDS below, so funannotate's
    // --left_norm/--right_norm/--trinity handling, which always resolves its
    // read/Trinity arguments to their realpath() before symlinking them into
    // its own output tree, silently produces a dangling symlink: the bind is
    // missing, so the container can create the symlink (doesn't require the
    // target to exist) but can never stat through it, and every
    // `funannotate train` invocation has been failing with "Read
    // normalization failed, .../left.norm.fq.gz does not exist" as a result
    // -- confirmed 2026-08-25 by reproducing train.py's exact symlink logic
    // both with and without this bind present.
    def rnaseqReadsDir = file("${launchDir}/rnaseq_reads").toRealPath()
    def rnaseqDataDir  = file("${params.rnaseq_data}").toRealPath()
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if the shared Trinity-GG assembly is too thin to train on ────────
    # A too-few-transcripts trinity_fa usually means it was assembled against the
    # wrong reference strain (same species name, divergent genome) rather than a
    # real expression signal -- PASA has nothing to build a training set from and
    # funannotate train either crashes or emits junk models. See
    # train_min_trinity_transcripts in conf/profile_funannotate.config.
    if [ -s "${trinity_fa}" ] && [ "${params.train_min_trinity_transcripts}" -gt 0 ]; then
        TRINITY_TX_COUNT=\$(grep -c '^>' "${trinity_fa}" || true)
        if [ "\$TRINITY_TX_COUNT" -lt "${params.train_min_trinity_transcripts}" ]; then
            echo "[WARN] ${out}: shared Trinity-GG assembly has only \$TRINITY_TX_COUNT transcripts (< ${params.train_min_trinity_transcripts}); likely assembled against the wrong reference strain. Skipping funannotate train -- rerun scripts/fix_low_trinity.py to pick a better reference." >&2
            exit 0
        fi
    fi

    # ── Skip if training is already resolved and evidence is not newer ────────────
    # "Resolved" is either a real published PASA GFF3, OR a durable marker recording
    # that a composite/composite_fallback-tier attempt legitimately couldn't train
    # from transcript evidence (see the TRAIN_STATUS handling near the end and
    # nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md). Without this marker, a
    # gracefully-degraded hybrid strain would never publish anything and would be
    # re-attempted on every single future run/-resume forever -- the marker makes
    # that outcome durable, the same way a real GFF3 makes success durable.
    # Accept a compressed prediction (.gbk.gz) as "done" so folders can be space-saved.
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    TRAIN_FAILED_MARKER="${params.training_target}/${out}/training/.composite_train_failed"
    PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk"
    [ -f "\$PREDICT_GBK" ] || PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk.gz"
    RESOLVED_MARKER=""
    if [ -f "\$TRAIN_GFF3" ]; then
        RESOLVED_MARKER="\$TRAIN_GFF3"
    elif [ -f "\$TRAIN_FAILED_MARKER" ]; then
        RESOLVED_MARKER="\$TRAIN_FAILED_MARKER"
    fi
    if [ -n "\$RESOLVED_MARKER" ]; then
        RETRAIN=0
        # Re-train if the genome assembly itself is newer than the resolved marker --
        # otherwise a swapped-in new assembly version silently keeps stale training
        # (or a stale "not trainable" verdict) forever, since the checks below only
        # ever look at rnaseq/trinity evidence.
        if [ -s "${genome_fa}" ] && [ "${genome_fa}" -nt "\$RESOLVED_MARKER" ]; then
            echo "[INFO] Genome assembly newer than training output for ${out}; retraining"
            rm -rf "${params.training_target}/${out}/training"
            RETRAIN=1
        fi
        # Re-train if the trinity/composite evidence itself is newer than the resolved
        # marker -- covers a rebuilt composite (new/changed parent set, retuned PASA
        # composite thresholds) for a strain previously marked not-trainable, which
        # r1/se newer-than-GBK below can never catch since hybrid r1/se are always
        # empty placeholders.
        if [ -s "${trinity_fa}" ] && [ "${trinity_fa}" -nt "\$RESOLVED_MARKER" ]; then
            echo "[INFO] Trinity/composite evidence newer than training output for ${out}; retraining"
            RETRAIN=1
        fi
        if [ -f "\$PREDICT_GBK" ]; then
            # Re-train if the rnaseq reads are newer than the existing prediction GBK.
            if [ -s "${r1}" ] && [ "${r1}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq R1 reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            elif [ -s "${se}" ] && [ "${se}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq SE reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            fi
        fi
        if [ \$RETRAIN -eq 0 ]; then
            if [ "\$RESOLVED_MARKER" = "\$TRAIN_GFF3" ]; then
                echo "[INFO] Training already complete for ${out}; skipping"
                # Still repoint any absolute cross-project symlinks left by an older run.
                # Called by absolute path: this process runs via `bash -l`, and the login
                # shell's profile/module init rebuilds PATH, dropping the bin/ dir nextflow
                # otherwise appends for us.
                ${workflow.projectDir}/bin/relink_training_symlinks.py --apply --quiet \\
                    "${params.training_target}/${out}/training" || true
            else
                echo "[INFO] ${out} previously determined not trainable from composite transcript evidence (pasa_tier=${pasa_tier}); skipping re-attempt"
            fi
            exit 0
        fi
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer
    # ── Containerized funannotate ──────────────────────────────────────────────
    # funannotate-live's Docker image (Rust-optimized PASA/EVM/Trinity + AVX2
    # bowtie2) run via `apptainer exec \${params.funannotate_sif}`. Kept
    # WITHOUT a `container =` directive on this process -- the mysqldb sidecar
    # below is started via `singularity instance start` from this same
    # host-executed script, and both .sifs stay siblings launched from the
    # host shell rather than nesting one Singularity container inside another.
    #
    # Confirmed about the mysql sidecar + this exec:
    #   - Network: Singularity shares the host network namespace by default (no
    #     --net in use anywhere here), so `singularity instance start` and the
    #     separate `\$SING` exec below reach each other via MYHOSTNAME:PORT
    #     exactly as the host-module funannotate did.
    #   - Env passthrough: `apptainer exec` passes through ALL host-shell env
    #     vars by default -- do not add `--cleanenv` without converting
    #     PASACONF/AUGUSTUS_CONFIG_PATH/FUNANNOTATE_DB/TMPDIR to explicit
    #     `--env` flags first, or PASA loses its mysql connection info silently.
    #   - `which` false alarm: `bash -c 'which <tool>'` inside this image can
    #     fail (Illegal option --) if the calling host shell exports a
    #     GNU-which-flavored `which()` bash function (Rocky's which2 package) --
    #     it leaks in via env passthrough as BASH_FUNC_which%%. This is a HOST
    #     artifact, not an image bug, and does NOT affect funannotate/PASA/
    #     Trinity's own Perl `system()`/backtick calls (those spawn `/bin/sh`
    #     = dash, which never imports BASH_FUNC_*%% vars at all -- confirmed
    #     `sh -c 'which perl'` works cleanly). Unset defensively below anyway
    #     so nothing that explicitly execs bash inside the container can hit it.
    #
    # All host paths funannotate/PASA touch that live outside Nextflow's
    # auto-bound work dir must be explicit -B binds -- Singularity does NOT
    # bind arbitrary host paths by default. See SING_BINDS below.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    # APPTAINERENV_ prefix is REQUIRED here, not a plain export: the
    # funannotate-live image bakes its own ENV defaults for these vars
    # (AUGUSTUS_CONFIG_PATH=/venv/config, FUNANNOTATE_DB=/opt/databases/...),
    # sourced from the image's /.singularity.d/env/ scripts AFTER the host
    # shell env is inherited, so a plain `export` is silently overwritten
    # inside the container (confirmed empirically 2026-08-24 against
    # FUNANNOTATE_PREDICT's identical setup). PASACONF is unaffected (the
    # image has no conflicting default) so stays a plain export.
    export APPTAINERENV_AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export APPTAINERENV_FUNANNOTATE_DB=${params.funannotate_db}
    export TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # \$PWD (the task workdir) holds genome_input.fa (inflated below) plus the
    # Nextflow-staged r1/r2/se/trinity_fa read files passed directly into
    # `funannotate train` -- none of those are under target/training_target,
    # so \$PWD needs its own explicit bind. Same missing-bind bug class as
    # GENEMARK_RUN/FUNANNOTATE_PREDICT (confirmed 2026-08-24).
    SING_BINDS="--bind \$PWD:\$PWD,${params.training_target}:${params.training_target},${params.target}:${params.target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},${rnaseqReadsDir}:${rnaseqReadsDir},${rnaseqDataDir}:${rnaseqDataDir},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        # Lives under \$TMPDIR (== node-local \$SCRATCH under SLURM, see the
        # `export TMPDIR` above), NOT under training_target on shared storage.
        # This datadir is pure per-task sidecar infra -- it's explicitly
        # excluded from the results published to TRAINDIR below and nothing
        # downstream ever reads it back -- so it never needs to survive past
        # this task. Previously it lived under training_target (persistent,
        # shared /bigdata) keyed only by species, so a prior attempt that
        # crashed mid-write (OOM, node failure, job kill) could leave an
        # orphaned InnoDB tablespace file behind; PASA's "-r"
        # (drop-then-recreate) DROP DATABASE couldn't rmdir the stale
        # directory (errno 39 "Directory not empty"), and the following
        # CREATE DATABASE then failed with "database exists". Confirmed
        # 2026-09-01 against Aspergillus_arachidicola_CBS_117612's
        # pasa-assembly.log (orphaned splice_variation.ibd from an Aug 10
        # attempt). Putting it on \$SCRATCH instead makes that bug class
        # impossible: each SLURM job (== each task attempt) gets its own
        # fresh, node-local scratch dir that SLURM tears down when the job
        # ends, so there is nothing left for a later attempt to inherit. The
        # \$TMPDIR:\$TMPDIR bind above already covers this path, so no extra
        # SING_BINDS entry is needed for it.
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
        # Point MariaDB's on-disk temp-table dir at this job's own node-local
        # \$TMPDIR (== \$SCRATCH, see the `export TMPDIR` above) instead of the
        # my.cnf default of /tmp. The mysqld sidecar is started below without
        # --containall, so it shares the host mount namespace and its "/tmp"
        # is literally the compute node's real, SHARED /tmp -- when several
        # FUNANNOTATE_TRAIN tasks land on the same node at once, their MariaDB
        # instances all write large PASA-alignment MyISAM temp tables into
        # that same shared /tmp concurrently, and once it fills, mysqld fails
        # mid-write on a temp table and then fails again trying to delete the
        # file it never finished creating -- surfaces in PASA as
        # "Thread N terminated abnormally: ... Error on delete of
        # '/tmp/#sql_....MAI' (Errcode: 2 'No such file or directory')".
        # Confirmed 2026-08-29 against Alternaria_alternata_Z14's
        # pasa-assembly.log. \$MYSQL_TMP is bound below to a fixed in-container
        # path so this rewrite doesn't depend on \$TMPDIR's host-side value.
        MYSQL_TMP="\$TMPDIR/pasa_mysql_tmp_${asmid}"
        mkdir -p "\$MYSQL_TMP"
        sed -i "s#^tmpdir[[:space:]]*=.*#tmpdir\t\t= /pasa_mysql_tmp#" \$MYSQL_SCRATCH/conf/my.cnf
        # No separate SING_BINDS entry needed for \$MYSQL_SCRATCH (which the
        # containerized funannotate exec needs to read \$PASACONF from) --
        # it's now a subdirectory of \$TMPDIR, already covered by the
        # \$TMPDIR:\$TMPDIR bind above.
        # >/dev/null too, not just stderr: apptainer/singularity prints its
        # "Usage: ... instance stop ..." help text to STDOUT (not stderr) when
        # the named instance is already gone, which happens on every run here
        # -- the EXIT trap below always re-fires this after the explicit
        # stop_mysqldb call a few lines from the end of the script already
        # stopped it, so the second, redundant call is expected and silenced.
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} >/dev/null 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        # No trailing "/usr/bin/mysqld_safe" here: `instance start` only runs
        # trailing args THROUGH the image's %startscript (see --help), and
        # Docker's ENTRYPOINT/CMD populate %runscript, not %startscript --
        # mariadb_sif is a custom build (nextflow/docs/mariadb-10.3.9.def)
        # whose %startscript itself execs mysqld_safe, so it needs no args.
        # Confirmed 2026-08-25: passing mysqld_safe as a trailing arg to a
        # plain docker://mariadb:10.3.9 (or an `apptainer pull`-cached sif of
        # it, both lacking a real %startscript) starts an empty instance that
        # never launches mysqld at all -- PASA then fails to connect.
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf,\$MYSQL_TMP:/pasa_mysql_tmp \\
            ${params.mariadb_sif} mysqldb_${asmid}
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    # Shared-Trinity rows only ever reach this process with pasa_tier
    # 'stringent', 'relaxed', 'composite', or 'composite_fallback' --
    # FUNANNOTATE_RNASEQ.nf filters 'skip'-tier rows (ANI < 90% to the
    # representative, or unmeasured) out before FUNANNOTATE_TRAIN is even
    # invoked, routing them to ab-initio-only like a genuinely RNA-seq-less
    # strain. 'relaxed' rows are non-representative strains whose ANI to the
    # representative (90-97%) says the shared Trinity-GG transcripts --
    # assembled from a DIFFERENT strain's own RNA-seq (see RNASEQ_PREPARE) --
    # are real but somewhat diverged, so PASA's stringent defaults (95%
    # identity / 90% aligned) would reject most of it; PASA_TIER_ARGS relaxes
    # those specifically for this run. 'stringent' rows (ANI >= 97%, or the
    # representative itself at ANI=100) keep PASA's defaults untouched.
    # 'composite'/'composite_fallback' rows are hybrid-cross species whose
    # trinity_fa is a concatenation of Trinity-GG assemblies from OTHER species
    # (see nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md) -- two DIFFERENT
    # tiers, not one, per bioinformatics review 2026-08-28: 'composite' is
    # parent-matched evidence, which should align to ITS OWN subgenome copy at
    # ~99-100% identity (same lineage) -- a loose threshold there risks the
    # WRONG parent's transcript cross-mapping onto a subgenome it doesn't
    # belong to, merging homeologs into chimeric models, so this tier is
    # actually STRICTER than 'relaxed', not looser. 'composite_fallback' (no
    # parent-specific data, genus-wide instead) keeps the original loose
    # thresholds, since that path's divergence assumption is genuinely wide.
    PASA_TIER_ARGS=""
    if [ "${pasa_tier}" = "relaxed" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_shared_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_shared_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_shared_num_bp_splice}"
    elif [ "${pasa_tier}" = "composite" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_composite_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_composite_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_composite_num_bp_splice}"
    elif [ "${pasa_tier}" = "composite_fallback" ]; then
        PASA_TIER_ARGS="--pasa_min_avg_per_id ${params.pasa_composite_fallback_min_avg_per_id} --pasa_min_pct_aligned ${params.pasa_composite_fallback_min_pct_aligned} --pasa_num_bp_splice ${params.pasa_composite_fallback_num_bp_splice}"
    fi

    # ── Local, per-attempt working directory for the funannotate run itself ────
    # `-o` points at a fresh local dir instead of the persistent, retry-reused
    # training_target path -- funannotate's own --left_norm/--right_norm symlink
    # logic (SafeRemove + os.symlink in funannotate/train.py) collides with a
    # leftover normalize/*.norm.fq.gz from a prior attempt when `-o` is the
    # reused persistent dir (FileExistsError), and separately its dirname-based
    # skip-check can leave that symlink never created at all ("does not exist")
    # depending on how tmpdir's path relates to the read path. A fresh local dir
    # avoids both. Published into training_target below via rsync, which
    # preserves symlinks as-is rather than copying the data they point to.
    LOCAL_TRAIN="\$PWD/train_local"
    mkdir -p "\$LOCAL_TRAIN"

    # Resolve the read paths to their canonical location (rnaseq_reads/
    # rnaseq_data) before handing them to funannotate. Nextflow stages
    # r1/r2/se as symlinks directly in this task's own workdir (\$PWD) -- the
    # same directory \$LOCAL_TRAIN sits inside -- so passing them as-is makes
    # funannotate's dirname(tmpdir) != dirname(left_norm) check (train.py)
    # compare \$PWD to \$PWD and evaluate equal, which skips creating
    # normalize/*.norm.fq.gz entirely. Resolving first also means the symlink
    # funannotate does create (and rsync later publishes) points straight at
    # the canonical file, never at this ephemeral workdir.
    R1_REAL="${r1}"
    [ -e "\$R1_REAL" ] && R1_REAL="\$(readlink -f "\$R1_REAL")"
    R2_REAL="${r2}"
    [ -e "\$R2_REAL" ] && R2_REAL="\$(readlink -f "\$R2_REAL")"
    SE_REAL="${se}"
    [ -e "\$SE_REAL" ] && SE_REAL="\$(readlink -f "\$SE_REAL")"

    # Whole invocation wrapped in a group + tee so composite/composite_fallback-tier
    # failures can be inspected below without re-running anything -- \${PIPESTATUS[0]}
    # (not plain \$?, which after a pipe would report tee's exit code) captures the
    # group's real exit status, i.e. whichever funannotate train branch ran last.
    {
    if [ -s "${trinity_fa}" ]; then
        if [ -s "${r1}" ]; then
            echo "[INFO] Running funannotate train (PASA+PE) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            \$SING funannotate train -i "\$GENOME_IN" -o "\$LOCAL_TRAIN" \\
                --trinity ${trinity_fa} --left_norm "\$R1_REAL" --right_norm "\$R2_REAL" \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        elif [ -s "${se}" ]; then
            echo "[INFO] Running funannotate train (PASA+SE) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            \$SING funannotate train -i "\$GENOME_IN" -o "\$LOCAL_TRAIN" \\
                --trinity ${trinity_fa} --single_norm "\$SE_REAL" \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        else
            # No reads at all -- r1/se are present-but-empty (0-byte) placeholders,
            # not missing paths (Nextflow path inputs can't be null). Passing
            # --left_norm/--right_norm pointed at those empty files (the previous
            # behavior here) handed funannotate empty FASTQs to normalize instead
            # of omitting the flags entirely -- latent bug, fixed alongside the
            # composite-tier work that made this branch a live, common path
            # (hybrid-cross species' composite parent-transcript evidence: see
            # nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md) rather than a rare
            # edge case.
            echo "[INFO] Running funannotate train (PASA only, no reads) for ${out} using shared Trinity (pasa_tier=${pasa_tier})"
            \$SING funannotate train -i "\$GENOME_IN" -o "\$LOCAL_TRAIN" \\
                --trinity ${trinity_fa} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$PASA_TIER_ARGS \\
                \$pasa_db_arg
        fi
    elif [ -s "${r1}" ]; then
        echo "[INFO] Running funannotate train (full PE, no shared Trinity) for ${out}"
        \$SING funannotate train -i "\$GENOME_IN" -o "\$LOCAL_TRAIN" \\
            --left_norm "\$R1_REAL" --right_norm "\$R2_REAL" --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (full SE, no shared Trinity) for ${out}"
        \$SING funannotate train -i "\$GENOME_IN" -o "\$LOCAL_TRAIN" \\
            --single_norm "\$SE_REAL" --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi
    } 2>&1 | tee funannotate_train_capture.log
    TRAIN_STATUS=\${PIPESTATUS[0]}
    if [ "\$TRAIN_STATUS" -ne 0 ]; then
        # 'composite'/'composite_fallback' failures get one chance to degrade
        # gracefully instead of the usual hard-fail+retry: only when PASA itself
        # actually completed its alignment/assignment step (its own
        # "PASA assigned N transcripts to M loci" summary line appears in the
        # captured output -- same signal DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md's
        # still-unimplemented Option 1 proposes generally, applied here just for
        # this tier). Absence of that line means PASA (or something before it --
        # MySQL, OOM, a crash) never finished, which is exactly the kind of infra
        # failure that SHOULD keep retrying with more memory, not get silently
        # absorbed into "this hybrid isn't trainable" -- conflating the two would
        # mask real infra bugs across every hybrid strain. See "Should we degrade
        # gracefully on composite-tier failure?" in
        # nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md.
        if { [ "${pasa_tier}" = "composite" ] || [ "${pasa_tier}" = "composite_fallback" ]; } && \\
           grep -qE 'PASA assigned [0-9]+ transcripts to [0-9]+ loci' funannotate_train_capture.log; then
            echo "[WARN] ${out}: funannotate train failed (exit \$TRAIN_STATUS) but PASA completed alignment/assignment for this composite-tier (pasa_tier=${pasa_tier}) run -- treating as 'not enough usable transcript evidence' rather than an infra failure. Degrading to ab-initio-only; predict will proceed without PASA evidence for this strain." >&2
            mkdir -p "${params.training_target}/${out}/training"
            : > "${params.training_target}/${out}/training/.composite_train_failed"
            printf "out\\tspecies\\tpasa_tier\\texit_code\\ttimestamp\\n%s\\t%s\\t%s\\t%s\\t%s\\n" \\
                "${out}" "${species}" "${pasa_tier}" "\$TRAIN_STATUS" "\$(date -Iseconds)" \\
                > "${out}.composite_train_failed.tsv"
            exit 0
        fi
        echo "[ERROR] funannotate train failed for ${out} (exit \$TRAIN_STATUS); not publishing \$LOCAL_TRAIN over ${params.training_target}/${out}/training" >&2
        exit "\$TRAIN_STATUS"
    fi

    # ── Publish the local run into the persistent, retry-reused training dir ───
    # rsync preserves symlinks as-is rather than following them, so
    # normalize/{left,right,single}.norm.fq.gz -- which funannotate wrote as
    # symlinks to the realpath()'d canonical rnaseq_reads/rnaseq_data fastq --
    # arrive in TRAINDIR still pointing straight at rnaseq_reads/rnaseq_data,
    # never at \$LOCAL_TRAIN/work. trinity_gg/ and hisat2/ are large per-run
    # intermediates nothing downstream needs (predict reads the *.bam/*.gff3
    # files funannotate leaves alongside them, not these subdirs), so they're
    # excluded from the copy entirely rather than copied then deleted.
    # mysql_db/ is excluded too -- it was never written under \$LOCAL_TRAIN (the
    # sidecar above binds it straight from TRAINDIR) -- and --delete would
    # otherwise wipe it from TRAINDIR since it's absent from the local source.
    TRAINDIR="${params.training_target}/${out}/training"
    mkdir -p "${params.training_target}/${out}"
    echo "[INFO] Publishing local training run for ${out} to \$TRAINDIR"
    rsync -a --delete \\
        --exclude 'trinity_gg/' --exclude 'hisat2/' --exclude 'mysql_db/' \\
        "\$LOCAL_TRAIN/training/" "\$TRAINDIR/"
    rm -rf "\$LOCAL_TRAIN"

    # ── Make funannotate's convenience symlinks local ──────────────────────────
    # funannotate writes funannotate_train.{coordSorted.bam,transcripts.gff3,
    # trinity-GG.fasta} and transcript.alignments.bam as absolute paths into
    # whatever tree it ran under -- now the just-deleted \$LOCAL_TRAIN -- so every
    # publish leaves these needing a repoint, not just a moved/retired project.
    ${workflow.projectDir}/bin/relink_training_symlinks.py --apply --quiet "\$TRAINDIR" || true

    echo "[INFO] Training cleanup complete for ${out}"
    echo "mysql is ${params.pasa_mysql}"
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}
