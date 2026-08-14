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
          val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    script:
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if training output already present and rnaseq is not newer than GBK ──
    # Accept a compressed prediction (.gbk.gz) as "done" so folders can be space-saved.
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk"
    [ -f "\$PREDICT_GBK" ] || PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk.gz"
    if [ -f "\$TRAIN_GFF3" ]; then
        RETRAIN=0
        # Re-train if the genome assembly itself is newer than the existing PASA training
        # output -- otherwise a swapped-in new assembly version silently keeps stale
        # training coordinates forever, since RETRAIN below only ever looks at rnaseq reads.
        if [ -s "${genome_fa}" ] && [ "${genome_fa}" -nt "\$TRAIN_GFF3" ]; then
            echo "[INFO] Genome assembly newer than training GFF3 for ${out}; retraining"
            rm -rf "${params.training_target}/${out}/training"
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
            echo "[INFO] Training already complete for ${out}; skipping"
            # Still repoint any absolute cross-project symlinks left by an older run.
            # Called by absolute path: this process runs via `bash -l`, and the login
            # shell's profile/module init rebuilds PATH, dropping the bin/ dir nextflow
            # otherwise appends for us.
            ${workflow.projectDir}/bin/relink_training_symlinks.py --apply --quiet \\
                "${params.training_target}/${out}/training" || true
            exit 0
        fi
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity
    # ── Containerized funannotate ──────────────────────────────────────────────
    # funannotate-live's Docker image (Rust-optimized PASA/EVM/Trinity + AVX2
    # bowtie2) run via `singularity exec \${params.funannotate_sif}`. Kept
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
    #   - Env passthrough: `singularity exec` passes through ALL host-shell env
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

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    export TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
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
        # Bind the WHOLE MYSQL_SCRATCH dir, not a subpath -- \$PASACONF
        # (MYSQL_SCRATCH/conf/pasa-local-*.config.txt) is what the containerized
        # funannotate exec needs to read, and MYSQL_SCRATCH/mysql_db (a prior
        # value here) was never created by anything above (only
        # MYSQL_SCRATCH/db and MYSQL_SCRATCH/conf are), so it was a dead bind
        # source pointing nowhere. Appended to SING_BINDS (not a separate
        # SINGULARITY_BINDPATH env var) so it's visible in the one \$SING exec
        # command below, matching this project's SING_BINDS/SING convention.
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

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    if [ -s "${trinity_fa}" ]; then
        if [ -s "${r1}" ]; then
            echo "[INFO] Running funannotate train (PASA+PE) for ${out} using shared Trinity"
            \$SING funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        elif [ -s "${se}" ]; then
            echo "[INFO] Running funannotate train (PASA+SE) for ${out} using shared Trinity"
            \$SING funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --single_norm ${se} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        else
            echo "[INFO] Running funannotate train (PASA only, no reads) for ${out} using shared Trinity"
            \$SING funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        fi
    elif [ -s "${r1}" ]; then
        echo "[INFO] Running funannotate train (full PE, no shared Trinity) for ${out}"
        \$SING funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (full SE, no shared Trinity) for ${out}"
        \$SING funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi

    # ── Remove large intermediates not needed for predict or update ─────────────
    # Keeps: *.bam, *.bai, *.pasa.gff3, *.stringtie.gtf, *.transcripts.gff3
    TRAINDIR="${params.training_target}/${out}/training"
    echo "[INFO] Removing large training intermediates in \$TRAINDIR"
    rm -rf "\$TRAINDIR/hisat2"
    rm -rf "\$TRAINDIR/trinity_gg"

    # ── Make funannotate's convenience symlinks local ──────────────────────────
    # funannotate writes funannotate_train.{coordSorted.bam,transcripts.gff3,
    # trinity-GG.fasta} and transcript.alignments.bam as absolute paths into
    # whatever tree it ran under, which breaks the moment the directory is moved
    # or the source project is retired. Repoint them at the sibling files.
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
