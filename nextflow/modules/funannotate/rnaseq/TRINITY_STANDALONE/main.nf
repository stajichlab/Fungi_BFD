// Fallback for species whose genome-guided Trinity (RNASEQ_PREPARE) collapsed to too
// few transcripts (train_min_trinity_transcripts guard, checked by
// COUNT_TRINITY_TRANSCRIPTS in FUNANNOTATE_RNASEQ) despite having substantial
// normalized RNA-seq -- usually because the representative genome is a poor match for
// the actual reads. Runs a de novo (non-genome-guided) Trinity assembly on the same
// already-normalized reads instead of trying a different reference genome (that fix is
// rnaseq_representative_override.csv / scripts/pick_rnaseq_representative_override.py).
//
// Publishes BOTH the real output (*.trinity-denovo.fasta) and a relative symlink named
// *.trinity-GG.fasta pointing at it, so every downstream consumer (FUNANNOTATE_TRAIN,
// BUILD_HYBRID_COMPOSITE_TRINITY) keeps reading the same trinity-GG.fasta path with no
// changes needed -- they just get a better assembly at it. Hybrid crosses need no
// special-casing: BUILD_HYBRID_COMPOSITE_TRINITY reads each parent's own trinity-GG.fasta,
// so fixing a parent species here fixes any hybrid composite built from it once rebuilt.
process TRINITY_STANDALONE {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    input:
    tuple val(species_tag), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared
    path("${species_tag}.trinity-denovo.fasta"), emit: denovo
    path("${species_tag}.trinity-standalone.log"), optional: true, emit: log

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer

    # Containerized Trinity: same funannotate-live image (rust-optimized Trinity
    # v2.16.1) that FUNANNOTATE_TRAIN runs PASA/EVM through, instead of the
    # host `module load funannotate/dev-1.9` build -- keeps every Trinity
    # invocation in the pipeline pinned to one image/version.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    export TMPDIR=\${SCRATCH:-/tmp}
    OUTDIR="\${SCRATCH:?}/trinity_${species_tag}"

    # r1/r2/se are staged by Nextflow as symlinks into a storeDir process's cache
    # (rnaseq_reads/, see SRA_FETCH*/WRITE_EMPTY_READS), i.e. OUTSIDE this task's own
    # workDir -- \$PWD:\$PWD alone does not make the symlink TARGETS visible inside the
    # container. Without this bind, Trinity's own create_full_path() can't stat the
    # (broken-from-its-perspective) symlink and dies with "cannot locate file:
    # <species>_norm_R1.fastq.gz" even though the file is right there on the host.
    # Confirmed 2026-09-05 against Microsporum_canis.
    #
    # ${launchDir}/rnaseq_reads alone is NOT enough: storeDir outputs are keyed by
    # the *launchDir of the run that produced them*, which need not match this run's
    # launchDir when Nextflow is invoked from a subdirectory (e.g. reads produced by
    # a run launched from the parent dir, this run launched from a child dir with its
    # own launchDir/rnaseq_reads being just a symlink back up -- the symlink resolves
    # fine on the host but the container never sees the real path unless it's bound
    # too). Confirmed 2026-09-06 against Albifimbria_verrucaria: r1/r2 resolved to
    # <parent>/rnaseq_reads/..., one level above this run's launchDir/rnaseq_reads.
    # So resolve each staged read file's real (symlink-followed) parent dir at
    # runtime and bind that too, in addition to the launchDir-relative guess above.
    EXTRA_BINDS=""
    for f in "${r1}" "${r2}" "${se}"; do
        if [ -e "\$f" ]; then
            realdir=\$(dirname "\$(readlink -f "\$f")")
            case ",\$EXTRA_BINDS," in
                *",\$realdir,"*) ;;
                *) EXTRA_BINDS="\${EXTRA_BINDS:+\$EXTRA_BINDS,}\$realdir:\$realdir" ;;
            esac
        fi
    done
    SING_BINDS="--bind \$PWD:\$PWD,${launchDir}/rnaseq_reads:${launchDir}/rnaseq_reads\${EXTRA_BINDS:+,\$EXTRA_BINDS},\$TMPDIR:\$TMPDIR"
    # --cleanenv: without it, apptainer passes the submitting shell's env straight
    # through, including whatever `module load java/...` set on the host (JAVA_HOME,
    # LD_LIBRARY_PATH pointing at that Java's lib/server/). That leaks the host JVM's
    # shared library over the container's own bundled JDK, producing a CDS
    # ("shared archive file version 0x12 does not match required version 0x13")
    # mismatch on EVERY java invocation -- and Trinity's Butterfly phase launches one
    # JVM per assembled component, so this isn't just log noise, it's real per-task
    # overhead multiplied hundreds of thousands of times. Confirmed 2026-09-07 against
    # Microsporum_canis (job 28174732: 19+ hours, 84% through Butterfly, 420MB log of
    # nothing but these warnings). --env TMPDIR=\$TMPDIR keeps the one host var Trinity
    # actually needs (java.io.tmpdir / scratch fallback) despite --cleanenv stripping it.
    SING="apptainer exec --cleanenv --env TMPDIR=\$TMPDIR \${SING_BINDS} ${params.funannotate_sif}"

    # Reads in rnaseq_reads/ are already normalized (in-silico read normalization by
    # RNASEQ_PREPARE's funannotate train --stop_after_trinity run), so skip Trinity's
    # own normalization step.
    if [ -s "${r1}" ]; then
        echo "[INFO] TRINITY_STANDALONE: paired-end de novo Trinity for ${species_tag}"
        \$SING Trinity --seqType fq --no_normalize_reads \\
            --left ${r1} --right ${r2} \\
            --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} \\
            --output "\$OUTDIR" --full_cleanup \\
            > ${species_tag}.trinity-standalone.log 2>&1
        TRINITY_EXIT=\$?
    else
        echo "[INFO] TRINITY_STANDALONE: single-end de novo Trinity for ${species_tag}"
        \$SING Trinity --seqType fq --no_normalize_reads \\
            --single ${se} \\
            --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} \\
            --output "\$OUTDIR" --full_cleanup \\
            > ${species_tag}.trinity-standalone.log 2>&1
        TRINITY_EXIT=\$?
    fi

    # Trinity itself crashing (bad bind, OOM, corrupt input, ...) must NOT be masked as
    # success: this process's outputs are storeDir-cached, so an exit-0 task with a
    # silently-empty trinity-denovo.fasta would be accepted as "done" forever and never
    # retried -- every downstream consumer (FUNANNOTATE_TRAIN, BUILD_HYBRID_COMPOSITE_TRINITY)
    # would then permanently fall back to their no-shared-Trinity path for this species.
    # Confirmed 2026-09-05: exactly this happened for Microsporum_canis via the missing
    # bind above. Only a genuine exit 0 with no assembly (Trinity ran fine but truly
    # assembled nothing) is treated as a legitimate empty result below.
    if [ "\$TRINITY_EXIT" -ne 0 ]; then
        echo "[ERROR] TRINITY_STANDALONE: Trinity exited \$TRINITY_EXIT for ${species_tag}; not caching an empty result -- see ${species_tag}.trinity-standalone.log" >&2
        tail -n 60 ${species_tag}.trinity-standalone.log >&2 || true
        exit "\$TRINITY_EXIT"
    fi

    # --full_cleanup removes \$OUTDIR itself and leaves the assembly as \$OUTDIR.Trinity.fasta
    TRINITY_FA="\${OUTDIR}.Trinity.fasta"
    if [ -s "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-denovo.fasta
    else
        echo "[WARN] TRINITY_STANDALONE: Trinity exited 0 but no Trinity.fasta produced for ${species_tag} (no assemblable transcripts)" >&2
        touch ${species_tag}.trinity-denovo.fasta
    fi

    ln -sf "${species_tag}.trinity-denovo.fasta" "${species_tag}.trinity-GG.fasta"

    rm -f "\$TRINITY_FA" "\${TRINITY_FA}.gene_trans_map"
    echo "[INFO] TRINITY_STANDALONE complete for ${species_tag}"
    """

    stub:
    """
    echo ">stub_denovo_${species_tag}" > ${species_tag}.trinity-denovo.fasta
    ln -sf ${species_tag}.trinity-denovo.fasta ${species_tag}.trinity-GG.fasta
    """
}
