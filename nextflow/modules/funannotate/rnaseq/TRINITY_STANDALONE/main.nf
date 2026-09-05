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

    SING_BINDS="--bind \$PWD:\$PWD,\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"

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
    else
        echo "[INFO] TRINITY_STANDALONE: single-end de novo Trinity for ${species_tag}"
        \$SING Trinity --seqType fq --no_normalize_reads \\
            --single ${se} \\
            --max_memory ${task.memory.toGiga()}G --CPU ${task.cpus} \\
            --output "\$OUTDIR" --full_cleanup \\
            > ${species_tag}.trinity-standalone.log 2>&1
    fi

    # --full_cleanup removes \$OUTDIR itself and leaves the assembly as \$OUTDIR.Trinity.fasta
    TRINITY_FA="\${OUTDIR}.Trinity.fasta"
    if [ -s "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-denovo.fasta
    else
        echo "[WARN] TRINITY_STANDALONE: no Trinity.fasta produced for ${species_tag}" >&2
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
