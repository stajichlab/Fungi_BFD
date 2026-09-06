// Run funannotate train on the representative (first) assembly of each species, then
// archive the Trinity-GG transcripts (normalized reads are in rnaseq_reads)
// reads into rnaseq_data/ so all other strains can skip those expensive steps.
// storeDir skips this process entirely if all five output files already exist.
process RNASEQ_PREPARE {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(species_tag), val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared
    path("${species_tag}.funannotate-trinity.log"), optional: true, emit: train_log

    script:
    // Real, symlink-resolved location of rnaseq_reads/ (itself a top-level
    // symlink into ../rnaseq_reads from launchDir) -- same missing-bind bug
    // class documented in FUNANNOTATE_TRAIN/main.nf: funannotate's
    // --left_norm/--right_norm handling resolves its read arguments to their
    // realpath() before symlinking them into its own scratch output tree, so
    // without this bind the container can create the symlink (doesn't
    // require the target to exist) but can never stat through it -- every
    // invocation fails with "Read normalization failed, .../left.norm.fq.gz
    // does not exist" even though the shell-level readlink -f below resolves
    // the path correctly. RNASEQ_PREPARE never got this fix when
    // FUNANNOTATE_TRAIN did (confirmed reproducing Albifimbria_verrucaria_Mv01
    // 2026-09-05).
    def rnaseqReadsDir = file("${launchDir}/rnaseq_reads").toRealPath()
    """
    # ── Empty-reads sentinel: no RNA-seq found by SRA_FETCH / SRA_FETCH_SE ──
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ]; then
        echo "[INFO] No RNAseq reads for ${species_tag}; writing empty shared markers"
        touch ${species_tag}.trinity-GG.fasta
        exit 0
    fi

    # ── If representative was already trained, just extract shared files ──────
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; extracting shared files to rnaseq_data"
        TRAINDIR="${params.training_target}/${out}/training"
        TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
        if [ -n "\$TRINITY_FA" ]; then
            cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
        else
            touch ${species_tag}.trinity-GG.fasta
        fi
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer

    # Containerized funannotate/Trinity: same funannotate-live image (rust-optimized
    # Trinity v2.16.1 + PASA/EVM/minimap2/hisat2/fastp) that FUNANNOTATE_TRAIN runs
    # through, instead of the host `module load funannotate/dev-1.9` build -- keeps
    # the genome-guided Trinity-GG assembly pinned to the same image/version as
    # every other Trinity invocation in the pipeline. See "which" note in
    # FUNANNOTATE_TRAIN/main.nf: harmless host-shell function leak, unset defensively.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    export APPTAINERENV_AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export APPTAINERENV_FUNANNOTATE_DB=${params.funannotate_db}
    export TMPDIR=\${SCRATCH:-/tmp}

    SING_BINDS="--bind \$PWD:\$PWD,${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},${rnaseqReadsDir}:${rnaseqReadsDir},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"

    # ── Run full funannotate train on the representative genome ───────────────
    # Use SCRATCH for the funannotate output dir so Trinity/HISAT2/normalize
    # intermediates land on fast local storage and don't consume project quota.
    echo "[INFO] RNASEQ_PREPARE: running funannotate train for representative ${out} (species: ${species_tag})"

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # Resolve r1/r2/se to their real, symlink-resolved location before handing
    # them to funannotate. Nextflow stages them as symlinks directly in this
    # task's own \$PWD -- passing them as-is can make funannotate's
    # dirname(tmpdir) != dirname(left_norm) check (train.py) compare \$PWD to
    # \$PWD and evaluate equal, which skips creating normalize/*.norm.fq.gz
    # entirely ("Read normalization failed, .../left.norm.fq.gz does not
    # exist") and silently produces a 0-transcript Trinity-GG assembly. Same
    # fix as FUNANNOTATE_TRAIN/main.nf; whether this bites depends on how the
    # task's \$PWD relates to the read path, so it doesn't reproduce every run.
    R1_REAL="${r1}"
    [ -e "\$R1_REAL" ] && R1_REAL="\$(readlink -f "\$R1_REAL")"
    R2_REAL="${r2}"
    [ -e "\$R2_REAL" ] && R2_REAL="\$(readlink -f "\$R2_REAL")"
    SE_REAL="${se}"
    [ -e "\$SE_REAL" ] && SE_REAL="\$(readlink -f "\$SE_REAL")"

    if [ -s "${r1}" ]; then
        \$SING funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --left_norm "\$R1_REAL" --right_norm "\$R2_REAL" --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    else
        echo "[INFO] RNASEQ_PREPARE: using single-end reads for ${out}"
        \$SING funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --single_norm "\$SE_REAL" --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    fi

    # ── Copy shared outputs to rnaseq_data/ ──────────────────────────────────
    TRAINDIR="\$SCRATCH/${out}/training"
    TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
    if [ -n "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
    else
        echo "[WARN] No trinity.fasta found under \$TRAINDIR for ${out}"
        touch ${species_tag}.trinity-GG.fasta
    fi

    # ── Preserve the funannotate train log before scratch is wiped ────────────
    # funannotate writes logfiles/funannotate-trinity.log under the scratch
    # output dir; that whole dir is removed below. Copy it into this task's
    # work dir, and into training_target/${out}/logfiles/ so it lands in the
    # same genome_annotation_training/<out>/logfiles/ schema that a normal
    # (non-representative) funannotate train run under FUNANNOTATE_TRAIN uses.
    TRAIN_LOG="\$SCRATCH/${out}/logfiles/funannotate-trinity.log"
    if [ -f "\$TRAIN_LOG" ]; then
        cp "\$TRAIN_LOG" ${species_tag}.funannotate-trinity.log
        mkdir -p "${params.training_target}/${out}/logfiles"
        cp "\$TRAIN_LOG" "${params.training_target}/${out}/logfiles/funannotate-trinity.log"
    else
        echo "[WARN] No funannotate-trinity.log found at \$TRAIN_LOG for ${out}"
    fi

    # ── Clean up scratch output dir (all intermediates were temporary) ────────
    rm -rf "\$SCRATCH/${out}"
    echo "[INFO] RNASEQ_PREPARE complete for ${species_tag}"
    """

    stub:
    """
    echo ">stub_trinity_${species_tag}" > ${species_tag}.trinity-GG.fasta
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}
