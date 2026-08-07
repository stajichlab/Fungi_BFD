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
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate/dev-1.9
    module load fastp

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    export TMPDIR=\${SCRATCH:-/tmp}

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

    if [ -s "${r1}" ]; then
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    else
        echo "[INFO] RNASEQ_PREPARE: using single-end reads for ${out}"
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --single_norm ${se} --aligners minimap2 \\
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
