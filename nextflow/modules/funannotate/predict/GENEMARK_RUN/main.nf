// Standalone GeneMark step, split out of funannotate predict's internal
// call so predict itself can move onto the rust container. GeneMark now runs
// from its own container (params.genemark_sif) instead of a host module.
// Default: docker://teambraker/braker3 (the Gaius-Augustus BRAKER team's
// public image), which happens to bundle GeneMark-ES/ET 4.72 -- same version
// this project's prior host module used -- alongside
// AUGUSTUS/bam2hints/join_mult_hints.pl. Confirmed working 2026-08-18 (all
// tools present on PATH, same license-key resolution as the earlier private
// hyphaltip/genemark-container build this replaced; see
// ENVIRONMENTS_INSTALLATIONS.md). GeneMark's own license still requires a
// user-obtained key (~/.gm_key) regardless of which image runs it -- the
// image being public just means the GeneMark *binary* ships in it; it does
// not grant a license. See nextflow/docs/GENEMARK_RUN_DESIGN.md for the full
// design and why this exists.
//
// Reuse mode (checked first, applies regardless of ES/ET -- gmes_petap.pl's
// --predict_with is its own mutually-exclusive run mode, not a variant of
// --ES/--ET): shared_mod set and force_independent != 'true' ->
// `gmes_petap.pl --predict_with <mod>` -- genuinely training-free prediction
// (verified in gmes_petap.pl source: mutually exclusive with --ES/--ET/--EP,
// no retraining at all). This is faster AND more correct than what
// funannotate predict's own internal reuse does today: `-p parameters.json`
// reuse only *seeds* --ini_mod into a full ES retrain (see RunGeneMarkES()
// in funannotate-live/funannotate/library.py) -- confirmed by
// analysis/funannotate_predict_stage_timing/ showing genemark_es_train_seconds
// stays 600-800s regardless of reuse eligibility. --predict_with was simply
// never wired up by funannotate's wrapper.
//
// Otherwise, fresh training, mode selected by `mode` + training_bam:
//   - mode == 'ET' and training_bam non-empty/exists: `gmes_petap.pl --ET`,
//     seeded by RNA-seq-informed intron hints derived from
//     training/transcript.alignments.bam (FUNANNOTATE_TRAIN's own PASA/
//     minimap2 transcript-to-genome alignment -- NOT the raw RNA-seq read
//     BAM). Real end-to-end validated 2026-08-12
//     (analysis/genemark_run_validation/et_eval2/, 10,776 genes): raw
//     bam2hints output is unstranded and must be run through
//     filterIntronsFindStrand.pl (vendored from BRAKER at
//     nextflow/bin/vendor/, Artistic License) to assign strand + drop
//     non-canonical splice sites, or GeneMark's branch-point training step
//     (bp_seq_select.pl) finds nothing usable and dies with "hash is empty".
//   - otherwise (mode == 'ES', or mode == 'ET' but no training_bam -- e.g. a
//     genome with no RNA-seq at all, see FUNANNOTATE_RNASEQ.nf's
//     predict_no_rnaseq branch): `gmes_petap.pl --ES` self-training,
//     identical to what predict's internal RunGeneMarkES() does today.
// Both fresh-training branches produce a new .mod (candidate for
// backfilling to the shared store by BACKFILL_ABINITIO_PARAMS -- this
// process never writes to the shared store itself; see design doc for why).
//
// --genemark_gtf expects GeneMark's raw native GTF (predict.py runs
// genemark_gtf2gff3.pl on it internally) -- no conversion needed here.
process GENEMARK_RUN {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '4h'

    input:
    tuple val(out), val(asmid), val(species), val(strain),
          val(genome_fa), val(transl_table), val(mode), val(training_bam),
          val(force_independent), val(shared_mod)

    output:
    tuple val(out), path("${out}.genemark.gtf"), emit: gtf
    // Keyed tuple, not a bare path: a batch can mix fresh-train rows (which
    // emit this) and fast-reuse rows (which don't, --predict_with produces
    // no new .mod) -- a bare `optional: true` path output would desync
    // positional alignment the moment that happens. See design doc point 1.
    tuple val(out), val(species), path("${out}.genemark.mod"), optional: true, emit: mod

    script:
    def filterIntronsFindStrand = "${workflow.projectDir}/bin/vendor/filterIntronsFindStrand.pl"
    """
    # ── Skip cleanly if a current GTF already exists (idempotent -resume) ───
    OUT_GTF="${out}.genemark.gtf"

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity
    # ── Containerized GeneMark ───────────────────────────────────────────────
    # gmes_petap.pl/bam2hints/join_mult_hints.pl are all already on \$PATH
    # inside the image, so no hardcoded GENEMARK_PATH/binary path is needed
    # (default image: docker://teambraker/braker3, whose PATH includes
    # /opt/ETP/bin/gmes and /opt/Augustus/bin -- see ENVIRONMENTS_INSTALLATIONS.md).
    # --bind /opt/linux:/opt/linux resolves the license key: gmes_petap.pl
    # reads ~/.gm_key, which on this host is \$HOME/.gm_key ->
    # /opt/linux/rocky/8.x/x86_64/pkgs/genemarkESET/4.72_lic/gm_key ->
    # gm_key_64 (a *relative* symlink resolved against /bigdata, which is
    # already auto-bound) -- without /opt/linux bound, that chain 404s inside
    # the container even though the final target under /bigdata is reachable.
    # \$HOME itself is auto-bound by singularity, so no separate bind is
    # needed for the symlink's starting point.
    SING="singularity exec --bind /opt/linux:/opt/linux,${workflow.projectDir}:${workflow.projectDir} ${params.genemark_sif}"

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) pigz -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac

    # ── Too-small/fragmented-genome pre-flight guard ─────────────────────────
    # Mirrors FUNANNOTATE_PREDICT's own guard (same params.predict_min_asm_bp/
    # predict_frag_max_n50/predict_frag_max_contigs, same bin/asm_preflight_stats.py
    # -- shared, not duplicated, see GENEMARK_RUN_DESIGN.md's "Known gap"
    # section) -- GENEMARK_RUN sits upstream of that guard in the DAG, so
    # without a check here, a genome predict would go on to skip anyway
    # (small AND fragmented) instead burns a full gmes_petap.pl --ES/--ET
    # attempt and hard-fails once GeneMark's own internal contig selection
    # (--min_contig 10000, after excluding soft-masked repeat sequence)
    # leaves too little usable training sequence, e.g. "error, input
    # sequence size is too small data/training.fna: 32259" for a
    # 1350-contig, 5.7Mb assembly. Skip cleanly here (empty GTF, no .mod) so
    # predict's own preflight guard is the one that actually flags/records
    # the skip in predict_skipped_too_small.tsv.
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" genome.fa \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs})
    if [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[WARN] GENEMARK_RUN ${out}: too small/fragmented (\$ASM_BP bp, \$ASM_CTG contigs, N50 \$ASM_N50); skipping GeneMark -- predict's own preflight guard will flag/report this genome" >&2
        touch "${out}.genemark.gtf"
        rm -f genome.fa
        exit 0
    fi

    # ── --gcode support probe (funannotate's _genemark_supports_gcode()) ────
    GCODE_ARGS=()
    if [ "${transl_table}" = "6" ] || [ "${transl_table}" = "26" ]; then
        if \$SING gmes_petap.pl 2>&1 | grep -q -- '--gcode'; then
            GCODE_ARGS=(--gcode "${transl_table}")
        else
            echo "[WARN] GeneMark does not support --gcode in this version; running with default code 1 -- gene calls in CUG/alt-table genomes may be unreliable" >&2
        fi
    elif [ "${transl_table}" != "1" ]; then
        echo "[WARN] GeneMark only supports --gcode 6 or 26; ignoring table ${transl_table} and running with default code 1" >&2
    fi

    # gmes_petap.pl's own "too small" abort (e.g. "error, input sequence size
    # is too small data/training.fna: 32259") prints to stdout, not a log
    # file, and gmes_petap.pl still exits 0 in that case -- only the missing
    # output/gmhmm.mod reveals the failure. Capture stdout alongside letting
    # it stream to the task log (tee), so the branches below can distinguish
    # "GeneMark declined because the input was too small/fragmented after its
    # own internal contig selection" (graceful skip, same outcome as the
    # pre-flight guard above) from a genuine unexpected failure (hard error).
    GMES_LOG="gmes_stdout.log"
    too_small_skip() {
        grep -qi "input sequence size is too small" "\$GMES_LOG" 2>/dev/null
    }

    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fast-reuse (--predict_with) against shared model ${shared_mod}"
        cp "${shared_mod}" genemark-shared.mod
        \$SING gmes_petap.pl --predict_with genemark-shared.mod \\
            --sequence genome.fa --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
    elif [ "${mode}" = "ET" ] && [ -n "${training_bam}" ] && [ -s "${training_bam}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fresh ET self-training seeded by RNA-seq intron hints from ${training_bam}"
        \$SING bam2hints --intronsonly --in="${training_bam}" --out=raw_hints.gff
        if [ ! -s raw_hints.gff ]; then
            echo "ERROR: bam2hints produced no intron hints from ${training_bam}" >&2
            exit 1
        fi
        \$SING perl "${filterIntronsFindStrand}" genome.fa raw_hints.gff --score > stranded_hints.gff
        if [ ! -s stranded_hints.gff ]; then
            echo "ERROR: filterIntronsFindStrand.pl dropped all hints (no canonical splice sites found) -- ${training_bam} may be too noisy for --ET; consider mode=ES" >&2
            exit 1
        fi
        sort -n -k4,4 stranded_hints.gff | sort -s -n -k5,5 | sort -s -k3,3 | sort -s -k1,1 \\
            | \$SING join_mult_hints.pl > genemark.intron-hints.gff
        \$SING gmes_petap.pl --ET genemark.intron-hints.gff --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        elif too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark-ET declined -- not enough usable (unmasked, >=10kb) training sequence after masking; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        else
            echo "ERROR: GeneMark-ET did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    else
        if [ "${mode}" = "ET" ]; then
            echo "[INFO] GENEMARK_RUN ${out}: mode=ET requested but no training_bam available (no RNA-seq for this genome) -- falling back to ES self-training"
        fi
        echo "[INFO] GENEMARK_RUN ${out}: fresh ES self-training (force_independent=${force_independent}, shared_mod='${shared_mod}')"
        \$SING gmes_petap.pl --ES --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        elif too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark-ES declined -- not enough usable (unmasked, >=10kb) training sequence after masking; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        else
            echo "ERROR: GeneMark-ES did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    fi

    if [ ! -s genemark.gtf ]; then
        if too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark declined -- not enough usable training sequence; skipping GeneMark for this genome" >&2
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        fi
        echo "ERROR: GeneMark did not produce genemark.gtf" >&2
        exit 1
    fi
    cp genemark.gtf "\$OUT_GTF"
    rm -f genome.fa

    # ── Also persist a fresh .mod to predict_misc/, like funannotate's own
    # internal GeneMark call used to (see GENEMARK_RUN_DESIGN.md's "Known
    # gap" section) -- this is Option B persistence, same pattern as
    # FUNANNOTATE_PREDICT writing directly into params.target/${out}, not a
    # Nextflow publishDir copy. Consumers of GENEMARK_RUN.out.mod (the wired
    # BACKFILL_ABINITIO_PARAMS join in FUNANNOTATE_PREDICTION.nf) are
    # unaffected and keep taking priority via their explicit channel value;
    # this only backstops consumers that read genemark.mod off disk instead
    # (workflows/backfill_abinitio.nf's standalone sweep,
    # species_reuse_clusters.py's own backfill) so they stop silently
    # missing genemark for anything predicted through this process.
    if [ -f "${out}.genemark.mod" ]; then
        LOWER_OUT=\$(printf '%s' "${out}" | tr '[:upper:]' '[:lower:]')
        AB_INITIO_DIR="${params.target}/${out}/predict_misc/ab_initio_parameters"
        mkdir -p "\$AB_INITIO_DIR"
        cp "${out}.genemark.mod" "\$AB_INITIO_DIR/\${LOWER_OUT}.genemark.mod"
    fi

    echo "[INFO] GENEMARK_RUN complete for ${out}"
    """

    stub:
    """
    touch "${out}.genemark.gtf"
    if [ -z "${shared_mod}" ] || [ "${force_independent}" = "true" ]; then
        touch "${out}.genemark.mod"
    fi
    """
}
