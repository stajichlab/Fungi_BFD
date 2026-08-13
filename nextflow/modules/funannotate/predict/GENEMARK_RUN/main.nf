// Standalone GeneMark step, split out of funannotate predict's internal
// call so predict itself can move onto the rust container (GeneMark's
// license forbids redistribution, so it can never be containerized and must
// stay on this host module). See nextflow/docs/GENEMARK_RUN_DESIGN.md for
// the full design and why this exists.
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
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate/dev-1.9

    if [ -z "\$GENEMARK_PATH" ] || [ ! -x "\$GENEMARK_PATH/gmes_petap.pl" ]; then
        echo "ERROR: GENEMARK_PATH not set or gmes_petap.pl not executable (\$GENEMARK_PATH)" >&2
        exit 1
    fi

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) pigz -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac

    # ── --gcode support probe (funannotate's _genemark_supports_gcode()) ────
    GCODE_ARGS=()
    if [ "${transl_table}" = "6" ] || [ "${transl_table}" = "26" ]; then
        if "\$GENEMARK_PATH/gmes_petap.pl" 2>&1 | grep -q -- '--gcode'; then
            GCODE_ARGS=(--gcode "${transl_table}")
        else
            echo "[WARN] GeneMark does not support --gcode in this version; running with default code 1 -- gene calls in CUG/alt-table genomes may be unreliable" >&2
        fi
    elif [ "${transl_table}" != "1" ]; then
        echo "[WARN] GeneMark only supports --gcode 6 or 26; ignoring table ${transl_table} and running with default code 1" >&2
    fi

    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fast-reuse (--predict_with) against shared model ${shared_mod}"
        cp "${shared_mod}" genemark-shared.mod
        "\$GENEMARK_PATH/gmes_petap.pl" --predict_with genemark-shared.mod \\
            --sequence genome.fa --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
    elif [ "${mode}" = "ET" ] && [ -n "${training_bam}" ] && [ -s "${training_bam}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fresh ET self-training seeded by RNA-seq intron hints from ${training_bam}"
        bam2hints --intronsonly --in="${training_bam}" --out=raw_hints.gff
        if [ ! -s raw_hints.gff ]; then
            echo "ERROR: bam2hints produced no intron hints from ${training_bam}" >&2
            exit 1
        fi
        perl "${filterIntronsFindStrand}" genome.fa raw_hints.gff --score > stranded_hints.gff
        if [ ! -s stranded_hints.gff ]; then
            echo "ERROR: filterIntronsFindStrand.pl dropped all hints (no canonical splice sites found) -- ${training_bam} may be too noisy for --ET; consider mode=ES" >&2
            exit 1
        fi
        sort -n -k4,4 stranded_hints.gff | sort -s -n -k5,5 | sort -s -k3,3 | sort -s -k1,1 \\
            | join_mult_hints.pl > genemark.intron-hints.gff
        "\$GENEMARK_PATH/gmes_petap.pl" --ET genemark.intron-hints.gff --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        else
            echo "ERROR: GeneMark-ET did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    else
        if [ "${mode}" = "ET" ]; then
            echo "[INFO] GENEMARK_RUN ${out}: mode=ET requested but no training_bam available (no RNA-seq for this genome) -- falling back to ES self-training"
        fi
        echo "[INFO] GENEMARK_RUN ${out}: fresh ES self-training (force_independent=${force_independent}, shared_mod='${shared_mod}')"
        "\$GENEMARK_PATH/gmes_petap.pl" --ES --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}"
        if [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        else
            echo "ERROR: GeneMark-ES did not produce output/gmhmm.mod" >&2
            exit 1
        fi
    fi

    if [ ! -s genemark.gtf ]; then
        echo "ERROR: GeneMark did not produce genemark.gtf" >&2
        exit 1
    fi
    cp genemark.gtf "\$OUT_GTF"
    rm -f genome.fa
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
