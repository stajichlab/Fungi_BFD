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
// Every empty-GTF skip (pre-flight small/fragmented, or GeneMark's own
// post-attempt "too small" decline in ES/ET/fast-reuse mode) is recorded to
// ${target}/genemark_skipped_too_small.tsv (added 2026-08-30) -- without
// this, the only trace of a skip was a per-task stderr [WARN], which
// `cleanup = true` deletes along with work/ after any successful pipeline
// run, leaving no durable way to notice a genome silently got no GeneMark
// evidence (FUNANNOTATE_PREDICT falls back to --auto-skip-genemark for an
// empty genemark_gtf with no further signal). Mirrors
// FUNANNOTATE_PREDICT's own predict_skipped_too_small.tsv pattern.
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
//
// ── DRAFT: scatter/gather route (train once, predict in chunks) ─────────────
// This is GENEMARK_RUN's old process body, renamed GENEMARK_TRAIN, with one added
// route. When params.genemark_chunked is true and the genome has at least
// params.genemark_chunk_min_contigs contigs (ASM_CTG from asm_preflight_stats.py),
// and the row is ES or fast-reuse (ET is NOT chunked yet: ET equivalence is untested):
//   fresh ES  : `gmes_petap.pl --ES --training_only` (writes ./gmhmm.mod, no GTF), then
//   reuse     : no training; the stored shared_mod is the model
// and the task then writes chunks/chunk_NN.fa.gz (bin/genemark_split_fasta.py) and stops.
// GENEMARK_PREDICT_CHUNK predicts each chunk with `--predict_with <model>`; GENEMARK_MERGE
// merges, validates and publishes. Every other row (small genomes, preflight skips,
// "too small" declines, Microsporidia, ET) takes exactly the old monolithic path and
// emits exactly the old outputs. In the chunked route this task publishes nothing and
// emits NO `${out}.genemark.gtf`, `${out}.genemark.mod` or `${out}.other.gff3` -- those
// names exist only on the monolithic route, so genemark_store never holds a partial result.
// Measured basis (2026-10-12): same-model chunked == monolithic byte for byte on 3 genomes;
// ES training itself is not deterministic, so the published .mod is the reference.
process GENEMARK_TRAIN {
    tag "$out"
    // Keep every result in genemark_store/<out>/ so later runs reuse it instead of
    // re-running GeneMark (utils.nf storedGenemarkFor; DECISIONS D135). Not under
    // predict_misc/, which FUNANNOTATE_PREDICT deletes on a stale re-predict.
    // Publishes only the monolithic-route results (GTF, other.gff3 and (fresh training) .mod);
    // the chunked route's files (chunks/, *.train.mod, *.reuse.mod, *.other.pending.gff3) are
    // never published here -- GENEMARK_MERGE publishes after validating the merged GTF.
    publishDir path: { "${params.genemark_store}/${out}" }, mode: 'copy', overwrite: true,
        saveAs: { fn -> (fn ==~ /.*\.genemark\.(gtf|mod)$/ || fn ==~ /.*\.other\.gff3$/) ? fn : null },
        enabled: params.genemark_store ? true : false

    cpus   16
    memory '32 GB'
    time   '4h'

    input:
    tuple val(out), val(asmid), val(species), val(strain),
          val(genome_fa), val(transl_table), val(mode), val(training_bam),
          val(force_independent), val(shared_mod), val(is_microsporidia)

    output:
    // Monolithic route (and every skip route): same files as the old GENEMARK_RUN, now keyed
    // optional outputs because the chunked route emits none of them. other.gff3 is touched
    // before every skip/exit, so the monolithic route always has gtf AND other_gff.
    tuple val(out), path("${out}.genemark.gtf"), optional: true, emit: gtf
    // Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    // empty file when is_microsporidia=false or this genome didn't qualify.
    tuple val(out), path("${out}.other.gff3"), optional: true, emit: other_gff
    // Keyed tuple, not a bare path: a batch can mix fresh-train rows (which
    // emit this) and fast-reuse rows (which don't, --predict_with produces
    // no new .mod) -- a bare `optional: true` path output would desync
    // positional alignment the moment that happens. See design doc point 1.
    tuple val(out), val(species), path("${out}.genemark.mod"), optional: true, emit: mod
    // Chunked route only.
    tuple val(out), val(species), path("${out}.train.mod"), optional: true, emit: train_mod     // fresh ES model
    tuple val(out), path("${out}.reuse.mod"), optional: true, emit: reuse_mod                   // stored shared model
    tuple val(out), path("${out}.other.pending.gff3"), optional: true, emit: other_pending
    tuple val(out), path("chunks/chunk_*.fa.gz"), optional: true, emit: chunks

    script:
    // gmes_petap.pl (CheckGcode) accepts only --gcode 1, 6 and 29; table 26 gives
    // "error, input genetic code 26 is not supported" (exit 1). Table 26 (CUG=Ala)
    // keeps the standard stop codons, so GeneMark runs with its default code 1.
    def gm_table = (transl_table.toString() == '26') ? '1' : transl_table
    // prodigal 2.6.3 accepts only tables 1-6, 9-16 and 21-25 (tested 0-31: others exit 15,
    // "Invalid translation table specified."). Any other table (e.g. 26) falls back to table 1.
    def prodigal_table = (transl_table.toString() in ['1','2','3','4','5','6','9','10','11','12','13','14','15','16','21','22','23','24','25']) ? transl_table.toString() : '1'
    def filterIntronsFindStrand = "${workflow.projectDir}/bin/vendor/filterIntronsFindStrand.pl"
    """
    # ── Skip cleanly if a current GTF already exists (idempotent -resume) ───
    OUT_GTF="${out}.genemark.gtf"

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer
    # ── Containerized GeneMark ───────────────────────────────────────────────
    # gmes_petap.pl/bam2hints/join_mult_hints.pl are all already on \$PATH
    # inside the image, so no hardcoded GENEMARK_PATH/binary path is needed
    # (default image: docker://teambraker/braker3, whose PATH includes
    # /opt/ETP/bin/gmes and /opt/Augustus/bin -- see ENVIRONMENTS_INSTALLATIONS.md).
    #
    # IMPORTANT: Nextflow's apptainer.autoMounts/runOptions ONLY apply when
    # Nextflow itself launches the container via a process `container =`
    # directive. This process builds \$SING manually in the script body, so
    # NONE of that automatic binding applies here -- every path the
    # container needs to read/write must be bound explicitly below.
    # Confirmed empirically (2026-08-24): without the task workdir bound,
    # genome.fa is invisible inside the container even though `apptainer
    # exec ... pwd` reports the correct path (the string matches, but the
    # filesystem isn't mounted). Built from Nextflow context vars and
    # dynamic resolution rather than a hardcoded filesystem prefix (e.g.
    # /bigdata) so this stays portable across HPC systems.
    export TMPDIR=\${SCRATCH:-/tmp}

    # GeneMark's license resolution (gmes_petap.pl reads ~/.gm_key) can hop
    # through host-specific symlinks that land outside \$HOME/\$PWD (e.g. a
    # module-installed license landing under this cluster's /opt/linux
    # package tree) -- resolve the real target dynamically and bind its
    # parent dir instead of hardcoding a host-specific prefix.
    GM_KEY_REAL=\$(readlink -f "\$HOME/.gm_key" 2>/dev/null || echo "\$HOME/.gm_key")
    GM_KEY_DIR=\$(dirname "\$GM_KEY_REAL")

    SING_BINDS="--bind \$PWD:\$PWD,${workflow.projectDir}:${workflow.projectDir},\$GM_KEY_DIR:\$GM_KEY_DIR,\$TMPDIR:\$TMPDIR"
    # training_bam (ET mode only) lives outside the task workdir -- it's
    # read directly by the containerized bam2hints call below, not staged
    # into \$PWD first, so its directory needs its own bind.
    TRAINING_BAM="${training_bam}"
    if [ -n "\$TRAINING_BAM" ]; then
        SING_BINDS="\$SING_BINDS,\$(dirname "\$TRAINING_BAM"):\$(dirname "\$TRAINING_BAM")"
    fi
    SING="apptainer exec \$SING_BINDS ${params.genemark_sif}"

    GENOME_GZ="${genome_fa}"
    case "\$GENOME_GZ" in
        *.gz) pigz -dc "\$GENOME_GZ" > genome.fa ;;
        *)    cp "\$GENOME_GZ" genome.fa ;;
    esac

    # `other_gff` is a non-optional output, so it must exist on every exit
    # path -- including the small_fragmented pre-flight skip below -- not
    # just the Microsporidia Prodigal supplement path further down. Touch it
    # unconditionally, right after genome.fa is inflated, before any guard
    # that can `exit 0`.
    touch "${out}.other.gff3"

    # The production shell is /bin/bash -l (nextflow.config process.shell: no `set -e`), and the gtf/other_gff/
    # chunks outputs are all optional, so a task that exits 0 with the wrong files would be COMPLETED and the
    # genome dropped silently. On any exit 0, require the outputs of exactly one route:
    #   monolithic / skip route : <out>.genemark.gtf AND <out>.other.gff3
    #   chunked route           : chunks/chunk_*.fa.gz AND <out>.other.pending.gff3 AND a non-empty model
    route_check() {
        local rc=\$?
        if [ "\$rc" -eq 0 ]; then
            if [ -e "${out}.genemark.gtf" ] && [ -e "${out}.other.gff3" ]; then return 0; fi
            if ls chunks/chunk_*.fa.gz > /dev/null 2>&1 && [ -e "${out}.other.pending.gff3" ] \\
                    && { [ -s "${out}.train.mod" ] || [ -s "${out}.reuse.mod" ]; }; then return 0; fi
            echo "ERROR: GENEMARK_TRAIN ${out}: task ended with exit 0 but its outputs match neither route" >&2
            exit 1
        fi
    }
    trap route_check EXIT

    # ── Persistent audit trail for empty/skipped GeneMark GTFs ──────────────
    # GeneMark's own "too small" decline (both the pre-flight guard below and
    # the post-attempt too_small_skip() cases further down) previously only
    # logged a per-task stderr [WARN] -- with `cleanup = true` deleting work/
    # after any successful pipeline run, that WARN and the empty GTF both
    # vanish permanently, leaving no way to notice or audit which genomes
    # silently got no GeneMark evidence (FUNANNOTATE_PREDICT falls back to
    # --auto-skip-genemark for an empty genemark_gtf with no further signal).
    # Mirrors FUNANNOTATE_PREDICT's own predict_skipped_too_small.tsv pattern
    # (same ${params.target} dir, distinct filename since this fires earlier
    # in the DAG and for a different reason: GeneMark's own internal contig
    # selection found too little usable sequence, not this module's own
    # pre-flight stats necessarily agreeing).
    GENEMARK_SKIP_REPORT="${params.target}/genemark_skipped_too_small.tsv"
    record_genemark_skip() {
        local reason="\$1" asm_bp="\${2:-NA}" asm_ctg="\${3:-NA}" asm_n50="\${4:-NA}"
        mkdir -p "${params.target}"
        [ -s "\$GENEMARK_SKIP_REPORT" ] || printf 'out\tasmid\treason\ttotal_bp\tcontigs\tN50\n' > "\$GENEMARK_SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "\$reason" "\$asm_bp" "\$asm_ctg" "\$asm_n50" >> "\$GENEMARK_SKIP_REPORT"
    }

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
            --max-contigs ${params.predict_frag_max_contigs} \\
            --min-contig-len ${params.predict_min_training_contig_len} \\
            --min-training-contigs ${params.predict_min_training_contigs} \\
            --abs-min-bp ${params.predict_abs_min_asm_bp})
    if [ "\$ASM_VERDICT" != "ok" ]; then
        echo "[WARN] GENEMARK_RUN ${out}: preflight verdict '\$ASM_VERDICT' (\$ASM_BP bp, \$ASM_CTG contigs, N50 \$ASM_N50); skipping GeneMark -- predict's own preflight guard will flag/report this genome" >&2
        record_genemark_skip "preflight_\$ASM_VERDICT" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50"
        touch "${out}.genemark.gtf"
        rm -f genome.fa
        exit 0
    fi

    # ── Microsporidia Prodigal supplement ────────────────────────────────────
    # Gated on PHYLUM (is_microsporidia, computed once in Groovy from
    # loadMicrosporidiaOutSet()) AND total assembled bp -- NOT the
    # small_fragmented composite verdict above, which also requires
    # fragmentation and would miss the genome this recipe was validated on
    # (Ordospora colligata OC4: 2.29 Mb, 15 contigs, N50 228,601 -- small but
    # not fragmented). Runs alongside, not instead of, the GeneMark attempt
    # below: the validated recipe used genemark:1 AND prodigal:5 together --
    # if GeneMark independently fails on this genome, the existing empty-GTF
    # degradation further down already absorbs that; this block does not
    # change GeneMark's own success/failure handling at all.
    # (${out}.other.gff3 is already touched above, right after genome.fa was
    # inflated, so every exit path -- including the small_fragmented
    # pre-flight guard -- has it.)
    if [ "${is_microsporidia}" = "true" ] && [ "\$ASM_BP" -lt "${params.predict_min_asm_bp}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: microsporidia + small (\$ASM_BP bp < ${params.predict_min_asm_bp}); running Prodigal supplement"
        module load prodigal
        if [ "${prodigal_table}" != "${transl_table}" ]; then
            echo "[WARN] GENEMARK_RUN ${out}: prodigal does not support translation table ${transl_table}; using table ${prodigal_table}" >&2
        fi
        prodigal -i genome.fa -o "${out}.prodigal.raw.gff3" -f gff -p single -g "${prodigal_table}"
        python "${workflow.projectDir}/bin/prodigal_hierarchy.py" \\
            "${out}.prodigal.raw.gff3" "${out}.other.gff3"
        PRODIGAL_REPORT="${params.target}/prodigal_selected_microsporidia.tsv"
        mkdir -p "${params.target}"
        [ -s "\$PRODIGAL_REPORT" ] || printf 'out\tasmid\treason\ttotal_bp\n' > "\$PRODIGAL_REPORT"
        printf '%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "microsporidia_phylum_and_small" "\$ASM_BP" >> "\$PRODIGAL_REPORT"
    fi

    # ── --gcode support probe (funannotate's _genemark_supports_gcode()) ────
    GCODE_ARGS=()
    if [ "${gm_table}" = "6" ] || [ "${gm_table}" = "26" ]; then
        if \$SING gmes_petap.pl 2>&1 | grep -q -- '--gcode'; then
            GCODE_ARGS=(--gcode "${gm_table}")
        else
            echo "[WARN] GeneMark does not support --gcode in this version; running with default code 1 -- gene calls in CUG/alt-table genomes may be unreliable" >&2
        fi
    elif [ "${gm_table}" != "1" ]; then
        echo "[WARN] GeneMark only supports --gcode 1, 6 and 29; ignoring table ${transl_table} and running with default code 1" >&2
    fi

    # gmes_petap.pl's own "too small" abort (e.g. "error, input sequence size
    # is too small data/training.fna: 32259") prints to stdout, not a log
    # file; gmes_petap.pl exits 1 in that case (gmes_petap.pl:1604) but the pipe to tee hides the status -- only the missing
    # output/gmhmm.mod reveals the failure. Capture stdout alongside letting
    # it stream to the task log (tee), so the branches below can distinguish
    # "GeneMark declined because the input was too small/fragmented after its
    # own internal contig selection" (graceful skip, same outcome as the
    # pre-flight guard above) from a genuine unexpected failure (hard error).
    GMES_LOG="gmes_stdout.log"
    too_small_skip() {
        grep -qi "input sequence size is too small" "\$GMES_LOG" 2>/dev/null
    }

    # ── Scatter/gather decision (DRAFT; off unless params.genemark_chunked) ──────
    # Chunk only large, fragmented genomes (contig count from the pre-flight stats above) and
    # only ES / fast-reuse rows. ET with a training BAM stays monolithic: the ET final
    # prediction should not use the hints, so chunking is expected to be safe, but that is untested.
    CHUNKED=false
    if [ "${params.genemark_chunked}" = "true" ] && [ "\$ASM_CTG" -ge ${params.genemark_chunk_min_contigs} ]; then
        if [ "${mode}" = "ET" ] && [ -n "${training_bam}" ] && [ -s "${training_bam}" ] \\
                && { [ -z "${shared_mod}" ] || [ "${force_independent}" = "true" ]; }; then
            echo "[INFO] GENEMARK_RUN ${out}: \$ASM_CTG contigs but fresh ET mode -- not chunked (ET equivalence untested)"
        else
            CHUNKED=true
            echo "[INFO] GENEMARK_RUN ${out}: \$ASM_CTG contigs >= ${params.genemark_chunk_min_contigs} -- scatter/gather route"
        fi
    fi
    # Durable record of which genomes took the chunked route (the work dir is deleted by cleanup = true).
    CHUNKED_REPORT="${params.target}/genemark_chunked.tsv"
    # One row per chunk: contigs and bases in the chunk, so a slow chunk can be matched to its size
    # (trace tag is <out>:<chunk>; chunks are cut by contig count only).
    record_genemark_chunked() {
        local route="\$1" n_chunks="\$2" genome_contigs="\$3"
        mkdir -p "${params.target}"
        [ -s "\$CHUNKED_REPORT" ] || printf 'out\tasmid\troute\tn_chunks\tchunk\tcontigs\tbp\tgenome_contigs\n' > "\$CHUNKED_REPORT"
        awk -F'\\t' -v out="${out}" -v asm="${asmid}" -v route="\$route" -v n="\$n_chunks" -v g="\$genome_contigs" \\
            '{ idx = \$1; sub(/^chunk_/, "", idx); c = \$2; sub(/^contigs=/, "", c); b = \$3; sub(/^bp=/, "", b);
               printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", out, asm, route, n, idx, c, b, g }' \\
            chunks/chunks.manifest.txt >> "\$CHUNKED_REPORT"
    }
    # finish_chunked <model file> <fresh|reuse>
    # Split genome.fa into chunks/chunk_NNNN.fa.gz and leave the monolithic output names unset, so only
    # GENEMARK_MERGE publishes. Every step is checked explicitly (no `set -e` here, and a pipe to tee would hide
    # the split script's exit status): the chunks must hold exactly the contigs of genome.fa.
    finish_chunked() {
        local model="\$1" route="\$2"
        [ -s "\$model" ] || { echo "ERROR: GENEMARK_TRAIN ${out}: model \$model is missing or empty" >&2; exit 1; }
        local nctg; nctg=\$(grep -c '^>' genome.fa)
        [ "\$nctg" -gt 0 ] || { echo "ERROR: GENEMARK_TRAIN ${out}: genome.fa has no contigs" >&2; exit 1; }
        mkdir -p chunks
        python "${workflow.projectDir}/bin/genemark_split_fasta.py" --input genome.fa --outdir chunks \\
            --contigs-per-chunk ${params.genemark_chunk_contigs} --min-chunk-bp ${params.genemark_chunk_min_bp} \\
            > chunks/chunks.manifest.txt || { echo "ERROR: GENEMARK_TRAIN ${out}: genome split failed" >&2; exit 1; }
        cat chunks/chunks.manifest.txt
        local nchunk; nchunk=\$(ls chunks/chunk_*.fa.gz | wc -l)
        [ "\$nchunk" -ge 1 ] || { echo "ERROR: GENEMARK_TRAIN ${out}: no chunk files written" >&2; exit 1; }
        local wrote; wrote=\$(for f in chunks/chunk_*.fa.gz; do pigz -dc "\$f"; done | grep -c '^>')
        [ "\$wrote" -eq "\$nctg" ] || { echo "ERROR: GENEMARK_TRAIN ${out}: chunks hold \$wrote of \$nctg contigs" >&2; exit 1; }
        mv "${out}.other.gff3" "${out}.other.pending.gff3" || exit 1
        record_genemark_chunked "\$route" "\$nchunk" "\$nctg"
        rm -f genome.fa
        echo "[INFO] GENEMARK_RUN ${out}: chunked route ready (\$nchunk chunks, \$nctg contigs)"
    }

    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: fast-reuse (--predict_with) against model ${shared_mod}"
        if [ "\$CHUNKED" = "true" ]; then
            cp "${shared_mod}" "${out}.reuse.mod" || { echo "ERROR: cannot copy ${shared_mod}" >&2; exit 1; }
            finish_chunked "${out}.reuse.mod" reuse
            exit 0
        fi
        cp "${shared_mod}" genemark-shared.mod
        # Same --max_intron/--soft_mask as fresh --ES: with them, --predict_with the
        # genome's own ES model gave 100% identical CDS chains to the ES run on 6 wave 1
        # genomes; without them, 99.0-100% (DECISIONS D135).
        \$SING gmes_petap.pl --predict_with genemark-shared.mod \\
            --sequence genome.fa --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
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
            record_genemark_skip "ET_post_attempt_too_small"
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
        # Chunked route: `--training_only` stops after training and leaves the model at
        # ./gmhmm.mod (NOT output/gmhmm.mod, which only PredictGenes creates); no GTF is written.
        # (Do not use --prediction_only: it skips the masking step, so predictions would differ.)
        TRAIN_ONLY=()
        [ "\$CHUNKED" = "true" ] && TRAIN_ONLY=(--training_only)
        \$SING gmes_petap.pl --ES "\${TRAIN_ONLY[@]}" --sequence genome.fa \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 \\
            --cores ${task.cpus} --fungus "\${GCODE_ARGS[@]}" 2>&1 | tee "\$GMES_LOG"
        if [ "\$CHUNKED" = "true" ] && [ -s gmhmm.mod ]; then
            GMHMM_FAILS=\$(grep -c "warning on: .*gmhmme3" gmes.log 2>/dev/null || true)
            if [ "\${GMHMM_FAILS:-0}" -gt 0 ]; then
                echo "ERROR: GENEMARK_RUN ${out}: gmhmme3 failed \$GMHMM_FAILS time(s) during training; likely out of memory at ${task.memory}" >&2
                exit 1
            fi
            cp gmhmm.mod "${out}.train.mod" || exit 1
            finish_chunked "${out}.train.mod" fresh
            exit 0
        elif [ "\$CHUNKED" != "true" ] && [ -f output/gmhmm.mod ]; then
            cp output/gmhmm.mod "${out}.genemark.mod"
        elif too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark-ES declined -- not enough usable (unmasked, >=10kb) training sequence after masking; skipping GeneMark for this genome" >&2
            record_genemark_skip "ES_post_attempt_too_small"
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        else
            echo "ERROR: GeneMark-ES did not produce its model (output/gmhmm.mod, or ./gmhmm.mod with --training_only)" >&2
            exit 1
        fi
    fi

    if [ ! -s genemark.gtf ]; then
        if too_small_skip; then
            echo "[WARN] GENEMARK_RUN ${out}: GeneMark declined -- not enough usable training sequence; skipping GeneMark for this genome" >&2
            record_genemark_skip "final_check_too_small"
            touch "${out}.genemark.gtf"
            rm -f genome.fa
            exit 0
        fi
        echo "ERROR: GeneMark did not produce genemark.gtf" >&2
        exit 1
    fi
    # GeneMark can stop part way (e.g. gmhmme3 children OOM-killed: 32 kills in an 8 GB,
    # 13-core task) and leave a partial genemark.gtf that still uses its internal
    # sequence names ("1_dna"); the pipe to tee hides its exit status, and predict
    # then fails in EVM on every retry (Fungi_BFD_runs DECISIONS D141). Require every
    # sequence name in the GTF to be a genome contig; failing here retries the task
    # with more memory instead of publishing a broken GTF.
    # gmhmme3 child failures also hit training and can leave a degraded model with a
    # well-formed GTF. In 280 sampled wave 1 runs they appeared only together with OOM
    # kills (D141). Fail so the retry gets more memory.
    GMHMM_FAILS=\$(grep -c "warning on: .*gmhmme3" gmes.log 2>/dev/null || true)
    if [ "\${GMHMM_FAILS:-0}" -gt 0 ]; then
        echo "ERROR: GENEMARK_RUN ${out}: gmhmme3 failed \$GMHMM_FAILS time(s) (gmes.log 'warning on:'); likely out of memory at ${task.memory}" >&2
        exit 1
    fi
    if ! awk 'NR == FNR { if (/^>/) { id = substr(\$1, 2); ids[id] = 1 } next }
              !/^#/ && NF >= 9 && !(\$1 in ids) { print "unknown sequence " \$1 > "/dev/stderr"; bad = 1; exit 1 }
              END { exit bad }' genome.fa genemark.gtf; then
        echo "ERROR: GENEMARK_RUN ${out}: genemark.gtf names sequences that are not in the genome (GeneMark stopped part way; see gmes.log)" >&2
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
    touch "${out}.other.gff3"
    if [ -z "${shared_mod}" ] || [ "${force_independent}" = "true" ]; then
        touch "${out}.genemark.mod"
    fi
    """
}
