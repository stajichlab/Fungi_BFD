// Option B persistence model: funannotate predict computes DIRECTLY into the persistent
// per-genome dir (${params.target}/${out}), symmetric with FUNANNOTATE_TRAIN writing to
// training_target. funannotate checkpoints into predict_misc/, so a restart after an
// OOM/timeout/orchestrator death resumes completed steps in place rather than starting
// over. There is no publishDir copy and no work-dir<->target rsync: the durable output is
// written where downstream steps already read it. Large intermediates still go to the
// node-local --tmpdir. The Nextflow output is a small marker file (nothing consumes the
// predict dir as a channel; postpredict rebuilds metadata from the CSV and gates on the
// on-disk GBK), so emitting a marker keeps the DAG edge without copying the result tree.
process FUNANNOTATE_PREDICT {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '32h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(shared_params_json), val(genemark_gtf), val(other_gff)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table), emit: metadata
    path("${out}.predict.done"), emit: done

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer
    # ── Containerized funannotate ──────────────────────────────────────────────
    # Same swap/rationale as FUNANNOTATE_TRAIN/main.nf: funannotate predict runs
    # via \$SING (apptainer exec \${params.funannotate_sif}) instead of
    # `module load funannotate`. GeneMark is never invoked from inside this
    # container -- it arrives pre-computed via --genemark_gtf from the
    # upstream GENEMARK_RUN process, which runs it from its own separate,
    # privately-built (non-redistributable) container (see the weight
    # comment below). No mysql/mariadb sidecar here (PASA-mysql is
    # train/update-specific); binds cover target/training_target (predict
    # writes directly into target, reads the training/ symlink target),
    # augustus_config, funannotate_db, params.proteins (the --protein_evidence
    # file), and \$TMPDIR. Same `which` false-alarm note as FUNANNOTATE_TRAIN
    # applies here -- unset defensively anyway.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    # APPTAINERENV_ prefix is REQUIRED here, not a plain export: the
    # funannotate-live image bakes its own ENV defaults for these vars
    # (AUGUSTUS_CONFIG_PATH=/venv/config, FUNANNOTATE_DB=/opt/databases/...),
    # which apptainer sources from the image's /.singularity.d/env/ scripts
    # AFTER inheriting the host shell env -- so a plain `export` alone is
    # silently overwritten inside the container (confirmed empirically
    # 2026-08-24: host FUNANNOTATE_DB was ignored, predict looked for the
    # repeat DB under /opt/databases instead). APPTAINERENV_/SINGULARITYENV_-
    # prefixed vars are applied last and always win. Also keep the plain
    # (unprefixed) AUGUSTUS_CONFIG_PATH export -- this script's own
    # --AUGUSTUS_CONFIG_PATH CLI arg below reads it as a HOST shell
    # variable, not through the container, so it needs both.
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export APPTAINERENV_AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export APPTAINERENV_FUNANNOTATE_DB=${params.funannotate_db}
    export TMPDIR=\${SCRATCH:-/tmp}
    # \$PWD (the task workdir, where genome_input.fa is inflated below) is NOT
    # covered by any of the other binds -- confirmed empirically (2026-08-24):
    # predict got through startup/training-file parsing (those paths are
    # bound via target/training_target) but then failed with "genome_input.fa
    # is not a valid file" once it actually needed the raw genome FASTA. Same
    # missing-bind bug class as GENEMARK_RUN; manually-built \$SING commands
    # get none of Nextflow's automatic task-workdir binding.
    SING_BINDS="--bind \$PWD:\$PWD,${params.target}:${params.target},${params.training_target}:${params.training_target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},${params.proteins}:${params.proteins},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"

    PREDICTDIR="${params.target}/${out}"
    PREDICT_GBK="\$PREDICTDIR/predict_results/${out}.gbk"

    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] out=${out} asmid=${asmid} species=${species} strain=${strain}"
        echo "[DEBUG] locustag=${locustag} busco=${busco_lineage} transl_table=${transl_table}"
        echo "[DEBUG] proteins=${params.proteins} genome_fa=${genome_fa}"
        echo "[DEBUG] PREDICTDIR=\$PREDICTDIR TMPDIR=\$TMPDIR pwd=\$(pwd)"
    fi

    # ── Skip vs. refresh decision ─────────────────────────────────────────────
    # The workflow schedules this process when the GBK is missing OR stale (rnaseq/trinity/
    # genome assembly newer than the GBK, per staleRnaseq()/staleGenome()). Re-derive
    # staleness here from the same on-disk timestamps so a current GBK short-circuits, but
    # a stale one forces a clean re-predict.
    # Accept a compressed prediction (.gbk.gz) as "done" so folders can be space-saved.
    SKIP_GBK="\$PREDICT_GBK"
    [ -s "\$SKIP_GBK" ] || SKIP_GBK="\$PREDICTDIR/predict_results/${out}.gbk.gz"
    if [ -s "\$SKIP_GBK" ]; then
        SPECIES_TAG=\$(printf '%s' "${species}" | sed -E 's/[[:space:]]+/_/g')
        STALE=0
        STALE_REASON=""
        for f in "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_R1.fastq.gz" \\
                 "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_SE.fastq.gz" \\
                 "${launchDir}/rnaseq_data/\${SPECIES_TAG}.trinity-GG.fasta" \\
                 "${shared_params_json}"; do
            if [ -n "\$f" ] && [ -s "\$f" ] && [ "\$f" -nt "\$SKIP_GBK" ]; then STALE=1; STALE_REASON="rnaseq/trinity/shared-params"; fi
        done
        # Genome assembly itself can be swapped/updated for the same asmid path without
        # any rnaseq change; catch that here too so a re-assembled genome doesn't keep an
        # annotation predicted against the old coordinates. See FUNANNOTATE_TRAIN's
        # analogous check for the same underlying gap.
        if [ -s "${genome_fa}" ] && [ "${genome_fa}" -nt "\$SKIP_GBK" ]; then
            STALE=1; STALE_REASON="genome assembly"
        fi
        if [ "\$STALE" -eq 0 ]; then
            echo "[INFO] Prediction already complete and current for ${out}; nothing to do"
            touch ${out}.predict.done
            exit 0
        fi
        echo "[INFO] Stale prediction for ${out}: \$STALE_REASON newer than GBK — clearing predict outputs for a fresh run"
        rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
    fi

    mkdir -p "\$PREDICTDIR"

    # ── Guard against a corrupt partial from a previous attempt ───────────────
    # funannotate resumes from predict_misc/. If predict_results/ exists without a
    # predict_misc/ (a half-written tree with no checkpoints and no GBK), clear it so
    # predict starts the consensus/output step from a clean state instead of choking on it.
    if [ ! -d "\$PREDICTDIR/predict_misc" ] && [ -d "\$PREDICTDIR/predict_results" ]; then
        echo "[WARN] predict_results/ present without predict_misc/ for ${out}; clearing stale partial"
        rm -rf "\$PREDICTDIR/predict_results"
    fi

    # funannotate predict expects training data at <outdir>/training; point it at the
    # persistent training dir. The symlink lives in the persistent project tree (no
    # publishDir to recursively copy the target), so it is left in place.
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$PREDICTDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy; funannotate
    # cannot read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Too-small-genome pre-flight guard ────────────────────────────────────
    # Assemblies that are both small AND fragmented cannot yield funannotate's
    # required 30 training models; predict would run for hours then abort with
    # "Not enough gene models N to train Augustus (30 required), exiting". Detect
    # that up front from cheap contig stats (bin/asm_preflight_stats.py -- shared
    # with GENEMARK_RUN, which needs the identical policy upstream of this
    # process; see GENEMARK_RUN_DESIGN.md) and skip cleanly (flag, no crash).
    # Requires BOTH gates so complete small genomes (e.g. Malassezia) are unaffected.
    # See analysis/funannotate_model_failures/. Disabled when predict_min_asm_bp=0.
    SKIP_REPORT="${params.target}/predict_skipped_too_small.tsv"
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" "\$GENOME_IN" \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs})
    echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}"
    if [ "\$ASM_VERDICT" = "small_fragmented" ] && [ ! -s "${other_gff}" ]; then
        echo "[WARN] ${out} is too small/fragmented for funannotate training (\${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
        mkdir -p "${params.target}"
        [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_small_fragmented" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
        touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
        touch ${out}.predict.done
        exit 0
    elif [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[INFO] ${out} is small/fragmented but has Prodigal evidence (${other_gff}); proceeding instead of skipping" >&2
    fi

    # Species-level ab-initio parameter reuse (todo/species_level_abinitio_reuse.md):
    # -p <parameters.json> tells funannotate predict to reuse pre-trained AUGUSTUS/
    # SNAP/GeneMark-ES parameters instead of (re)training them from scratch. Empty
    # shared_params_json means today's independent-training behavior, unchanged.
    ABINITIO_REUSE_FLAG=()
    if [ -n "${shared_params_json}" ]; then
        echo "[INFO] ${out}: reusing shared ab-initio parameters from ${shared_params_json}"
        ABINITIO_REUSE_FLAG=(-p "${shared_params_json}")
    fi

    # GeneMark now runs as its own upstream GENEMARK_RUN process (own
    # privately-built, non-redistributable container -- see nextflow/docs/
    # GENEMARK_RUN_DESIGN.md), which hands predict a pre-computed GTF here.
    # --genemark_gtf makes predict skip its internal GeneMark call entirely
    # (verified in funannotate-live/funannotate/predict.py: the --genemark_gtf
    # branch short-circuits before RunGeneMarkES/RunGeneMarkET ever run), so
    # -p's own genemark_mod reuse never gets a chance to redundantly re-run
    # GeneMark -- Augustus/SNAP reuse from -p is unaffected. Empty
    # genemark_gtf (run_genemark=false, or GENEMARK_RUN unavailable) falls
    # back to --auto-skip-genemark's existing graceful degradation, unchanged.
    #
    # -w genemark:1 MUST be passed explicitly whenever genemark_gtf is set --
    # predict.py:567 unconditionally zeroes StartWeights["genemark"] when
    # gmes_petap.pl isn't on the host running predict (`if not genemarkcheck:
    # StartWeights["genemark"] = 0`), with NO check for whether --genemark_gtf
    # was supplied as an alternative. This process already runs funannotate
    # predict from the funannotate_sif container, which has no gmes_petap.pl
    # at all, so genemarkcheck is always False here -- without this override,
    # a correctly-supplied --genemark_gtf would silently get EVM weight 0
    # (the precomputed evidence consumed but then discarded). Confirmed by
    # reading predict.py and validated end-to-end (Test 3,
    # GENEMARK_RUN_DESIGN.md: final gene count 11,202 vs. baseline 11,198).
    #
    # All -w values MUST go in a single -w group, not two separate -w flags:
    # funannotate's argparse `-w/--weights` (nargs='+', no action='append')
    # replaces the whole list on a second -w occurrence rather than merging
    # (verified directly against argparse) -- a second `-w genemark:1` would
    # silently drop codingquarry:0/glimmerhmm:0 entirely.
    # Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    # other_gff non-empty means GENEMARK_RUN ran the Prodigal supplement for
    # this genome (is_microsporidia=true AND below predict_min_asm_bp).
    # AUGUSTUS/SNAP forced off here specifically -- both need >=200 BUSCO
    # training models these genomes don't have (Microsporidia_predict
    # PLAN.md 9.15) -- this does NOT change augustus/snap weighting for any
    # genome where other_gff is empty, which keeps today's default-on
    # behavior. --no-evm-partitions and --min_protlen 30 mirror the
    # validated recipe (microsporidia-default.json) for these near-zero-
    # intergenic compact genomes; both are new flags this pipeline didn't
    # pass before, applied ONLY on this branch.
    #
    # Why zeroing the weight is enough, and no further BUSCO change is
    # needed: funannotate's predict.py only populates RunModes["augustus"]/
    # ["snap"] when their StartWeight > 0; RunBusco is only set True when
    # some RunModes value == "busco". With augustus/snap/glimmerhmm all at
    # weight 0, RunBusco never becomes True and BUSCO is never invoked at
    # all -- not "invoked and tolerated despite too few models". The
    # --min_training_models 30 check a few lines below is itself nested
    # inside `if "augustus" in RunModes:`, so it stays unreachable dead code
    # on this branch. This was confirmed by tracing predict.py directly (not
    # just inferred) and matches the failure already documented in
    # Microsporidia_predict/STATUS.md:73-82 (passing --busco_db
    # microsporidia_odb10 with augustus/snap weight > 0 found only 36
    # complete BUSCO models against a >=200 requirement and failed; the fix
    # there was this same weight-zeroing, not any BUSCO-side workaround).
    # Do NOT also add --min_training_models 0 or drop --busco_db here --
    # neither is necessary and both would be pure noise.
    GENEMARK_GTF_FLAG=()
    OTHER_GFF_FLAG=()
    EXTRA_PREDICT_ARGS=()
    WEIGHT_ARGS=(codingquarry:0 glimmerhmm:0)
    if [ -s "${other_gff}" ]; then
        echo "[INFO] ${out}: using Prodigal evidence from ${other_gff} (--other_gff weight 5)"
        OTHER_GFF_FLAG=(--other_gff "${other_gff}:5")
        WEIGHT_ARGS+=(augustus:0 snap:0)
        # --busco_seed_species microsporidia is REQUIRED here, not cosmetic:
        # without it funannotate predict hard-aborts immediately with
        # "ERROR: --busco_seed_species {} is not valid as it is not in
        # database" regardless of augustus/snap weight (confirmed by a real,
        # non-stub run against Ordospora colligata OC4 -- omitting this flag
        # was a real gap this plan's original Task 4 Step 5 missed relative
        # to the validated reference implementation, which always passes
        # --busco_seed_species microsporidia; see
        # ../../Microsporidia_predict/scripts/run_microsporidia_predict.py).
        EXTRA_PREDICT_ARGS+=(--no-evm-partitions --min_protlen 30 --busco_seed_species microsporidia)
    fi
    # -s not -n: GENEMARK_RUN's too-small-genome skip path emits a real but
    # deliberately empty ${out}.genemark.gtf (see GENEMARK_RUN/main.nf).
    if [ -s "${genemark_gtf}" ]; then
        echo "[INFO] ${out}: using pre-computed GeneMark GTF from ${genemark_gtf}"
        GENEMARK_GTF_FLAG=(--genemark_gtf "${genemark_gtf}")
        WEIGHT_ARGS+=(genemark:1)
    fi

    \$SING funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w "\${WEIGHT_ARGS[@]}" \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} --auto-skip-genemark \\
        "\${ABINITIO_REUSE_FLAG[@]}" "\${GENEMARK_GTF_FLAG[@]}" "\${OTHER_GFF_FLAG[@]}" "\${EXTRA_PREDICT_ARGS[@]}" || true

    # ── Post-predict catch ────────────────────────────────────────────────────
    # If predict produced no GBK, distinguish the known "too few training models"
    # outcome (an unfixable property of the assembly) from a genuine error. The
    # former is flagged and skipped so it does not abort the batch; anything else
    # still hard-fails so real problems surface.
    if [ ! -s "\$PREDICT_GBK" ]; then
        PLOG="\$PREDICTDIR/logfiles/funannotate-predict.log"
        if [ -f "\$PLOG" ] && grep -q "Not enough gene models .* to train Augustus" "\$PLOG"; then
            NMODELS=\$(grep -oE "Not enough gene models [0-9]+" "\$PLOG" | grep -oE "[0-9]+" | tail -1)
            echo "[WARN] ${out}: funannotate found only \${NMODELS:-<min} training models (needs 30); too small/fragmented to annotate — skipping" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "funannotate_too_few_models:\${NMODELS:-NA}" "" "" "" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
        echo "ERROR: funannotate predict did not produce expected GBK: \$PREDICT_GBK" >&2
        exit 1
    fi
    if [ -d "\$PREDICTDIR/predict_misc/ab_initio_parameters" ]; then
        mv "\$PREDICTDIR/predict_misc/ab_initio_parameters" "\$PREDICTDIR"
        mv "\$PREDICTDIR/predict_misc/trnascan.no-overlaps.gff3" "\$PREDICTDIR"
        rm -rf "\$PREDICTDIR/predict_misc"
        mkdir -p "\$PREDICTDIR/predict_misc"
        mv "\$PREDICTDIR/ab_initio_parameters" "\$PREDICTDIR/trnascan.no-overlaps.gff3" "\$PREDICTDIR/predict_misc"
    fi
    find "\$PREDICTDIR/predict_results/" -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    # Provenance marker for ab-initio parameter reuse (decision 7,
    # todo/species_level_abinitio_reuse.md S6) -- lets a later BUSCO/QC regression be
    # traced back to "shared params" vs "independent training" without re-parsing logs.
    rm -f "\$PREDICTDIR/${out}.predict.abinitio_reused"
    if [ -n "${shared_params_json}" ]; then
        printf '%s\\n' "${shared_params_json}" > "\$PREDICTDIR/${out}.predict.abinitio_reused"
    fi
    sync
    touch ${out}.predict.done
    echo "[INFO] Prediction complete for ${out} at \$PREDICTDIR"
    """

    stub:
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || [ -f "${genome_fa}.gz" ] || { echo "ERROR: genome not found at ${genome_fa}[.gz]" >&2; exit 1; }
    mkdir -p ${params.target}/${out}/predict_results ${params.target}/${out}/predict_misc
    touch ${params.target}/${out}/predict_results/${out}.gbk ${params.target}/${out}/predict_results/${out}.proteins.fa
    touch ${out}.predict.done
    """
}
