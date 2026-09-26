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
    # \$PWD (the task workdir) is NOT covered by any of the other binds -- confirmed
    # empirically (2026-08-24): predict got through startup/training-file parsing
    # (those paths are bound via target/training_target) but then failed with
    # "genome_input.fa is not a valid file" once it actually needed the raw genome
    # FASTA (at the time, inflated into \$PWD; it now inflates into \$TMPDIR instead --
    # see below -- but \$PWD still holds this task's other relative-path outputs). Same
    # missing-bind bug class as GENEMARK_RUN; manually-built \$SING commands
    # get none of Nextflow's automatic task-workdir binding.
    #
    # gene_prediction_shared_abinitio and workflow.workDir are ALSO required here,
    # not covered above: shared_params_json and genemark_gtf are declared as `val`
    # (not `path`) inputs, so Nextflow never stages/symlinks them into this task's
    # own \$PWD -- they are raw absolute host paths pointing at another process's
    # persistent store (gene_prediction_shared_abinitio/<species>/parameters.json)
    # or a sibling GENEMARK_RUN task's own work subdir (workDir/funannotate/xx/...).
    # Without these binds, funannotate's `-p <parameters.json>` open() call fails
    # with FileNotFoundError even though the file exists and is readable on the
    # host -- confirmed 2026-09-04 against Neurospora_tetrasperma_FGSC_2509.
    # Root cause of a prior OOM + Pool.join() deadlock in soft-mask parsing
    # (confirmed 2026-09-06 against Austropuccinia_psidii_Au3-C622-A115012,
    # reproduced on 4 separate SLURM attempts via `sacct -j <id>.batch
    # --format=State,MaxRSS`, every one OUT_OF_MEMORY with MaxRSS pinned at the task's
    # memory limit, plus "Detected N oom_kill event(s)" in .command.log): library.py's
    # checkMasklowMem() spawned a multiprocessing.Pool to compute masking stats per
    # scaffold; its worker function (maskingstats2bed) built a Python list with one
    # int per masked base (~34 bytes/base) before grouping into BED runs -- on this
    # genome's largest scaffolds (up to 86 Mbp, >90% soft-masked) that is several GB
    # PER WORKER, OOM-killing the cgroup. Once SIGKILL took out a worker mid-task,
    # Pool.join() blocked forever waiting on that worker's never-arriving result
    # (CPython bpo-22393) -- explaining the observed 12+ hour near-zero-CPU stall on
    # every attempt: not an I/O hang, just OOM-then-deadlock.
    #
    # Fixed upstream in funannotate-live commit 17c58f5 (rewrites maskingstats2bed
    # to scan runs via regex instead of a per-base list, and checkMasklowMem to use
    # ProcessPoolExecutor so a killed worker raises instead of hanging). This was
    # initially deployed as a local bind-mount overlay of predict.py/library.py
    # (nextflow/patches/funannotate/) while waiting for a container rebuild; that
    # overlay is now retired (2026-09-17) because funannotate-1.9.0-beta.11.sif
    # (built 2026-09-09, after 17c58f5) already has the fix baked in at its native
    # path -- verified byte-identical against the patch files via
    # `apptainer exec ... cat /pixi/.../funannotate/{predict,library}.py`. Still
    # byte-identical in funannotate-1.9.0-beta.12.sif (re-verified 2026-09-19).
    # augustus_parallel.py: beta.11's aux_scripts/augustus_parallel.py unconditionally
    # built `hints_input = '--hintsfile='+args.hints` even when --hints was never
    # passed to it (args.hints is None). A predict run with zero protein AND zero
    # RNA-seq evidence omits --hints entirely, so that line raised TypeError inside
    # every one of the multiprocessing Augustus workers -- instantly, before augustus
    # itself ever ran (confirmed 2026-09-18 against Neopereziidae_sp._gmOTU29_Ox002475:
    # 0/463 chunks completed, 463/463 failed in ~1s). Carried here as a bind-mount
    # overlay (nextflow/patches/funannotate/aux_scripts/augustus_parallel.py), same
    # pattern as the predict.py/library.py patches above; that overlay is now retired
    # (2026-09-19) because funannotate-1.9.0-beta.12.sif ships the guarded form
    # (`'--hintsfile='+args.hints if args.hints else ''`) at its native path --
    # verified code-identical to the patch file (comments aside) via
    # `apptainer exec ... cat /pixi/.../funannotate/aux_scripts/augustus_parallel.py`.
    SING_BINDS="--bind \$PWD:\$PWD,${params.target}:${params.target},${params.training_target}:${params.training_target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},${params.proteins}:${params.proteins},${params.proteins_microsporidia}:${params.proteins_microsporidia},${params.gene_prediction_shared_abinitio}:${params.gene_prediction_shared_abinitio},${workflow.workDir}:${workflow.workDir},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"

    # Microsporidia protein evidence override (see profile_funannotate.config's
    # proteins_microsporidia comment for the exonerate-alignment rationale).
    # busco_lineage=="microsporidia" is a verified 1:1 proxy for samples.csv
    # PHYLUM=Microsporidia (163/163 rows, 2026-09-18) and is already available in
    # this process's input tuple, so no channel/tuple-arity changes are needed.
    PROTEIN_EVIDENCE="${params.proteins}"
    if [ "${busco_lineage}" = "microsporidia" ]; then
        echo "[INFO] ${out}: busco_lineage=microsporidia; using full swissprot (${params.proteins_microsporidia}) as protein evidence instead of the fungal-only subset"
        PROTEIN_EVIDENCE="${params.proteins_microsporidia}"
    fi

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
    # genome assembly/training output newer than the GBK, per staleRnaseq()/staleGenome()/
    # staleTraining()). Re-derive staleness here from the same on-disk timestamps so a
    # current GBK short-circuits, but a stale one forces a clean re-predict.
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
        # FUNANNOTATE_TRAIN can produce a newer/better training set (or flip a prior
        # .pasa_train_failed verdict to trainable) without any rnaseq/trinity/genome file
        # itself changing -- e.g. a retrain that only succeeded this time because an
        # earlier attempt hit infra flakiness (SLURM preemption, an NFS exit-status-read
        # timeout). Mirrors the Groovy-level staleTraining() gate that schedules this task
        # in the first place; without this check here too, the task would reach this point,
        # find its own GBK "current" by every OTHER measure, and silently no-op forever
        # (the GBK mtime never advances past the training GFF3, so staleTraining() keeps
        # rescheduling it on every future -resume). .pasa_train_failed is deliberately
        # 0-byte (`: > marker`), so it's checked with `-e`, not `-s` like the files above.
        TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
        TRAIN_FAILED_MARKER="${params.training_target}/${out}/training/.pasa_train_failed"
        if [ -s "\$TRAIN_GFF3" ] && [ "\$TRAIN_GFF3" -nt "\$SKIP_GBK" ]; then
            STALE=1; STALE_REASON="training output"
        elif [ -e "\$TRAIN_FAILED_MARKER" ] && [ "\$TRAIN_FAILED_MARKER" -nt "\$SKIP_GBK" ]; then
            STALE=1; STALE_REASON="training output"
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
    #
    # \$PREDICTDIR/training must end up as EITHER a symlink to the canonical training dir
    # OR absent -- never a real directory. `ln -sfn` does NOT replace an existing real
    # (non-symlink) directory: confirmed empirically that it silently exits 0 and nests a
    # stray symlink *inside* it instead, leaving the stale real directory (and whatever
    # training data is or isn't in it) untouched with no error and no warning. That is the
    # exact mechanism behind the "cryptic" 2026-09-10 training-symlink incident (broken/
    # missing symlinks plus 776 duplicated real training dirs found dataset-wide -- see
    # .living/learnings.md in ../Fungi_BFD and scripts/one-off/reconcile_training_duplicates.py).
    # Hard-fail instead of silently proceeding against stale/wrong data.
    if [ -e "\$PREDICTDIR/training" ] && [ ! -L "\$PREDICTDIR/training" ]; then
        echo "[ERROR] ${out}: \$PREDICTDIR/training exists as a real (non-symlink) directory, not a symlink to the canonical training dir at ${params.training_target}/${out}/training. This is data corruption, not a normal predict failure -- reconcile it (e.g. scripts/one-off/reconcile_training_duplicates.py) before predict can proceed for this species." >&2
        exit 1
    fi
    if [ -L "\$PREDICTDIR/training" ] && [ ! -e "\$PREDICTDIR/training" ]; then
        echo "[WARN] ${out}: \$PREDICTDIR/training was a broken symlink; removing before relinking" >&2
        rm -f "\$PREDICTDIR/training"
    fi
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$PREDICTDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy; funannotate
    # cannot read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through.
    # Inflated into \$TMPDIR (node-local \$SCRATCH, exported above), NOT \$PWD (this
    # task's workdir, which lives on /bigdata/NFS): predict rereads this multi-GB
    # genome file across its own steps, and node-local scratch is the right place for
    # it regardless. NOTE (2026-09-06): this was originally written up as fixing an
    # "NFS-latency-bound genome parse" for Austropuccinia_psidii_Au3-C622-A115012 --
    # that theory was wrong (the actual 20-scaffold, ~1.1 GB genome parses to
    # per-scaffold FASTA files in ~10s from NFS; the real 12+ hour stall was an
    # OOM-kill + multiprocessing.Pool.join() deadlock inside checkMasklowMem, now
    # fixed upstream -- see SING_BINDS comment above). Kept as a legitimate hygiene
    # improvement, not because it was load-bearing for that failure.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA to \$TMPDIR"; pigz -dc "\$GENOME_FA" > "\$TMPDIR/genome_input.fa"; GENOME_IN="\$TMPDIR/genome_input.fa" ;;
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
    read ASM_BP ASM_CTG ASM_N50 ASM_VERDICT ASM_REPEAT_PCT < <(
        python "${workflow.projectDir}/bin/asm_preflight_stats.py" "\$GENOME_IN" \\
            --min-bp ${params.predict_min_asm_bp} --max-n50 ${params.predict_frag_max_n50} \\
            --max-contigs ${params.predict_frag_max_contigs} \\
            --min-contig-len ${params.predict_min_training_contig_len} \\
            --min-training-contigs ${params.predict_min_training_contigs} \\
            --abs-min-bp ${params.predict_abs_min_asm_bp} \\
            --report-repeat-pct)
    echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}, \${ASM_REPEAT_PCT}% repeat-masked"
    if [ "\$ASM_VERDICT" != "ok" ] && ! { [ "\$ASM_VERDICT" = "small_fragmented" ] && [ -s "${other_gff}" ]; }; then
        echo "[WARN] ${out} failed preflight ('\$ASM_VERDICT': \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
        mkdir -p "${params.target}"
        [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_\$ASM_VERDICT" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
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
    # STATUS.md:73-75) -- this does NOT change augustus/snap weighting for any
    # genome where other_gff is empty, which keeps today's default-on
    # behavior. --no-evm-partitions and --min_protlen 30 mirror the
    # validated recipe (microsporidia-default.json) for these near-zero-
    # intergenic compact genomes; both are new flags this pipeline didn't
    # pass before, applied ONLY on this branch.
    #
    # CORRECTION (2026-09-04, superseding the paragraph this replaces): a
    # real, non-stub run against Ordospora colligata OC4 proved the claim
    # below wrong. predict.py sets RunModes["augustus"] (to "pretrained",
    # "busco", or "pasa") UNCONDITIONALLY, independent of augustus/snap/
    # glimmerhmm's EVM weight -- weight only controls whether that
    # predictor's OUTPUT is included in the final EVM consensus, not
    # whether predict.py attempts to prepare/train it at all. This is why
    # --busco_seed_species (see below) is required regardless of weight.
    # What zeroing the weight DOES still guarantee, confirmed by the real
    # run's own gene count matching the standalone reference within ~1%
    # (see Task 5 Step 6 in nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    # augustus/snap/glimmerhmm's training/reuse still happens (fast, via a
    # pre-existing --busco_seed_species entry's info.json, not a fresh
    # from-scratch BUSCO run against this genome's own too-few complete
    # models -- Microsporidia_predict/STATUS.md:73-82's 36-vs-200 failure),
    # but their EVM weight of 0 excludes their output from the final gene
    # models regardless of what that training produced.
    # Do NOT also add --min_training_models 0 or drop --busco_db here --
    # neither is necessary and both would be pure noise.
    GENEMARK_GTF_FLAG=()
    OTHER_GFF_FLAG=()
    EXTRA_PREDICT_ARGS=()
    # PASA training-set gate: empty unless the param is set, because a
    # funannotate_sif built before the gate existed rejects the unknown flag;
    # a gate-aware sif applies its own default (500 complete-ORF models).
    if [ -n "${params.predict_min_pasa_complete_models != null ? params.predict_min_pasa_complete_models : ''}" ]; then
        EXTRA_PREDICT_ARGS+=(--min_pasa_complete_models ${params.predict_min_pasa_complete_models})
    fi
    WEIGHT_ARGS=(codingquarry:0 glimmerhmm:0)
    if [ -s "${other_gff}" ]; then
        echo "[INFO] ${out}: using Prodigal evidence from ${other_gff} (--other_gff weight 5)"
        OTHER_GFF_FLAG=(--other_gff "${other_gff}:5")
        WEIGHT_ARGS+=(augustus:0 snap:0)
        # --busco_seed_species is REQUIRED here, not cosmetic: predict.py
        # unconditionally sets RunModes["augustus"]="busco" and requires a
        # pre-existing trained_species entry to seed it, REGARDLESS of
        # augustus/snap EVM weight (confirmed directly: RunModes["augustus"]
        # population happens before/independent of the weight-gated EVM
        # combiner step -- our earlier belief that zeroing the weight alone
        # skips this entirely was wrong for THIS code path). Without a valid
        # entry, predict hard-aborts with "ERROR: --busco_seed_species {} is
        # not valid" (confirmed by a real, non-stub run against Ordospora
        # colligata OC4). The literal string "microsporidia" is NOT a valid
        # seed species -- confirmed neither this project nor the shared
        # production funannotate_db has ever registered an Augustus species
        # by that name (checked directly: only real microsporidia strain
        # entries like Encephalitozoon_cuniculi exist). Use an existing real
        # trained species as the seed instead -- its own weight is 0 so its
        # actual training output never influences the final EVM gene models,
        # it only exists to satisfy this startup precondition.
        EXTRA_PREDICT_ARGS+=(--no-evm-partitions --min_protlen 30 --busco_seed_species Encephalitozoon_cuniculi)
    fi
    # -s not -n: GENEMARK_RUN's too-small-genome skip path emits a real but
    # deliberately empty ${out}.genemark.gtf (see GENEMARK_RUN/main.nf).
    if [ -s "${genemark_gtf}" ]; then
        echo "[INFO] ${out}: using pre-computed GeneMark GTF from ${genemark_gtf}"
        GENEMARK_GTF_FLAG=(--genemark_gtf "${genemark_gtf}")
        WEIGHT_ARGS+=(genemark:1)
    fi

    # other_gff/genemark_gtf are `val`, not `path`, in this process's input
    # tuple -- Nextflow never stages/symlinks them into THIS task's own
    # \$PWD, so they're referenced by their original absolute path under
    # GENEMARK_RUN's own separate task work directory. Same missing-bind bug
    # class already documented above for \$PWD/genome_input.fa: a path
    # outside SING_BINDS's explicit list is invisible inside the container
    # even though the host shell can read it fine. Confirmed empirically
    # (2026-09-04, real non-stub run against Ordospora colligata OC4):
    # without this, funannotate predict fails with "<path>/Ordospora_
    # colligata_OC4.other.gff3 is not a valid file, exiting" despite the
    # host-side `[ -s "${other_gff}" ]` check above having already confirmed
    # the file exists and is non-empty. Mirrors GENEMARK_RUN's own
    # training_bam dirname-binding pattern (GENEMARK_RUN/main.nf).
    # -s not -n (matching the flag-construction checks above): a set-but-empty
    # or set-but-missing path must not add a bind whose source doesn't exist --
    # apptainer refuses to start if it does.
    if [ -s "${other_gff}" ]; then
        SING_BINDS="\$SING_BINDS,\$(dirname "${other_gff}"):\$(dirname "${other_gff}")"
    fi
    if [ -s "${genemark_gtf}" ]; then
        SING_BINDS="\$SING_BINDS,\$(dirname "${genemark_gtf}"):\$(dirname "${genemark_gtf}")"
    fi
    SING="apptainer exec \${SING_BINDS} ${params.funannotate_sif}"

    # ── Repeat-aware EVM mode ─────────────────────────────────────────────────
    # See profile_funannotate.config's predict_evm_repeat_pct_threshold comment
    # for the Austropuccinia_psidii/GCA_003724095.1 "Evidence modeler has
    # failed" failures this addresses. \$ASM_REPEAT_PCT comes from the same
    # preflight FASTA pass as the small-genome guard above (no second scan).
    EVM_REPEAT_FLAGS=()
    THRESHOLD=${params.predict_evm_repeat_pct_threshold}
    if [ "\$THRESHOLD" != "0" ] && awk -v p="\$ASM_REPEAT_PCT" -v t="\$THRESHOLD" 'BEGIN{exit !(p>=t)}'; then
        echo "[INFO] ${out}: \${ASM_REPEAT_PCT}% repeat-masked >= \${THRESHOLD}% threshold -- enabling repeat-aware EVM mode (--repeats2evm, --evm-partition-interval ${params.predict_evm_repeat_aware_interval})"
        EVM_REPEAT_FLAGS=(--repeats2evm --evm-partition-interval ${params.predict_evm_repeat_aware_interval})
        if [ "${params.predict_evm_repeat_aware_drop_snap}" = "true" ]; then
            echo "[INFO] ${out}: repeat-aware mode also zeroing SNAP's EVM weight (predict_evm_repeat_aware_drop_snap=true)"
            WEIGHT_ARGS+=(snap:0)
        fi
    fi

    \$SING funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w "\${WEIGHT_ARGS[@]}" \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence "\$PROTEIN_EVIDENCE" \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} --auto-skip-genemark \\
        "\${ABINITIO_REUSE_FLAG[@]}" "\${GENEMARK_GTF_FLAG[@]}" "\${OTHER_GFF_FLAG[@]}" "\${EXTRA_PREDICT_ARGS[@]}" "\${EVM_REPEAT_FLAGS[@]}" || true

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
        # Besides the ab-initio parameters and tRNAs, keep the small files needed
        # to diagnose a gene-count drop after the fact (2026-09-25: the pilot's
        # EVM undercall could not be traced because these were deleted):
        #   weights.evm.txt           -- EVM weights actually used
        #   final_training_models.gff3 -- models Augustus/SNAP were trained on
        # Each is a few hundred KB at most; missing ones are skipped.
        KEEP_DIR="\$PREDICTDIR/.predict_misc_keep"
        rm -rf "\$KEEP_DIR"; mkdir -p "\$KEEP_DIR"
        for f in ab_initio_parameters trnascan.no-overlaps.gff3 weights.evm.txt final_training_models.gff3; do
            [ -e "\$PREDICTDIR/predict_misc/\$f" ] && mv "\$PREDICTDIR/predict_misc/\$f" "\$KEEP_DIR/"
        done
        [ -f "\$KEEP_DIR/final_training_models.gff3" ] && pigz "\$KEEP_DIR/final_training_models.gff3"
        rm -rf "\$PREDICTDIR/predict_misc"
        mv "\$KEEP_DIR" "\$PREDICTDIR/predict_misc"
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
