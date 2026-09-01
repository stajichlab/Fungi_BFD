//
// FUNANNOTATE_RNASEQ — RNA-seq acquisition and funannotate train.
//
// With --run_sra_fetch on: SRA_QUERY_BATCH resolves accessions per species,
// SRA_FETCH/SRA_FETCH_SE download them once per species, RNASEQ_PREPARE runs
// funannotate train on the representative assembly and archives Trinity-GG plus
// trimmed/normalised reads to rnaseq_data/, and every other strain of that
// species runs FUNANNOTATE_TRAIN --trinity against the archived assembly.
// With it off, samples pass straight through untrained.
//
// "Representative assembly" here is the same strain FUNANNOTATE_PREDICTION
// uses for ab-initio training reuse (abinitioReuseMap, loaded once in
// funannotate.nf and passed to both subworkflows) -- picked by ANI+BUSCO
// completeness/N50 in PICK_REPRESENTATIVE_STRAIN, not an arbitrary species-
// group ordering. When abinitioReuseMap has no assignment for a species
// (--run_ani_reuse false, or share_abinitio_params off, or the CSV doesn't
// exist yet), repr_ch falls back to the first assembly in the group, exactly
// as before this was wired up.
//
// rnaseqRepOverride (loaded once in funannotate.nf, same as abinitioReuseMap) lets
// scripts/pick_rnaseq_representative_override.py force a *different* strain for this
// species-level pick specifically -- ANI+BUSCO optimizes for assembly quality, which
// is not always the genome the species' actual RNA-seq reads align to well. See that
// script and loadRnaseqRepresentativeOverride() in modules/funannotate/utils.nf.
//

include { SRA_QUERY_BATCH   } from '../../modules/funannotate/rnaseq/SRA_QUERY_BATCH/main.nf'
include { COLLECT_SRA_QUERY } from '../../modules/funannotate/rnaseq/COLLECT_SRA_QUERY/main.nf'
include { WRITE_EMPTY_READS } from '../../modules/funannotate/rnaseq/WRITE_EMPTY_READS/main.nf'
include { SRA_FETCH         } from '../../modules/funannotate/rnaseq/SRA_FETCH/main.nf'
include { SRA_FETCH_SE      } from '../../modules/funannotate/rnaseq/SRA_FETCH_SE/main.nf'
include { RNASEQ_PREPARE    } from '../../modules/funannotate/rnaseq/RNASEQ_PREPARE/main.nf'
include { COUNT_TRINITY_TRANSCRIPTS } from '../../modules/funannotate/rnaseq/COUNT_TRINITY_TRANSCRIPTS/main.nf'
include { TRINITY_STANDALONE } from '../../modules/funannotate/rnaseq/TRINITY_STANDALONE/main.nf'
include { FUNANNOTATE_TRAIN } from '../../modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf'
// A dedicated process/storeDir, NOT a WRITE_EMPTY_READS alias -- see that module's
// header comment for why aliasing it (same storeDir + filename as SRA_FETCH's real
// output) silently resurrected stale real reads for hybrid species that had
// accumulated them before this feature existed (confirmed 2026-08-28,
// Saccharomyces_x_bayanus_NBRC1948).
include { WRITE_EMPTY_HYBRID_READS } from '../../modules/funannotate/rnaseq/WRITE_EMPTY_HYBRID_READS/main.nf'
include { BUILD_HYBRID_COMPOSITE_TRINITY } from '../../modules/funannotate/rnaseq/BUILD_HYBRID_COMPOSITE_TRINITY/main.nf'

include { gbkResult; staleRnaseq } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_RNASEQ {
    take:
    predict_genome_ch   // tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid)
    abinitioReuseMap     // out -> [species, reuse_eligible, is_representative], from loadAbinitioReuseMap()
    rnaseqRepOverride     // species_tag -> out, from loadRnaseqRepresentativeOverride()
    hybridParentage       // hybrid_species_tag -> [parent_species_tag, ...], from loadHybridParentage()
                           // see nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md

    main:
    // When SRA is enabled: SRA_FETCH fetches reads once per species; RNASEQ_PREPARE runs
    // funannotate train on the representative assembly and archives Trinity-GG, trimmed, and
    // normalized reads to rnaseq_data/; all other strains run FUNANNOTATE_TRAIN --trinity.
    def predict_input_ch
    def reads_ch = Channel.empty()
    if (params.run_sra_fetch.toBoolean()) {
        // ── Hybrid-cross species never go through SRA acquisition ─────────────
        // (nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md). Querying SRA per
        // hybrid-cross taxid is fragile/sparse (narrow nothospecies taxids
        // routinely return next to nothing -- see Saccharomyces_x_bayanus_FM677,
        // 2026-08-28), and sharing one hybrid strain's own admixed Trinity
        // across differently-recombined siblings assumes similarity the
        // biology doesn't support. Instead they get composite parent-transcript
        // evidence built further down. hybridParentage is a plain offline map
        // (loaded once in funannotate.nf), so this split is a synchronous
        // per-item test, not a new channel dependency.
        def genome_branched = predict_genome_ch.branch { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
            hybrid: hybridParentage.containsKey(species.replaceAll(/\s+/, '_'))
            normal: true
        }
        def normalGenomeCh = genome_branched.normal
        def hybridGenomeCh = genome_branched.hybrid

        // Build per-species input: group assemblies, keep first taxonid per species.
        def sra_input = normalGenomeCh
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                def species_tag = species.replaceAll(/\s+/, '_')
                tuple(species_tag, taxonid)
            }
            .groupTuple(by: 0)
            .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

        // Step 1: query or reuse cached per-species SRA query results.
        // skip_sra_query=true reads existing CSVs from rnaseq_reads/sra_query/ directly,
        // bypassing SRA_QUERY_BATCH entirely (no SLURM jobs submitted).
        def sra_query_results
        if (params.skip_sra_query.toBoolean()) {
            sra_query_results = sra_input
                .map { species_tag, _taxonid ->
                    def csv = file("${launchDir}/rnaseq_reads/sra_query/${species_tag}.sra_query.csv")
                    if (!csv.exists()) {
                        log.warn "skip_sra_query: no cached CSV for ${species_tag} — skipping this species"
                        return null
                    }
                    tuple(species_tag, csv)
                }
                .filter { it != null }
        } else {
            // Filter out species whose CSV is already cached on disk *before* batching,
            // so a batch made up entirely of cached species never gets submitted to
            // SLURM at all. SRA_QUERY_BATCH's own inline `[ -s "$cached" ]` check still
            // covers a batch with a MIX of cached/uncached species (that job still runs,
            // just does no esearch/efetch work for the cached ones), but this branch is
            // what stops a fully-cached batch from launching a short-queue job in the
            // first place -- mirrors what SRA_FETCH gets for free from storeDir.
            def sra_branched = sra_input
                .branch { species_tag, _taxonid ->
                    def csv = file("${launchDir}/rnaseq_reads/sra_query/${species_tag}.sra_query.csv")
                    cached:     csv.exists() && csv.size() > 0
                    needs_query: true
                }
            def cached_results = sra_branched.cached
                .map { species_tag, _taxonid ->
                    tuple(species_tag, file("${launchDir}/rnaseq_reads/sra_query/${species_tag}.sra_query.csv"))
                }
            def sra_batched = sra_branched.needs_query
                .collate(params.sra_query_batch_size)
                .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
            SRA_QUERY_BATCH(sra_batched)
            def queried_results = SRA_QUERY_BATCH.out.query_results
                .flatten()
                .map { csv -> tuple(csv.baseName.replaceAll(/\.sra_query$/, ''), csv) }
            sra_query_results = cached_results.mix(queried_results)
        }

        // Step 2: Collect all per-species results into {stem}.rnaseq_sra.csv
        def stem = file(params.samples).baseName
        COLLECT_SRA_QUERY(
            sra_query_results.map { _stag, csv -> csv }.collect(),
            stem
        )

        if (!params.stop_after_sra_query.toBoolean()) {
        // Step 3: Classify each species CSV for routing.
        // Read blacklist once so closures below can check SE_trinity accessions.
        // Uses a Map<accession, action> for O(1) lookup.
        def blPath = file("${launchDir}/rnaseq_blacklist.csv")
        def blMap = blPath.exists()
            ? blPath.readLines().drop(1)
                  .findAll { it.trim() && !it.startsWith('#') }
                  .collectEntries { line ->
                      def cols = line.split(',')
                      cols.size() >= 4 ? [(cols[0].trim()): cols[3].trim()] : [:]
                  }
            : [:]

        // csvHasPE: CSV has at least one PAIRED accession not blocked or overridden to SE.
        def csvHasPE = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                if (cols.size() < 3) return false
                def layout = cols.size() > 5 ? cols[5].trim() : 'PAIRED'
                def action = blMap.get(cols[2].trim(), '')
                layout == 'PAIRED' && action != 'skip' && action != 'SE_trinity'
            }
        }

        // csvHasSEtrinity: CSV has PAIRED accessions overridden to SE via SE_trinity blacklist.
        // These bypass the enable_single_end gate — they are a manual per-accession override.
        def csvHasSEtrinity = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                cols.size() >= 3 && blMap.get(cols[2].trim(), '') == 'SE_trinity'
            }
        }

        // csvHasSingleLayout: CSV has at least one genuine SINGLE-layout accession.
        // Only active when enable_single_end=true.
        def csvHasSingleLayout = { csv ->
            csv.readLines().drop(1).findAll { it.trim() }.any { line ->
                def cols = line.split(',')
                cols.size() > 5 && cols[5].trim() == 'SINGLE' && blMap.get(cols[2].trim(), '') != 'skip'
            }
        }

        // Three-way branch:
        //   has_pe  → SRA_FETCH  (PE wins; SE_trinity entries ignored here, handled by SRA_FETCH)
        //   has_se  → SRA_FETCH_SE (SE_trinity always; SINGLE layout only if enable_single_end)
        //   no_data → WRITE_EMPTY_READS
        def branched_sra = sra_query_results
            .branch {
                has_pe: csvHasPE.call(it[1])
                has_se: csvHasSEtrinity.call(it[1]) ||
                        (params.enable_single_end.toBoolean() && csvHasSingleLayout.call(it[1]))
                no_data: true
            }

        SRA_FETCH(branched_sra.has_pe)
        SRA_FETCH_SE(branched_sra.has_se)
        WRITE_EMPTY_READS(branched_sra.no_data.map { stag, _csv -> stag })

        // Accessions found to be single-end (non-empty _1, no _2) during the PE fetch are
        // recorded as blacklist-ready rows; merge all per-task notes into one reviewable
        // file at the project root. Add these to rnaseq_blacklist.csv as SE_trinity and
        // rerun to route them through SRA_FETCH_SE.
        SRA_FETCH.out.se_candidates
            .collectFile(name: 'rnaseq_se_candidates.csv', storeDir: launchDir, newLine: false)

        // Accessions whose download failed outright (pfd + EBI FTP both produced nothing)
        // are recorded in rnaseq_blacklist.csv column order so they can be reviewed and
        // pasted straight into the blacklist as skip entries; merged into one file at root.
        SRA_FETCH.out.blacklist_candidates
            .collectFile(name: 'rnaseq_blacklist_candidates.csv', storeDir: launchDir, newLine: false)

        reads_ch = SRA_FETCH.out.reads
            .mix(SRA_FETCH_SE.out.reads)
            .mix(WRITE_EMPTY_READS.out.reads)

        if (!params.stop_after_sra_fetch.toBoolean()) {
        // Build per-assembly channel keyed by species_tag with SRA reads joined.
        // reads_ch is now a 4-tuple: (species_tag, r1, r2, se)
        def assembly_with_reads = normalGenomeCh
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                def species_tag = species.replaceAll(/\s+/, '_')
                tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
            }
            .combine(reads_ch, by: 0)
        // assembly_with_reads tuple: (species_tag, out, asmid, species, strain, locustag,
        //                             busco, hlen, ttable, genome_fa, r1, r2, se)

        // RNASEQ_PREPARE: run funannotate train --stop_after_trinity once per species on
        // the representative assembly, then cache the Trinity-GG FASTA in rnaseq_data/
        // so all other strains share it. Normalized reads stay in rnaseq_reads/ (SRA_FETCH storeDir).
        // pasa.gff3 is NOT produced here (--stop_after_trinity stops before PASA);
        // it is produced by FUNANNOTATE_TRAIN for every strain including the representative.
        // Species whose representative r1 and se are both zero-length skip RNASEQ_PREPARE
        // entirely; an empty trinity FASTA is written locally without submitting a SLURM job.
        //
        // "Representative" here is abinitioReuseMap's is_representative pick (ANI+BUSCO),
        // the same strain FUNANNOTATE_PREDICTION uses for ab-initio reuse -- not an
        // arbitrary first-in-group assembly. Falls back to index 0 when the map has no
        // assignment for this species (see module header comment).
        //
        // rnaseqRepOverride (species_tag -> out) takes priority over that ANI+BUSCO pick
        // when present: ANI+BUSCO optimizes for assembly quality, not for which genome
        // the species' actual RNA-seq reads align to, and the two occasionally disagree
        // badly (Ascochyta_rabiei's ANI pick assembles fine but its real RNA-seq
        // genome-guided-Trinity's down to 7 transcripts against it; GCF_004011695.2/Me14
        // works). See loadRnaseqRepresentativeOverride() and
        // scripts/pick_rnaseq_representative_override.py.
        def repr_ch = assembly_with_reads
            .groupTuple(by: 0)
            .map { species_tag, outs, asmids, species_list, strains, locustags,
                   buscos, hlens, ttables, genomes, r1s, r2s, ses ->
                def overrideOut = rnaseqRepOverride[species_tag]
                def repIdx = overrideOut
                    ? outs.findIndexOf { out -> out == overrideOut }
                    : outs.findIndexOf { out -> abinitioReuseMap[out]?.is_representative }
                if (overrideOut && repIdx < 0) {
                    log.warn "rnaseq_representative_override names out '${overrideOut}' for " +
                              "${species_tag} but no such assembly is in this run -- falling " +
                              "back to the ANI+BUSCO pick"
                }
                def i = repIdx >= 0 ? repIdx : 0
                tuple(species_tag, outs[i], asmids[i], species_list[i], strains[i],
                      locustags[i], buscos[i], hlens[i], ttables[i], genomes[i], r1s[i], r2s[i], ses[i])
            }

        def repr_branched = repr_ch.branch {
            has_reads: it[10].size() > 0 || it[12].size() > 0  // r1=[10] or se=[12]
            no_reads:  true
        }

        RNASEQ_PREPARE(repr_branched.has_reads)

        // For species with no RNA-seq reads, write an empty trinity FASTA to rnaseq_data/
        // in the driver process (no SLURM job) and emit it directly as a shared channel item.
        def empty_shared_ch = repr_branched.no_reads
            .map { species_tag, _out, _asmid, _sp, _st, _lt, _bl, _hl, _tt, _gfa, _r1, _r2, _se ->
                def empty_fa = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
                if (!empty_fa.exists()) empty_fa.text = ''
                tuple(species_tag, empty_fa)
            }

        // ── Genome-guided Trinity fallback: standalone (non-GG) Trinity ──────────
        // Genome-guided Trinity occasionally collapses to a near-empty assembly when
        // the representative genome is a poor match for the actual RNA-seq (see
        // rnaseq_representative_override.csv for the "pick a better genome" fix; this
        // is the "just run Trinity without the genome" fix). Reuses
        // train_min_trinity_transcripts -- the same guard FUNANNOTATE_TRAIN already
        // applies -- so a single threshold governs both checks, and setting it to 0
        // disables this fallback too.
        COUNT_TRINITY_TRANSCRIPTS(RNASEQ_PREPARE.out.shared)

        def prepare_reads_ch = repr_branched.has_reads
            .map { species_tag, _out, _asmid, _sp, _st, _lt, _bl, _hl, _tt, _gfa, r1, r2, se ->
                tuple(species_tag, r1, r2, se)
            }

        // (species_tag, r1, r2, se, trinity_fa, n_transcripts)
        def gg_branched = prepare_reads_ch
            .combine(COUNT_TRINITY_TRANSCRIPTS.out.counted, by: 0)
            .map { species_tag, r1, r2, se, trinity_fa, count_file ->
                tuple(species_tag, r1, r2, se, trinity_fa, count_file.text.trim() as int)
            }
            .branch {
                low: (it[5] as int) < (params.train_min_trinity_transcripts as int)
                ok:  true
            }

        TRINITY_STANDALONE(gg_branched.low.map { species_tag, r1, r2, se, _tf, _n ->
            tuple(species_tag, r1, r2, se)
        })

        def gg_ok_shared_ch = gg_branched.ok
            .map { species_tag, _r1, _r2, _se, trinity_fa, _n -> tuple(species_tag, trinity_fa) }

        def shared_ch = gg_ok_shared_ch
            .mix(TRINITY_STANDALONE.out.shared)
            .mix(empty_shared_ch)

        // ── ANI-tiered gating for shared-Trinity training ─────────────────────
        // Non-representative strains reuse the SAME Trinity-GG assembly built
        // from the representative strain's own reads/genome (RNASEQ_PREPARE), so
        // PASA validation identity on re-alignment reflects real strain-to-strain
        // divergence, not assembly error (see FUNANNOTATE_TRAIN module header and
        // .living/learnings.md 2026-08-17 for the Hansenula_anomala_K955 case
        // that surfaced this). Tier by the ANI already computed in
        // PICK_REPRESENTATIVE_STRAIN (abinitioReuseMap[out].ani_to_representative,
        // 100.0 for the representative itself):
        //   ANI >= 97%          -> 'stringent' (PASA defaults; near-identical strains)
        //   90% <= ANI < 97%    -> 'relaxed'   (pasa_shared_* thresholds)
        //   ANI < 90% / missing -> 'skip'      (skani/mash/sourmash returning no
        //                          value is itself a divergence signal, not "not
        //                          yet computed" -- route to ab-initio-only)
        // Tiering only applies when abinitioReuseMap actually has data for this
        // run (share_abinitio_params=true and the CSV exists); when the ANI-reuse
        // feature is off entirely there is no per-strain signal to tier on, so
        // every shared-Trinity strain falls back to 'relaxed' unconditionally.
        def aniSystemActive = !abinitioReuseMap.isEmpty()
        def pasaTierFor = { String out ->
            if (!aniSystemActive) return 'relaxed'
            def ani = abinitioReuseMap[out]?.ani_to_representative
            if (ani == null) return 'skip'
            if (ani >= 97.0) return 'stringent'
            if (ani >= 90.0) return 'relaxed'
            return 'skip'
        }

        // Join shared Trinity from rnaseq_data back to every assembly for FUNANNOTATE_TRAIN.
        // Normalized reads (r1/r2/se) come from SRA_FETCH/SRA_FETCH_SE via assembly_with_reads.
        def train_input = assembly_with_reads
            .combine(shared_ch, by: 0)
            .map { species_tag, out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, se, trinity_fa ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, se, trinity_fa, pasaTierFor.call(out as String))
            }
        // train_input tuple indices: out=0,asmid=1,sp=2,st=3,lt=4,bl=5,hl=6,tt=7,
        //                            genome_fa=8, r1=9, r2=10, se=11, trinity_fa=12, pasa_tier=13

        // ── Hybrid-cross species: composite parent-transcript evidence ────────
        // (nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md). Builds one composite
        // Trinity FASTA per hybrid species_tag (never per strain) by concatenating
        // its PARENT species' own Trinity-GG assemblies -- free, since every parent
        // species already builds its own via the normal RNASEQ_PREPARE path above.
        // Mixed into train_input below (same tuple shape) so the existing
        // ani_skip/has_rnaseq/no_rnaseq branching, staleness filtering, and
        // FUNANNOTATE_TRAIN invocation below need zero duplication for hybrids.
        def hybrid_strains_ch = hybridGenomeCh
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, _taxonid ->
                def species_tag = species.replaceAll(/\s+/, '_')
                tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
            }

        // One composite build per hybrid species_tag: dedupe by species_tag, keep
        // the species string (needed for the genus-fallback below -- GENUS itself
        // isn't threaded through predict_genome_ch's tuple shape, but every hybrid
        // in this dataset has GENUS == the first token of its SPECIES string, so
        // this is a safe proxy rather than a wider tuple-shape change).
        def hybrid_species_tags_ch = hybrid_strains_ch
            .map { species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa ->
                tuple(species_tag, species)
            }
            .unique { it[0] }

        // shared_ch (this run's freshly-built + already-cached non-hybrid Trinity
        // outputs) collected into a single offline map -- a deliberate
        // synchronization point, so every hybrid composite build waits for all of
        // THIS run's normal-species Trinity builds to finish first. Hybrids are
        // lower priority, so running last is acceptable, and this sidesteps
        // fragile per-parent channel joins entirely.
        def parentTrinityMap_ch = shared_ch
            .map { species_tag, trinity_fa -> [species_tag, trinity_fa] }
            .toList()
            .map { pairs -> pairs.collectEntries { st, fa -> [(st): fa] } }

        // evidence_source tracks WHICH path resolved the composite -- 'parents' (a
        // real, species-specific match via hybrid_parentage.csv) vs 'genus_fallback'
        // (no parent-specific data, fell back to any genus-mate's Trinity). Carried
        // through to FUNANNOTATE_TRAIN so it can pick the right PASA identity tier:
        // parent-matched evidence aligns near-identically to its own subgenome copy
        // (~99-100%, same lineage), so it needs a STRICTER tier than genus-fallback
        // evidence, whose divergence is genuinely wide open -- see "New PASA tier:
        // composite" in nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md (reversed
        // from the original single loose tier per bioinformatics review 2026-08-28:
        // a loose threshold on parent-matched evidence risks the WRONG parent's
        // transcript cross-mapping onto a subgenome copy it doesn't belong to,
        // merging homeologs into chimeric gene models).
        def hybrid_parent_lookup_ch = hybrid_species_tags_ch
            .combine(parentTrinityMap_ch)
            .map { species_tag, species, thisRunMap ->
                def parentTags = hybridParentage[species_tag] ?: []
                def parentFastas = []
                parentTags.each { pt ->
                    def fa = thisRunMap[pt]
                    if (!fa) {
                        def onDisk = file("${launchDir}/rnaseq_data/${pt}.trinity-GG.fasta")
                        if (onDisk.exists() && onDisk.size() > 0) fa = onDisk
                    }
                    if (fa && fa.size() > 0) parentFastas << fa
                }
                def evidence_source = 'parents'
                if (parentFastas.isEmpty()) {
                    // Genus-wide fallback: every non-hybrid Trinity under the same
                    // GENUS (approximated from `species`'s first token -- see above).
                    // Excludes every OTHER hybrid species_tag by checking
                    // hybridParentage.containsKey() on the candidate's own tag (NOT a
                    // "_x_" filename-substring heuristic, which would incorrectly
                    // admit a hypothetical hybrid named without an "x" token, e.g. a
                    // formal "Genus ×species" nothospecies with no ASCII x at all).
                    evidence_source = 'genus_fallback'
                    def genus = species.split(/\s+/)[0]
                    def dataDir = file("${launchDir}/rnaseq_data")
                    if (dataDir.exists()) {
                        dataDir.toFile().listFiles()?.each { f ->
                            def name = f.getName()
                            if (name.endsWith('.trinity-GG.fasta') && f.length() > 0) {
                                def candidateTag = name - '.trinity-GG.fasta'
                                if (candidateTag.startsWith("${genus}_") && !hybridParentage.containsKey(candidateTag)) {
                                    parentFastas << file(f.toString())
                                }
                            }
                        }
                    }
                }
                tuple(species_tag, parentFastas, evidence_source)
            }

        def hybrid_parent_branched = hybrid_parent_lookup_ch.branch {
            has_parents: it[1].size() > 0
            no_parents:  true
        }
        BUILD_HYBRID_COMPOSITE_TRINITY(hybrid_parent_branched.has_parents)

        // No parents and no genus-mates found at all -- last-resort true skip
        // (ab-initio-only), reached via the existing no_rnaseq branch below (empty
        // trinity_fa + empty r1/se already routes there, no special-casing needed).
        // evidence_source is moot here (tier resolves to 'skip' regardless of it,
        // since trinity_fa is empty) but kept for tuple-shape consistency with the
        // has_parents branch.
        def hybrid_empty_shared_ch = hybrid_parent_branched.no_parents
            .map { species_tag, _fastas, evidence_source ->
                def empty_fa = file("${launchDir}/rnaseq_data/${species_tag}.composite-parents.trinity-GG.fasta")
                if (!empty_fa.exists()) empty_fa.text = ''
                tuple(species_tag, empty_fa, evidence_source)
            }
        def hybrid_shared_ch = BUILD_HYBRID_COMPOSITE_TRINITY.out.shared.mix(hybrid_empty_shared_ch)

        // Hybrid strains never queried/fetched their own reads (see hybridGenomeCh
        // above) -- 0-byte placeholders via a dedicated storeDir (see
        // WRITE_EMPTY_HYBRID_READS module header for why it's not just an alias of
        // WRITE_EMPTY_READS), so hybrid rows slot into the identical train_input
        // tuple shape as everything else.
        WRITE_EMPTY_HYBRID_READS(hybrid_species_tags_ch.map { it[0] })
        def hybrid_assembly_with_reads = hybrid_strains_ch
            .combine(WRITE_EMPTY_HYBRID_READS.out.reads, by: 0)

        def hybrid_train_input = hybrid_assembly_with_reads
            .combine(hybrid_shared_ch, by: 0)
            .map { species_tag, out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, se, trinity_fa, evidence_source ->
                def tier = trinity_fa.size() == 0 ? 'skip' : (evidence_source == 'genus_fallback' ? 'composite_fallback' : 'composite')
                tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, se, trinity_fa, tier)
            }

        train_input = train_input.mix(hybrid_train_input)

        // Branch on r1 (idx 9), se (idx 11), or trinity_fa (idx 12) sizes; ani_skip
        // (checked first) catches shared-Trinity rows whose ANI tier says the
        // representative's transcriptome is too divergent (or unmeasured) to trust --
        // these bypass FUNANNOTATE_TRAIN exactly like genuinely RNA-seq-less strains.
        def branched = train_input.branch {
            ani_skip:   it[13] == 'skip' && it[12].size() > 0
            has_rnaseq: it[9].size() > 0 || it[11].size() > 0 || it[12].size() > 0
            no_rnaseq:  true
        }
        def predict_no_rnaseq = branched.no_rnaseq.mix(branched.ani_skip)
            .map { out, asmid, sp, st, lt, bl, hl, tt, genome_fa, _r1, _r2, _se, _tf, _tier ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
            }

        // Skip TRAIN at the channel level when pasa.gff3 already exists and is non-empty,
        // UNLESS the rnaseq reads or trinity FASTA is newer than the existing prediction GBK
        // (staleRnaseq), in which case we re-run training so predict can be refreshed too.
        def train_todo = branched.has_rnaseq.filter { out, _a, sp, _st, _lt, _bl, _hl, _tt, _gfa, _r1, _r2, _se, _tf, _tier ->
            def gff3 = file("${params.training_target}/${out}/training/funannotate_train.pasa.gff3")
            !gff3.exists() || gff3.size() == 0 || staleRnaseq(out as String, sp as String)
        }
        def train_done = branched.has_rnaseq
            .filter { out, _a, sp, _st, _lt, _bl, _hl, _tt, _gfa, _r1, _r2, _se, _tf, _tier ->
                def gff3 = file("${params.training_target}/${out}/training/funannotate_train.pasa.gff3")
                gff3.exists() && gff3.size() > 0 && !staleRnaseq(out as String, sp as String)
            }
            .map { out, asmid, sp, st, lt, bl, hl, tt, genome_fa, _r1, _r2, _se, _tf, _tier ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
            }
        FUNANNOTATE_TRAIN(train_todo)
        predict_input_ch = FUNANNOTATE_TRAIN.out.predict_input.mix(train_done).mix(predict_no_rnaseq)

        // Audit trail for composite/composite_fallback-tier graceful degrades --
        // one reviewable file, same convention as rnaseq_se_candidates.csv/
        // rnaseq_blacklist_candidates.csv above. Empty when nothing degraded.
        FUNANNOTATE_TRAIN.out.composite_failed
            .collectFile(name: 'composite_train_failed.tsv', storeDir: launchDir,
                         keepHeader: true, skip: 1)
        } // end if (!params.stop_after_sra_fetch)
        } // end if (!params.stop_after_sra_query)
    } else {
        predict_input_ch = predict_genome_ch
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, _taxonid ->
                tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
            }
    }

    emit:
    predict_input = predict_input_ch
    reads = reads_ch
}
