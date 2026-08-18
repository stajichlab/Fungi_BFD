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

include { SRA_QUERY_BATCH   } from '../../modules/funannotate/rnaseq/SRA_QUERY_BATCH/main.nf'
include { COLLECT_SRA_QUERY } from '../../modules/funannotate/rnaseq/COLLECT_SRA_QUERY/main.nf'
include { WRITE_EMPTY_READS } from '../../modules/funannotate/rnaseq/WRITE_EMPTY_READS/main.nf'
include { SRA_FETCH         } from '../../modules/funannotate/rnaseq/SRA_FETCH/main.nf'
include { SRA_FETCH_SE      } from '../../modules/funannotate/rnaseq/SRA_FETCH_SE/main.nf'
include { RNASEQ_PREPARE    } from '../../modules/funannotate/rnaseq/RNASEQ_PREPARE/main.nf'
include { FUNANNOTATE_TRAIN } from '../../modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf'

include { gbkResult; staleRnaseq } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_RNASEQ {
    take:
    predict_genome_ch   // tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid)
    abinitioReuseMap     // out -> [species, reuse_eligible, is_representative], from loadAbinitioReuseMap()

    main:
    // When SRA is enabled: SRA_FETCH fetches reads once per species; RNASEQ_PREPARE runs
    // funannotate train on the representative assembly and archives Trinity-GG, trimmed, and
    // normalized reads to rnaseq_data/; all other strains run FUNANNOTATE_TRAIN --trinity.
    def predict_input_ch
    def reads_ch = Channel.empty()
    if (params.run_sra_fetch.toBoolean()) {
        // Build per-species input: group assemblies, keep first taxonid per species.
        def sra_input = predict_genome_ch
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
            def sra_batched = sra_input
                .collate(params.sra_query_batch_size)
                .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }) }
            SRA_QUERY_BATCH(sra_batched)
            sra_query_results = SRA_QUERY_BATCH.out.query_results
                .flatten()
                .map { csv -> tuple(csv.baseName.replaceAll(/\.sra_query$/, ''), csv) }
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
        def assembly_with_reads = predict_genome_ch
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
        def repr_ch = assembly_with_reads
            .groupTuple(by: 0)
            .map { species_tag, outs, asmids, species_list, strains, locustags,
                   buscos, hlens, ttables, genomes, r1s, r2s, ses ->
                def repIdx = outs.findIndexOf { out -> abinitioReuseMap[out]?.is_representative }
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

        def shared_ch = RNASEQ_PREPARE.out.shared.mix(empty_shared_ch)

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
        predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
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
