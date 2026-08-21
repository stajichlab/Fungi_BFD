//
// FUNANNOTATE_ANNOTATION — optional pre-annotate steps, annotate and update.
//
// Receives the metadata channel from FUNANNOTATE_PREDICTION (genomes
// predicted in THIS run) and builds a postpredict channel from the
// samples CSV (prior-run genomes with current predictions). The two
// sets are disjoint by the predict/postpredict filters, so a plain mix
// needs no dedup.
//
// Optional pre-annotate steps (ANTISMASH, INTERPRO, SIGNALP) split the
// channel into "needs to run" vs "already done". UPDATE joins SRA reads
// from FUNANNOTATE_RNASEQ. ANNOTATE fires only on genomes missing an
// annotate_results GBK.
//

include { ANTISMASH_RUN        } from '../../modules/funannotate/function/ANTISMASH_RUN/main.nf'
include { SETUP_ANTISMASH_DB   } from '../../modules/funannotate/setup/SETUP_ANTISMASH_DB/main.nf'
include { INTERPROSCAN_RUN     } from '../../modules/funannotate/function/INTERPROSCAN_RUN/main.nf'
include { SIGNALP_RUN          } from '../../modules/funannotate/function/SIGNALP_RUN/main.nf'
include { FUNANNOTATE_ANNOTATE } from '../../modules/funannotate/function/FUNANNOTATE_ANNOTATE/main.nf'
include { FUNANNOTATE_UPDATE   } from '../../modules/funannotate/predict/FUNANNOTATE_UPDATE/main.nf'

include { makeSampleTag } from '../../modules/common/utils.nf'
include { gbkResult; staleRnaseq; staleGenome } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_ANNOTATION {
    take:
    metadata_ch            // tuple(out, asmid, sp, st, lt, bl, hl, tt) from FUNANNOTATE_PREDICTION
    reads_ch               // tuple(species_tag, r1, r2, se) from FUNANNOTATE_RNASEQ
    taxonFilter            // row closures, shared with the caller's sample parsing
    asmidFilter
    suppressFilter         // suppressRowFilter(loadSuppressSet()), same instance as the caller's

    main:
    // postpredict: all samples with a completed predict_results/*.gbk, whether
    // produced in this run or a prior one. This is the source for all optional
    // pre-annotate steps and for FUNANNOTATE_ANNOTATE itself.
    def postpredict = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .filter(suppressFilter)
        .map { row ->
            def species       = (row.SPECIES?.trim() ?: '').replaceAll(/['"]/, '')
            def strain        = (row.STRAIN?.trim() ?: '').replaceAll(/['"]/, '').replaceAll(/;.*$/, '').trim().replace(':', ' ')
            def out           = makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def asmid         = row.ASMID?.trim()
            def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco         = row.BUSCO_LINEAGE?.trim()
            def header_length = 24
            def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        // Only genomes whose prediction was already complete AND current in a PRIOR run.
        // This is the exact logical complement of the predict_ch filter, so this set is
        // disjoint from the genomes (re)predicted in THIS run (which arrive via
        // metadata_ch below). Keeping them disjoint means no genome is fed downstream
        // twice and stale genomes correctly wait for the fresh predict.
        .filter { out, asmid, sp, _st, _lt, _bl, _hl, _tt ->
            gbkResult("${params.target}/${out}/predict_results", out as String) != null &&
                !staleRnaseq(out as String, sp as String) &&
                !staleGenome(out as String, asmid as String)
        }

    // annotate_ready_ch threads through optional pre-annotate steps. Each optional
    // step splits the channel into "needs to run" vs "already done", processes the
    // former, then mixes the freshly-completed items back. FUNANNOTATE_ANNOTATE only
    // fires once all requested optional steps are complete for a given sample.
    // Joining ANTISMASH/INTERPRO/SIGNALP output back through postpredict reconstructs
    // the metadata tuple while encoding the dependency edge in the channel DAG.
    //
    // Same-run completion gate: genomes predicted in THIS run flow in via
    // metadata_ch (a real channel edge, so downstream waits for predict to finish),
    // while prior-run genomes flow in via postpredict (available immediately). The
    // two sets are disjoint by the filter above, so a plain mix needs no dedup.
    // NOTE: the optional steps below are still each gated behind their params
    // (run_antismash/interpro/signalp/update/annotate), all of which default to
    // false -- so by default nothing downstream of predict runs.
    //
    // Combined metadata for every genome with a current prediction: prior-run genomes
    // (postpredict) plus genomes predicted in THIS run (metadata_ch). Used both as the
    // annotate source and as the right-hand side of the metadata-reconstruction joins
    // below, so this-run genomes are not dropped when an optional step is enabled.
    // (Reused multiple times, exactly as postpredict was before.)
    def predict_meta = postpredict.mix(metadata_ch)
    def annotate_ready_ch = predict_meta

    if (params.run_antismash.toBoolean()) {
        def as_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            def asDir = file("${params.target}/${out}/antismash_local")
            !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
        }
        def as_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            def asDir = file("${params.target}/${out}/antismash_local")
            asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') }
        }
        // SETUP_ANTISMASH_DB uses storeDir so it runs at most once across all pipeline
        // runs (skipped entirely if params.antismash_databases already holds a
        // populated DB -- see that module's header comment for the version-match
        // assumption when reusing a pre-existing directory).
        SETUP_ANTISMASH_DB()
        def antismash_db_ch = SETUP_ANTISMASH_DB.out.ready.map { params.antismash_databases }
        ANTISMASH_RUN(as_todo.combine(antismash_db_ch))
        def as_completed = ANTISMASH_RUN.out
            .map { out, _files -> tuple(out, 'done') }
            .join(predict_meta)
            .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
        annotate_ready_ch = as_completed.mix(as_done)
    }

    if (params.run_interpro.toBoolean()) {
        def ipr_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            !file("${params.target}/${out}/annotate_misc/iprscan.xml").exists()
        }
        def ipr_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            file("${params.target}/${out}/annotate_misc/iprscan.xml").exists()
        }
        INTERPROSCAN_RUN(ipr_todo)
        def ipr_completed = INTERPROSCAN_RUN.out
            .map { out, _xml -> tuple(out, 'done') }
            .join(predict_meta)
            .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
        annotate_ready_ch = ipr_completed.mix(ipr_done)
    }

    if (params.run_signalp.toBoolean()) {
        def sp_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            !file("${params.target}/${out}/annotate_misc/signalp.results.txt").exists()
        }
        def sp_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            file("${params.target}/${out}/annotate_misc/signalp.results.txt").exists()
        }
        SIGNALP_RUN(sp_todo)
        def sp_completed = SIGNALP_RUN.out
            .map { out, _txt -> tuple(out, 'done') }
            .join(predict_meta)
            .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
        annotate_ready_ch = sp_completed.mix(sp_done)
    }

    if (params.run_update.toBoolean()) {
        if (params.run_sra_fetch.toBoolean()) {
            // UPDATE runs from predict results in parallel with antismash/interpro/signalp.
            // Reads are joined from SRA_FETCH (storeDir-cached, so prior-run reads are reused).
            // The join on upd_signal gates annotate_ready_ch so ANNOTATE waits for UPDATE.
            //
            // remainder: true (was a plain combine(by:0), an implicit inner join) -- a
            // species with no reads_ch entry (RNA-seq acquisition failed/was skipped
            // upstream for that species_tag) must NOT silently vanish from this join.
            // combine(by:0) previously dropped such a species from upd_input entirely,
            // which meant it never got an upd_signal entry either, which meant the
            // second join below (annotate_ready_ch.join(upd_signal)) ALSO dropped it --
            // silently removing an otherwise-fully-predicted genome from ANNOTATE too,
            // not just UPDATE. Found via design review, 2026-08-14.
            // join(by:0, remainder:true) turned out to be the wrong operator here, not
            // just awkward to destructure: Nextflow's join pairs by key like a 1:1
            // zipper, not a broadcast. When multiple strains of the same species (a
            // normal case -- RNA-seq is fetched once per species, shared across
            // strains) hit this join against reads_ch's single per-species entry, only
            // the FIRST matching predict_meta row actually paired with the real reads;
            // the rest silently got remainder-null padding and were logged as "no
            // RNA-seq reads available" even though the species' reads genuinely exist.
            // Caught for real: NRRL_1841 got real reads (and a real UPDATE) while
            // sibling strains B8014/NRRL_756 of the SAME species did not, on a run
            // where rnaseq_reads/Penicillium_citrinum_norm_R1.fastq.gz was present and
            // non-empty for all three. combine(by:0) is the correct operator for
            // fan-out (broadcasts the single reads_ch item to every matching left-side
            // row) but is an inner join -- it drops a species with no reads_ch entry at
            // all, which is exactly the silent-vanishing bug the remainder:true was
            // added for in the first place (see note above). So: use combine(by:0) for
            // the broadcast-matched case, and separately preserve genuinely-unmatched
            // species via an explicit collected-keys membership check rather than
            // relying on join's remainder semantics. Found via real run, 2026-08-18.
            def reads_species_ch = reads_ch.map { st, _r1, _r2, _se -> st }.toList()
            def upd_input_base = predict_meta
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable ->
                    def species_tag = species.replaceAll(/\s+/, '_')
                    tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable)
                }
            def upd_input_matched = upd_input_base
                .combine(reads_ch, by: 0)
                .map { _st, out, asmid, species, strain, locustag, busco, hlen, ttable, r1, r2, _se ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, r1, r2)
                }
            def upd_input_nomatch = upd_input_base
                .combine(reads_species_ch)
                .filter { species_tag, _out, _a, _sp, _st, _lt, _bl, _hl, _tt, matchedKeys -> !(species_tag in matchedKeys) }
                .map { _st, out, asmid, species, strain, locustag, busco, hlen, ttable, _keys ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, null, null)
                }
            def upd_input = upd_input_matched.mix(upd_input_nomatch)
            def hasReads = { r1 -> r1 != null && file(r1 as String).exists() && file(r1 as String).size() > 0 }
            def upd_todo = upd_input.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, r1, _r2 ->
                hasReads.call(r1) && gbkResult("${params.target}/${out}/update_results", out as String) == null
            }
            def upd_done_signal = upd_input
                .filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, r1, _r2 ->
                    hasReads.call(r1) && gbkResult("${params.target}/${out}/update_results", out as String) != null
                }
                .map { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _r1, _r2 -> tuple(out, 'upd') }
            // No reads (never fetched for this species, or genuinely RNA-seq-less):
            // skip UPDATE -- routine, not a fault, so log.warn only (not a report file) --
            // but still emit an upd_signal entry so ANNOTATE isn't blocked on it.
            def upd_skip_signal = upd_input
                .filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, r1, _r2 -> !(hasReads.call(r1)) }
                .map { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _r1, _r2 ->
                    log.warn "FUNANNOTATE_UPDATE skipped for ${out}: no RNA-seq reads available"
                    tuple(out, 'upd')
                }
            FUNANNOTATE_UPDATE(upd_todo)
            def upd_signal = FUNANNOTATE_UPDATE.out
                .map { out, _a, _sp, _st, _lt, _bl, _hl, _tt -> tuple(out, 'upd') }
                .mix(upd_done_signal)
                .mix(upd_skip_signal)
            annotate_ready_ch = annotate_ready_ch
                .join(upd_signal)
                .map { out, asmid, sp, st, lt, bl, hl, tt, _flag -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
        } else {
            log.warn "run_update=true but run_sra_fetch=false; funannotate update skipped (no reads available)"
        }
    }

    if (params.run_annotate.toBoolean()) {
        FUNANNOTATE_ANNOTATE(annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
            gbkResult("${params.target}/${out}/annotate_results", out as String) == null
        })
    }
}
