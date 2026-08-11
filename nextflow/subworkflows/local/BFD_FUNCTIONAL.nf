//
// BFD_FUNCTIONAL — the nine per-protein functional-annotation tools and their merges.
//
// Each tool is independently gated by params.run_<tool>. Every RUN_* process is
// storeDir-cached, so a genome whose outputs already exist is skipped without
// -resume; clearIfStale() first deletes any cached output older than the input
// proteins file, so a re-annotated genome forces a re-run of every tool.
//
// Merge scope (see MERGE_MODE below):
//   glob mode  — merge every result file in the output dir, including genomes from
//                previous runs. A sync barrier on this run's outputs ensures newly
//                written files are on disk before globbing; if a tool was not run
//                this session the barrier is an immediate Channel.of(true).
//   run mode   — merge only what this run produced.
//

include { RUN_PFAM      } from '../../modules/BFD/PFAM/main.nf'
include { RUN_CAZY      } from '../../modules/BFD/CAZY/main.nf'
include { RUN_MEROPS    } from '../../modules/BFD/MEROPS/main.nf'
include { RUN_SIGNALP   } from '../../modules/BFD/SIGNALP/main.nf'
include { RUN_TMHMM     } from '../../modules/BFD/TMHMM/main.nf'
include { RUN_TARGETP   } from '../../modules/BFD/TARGETP/main.nf'
include { RUN_IDP       } from '../../modules/BFD/IDP/main.nf'
include { RUN_WOLFPSORT } from '../../modules/BFD/WOLFPSORT/main.nf'
include { RUN_PREDGPI   } from '../../modules/BFD/PREDGPI/main.nf'
include { RUN_SWISSPROT } from '../../modules/BFD/SWISSPROT/main.nf'

include { MERGE_PFAM      } from '../../modules/BFD/MERGE_PFAM/main.nf'
include { MERGE_CAZY      } from '../../modules/BFD/MERGE_CAZY/main.nf'
include { MERGE_MEROPS    } from '../../modules/BFD/MERGE_MEROPS/main.nf'
include { MERGE_SWISSPROT } from '../../modules/BFD/MERGE_SWISSPROT/main.nf'
include { BUILD_SWISSPROT_ANNOT } from '../../modules/BFD/BUILD_SWISSPROT_ANNOT/main.nf'
include { MERGE_SIGNALP   } from '../../modules/BFD/MERGE_SIGNALP/main.nf'
include { MERGE_TMHMM     } from '../../modules/BFD/MERGE_TMHMM/main.nf'
include { MERGE_TARGETP   } from '../../modules/BFD/MERGE_TARGETP/main.nf'
include { MERGE_IDP       } from '../../modules/BFD/MERGE_IDP/main.nf'
include { MERGE_WOLFPSORT } from '../../modules/BFD/MERGE_WOLFPSORT/main.nf'
include { MERGE_PREDGPI   } from '../../modules/BFD/MERGE_PREDGPI/main.nf'

include { clearIfStale; gatedGlobIn; hashBucketForType } from '../../modules/common/utils.nf'

// Attach a staleness check to the protein channel for CAZY specifically -- CAZY
// uses a different pre-existing layout (one subdirectory per genome, not
// hash-bucketed) and is out of scope for the genome_stats/function hash-bucket
// reorg (see KNOWN_SUBDIRECTORY_LAYOUT_TYPES in
// scripts/one-off/reorg_genome_stats_hash_buckets.py). Left exactly as before.
def staleGuard(ch, List relOuts) {
    ch.map { meta, prot ->
        clearIfStale(prot, relOuts.collect { rel ->
            file("${params.outdir}/${rel.replace('@', meta.id)}")
        })
        tuple(meta, prot)
    }
}

// Same idea, but for the hash-bucketed function/* tools (everything except
// CAZY): relOuts templates are filename-only (no type prefix, no bucket) --
// both are computed here from `type` and meta.locustag, and must match exactly
// what the real RUN_* module actually writes (see the corresponding module
// under nextflow/modules/BFD/), or clearIfStale()'s check -- and therefore
// storeDir's own existence check -- looks in the wrong place and every genome
// gets needlessly recomputed.
def staleGuardBucketed(ch, String type, List relOuts) {
    ch.map { meta, prot ->
        def bucket = hashBucketForType(type, meta.locustag)
        clearIfStale(prot, relOuts.collect { rel ->
            file("${params.outdir}/${type}/${bucket}/${rel.replace('@', meta.locustag)}")
        })
        tuple(meta, prot)
    }
}

// Root a gated glob in params.outdir. Globs for the hash-bucketed tools use
// "**/*.ext", not "*.ext" -- these directories are hash-bucketed one level
// deep (T-014); a flat glob would silently match nothing for any bucketed
// genome, the same silent-loss failure class already documented for
// gated-glob-vs-live-channel elsewhere in this repo (.living/learnings.md).
// CAZY is the one exception -- its own pre-existing per-genome-subdirectory
// layout (not a hash bucket) already uses "*/*.ext" and is unaffected.
def gatedGlob(sync_ch, String glob) {
    gatedGlobIn(sync_ch, params.outdir, glob)
}

// Pick the merge input for a tool: this run's collected outputs, or a glob over
// everything ever produced, depending on mode. `ran` says whether the RUN step
// executed this session — when it did, its output channel is the sync barrier.
def mergeInput(boolean useGlob, boolean ran, out_ch, String glob) {
    if (!useGlob) {
        return ran ? out_ch.collect() : null
    }
    gatedGlob(ran ? out_ch.collect() : Channel.of(true), glob)
}

workflow BFD_FUNCTIONAL {
    take:
    proteins_ch   // tuple(meta, proteins)
    use_glob      // merge_all is set and no --taxon narrowing is active
    skip_merge

    main:
    def run_pfam      = params.run_pfam.toBoolean()
    def run_cazy      = params.run_cazy.toBoolean()
    def run_merops    = params.run_merops.toBoolean()
    def run_signalp   = params.run_signalp.toBoolean()
    def run_tmhmm     = params.run_tmhmm.toBoolean()
    def run_targetp   = params.run_targetp.toBoolean()
    def run_idp       = params.run_idp.toBoolean()
    def run_wolfpsort = params.run_wolfpsort.toBoolean()
    def run_predgpi   = params.run_predgpi.toBoolean()
    def run_swissprot = params.run_swissprot.toBoolean()
    def run_swissprot_annot = params.run_swissprot_annot.toBoolean()

    // ── Per-species RUN steps ────────────────────────────────────────────────
    if (run_pfam)
        RUN_PFAM(staleGuardBucketed(proteins_ch, 'pfam_hmmscan', ['@.domtblout.gz', '@.tblout.gz']))
    if (run_cazy)
        RUN_CAZY(staleGuard(proteins_ch, ['cazy/@/@.overview.tsv.gz', 'cazy/@/@.cazymes.tsv.gz', 'cazy/@/@.substrates.tsv.gz']))
    if (run_merops)
        RUN_MEROPS(staleGuardBucketed(proteins_ch, 'merops', ['@.blasttab.gz']))
    if (run_swissprot)
        RUN_SWISSPROT(staleGuardBucketed(proteins_ch, 'swissprot', ['@.blasttab.gz']))
    if (run_signalp)
        RUN_SIGNALP(staleGuardBucketed(proteins_ch, 'signalp', ['@.signalp.gff3.gz', '@.signalp.results.txt.gz']))
    if (run_tmhmm)
        RUN_TMHMM(staleGuardBucketed(proteins_ch, 'tmhmm', ['@.tmhmm_short.tsv.gz', '@.tmhmm_results.tsv.gz']))
    if (run_targetp)
        RUN_TARGETP(staleGuardBucketed(proteins_ch, 'targetP', ['@_summary.targetp2.gz']))
    if (run_idp)
        RUN_IDP(staleGuardBucketed(proteins_ch, 'aiupred', ['@.aiupred.txt.gz', '@.idp.csv.gz', '@.idp_summary.csv.gz']))
    if (run_wolfpsort)
        RUN_WOLFPSORT(staleGuardBucketed(proteins_ch, 'wolfpsort', ['@.wolfpsort.results.txt.gz']))
    if (run_predgpi)
        RUN_PREDGPI(staleGuardBucketed(proteins_ch, 'predgpi', ['@.predgpi.gff3.gz']))

    // ── MERGE steps ──────────────────────────────────────────────────────────
    // In glob mode a tool is merged even when it did not run this session, so
    // previously-computed genomes still make it into the table. In run mode a
    // tool that did not run has nothing to merge and is skipped.
    if (!skip_merge) {
        def pfam_in = mergeInput(use_glob, run_pfam, run_pfam ? RUN_PFAM.out.domtbl : null, "pfam_hmmscan/**/*.domtblout.gz")
        if (pfam_in) MERGE_PFAM(pfam_in)

        if (use_glob || run_cazy) {
            def ov_sync = run_cazy ? RUN_CAZY.out.overview.collect() : Channel.of(true)
            def ca_sync = run_cazy ? RUN_CAZY.out.cazymes.collect()  : Channel.of(true)
            MERGE_CAZY(
                use_glob ? gatedGlob(ov_sync, "cazy/*/*.overview.tsv.gz") : ov_sync,
                use_glob ? gatedGlob(ca_sync, "cazy/*/*.cazymes.tsv.gz")  : ca_sync
            )
        }

        def merops_in = mergeInput(use_glob, run_merops, run_merops ? RUN_MEROPS.out.blasttab : null, "merops/**/*.blasttab.gz")
        if (merops_in) MERGE_MEROPS(merops_in)

        def swissprot_in = mergeInput(use_glob, run_swissprot, run_swissprot ? RUN_SWISSPROT.out.blasttab : null, "swissprot/**/*.blasttab.gz")
        if (swissprot_in) MERGE_SWISSPROT(swissprot_in)

        // One-shot build of the SwissProt accession -> annotation table,
        // independent of the per-genome search engine used above (or of
        // run_swissprot itself, so it can be refreshed standalone).
        if (run_swissprot_annot)
            BUILD_SWISSPROT_ANNOT(Channel.of(file(params.swissprot_dat)))

        def signalp_in = mergeInput(use_glob, run_signalp, run_signalp ? RUN_SIGNALP.out.gff3 : null, "signalp/**/*.signalp.gff3.gz")
        if (signalp_in) MERGE_SIGNALP(signalp_in)

        def tmhmm_in = mergeInput(use_glob, run_tmhmm, run_tmhmm ? RUN_TMHMM.out.short_tsv : null, "tmhmm/**/*.tmhmm_short.tsv.gz")
        if (tmhmm_in) MERGE_TMHMM(tmhmm_in)

        def targetp_in = mergeInput(use_glob, run_targetp, run_targetp ? RUN_TARGETP.out.summary : null, "targetP/**/*_summary.targetp2.gz")
        if (targetp_in) MERGE_TARGETP(targetp_in)

        if (use_glob || run_idp) {
            def idp_sync = run_idp ? RUN_IDP.out.idp_csv.collect()         : Channel.of(true)
            def sum_sync = run_idp ? RUN_IDP.out.idp_summary_csv.collect() : Channel.of(true)
            MERGE_IDP(
                use_glob ? gatedGlob(idp_sync, "aiupred/**/*.idp.csv.gz")         : idp_sync,
                use_glob ? gatedGlob(sum_sync, "aiupred/**/*.idp_summary.csv.gz") : sum_sync
            )
        }

        def wolf_in = mergeInput(use_glob, run_wolfpsort, run_wolfpsort ? RUN_WOLFPSORT.out.results : null, "wolfpsort/**/*.wolfpsort.results.txt.gz")
        if (wolf_in) MERGE_WOLFPSORT(wolf_in)

        def predgpi_in = mergeInput(use_glob, run_predgpi, run_predgpi ? RUN_PREDGPI.out.gff3 : null, "predgpi/**/*.predgpi.gff3.gz")
        if (predgpi_in) MERGE_PREDGPI(predgpi_in)
    }
}
