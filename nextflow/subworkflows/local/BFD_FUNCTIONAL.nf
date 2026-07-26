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

include { MERGE_PFAM      } from '../../modules/BFD/MERGE_PFAM/main.nf'
include { MERGE_CAZY      } from '../../modules/BFD/MERGE_CAZY/main.nf'
include { MERGE_MEROPS    } from '../../modules/BFD/MERGE_MEROPS/main.nf'
include { MERGE_SIGNALP   } from '../../modules/BFD/MERGE_SIGNALP/main.nf'
include { MERGE_TMHMM     } from '../../modules/BFD/MERGE_TMHMM/main.nf'
include { MERGE_TARGETP   } from '../../modules/BFD/MERGE_TARGETP/main.nf'
include { MERGE_IDP       } from '../../modules/BFD/MERGE_IDP/main.nf'
include { MERGE_WOLFPSORT } from '../../modules/BFD/MERGE_WOLFPSORT/main.nf'
include { MERGE_PREDGPI   } from '../../modules/BFD/MERGE_PREDGPI/main.nf'

include { clearIfStale } from '../../modules/common/utils.nf'

// Attach a staleness check to the protein channel for one tool: delete any cached
// output older than the proteins file, then pass the tuple through unchanged.
def staleGuard(ch, List relOuts) {
    ch.map { meta, prot ->
        clearIfStale(prot, relOuts.collect { rel ->
            file("${params.outdir}/${rel.replace('@', meta.id)}")
        })
        tuple(meta, prot)
    }
}

// Build a gated glob channel rooted in params.outdir.
//   sync_ch — emits once the RUN step is done (or Channel.of(true))
//   glob    — shell-style glob relative to params.outdir
// Returns a channel of matching non-empty Paths, or empty if none found.
def gatedGlob(sync_ch, String glob) {
    sync_ch
        .flatMap { files("${params.outdir}/${glob}") }
        .filter  { it.size() > 0 }
        .collect()
        .filter  { !it.isEmpty() }
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

    // ── Per-species RUN steps ────────────────────────────────────────────────
    if (run_pfam)
        RUN_PFAM(staleGuard(proteins_ch, ['pfam_hmmscan/@.domtblout.gz', 'pfam_hmmscan/@.tblout.gz']))
    if (run_cazy)
        RUN_CAZY(staleGuard(proteins_ch, ['cazy/@/@.overview.tsv.gz', 'cazy/@/@.cazymes.tsv.gz', 'cazy/@/@.substrates.tsv.gz']))
    if (run_merops)
        RUN_MEROPS(staleGuard(proteins_ch, ['merops/@.blasttab.gz']))
    if (run_signalp)
        RUN_SIGNALP(staleGuard(proteins_ch, ['signalp/@.signalp.gff3.gz', 'signalp/@.signalp.results.txt.gz']))
    if (run_tmhmm)
        RUN_TMHMM(staleGuard(proteins_ch, ['tmhmm/@.tmhmm_short.tsv.gz', 'tmhmm/@.tmhmm_results.tsv.gz']))
    if (run_targetp)
        RUN_TARGETP(staleGuard(proteins_ch, ['targetP/@_summary.targetp2.gz']))
    if (run_idp)
        RUN_IDP(staleGuard(proteins_ch, ['aiupred/@.aiupred.txt.gz', 'aiupred/@.idp.csv.gz', 'aiupred/@.idp_summary.csv.gz']))
    if (run_wolfpsort)
        RUN_WOLFPSORT(staleGuard(proteins_ch, ['wolfpsort/@.wolfpsort.results.txt.gz']))
    if (run_predgpi)
        RUN_PREDGPI(staleGuard(proteins_ch, ['predgpi/@.predgpi.gff3.gz']))

    // ── MERGE steps ──────────────────────────────────────────────────────────
    // In glob mode a tool is merged even when it did not run this session, so
    // previously-computed genomes still make it into the table. In run mode a
    // tool that did not run has nothing to merge and is skipped.
    if (!skip_merge) {
        def pfam_in = mergeInput(use_glob, run_pfam, run_pfam ? RUN_PFAM.out.domtbl : null, "pfam_hmmscan/*.domtblout.gz")
        if (pfam_in) MERGE_PFAM(pfam_in)

        if (use_glob || run_cazy) {
            def ov_sync = run_cazy ? RUN_CAZY.out.overview.collect() : Channel.of(true)
            def ca_sync = run_cazy ? RUN_CAZY.out.cazymes.collect()  : Channel.of(true)
            MERGE_CAZY(
                use_glob ? gatedGlob(ov_sync, "cazy/*/*.overview.tsv.gz") : ov_sync,
                use_glob ? gatedGlob(ca_sync, "cazy/*/*.cazymes.tsv.gz")  : ca_sync
            )
        }

        def merops_in = mergeInput(use_glob, run_merops, run_merops ? RUN_MEROPS.out.blasttab : null, "merops/*.blasttab.gz")
        if (merops_in) MERGE_MEROPS(merops_in)

        def signalp_in = mergeInput(use_glob, run_signalp, run_signalp ? RUN_SIGNALP.out.gff3 : null, "signalp/*.signalp.gff3.gz")
        if (signalp_in) MERGE_SIGNALP(signalp_in)

        def tmhmm_in = mergeInput(use_glob, run_tmhmm, run_tmhmm ? RUN_TMHMM.out.short_tsv : null, "tmhmm/*.tmhmm_short.tsv.gz")
        if (tmhmm_in) MERGE_TMHMM(tmhmm_in)

        def targetp_in = mergeInput(use_glob, run_targetp, run_targetp ? RUN_TARGETP.out.summary : null, "targetP/*_summary.targetp2.gz")
        if (targetp_in) MERGE_TARGETP(targetp_in)

        if (use_glob || run_idp) {
            def idp_sync = run_idp ? RUN_IDP.out.idp_csv.collect()         : Channel.of(true)
            def sum_sync = run_idp ? RUN_IDP.out.idp_summary_csv.collect() : Channel.of(true)
            MERGE_IDP(
                use_glob ? gatedGlob(idp_sync, "aiupred/*.idp.csv.gz")         : idp_sync,
                use_glob ? gatedGlob(sum_sync, "aiupred/*.idp_summary.csv.gz") : sum_sync
            )
        }

        def wolf_in = mergeInput(use_glob, run_wolfpsort, run_wolfpsort ? RUN_WOLFPSORT.out.results : null, "wolfpsort/*.wolfpsort.results.txt.gz")
        if (wolf_in) MERGE_WOLFPSORT(wolf_in)

        def predgpi_in = mergeInput(use_glob, run_predgpi, run_predgpi ? RUN_PREDGPI.out.gff3 : null, "predgpi/*.predgpi.gff3.gz")
        if (predgpi_in) MERGE_PREDGPI(predgpi_in)
    }
}
