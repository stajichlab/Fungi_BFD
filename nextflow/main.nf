#!/usr/bin/env nextflow

//
// BigFungiData pipeline entry point.
//
// Process definitions live in modules/; orchestration lives in workflows/ and
// subworkflows/ (see REFACTOR_NEXTFLOW_PLAN.md).
//
// Pipeline selection is by --pipeline, NOT by -entry: Nextflow's strict parser
// rejects -entry outright ("the `-entry` option is not supported with the strict
// parser -- use a param to run a named workflow from the entry workflow"), so a
// single entry workflow dispatching on a param is the forward-compatible form.
//
// Usage:
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD       -resume
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile ani       --pipeline compare_ani -resume
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile ani_query --pipeline query_ani   -resume
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile funannotate --pipeline funannotate -resume
//
// --pipeline defaults to BFD so existing `-profile BFD` invocations keep working.
//

include { BFD }         from './workflows/BFD.nf'
include { COMPARE_ANI } from './workflows/compare_ANI.nf'
include { QUERY_ANI }   from './workflows/query_ANI.nf'
include { FUNANNOTATE } from './workflows/funannotate.nf'
include { EARLGREY_MASK } from './workflows/earlgrey_mask.nf'
include { COMPARATIVE }  from './workflows/comparative_genomics.nf'
include { PHYling }      from './subworkflows/local/PHYLING_ALIGN.nf'
include { BACKFILL_ABINITIO } from './workflows/backfill_abinitio.nf'
include { SETUP_SYMLINKS_ONLY } from './workflows/setup_symlinks.nf'
include { PARALOGOSCOPE_RUN } from './workflows/paralogoscope.nf'

// No default: an omitted --pipeline used to silently fall back to BFD (via the
// line this replaced), which runs BFD() against whatever profile/params are
// actually active. If that profile doesn't define BFD.nf's own params (e.g.
// -profile funannotate lacks run_setup), it doesn't fail here -- it fails
// deep inside BFD.nf with a "Cannot invoke method toBoolean() on null object"
// NullPointerException that gives no hint --pipeline was ever the problem.
// Requiring --pipeline explicit turns that into the clear error below instead.

// ANI_REPRESENTATIVE_SELECT (funannotate + ani profiles) needs a pip cache
// dir for DuckDB install. Singularity hard-FATALs ("mount source ... doesn't
// exist") if a bind source path is missing, and the bind is applied when the
// first task's container launches -- before any process body can run.
// Called at the top of workflow{} below, before any pipeline is dispatched,
// so the dir exists before any container starts. A function declaration is
// allowed at DSL2 top level (a bare statement there is not), and it's plain
// Groovy in the .nf script rather than nextflow.config, so it isn't a
// params/config key and neither ConfigValidator nor nf-schema's
// validateParameters flags it. Idempotent, so -resume is safe.
// (pip's --target install dir is node-local \$SCRATCH now, not a bigdata dir
// -- see profile_ANI.config -- so it needs no pre-created bind source here.)
def setupBindSourceDirs() {
    ['work/ANI/pip_cache'].each { new File("${launchDir}/${it}").mkdirs() }
}

workflow {
    setupBindSourceDirs()

    if (!params.pipeline) {
        error "--pipeline is required (no default -- an omitted --pipeline used to " +
              "silently fall back to BFD and fail confusingly deep inside BFD.nf).\n" +
              "  --pipeline must be one of: BFD, compare_ani, query_ani, funannotate, " +
              "earlgrey_mask, comparative, phyling, backfill_abinitio, setup_symlinks, paralogoscope"
    }

    // if/else rather than switch: the strict parser rejects switch statements.
    def pipeline = (params.pipeline as String).toLowerCase()

    if (pipeline == 'bfd') {
        BFD()
    }
    else if (pipeline == 'compare_ani') {
        COMPARE_ANI()
    }
    else if (pipeline == 'query_ani') {
        QUERY_ANI()
    }
    else if (pipeline == 'funannotate') {
        FUNANNOTATE()
    }
    else if (pipeline == 'earlgrey_mask') {
        EARLGREY_MASK()
    }
    else if (pipeline == 'comparative') {
        COMPARATIVE()
    }
    else if (pipeline == 'phyling') {
        PHYling()
    }
    else if (pipeline == 'backfill_abinitio') {
        BACKFILL_ABINITIO()
    }
    else if (pipeline == 'setup_symlinks') {
        SETUP_SYMLINKS_ONLY()
    }
    else if (pipeline == 'paralogoscope') {
        PARALOGOSCOPE_RUN()
    }
    else {
        error "--pipeline must be one of: BFD, compare_ani, query_ani, funannotate, earlgrey_mask, comparative, phyling, backfill_abinitio, setup_symlinks, paralogoscope (got '${params.pipeline}')\n" +
              "  BFD           — functional annotation + genome stats\n" +
              "  funannotate    — gene prediction + annotation (modular: modules/funannotate/)\n" +
              "  compare_ani    — all-vs-all ANI clustering\n" +
              "  query_ani      — ANI query against existing sketches\n" +
              "  earlgrey_mask  — EarlGrey repeat masking\n" +
              "  comparative    — comparative genomics clustering\n" +
              "  phyling        — PHYling phylogenomics\n" +
              "  backfill_abinitio — sweep: backfill shared ab-initio params for already-predicted representatives\n" +
              "  setup_symlinks — only run SETUP_SYMLINKS (input/ symlinking), skip all downstream steps\n" +
              "  paralogoscope — per-species WGD duplication dating via wgd (dmd -> ksd [+ syn])"
    }
}
