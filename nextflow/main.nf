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

params.pipeline = 'BFD'

// ANI_REPRESENTATIVE_SELECT (funannotate + ani profiles) needs pip cache +
// python_packages dirs for DuckDB install. Singularity hard-FATALs ("mount
// source ... doesn't exist") if a bind source path is missing, and the bind
// is applied when the first task's container launches -- before any process
// body can run. Called at the top of workflow{} below, before any pipeline
// is dispatched, so the dirs exist before any container starts. A function
// declaration is allowed at DSL2 top level (a bare statement there is not),
// and it's plain Groovy in the .nf script rather than nextflow.config, so it
// isn't a params/config key and neither ConfigValidator nor nf-schema's
// validateParameters flags it. Idempotent, so -resume is safe.
def setupBindSourceDirs() {
    ['work/ANI/pip_cache', 'work/ANI/python_packages'].each { new File("${launchDir}/${it}").mkdirs() }
}

workflow {
    setupBindSourceDirs()

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
    else {
        error "--pipeline must be one of: BFD, compare_ani, query_ani, funannotate, earlgrey_mask, comparative, phyling, backfill_abinitio, setup_symlinks (got '${params.pipeline}')\n" +
              "  BFD           — functional annotation + genome stats\n" +
              "  funannotate    — gene prediction + annotation (modular: modules/funannotate/)\n" +
              "  compare_ani    — all-vs-all ANI clustering\n" +
              "  query_ani      — ANI query against existing sketches\n" +
              "  earlgrey_mask  — EarlGrey repeat masking\n" +
              "  comparative    — comparative genomics clustering\n" +
              "  phyling        — PHYling phylogenomics\n" +
              "  backfill_abinitio — sweep: backfill shared ab-initio params for already-predicted representatives\n" +
              "  setup_symlinks — only run SETUP_SYMLINKS (input/ symlinking), skip all downstream steps"
    }
}
