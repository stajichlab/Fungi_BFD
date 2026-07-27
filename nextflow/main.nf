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

params.pipeline = 'BFD'

workflow {
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
    else {
        error "--pipeline must be one of: BFD, compare_ani, query_ani, funannotate (got '${params.pipeline}')"
    }
}
