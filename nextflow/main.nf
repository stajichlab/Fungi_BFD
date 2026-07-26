#!/usr/bin/env nextflow

//
// BigFungiData pipeline entry point.
//
// Routes -profile to the appropriate workflow in workflows/.
// All process definitions live in modules/; see workflows/ for orchestration.
//
// Usage:
//   nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD -resume
//

include { BFD } from './workflows/BFD.nf'

workflow {
    BFD()
}
