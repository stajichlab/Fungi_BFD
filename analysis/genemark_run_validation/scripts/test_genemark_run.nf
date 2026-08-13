#!/usr/bin/env nextflow
//
// Standalone real (non-stub) smoke test of the GENEMARK_RUN module wired
// into FUNANNOTATE_PREDICTION.nf -- invokes the SAME module file, not a
// hand-copied duplicate, so this tests the actual production code path.
// Self-contained: no profile_funannotate.config include, just the params
// GENEMARK_RUN's script block actually reads (params.max_intronlen).
//
// Usage:
//   nextflow run analysis/genemark_run_validation/scripts/test_genemark_run.nf \
//       --genome_fa <path> --out <name> --species "..." --strain "..." \
//       --transl_table 1 --shared_mod '' --force_independent false \
//       -process.executor local
//
nextflow.enable.dsl = 2

params.max_intronlen    = 3000
params.genome_fa        = null
params.out               = null
params.asmid             = ''
params.species           = ''
params.strain             = ''
params.transl_table       = '1'
params.shared_mod         = ''
params.force_independent  = 'false'

include { GENEMARK_RUN } from '../../../nextflow/modules/funannotate/predict/GENEMARK_RUN/main.nf'

workflow {
    if (!params.genome_fa || !params.out) {
        error "Usage: --genome_fa <path> --out <name> [--species ... --strain ... --transl_table ... --shared_mod ... --force_independent true|false]"
    }
    def genemark_in = Channel.of(
        tuple(params.out, params.asmid, params.species, params.strain,
              file(params.genome_fa), params.transl_table,
              params.force_independent, params.shared_mod)
    )
    GENEMARK_RUN(genemark_in)
    def outDir = file("${projectDir}/../outputs")
    outDir.mkdirs()
    GENEMARK_RUN.out.gtf.subscribe { out, gtf ->
        def dest = outDir.resolve("${out}.genemark.gtf")
        gtf.copyTo(dest)
        println "[RESULT] gtf: out=${out} lines=${gtf.readLines().size()} -> ${dest}"
    }
    GENEMARK_RUN.out.mod.subscribe { out, sp, mod ->
        def dest = outDir.resolve("${out}.genemark.mod")
        mod.copyTo(dest)
        println "[RESULT] mod: out=${out} species=${sp} -> ${dest}"
    }
}
