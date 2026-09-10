#!/usr/bin/env nextflow
//
// Standalone real (non-stub) smoke test of the GENEMARK_RUN module wired
// into FUNANNOTATE_PREDICTION.nf -- invokes the SAME module file, not a
// hand-copied duplicate. Lives at nextflow/ TOP LEVEL, sibling to main.nf --
// NOT under a subdirectory (tests/manual/, analysis/, etc.) -- specifically
// so `workflow.projectDir` resolves to nextflow/ exactly like the real
// entry point (nextflow/main.nf) does: it resolves to whichever script was
// originally launched's own containing directory, not a fixed repo root.
// GENEMARK_RUN's ET branch resolves the vendored filterIntronsFindStrand.pl
// via "${workflow.projectDir}/bin/vendor/...", and a copy of this test
// script living in a subdirectory would silently resolve that to the wrong
// path -- hit exactly this on the first two attempts (analysis/genemark_run_validation/
// and nextflow/tests/manual/), see analysis/genemark_run_validation/GENEMARK_RUN_VALIDATION.md
// and .living/learnings.md (2026-08-12, "workflow.projectDir resolves
// per-launched-script" entry).
//
// Usage:
//   nextflow run nextflow/genemark_run_smoke.nf \
//       -c <a process-resource-override config> \
//       --genome_fa <path> --out <name> --species "..." --strain "..." \
//       --transl_table 1 --mode ES|ET --training_bam <path, ET only> \
//       --force_independent false [--shared_mod <path>] [--out_dir <path>]
//
nextflow.enable.dsl = 2

params.max_intronlen    = 3000
params.genome_fa        = null
params.out               = null
params.asmid             = ''
params.species           = ''
params.strain             = ''
params.transl_table       = '1'
params.mode               = 'ES'
params.training_bam       = ''
params.shared_mod         = ''
params.force_independent  = 'false'
// Where to copy results -- defaults next to this script for convenience,
// override with --out_dir to land them in analysis/genemark_run_validation/outputs/.
params.out_dir            = "${projectDir}/manual_outputs"

include { GENEMARK_RUN } from './modules/funannotate/predict/GENEMARK_RUN/main.nf'

workflow {
    if (!params.genome_fa || !params.out) {
        error "Usage: --genome_fa <path> --out <name> [--species ... --strain ... --transl_table ... --mode ES|ET --training_bam ... --shared_mod ... --force_independent true|false --out_dir ...]"
    }
    def genemark_in = Channel.of(
        tuple(params.out, params.asmid, params.species, params.strain,
              file(params.genome_fa), params.transl_table,
              params.mode, params.training_bam,
              params.force_independent, params.shared_mod, 'false')
    )
    GENEMARK_RUN(genemark_in)
    def outDir = file(params.out_dir)
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
