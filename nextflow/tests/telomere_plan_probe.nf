#!/usr/bin/env nextflow
//
// Throwaway probe for nextflow/tests/test_telomere_batch_plan.py. Not part of
// any production workflow. It runs the real planFreq() and batchFreq() from
// BFD_GENOME_STATS.nf on the telomere step's settings (subdir telomeres_v2,
// suffix telomeres.tsv.gz, key asmid) over a folder of fake genome files, and
// prints one marker line per result so the Python test can read them:
//   CACHED:<asmid>\t<cached output path>   genome skipped, passed to the merge
//   BATCH:<asmid>,<asmid>,...              one BATCH_TELOMERES task
//
nextflow.enable.dsl = 2

include { planFreq; batchFreq } from '../subworkflows/local/BFD_GENOME_STATS.nf'

params.genome_dir          = null
params.genome_stats_outdir = null
params.telomere_batch_size = 100

workflow {
    if (!params.genome_dir || !params.genome_stats_outdir) {
        error "--genome_dir and --genome_stats_outdir are required"
    }
    def genomes = Channel.fromPath("${params.genome_dir}/*.fa")
        .map { f -> tuple([asmid: f.baseName], f) }
    def plan = planFreq(genomes, 'telomeres_v2', 'telomeres.tsv.gz', 'asmid')
    plan.cached.subscribe { meta, _g, out, _c -> println "CACHED:${meta.asmid}\t${out}" }
    batchFreq(plan.todo, params.telomere_batch_size as int)
        .subscribe { metas, _gs -> println "BATCH:${metas.collect { m -> m.asmid }.join(',')}" }
}
