#!/usr/bin/env nextflow
//
// Throwaway probe script for the Groovy/Python hash_bucket() parity test
// (nextflow/tests/test_hash_bucket_parity.py). Not part of any production
// workflow -- reads a TSV of `key<TAB>width` pairs and prints
// `BUCKET_RESULT:<key>\t<width>\t<bucket>` for each, so the Python test can
// grep the marker out of Nextflow's own stdout noise and compare against
// genome_stats_paths.py::hash_bucket().
//
nextflow.enable.dsl = 2

include { hashBucket } from '../modules/common/utils.nf'

params.keys_file = null

workflow {
    if (!params.keys_file) {
        error "--keys_file is required"
    }
    file(params.keys_file).readLines().each { line ->
        if (!line.trim()) return
        def (key, width) = line.split('\t')
        println "BUCKET_RESULT:${key}\t${width}\t${hashBucket(key, width as int)}"
    }
}
