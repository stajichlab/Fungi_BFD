#!/usr/bin/env nextflow
//
// Standalone real (non-preview, non-stub) smoke test of trainingSymlinkPreflight()
// (modules/funannotate/utils.nf) and the 3-way blocked/todo/done .branch{} pattern
// wired into FUNANNOTATE_RNASEQ.nf's train_branched, 2026-09-10. -preview does NOT
// execute operators (confirmed earlier in this repo's history -- see .living/learnings.md
// 2026-06-25 "combine operator" entry, ../Fungi_BFD), so a real run is required to
// validate the .branch{} clause syntax and the java.nio.file calls actually execute
// without exceptions. Lives at nextflow/ top level so workflow.projectDir resolves the
// same way the real entry point does (see genemark_run_smoke.nf's header comment for
// why this matters -- not relevant to this particular test, but kept for consistency).
//
// Usage:
//   nextflow run nextflow/training_preflight_smoke.nf --target <scratch dir>
//
// --target MUST be passed on the command line, not defaulted in-script: a bare
// top-of-script `params.target = ...` is shadowed by nextflow.config's own
// (production) default for this directory, and reassigning params.target from
// *inside* the workflow{} body is silently a no-op -- confirmed empirically
// (2026-09-10) that trainingSymlinkPreflight() kept reading the real production
// genome_annotation path in both cases, not the intended fixture root. Only a
// genuine CLI-level `--target` override actually takes effect.
nextflow.enable.dsl = 2

include { trainingSymlinkPreflight } from './modules/funannotate/utils.nf'

workflow {
    if (!params.target) {
        error "training_preflight_smoke.nf requires --target <scratch dir> on the command line"
    }
    def fixtureRoot = file(params.target as String)
    fixtureRoot.deleteDir()
    fixtureRoot.mkdirs()

    def canonicalDir = file("${fixtureRoot}/_canonical/REAL_OUT/training" as String)
    canonicalDir.mkdirs()
    file("${canonicalDir}/marker.txt" as String).text = "canonical training data\n"

    // Case 1: absent entirely -> 'ok'
    def caseAbsent = "ABSENT_OUT"

    // Case 2: valid symlink resolving to real data -> 'ok'
    def caseValid = "VALID_OUT"
    file("${fixtureRoot}/${caseValid}" as String).mkdirs()
    java.nio.file.Files.createSymbolicLink(
        java.nio.file.Paths.get("${fixtureRoot}/${caseValid}/training" as String),
        java.nio.file.Paths.get(canonicalDir.toString())
    )

    // Case 3: broken symlink -> 'broken_symlink', auto-removed
    def caseBroken = "BROKEN_OUT"
    file("${fixtureRoot}/${caseBroken}" as String).mkdirs()
    java.nio.file.Files.createSymbolicLink(
        java.nio.file.Paths.get("${fixtureRoot}/${caseBroken}/training" as String),
        java.nio.file.Paths.get("${fixtureRoot}/_nonexistent_target" as String)
    )

    // Case 4: real (non-symlink) directory -> 'real_directory', left untouched
    def caseRealDir = "REALDIR_OUT"
    file("${fixtureRoot}/${caseRealDir}/training" as String).mkdirs()
    file("${fixtureRoot}/${caseRealDir}/training/stray_data.txt" as String).text = "leftover real data\n"

    def cases = [caseAbsent, caseValid, caseBroken, caseRealDir]
    def results = cases.collectEntries { out -> [(out): trainingSymlinkPreflight(out)] }

    log.info "trainingSymlinkPreflight() results: ${results}"

    def expected = [
        (caseAbsent) : 'ok',
        (caseValid)  : 'ok',
        (caseBroken) : 'broken_symlink',
        (caseRealDir): 'real_directory',
    ]
    def mismatches = expected.findAll { k, v -> results[k] != v }
    if (mismatches) {
        error "trainingSymlinkPreflight() MISMATCH: expected ${expected}, got ${results}"
    }

    // Post-conditions: broken symlink actually removed; real dir left untouched.
    def brokenLinkPath = java.nio.file.Paths.get("${fixtureRoot}/${caseBroken}/training" as String)
    if (java.nio.file.Files.exists(brokenLinkPath, java.nio.file.LinkOption.NOFOLLOW_LINKS)) {
        error "broken symlink was NOT removed as expected: ${brokenLinkPath}"
    }
    def realDirPath = file("${fixtureRoot}/${caseRealDir}/training" as String)
    if (!realDirPath.isDirectory() || realDirPath.isLink()) {
        error "real directory case was unexpectedly modified: ${realDirPath}"
    }

    // Now exercise the exact 3-way blocked/todo/done .branch{} + collectFile pattern
    // used in FUNANNOTATE_RNASEQ.nf's train_branched, against these same 4 synthetic
    // rows (todo/done classification here is arbitrary/synthetic -- only the branch
    // mechanics + collectFile output are under test).
    def rows = Channel.fromList(cases)
    def branched = rows.branch { out ->
        def status = trainingSymlinkPreflight(out as String)
        blocked: status == 'real_directory'
        todo: out == caseBroken
        done: true
    }
    branched.blocked
        .map { out -> "${out}\treal_directory_at_training_symlink_path" }
        .collectFile(name: 'smoke_train_blocked.tsv', storeDir: params.target,
                     newLine: true, sort: true, seed: "out\treason")
    branched.todo.view { out -> "[todo] ${out}" }
    branched.done.view { out -> "[done] ${out}" }

    log.info "training_preflight_smoke: ALL CHECKS PASSED"
}
