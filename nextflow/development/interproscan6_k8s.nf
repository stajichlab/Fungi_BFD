#!/usr/bin/env nextflow

include { makeSampleTag } from '../modules/common/utils.nf'

/*
 * interproscan6_k8s.nf — Kubernetes-native IPS6 wrapper for Fungi_BFD.
 *
 * SOURCE: nextflow/interproscan6.nf (stajichlab/Fungi_BFD)
 * K8s adaptation: replaces nested `nextflow run ebi-pf-team/interproscan6`
 * (which relies on Singularity + module load) with a direct invocation of the
 * official interproscan Docker image (docker.io/interpro/interproscan:6.0.0).
 * The IPS6 data directory is expected on the shared PVC at params.ips6_datadir.
 *
 * Usage (from project root, with kubeconfig pointing at Nautilus):
 *   nextflow run nextflow/interproscan6_k8s.nf \
 *       -c nextflow/nextflow.config \
 *       -profile interproscan6_k8s \
 *       -resume
 *
 * Stub / dry-run test:
 *   nextflow run nextflow/interproscan6_k8s.nf \
 *       -c nextflow/nextflow.config \
 *       -profile interproscan6_k8s \
 *       -stub-run --n_test 2
 *
 * Prerequisites:
 *   1. A Nautilus PVC with ReadWriteMany access (see conf/profile_interproscan6_k8s.config).
 *   2. IPS6 data directory pre-staged to the PVC:
 *        kubectl cp /bigdata/stajichlab/shared/lib/interproscan6_data \
 *            <loader-pod>:/workspace/lib/interproscan6_data
 *   3. genome_annotation/ directory (with predict_results/<out>.proteins.fa files)
 *      pre-staged to the PVC at params.target.
 *   4. samples.csv copied to the PVC or passed via --samples.
 *
 * Key differences from the HPC version:
 *   - No `module load apptainer`; Docker is the container runtime.
 *   - interproscan.sh is called directly (it is the container ENTRYPOINT).
 *   - NXF_OFFLINE / nested nextflow run removed entirely.
 *   - SCRATCH replaced with a task-local tmpdir inside the pod workDir.
 *   - params.singularity_cache is unused on K8s (kept in config for HPC compat).
 */

// ─────────────────────────────────────────────────────────────────────────────
// INTERPROSCAN6_RUN
//   Runs interproscan.sh directly inside the official IPS6 Docker container.
//   One pod per fungal genome; cpus / memory / time come from the K8s profile.
// ─────────────────────────────────────────────────────────────────────────────
process INTERPROSCAN6_RUN {
    tag "$out"

    // Container is declared here as a fallback; the profile's process block
    // overrides it if you want to pin a different version.
    container 'docker.io/interpro/interproscan:6.0.0'

    // Resource defaults — overridden by withName in profile_interproscan6_k8s.config
    cpus   16
    memory '64 GB'
    time   '72h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(proteins_fa)

    output:
    tuple val(out), path("${out}/annotate_misc/iprscan.xml")

    script:
    // Build optional flag strings exactly as in the original
    def appsArg = params.ips6_applications ? "--applications '${params.ips6_applications}'" : ""
    def goArg   = params.ips6_goterms  ? "--goterms"  : ""
    def pathArg = params.ips6_pathways ? "--pathways" : ""
    def interpro = params.ips6_interpro ?: "latest"
    """
    set -euo pipefail

    # Validate input — fail fast before requesting the pod
    if [ ! -f "${proteins_fa}" ]; then
        echo "ERROR: proteins FASTA not found: ${proteins_fa}" >&2
        exit 1
    fi

    mkdir -p ${out}/annotate_misc

    # interproscan.sh is the container ENTRYPOINT; call it directly.
    # --disable-precalc avoids attempting the lookup service from inside a pod.
    interproscan.sh \
        --input          ${proteins_fa} \
        --output-dir     ${out}/annotate_misc \
        --outfilebase    iprscan \
        --formats        xml \
        --iprlookup      \
        --interproscan   ${interpro} \
        --cpu            ${task.cpus} \
        --tempdir        ./ips6_tmp_${out} \
        --datadir        ${params.ips6_datadir} \
        ${goArg} ${pathArg} ${appsArg}

    # Clean pod-local temp
    rm -rf ./ips6_tmp_${out}
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// WORKFLOW — identical sample selection logic as the original
// ─────────────────────────────────────────────────────────────────────────────
workflow {
    // Build suppress set from file (same as original)
    def suppressSet = (params.suppress && file(params.suppress).exists())
        ? file(params.suppress).readLines()
              .collect  { it.trim().split(',')[0].trim() }
              .findAll  { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)

    if (suppressSet) {
        log.info "Suppress list: ${suppressSet.size()} ASMIDs will be skipped"
    }

    channel.fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species = (row.SPECIES?.trim() ?: '').replaceAll(/['"]/, '')
            def strain  = (row.STRAIN?.trim() ?: '').replaceAll(/['"]/, '').replaceAll(/;.*$/, '').trim().replace(':', ' ')
            def out     = makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            def asmid  = row.ASMID?.trim()
            tuple(out, asmid)
        }
        .filter { out, asmid -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .filter { out, asmid -> !suppressSet.contains(asmid) }
        // Skip if proteins FASTA is absent
        .filter { out, _asmid ->
            def proteins = file("${params.target}/${out}/predict_results/${out}.proteins.fa")
            if (!proteins.exists()) {
                if (params.debug.toBoolean()) log.info "Skipping ${out}: no proteins FASTA"
                return false
            }
            return true
        }
        // Skip if iprscan.xml already exists and is newer than the proteins FASTA
        .filter { out, _asmid ->
            def xml  = file("${params.target}/${out}/annotate_misc/iprscan.xml")
            def prot = file("${params.target}/${out}/predict_results/${out}.proteins.fa")
            if (xml.exists() && xml.lastModified() >= prot.lastModified()) {
                if (params.debug.toBoolean()) log.info "Skipping ${out}: iprscan.xml current"
                return false
            }
            return true
        }
        .map { out, _asmid ->
            // Resolve the absolute path; on K8s this must be within the PVC mount
            tuple(out, "${params.target}/${out}/predict_results/${out}.proteins.fa")
        }
        | INTERPROSCAN6_RUN
}
