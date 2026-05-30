#!/usr/bin/env nextflow

/*
 * SOURCE: ../../../1KFG/common_annotate/pipeline/nextflow/interproscan6.nf
 * Last synced: 2026-05-23
 * Changes vs source: removed nextflow.enable.dsl=2; params block moved to
 *                    conf/profile_interproscan6.config.
 *
 * Usage (from project root):
 *   nextflow run nextflow/interproscan6.nf -c nextflow/nextflow.config -profile interproscan6 -resume
 *
 * interproscan6 requires Nextflow >= 25.04.6; run:
 *   nextflow pull ebi-pf-team/interproscan6 -r <ips6_version>
 * once on the head node so compute nodes can run with NXF_OFFLINE=true.
 */

process INTERPROSCAN6_RUN {
    tag "$out"

    cpus   16
    memory '64 GB'
    time   '72h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(proteins_fa)

    output:
    tuple val(out), path("${out}/annotate_misc/iprscan.xml")

    script:
    def appsArg = params.ips6_applications ? "--applications '${params.ips6_applications}'" : ""
    def goArg   = params.ips6_goterms  ? "--goterms"  : ""
    def pathArg = params.ips6_pathways ? "--pathways" : ""
    """
    if [ ! -f "${proteins_fa}" ]; then
        echo "ERROR: proteins FASTA not found: ${proteins_fa}" >&2
        exit 1
    fi

    module load singularity
    export NXF_SINGULARITY_CACHEDIR=${params.singularity_cache}
    export NXF_OFFLINE=true

    NF_WORK=\${SCRATCH:-/tmp}/nf_ips6_${out}
    mkdir -p \$NF_WORK ${out}/annotate_misc

    nextflow run ebi-pf-team/interproscan6 \\
        -r ${params.ips6_version} \\
        -profile singularity,local \\
        -w \$NF_WORK \\
        --input ${proteins_fa} \\
        --datadir ${params.ips6_datadir} \\
        --outprefix ${out}/annotate_misc/iprscan \\
        --formats xml \\
        --interpro ${params.ips6_interpro} \\
        --cpus ${task.cpus} \\
        ${goArg} ${pathArg} ${appsArg}

    rm -rf \$NF_WORK
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}

workflow {
    def suppressSet = (params.suppress && file(params.suppress).exists())
        ? file(params.suppress).readLines()
              .collect { it.trim().split(',')[0].trim() }
              .findAll { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }

    channel.fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
            def strain  = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
            strain = strain.replaceAll(/;.*$/, '').trim()
            def out     = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
            def asmid   = row.ASMID?.trim()
            tuple(out, asmid)
        }
        .filter { out, asmid -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .filter { out, asmid -> !suppressSet.contains(asmid) }
        .filter { out, _asmid ->
            def proteins = file("${params.target}/${out}/predict_results/${out}.proteins.fa")
            if (!proteins.exists()) {
                if (params.debug.toBoolean()) log.info "Skipping ${out}: no proteins FASTA found"
                return false
            }
            return true
        }
        .filter { out, _asmid ->
            def xml   = file("${params.target}/${out}/annotate_misc/iprscan.xml")
            def prot  = file("${params.target}/${out}/predict_results/${out}.proteins.fa")
            if (xml.exists() && xml.lastModified() >= prot.lastModified()) {
                if (params.debug.toBoolean()) log.info "Skipping ${out}: iprscan.xml exists and is current"
                return false
            }
            return true
        }
        .map { out, _asmid ->
            tuple(out, "${params.target}/${out}/predict_results/${out}.proteins.fa")
        }
        | INTERPROSCAN6_RUN
}
