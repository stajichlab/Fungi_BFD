#!/usr/bin/env nextflow

/*
 * compare_ANI — Average Nucleotide Identity comparison pipeline
 *
 * Groups genomes from samples.csv by a taxonomic rank (--compare GENUS by default),
 * runs fastANI all-vs-all within each group, then generates a clustering + outlier
 * report per group.
 *
 * Usage (from project root):
 *   nextflow run nextflow/compare_ANI.nf \
 *       -c nextflow/nextflow.config -profile ani -resume
 *
 *   nextflow run nextflow/compare_ANI.nf \
 *       -c nextflow/nextflow.config -profile ani \
 *       --taxon PHYLUM:Ascomycota --compare FAMILY -resume
 *
 *   nextflow run nextflow/compare_ANI.nf \
 *       -c nextflow/nextflow.config -profile ani \
 *       -params-file nextflow/params_ani.yaml -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/compare_ANI.nf \
 *       -c nextflow/nextflow.config -profile ani -stub-run --n_test 3
 */

// ════════════════════════════════════════════════════════════════════════════
// PARAMETER DEFAULTS
// (override via --param value on CLI or in a params-file)
// ════════════════════════════════════════════════════════════════════════════

params.genome_dir             = "${launchDir}/input/dna"
params.genome_suffix          = '.scaffolds.fa'   // alt: '.masked.fa' with input_clean_genomes/
// Controls how the genome filename stem is derived from samples.csv:
//   'species' (default) → ${SPECIES}_${STRAIN}${genome_suffix}  e.g. Fusarium_oxysporum_CBS123.scaffolds.fa
//   'asmid'             → ${ASMID}${genome_suffix}               e.g. GCF_010015735.1_Aaoar1.masked.fa
params.genome_name_style      = 'species'
params.compare                = 'GENUS'
params.ani_cluster_threshold  = 95.0
params.ani_outlier_threshold  = 90.0
params.min_group_size         = 2
params.ani_batch_size         = 0       // 0 = single job per group (recommended)
params.outdir                 = "${launchDir}/results/ANI"
params.fastani_fraglen        = 3000
params.fastani_kmer           = 16

// ════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ════════════════════════════════════════════════════════════════════════════

// ── ANI_COMPARE ──────────────────────────────────────────────────────────────
// Runs fastANI with a query list vs a reference list.
// For un-batched groups: query == ref == all genomes in the group.
// For batched groups:    query and ref are subset lists (batch i × batch j).
// storeDir means: if the output TSV already exists on disk, skip this job.
//
// group_size drives CPU/memory scaling (total genomes in the group, not batch size):
//   <= 200 genomes →  8 CPUs / 16 GB
//   <= 500 genomes → 24 CPUs / 48 GB
//    > 500 genomes → 64 CPUs / 128 GB

process ANI_COMPARE {
    tag   "${group_name} [${batch_tag}] n=${group_size}"
    label 'fastani'

    cpus   { group_size > 500 ? 64 : group_size > 200 ? 24 : 8 }
    memory { group_size > 500 ? '128 GB' : group_size > 200 ? '48 GB' : '16 GB' }

    storeDir "${params.outdir}/${params.compare}/${group_name}/batches"

    input:
        tuple val(group_name), path(query_genomes), path(ref_genomes), val(batch_tag), val(group_size)

    output:
        tuple val(group_name), path("${group_name}.${batch_tag}.ani.tsv")

    script:
    """
    ls ${query_genomes} > query_list.txt
    ls ${ref_genomes}   > ref_list.txt
    fastANI \\
        --ql query_list.txt \\
        --rl ref_list.txt \\
        -o ${group_name}.${batch_tag}.ani.tsv \\
        --fragLen ${params.fastani_fraglen} \\
        -k ${params.fastani_kmer} \\
        -t ${task.cpus}
    """

    stub:
    """
    printf 'query\\treference\\tANI\\tmapped\\ttotal\\n' > ${group_name}.${batch_tag}.ani.tsv
    """
}

// ── MERGE_ANI_BATCHES ────────────────────────────────────────────────────────
// Merges per-batch TSVs into a single file for the group.
// Only used when ani_batch_size > 0 and the group was split into batches.

process MERGE_ANI_BATCHES {
    tag   "${group_name}"
    label 'report'

    storeDir "${params.outdir}/${params.compare}/${group_name}"

    input:
        tuple val(group_name), path(batch_tsvs)

    output:
        tuple val(group_name), path("${group_name}.ani.tsv")

    script:
    """
    cat ${batch_tsvs} > ${group_name}.ani.tsv
    """

    stub:
    """
    printf 'query\\treference\\tANI\\tmapped\\ttotal\\n' > ${group_name}.ani.tsv
    """
}

// ── REPORT_ANI ───────────────────────────────────────────────────────────────
// Parses the merged ANI TSV and writes a human-readable report.
// report_ani.py is in bin/ which Nextflow adds to PATH automatically.

process REPORT_ANI {
    tag   "${group_name}"
    label 'report'

    publishDir { "${params.outdir}/${params.compare}/${group_name}" }, mode: 'copy'

    input:
        tuple val(group_name), path(ani_tsv), path(names_tsv)

    output:
        path("${group_name}_ANI_report.txt")

    script:
    """
    report_ani.py \\
        --input    ${ani_tsv} \\
        --names    ${names_tsv} \\
        --group    "${group_name}" \\
        --level    "${params.compare}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --outlier-threshold ${params.ani_outlier_threshold} \\
        --output   ${group_name}_ANI_report.txt
    """

    stub:
    """
    printf '=== ANI Report: ${group_name} (stub) ===\\n' > ${group_name}_ANI_report.txt
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Validate params ───────────────────────────────────────────────────────
    def validRanks = ['PHYLUM','SUBPHYLUM','CLASS','SUBCLASS','ORDER','FAMILY','GENUS']
    def compareRank = (params.compare as String).toUpperCase()
    if (!(compareRank in validRanks)) {
        error "--compare must be one of: ${validRanks.join(', ')}"
    }

    def nameStyle = (params.genome_name_style as String).toLowerCase()
    if (!(nameStyle in ['species', 'asmid'])) {
        error "--genome_name_style must be 'species' or 'asmid'"
    }
    log.info "Genome filename style: ${nameStyle} (suffix: ${params.genome_suffix})"

    // ── Taxonomy filter ───────────────────────────────────────────────────────
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be RANK:VALUE, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { _row -> true }
    }

    // ── Build sample channel ──────────────────────────────────────────────────
    // Emits: tuple(group_key, locustag, species, genome_path)
    def sample_ch = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def locustag  = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def species   = row.SPECIES?.trim() ?: ''
            def strain    = (row.STRAIN?.trim() ?: '').split(';')[0].trim()
                                .replace("'", '').replace(':', ' ')
            def groupKey  = row[compareRank]?.trim() ?: ''
            def stem      = nameStyle == 'asmid'
                                ? row.ASMID?.trim()
                                : [species, strain].findAll { s -> s }.join('_')
                                      .replaceAll(/[\s\/\#]+/, '_')
            def genome    = file("${params.genome_dir}/${stem}${params.genome_suffix}", glob: false)

            if (!groupKey) {
                log.warn "Skipping ${locustag}: empty ${compareRank} field"
                return null
            }
            if (!genome.exists()) {
                log.warn "Skipping ${locustag}: genome not found at ${genome}"
                return null
            }
            tuple(groupKey, locustag, species, genome)
        }
        .filter { item -> item != null }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // ── Group by taxonomic rank ───────────────────────────────────────────────
    // After groupTuple: tuple(group_name, [[locustag, species, genome], ...])
    def grouped_ch = sample_ch
        .map { groupKey, locustag, species, genome ->
            tuple(groupKey, tuple(locustag, species, genome))
        }
        .groupTuple()
        .filter { _gname, members -> members.size() >= params.min_group_size as int }

    // ── Prepare per-group genome list and names TSV ───────────────────────────
    // names_tsv: tab-separated filename\tspecies used by REPORT_ANI
    def prepared_ch = grouped_ch
        .map { group_name, members ->
            def genomes   = members.collect { m -> m[2] }
            def namesText = members.collect { m ->
                "${m[2].getName()}\t${m[1]}"
            }.join('\n')

            def nameFile = file("${workflow.workDir}/names_${group_name}.tsv")
            nameFile.text = "filename\tspecies\n" + namesText + "\n"

            tuple(group_name, genomes, nameFile)
        }

    // ── Route: batched vs single-job ─────────────────────────────────────────
    def batchSize = params.ani_batch_size as int

    def ani_out_ch

    if (batchSize > 0) {
        // Split into upper-triangle batch pairs; carry total group size for CPU scaling
        def batch_pairs_ch = prepared_ch
            .flatMap { group_name, genomes, _nameFile ->
                def N       = genomes.size()
                def nBatch  = (N + batchSize - 1).intdiv(batchSize)
                def batches = (0..<nBatch).collect { i ->
                    def end = (i + 1) * batchSize < N ? (i + 1) * batchSize : N
                    genomes.subList(i * batchSize, end)
                }
                // Upper triangle pairs (i <= j) — no for loops: use collectMany
                (0..<nBatch).collectMany { i ->
                    (i..<nBatch).collect { j ->
                        tuple(group_name, batches[i], batches[j], "b${i}_b${j}", N)
                    }
                }
            }

        def raw_ch     = ANI_COMPARE(batch_pairs_ch)
        def names_map  = prepared_ch.map { gn, _genomes, nf -> tuple(gn, nf) }

        ani_out_ch = raw_ch
            .groupTuple()
            .map { gn, tsv_list -> tuple(gn, tsv_list.flatten()) }
            | MERGE_ANI_BATCHES
            | combine(names_map, by: 0)

    } else {
        // Single job per group; pass total group size for CPU scaling
        def single_ch  = prepared_ch
            .map { group_name, genomes, _nameFile ->
                tuple(group_name, genomes, genomes, "full", genomes.size())
            }

        def raw_ch    = ANI_COMPARE(single_ch)
        def names_map = prepared_ch.map { gn, _genomes, nf -> tuple(gn, nf) }

        ani_out_ch = raw_ch
            | combine(names_map, by: 0)
    }

    // ── Generate reports ──────────────────────────────────────────────────────
    REPORT_ANI(ani_out_ch)
}
