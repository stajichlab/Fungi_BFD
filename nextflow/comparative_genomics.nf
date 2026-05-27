#!/usr/bin/env nextflow

/*
 * comparative_genomics.nf — comparative genomics clustering workflows.
 *
 * Usage (from project root):
 *   sbatch nextflow/run_comparative.sh --project MyProject --group group.csv
 *   sbatch nextflow/run_comparative.sh --project MyProject --taxon CLASS:Dothideomycetes
 *   sbatch nextflow/run_comparative.sh --project MyProject \
 *       --taxon PHYLUM:Ascomycota,CLASS:Sordariomycetes --ignore exclude.txt
 *
 * Organism selection (OR/UNION, then --ignore exclusion):
 *   --group  PATH  CSV with LOCUSTAG,GROUP columns
 *   --taxon  STR   Comma-separated RANK:VALUE pairs (OR logic)
 *   --ignore PATH  One LOCUSTAG per line; excluded from the analysis
 *
 * Clustering sub-workflows (toggled by run_* params):
 *   CLUSTER_MMSEQS2     MMseqs2 all-vs-all → sequence clusters
 *   CLUSTER_MCL         DIAMOND all-vs-all + MCL graph clustering
 *   CLUSTER_ORTHOFINDER OrthoFinder orthogroup inference
 *
 * Output root: {outdir}/{project}/
 */

// ─────────────────────────────────────────────────────────────────────────────
// PREPARE_COMPARATIVE — filter taxa, write manifest, symlink input files
// ─────────────────────────────────────────────────────────────────────────────

process STAGE_FILES {
    label    'setup'
    tag      params.project

    publishDir "${params.outdir}/${params.project}", mode: 'copy', pattern: '*.manifest.csv'

    input:
    path manifest
    val  pep_dir
    val  cds_dir

    output:
    val  "${params.outdir}/${params.project}/proteins", emit: proteins_dir
    val  "${params.outdir}/${params.project}/cds",      emit: cds_dir
    path manifest,                                      emit: manifest

    script:
    def proteins_out = "${params.outdir}/${params.project}/proteins"
    def cds_out      = "${params.outdir}/${params.project}/cds"
    """
    mkdir -p "${proteins_out}" "${cds_out}"
    tail -n +2 ${manifest} | while IFS=',' read -r locustag group; do
        [ -z "\$locustag" ] && continue
        src="${pep_dir}/\${locustag}.faa"
        dst="${proteins_out}/\${locustag}.faa"
        [ -f "\$src" ] && [ ! -e "\$dst" ] && ln -sf "\$src" "\$dst"
        src="${cds_dir}/\${locustag}.cds.fa"
        dst="${cds_out}/\${locustag}.cds.fa"
        [ -f "\$src" ] && [ ! -e "\$dst" ] && ln -sf "\$src" "\$dst"
    done
    echo "Staged \$(tail -n +2 ${manifest} | wc -l) species into ${proteins_out}"
    """
}

workflow PREPARE_COMPARATIVE {
    main:
    // Load group CSV: LOCUSTAG → GROUP (skip header row)
    def group_map = [:]
    if (params.group) {
        def gf = file(params.group)
        def first = true
        gf.eachLine { line ->
            if (first) { first = false; return }
            def cols = line.trim().split(',', 2)
            if (cols.size() == 2 && cols[0]) group_map[cols[0].trim()] = cols[1].trim()
        }
    }

    // Load ignore list
    def ignored = [] as Set
    if (params.ignore) {
        file(params.ignore).eachLine { line ->
            def lt = line.trim()
            if (lt) ignored << lt
        }
    }

    // Parse comma-separated RANK:VALUE taxon filters (OR logic)
    def taxon_filters = []
    if (params.taxon) {
        params.taxon.split(',').each { spec ->
            def kv = spec.trim().split(':', 2)
            if (kv.size() == 2) taxon_filters << [rank: kv[0].trim().toUpperCase(), value: kv[1].trim()]
        }
    }

    // Filter samples.csv; emit [LOCUSTAG, GROUP] tuples
    def species_ch = Channel.fromPath(params.samples)
        .splitCsv(header: true, strip: true)
        .filter { row ->
            def lt = row['LOCUSTAG']
            if (!lt || ignored.contains(lt)) return false
            def in_group = group_map.containsKey(lt)
            def in_taxon = taxon_filters.any { f -> row[f.rank]?.equalsIgnoreCase(f.value) }
            return in_group || in_taxon
        }
        .map { row ->
            def lt  = row['LOCUSTAG']
            def grp = group_map.containsKey(lt)
                        ? group_map[lt]
                        : (taxon_filters.find { f -> row[f.rank]?.equalsIgnoreCase(f.value) }?.value ?: 'default')
            [lt, grp]
        }

    // Assemble manifest CSV from the filtered channel
    def manifest_ch = species_ch
        .collectFile(
            name:    "${params.project}.manifest.csv",
            seed:    "LOCUSTAG,GROUP\n",
            newLine: true
        ) { lt, grp -> "${lt},${grp}" }

    STAGE_FILES(manifest_ch, params.pep_dir, params.cds_dir)

    emit:
    proteins_dir = STAGE_FILES.out.proteins_dir
    cds_dir      = STAGE_FILES.out.cds_dir
    manifest     = STAGE_FILES.out.manifest
}

// ─────────────────────────────────────────────────────────────────────────────
// CLUSTER_MMSEQS2 — MMseqs2 all-vs-all search + sequence clustering
// ─────────────────────────────────────────────────────────────────────────────

process MMSEQS_CREATEDB {
    label    'mmseqs'
    tag      "createdb"
    storeDir "${params.outdir}/${params.project}/mmseqs2/db"

    input:
    val proteins_dir

    output:
    path "combined.faa",    emit: combined_faa
    path "combined_db*",    emit: db_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    cat "${proteins_dir}"/*.faa > combined.faa
    mmseqs createdb combined.faa combined_db
    """

    stub:
    """
    touch combined.faa combined_db combined_db.index combined_db.dbtype
    touch combined_db_h combined_db_h.index combined_db_h.dbtype
    """
}

process MMSEQS_SEARCH {
    label    'mmseqs'
    tag      "search"
    storeDir "${params.outdir}/${params.project}/mmseqs2/search"

    input:
    path db_files

    output:
    path "search_result*", emit: result_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mkdir -p tmp
    mmseqs search combined_db combined_db search_result tmp \\
        --min-seq-id ${params.mmseqs_min_id} \\
        -c ${params.mmseqs_cov} --cov-mode 0 \\
        -s ${params.mmseqs_sensitivity} \\
        --threads ${task.cpus}
    rm -rf tmp
    """

    stub:
    """
    touch search_result search_result.index search_result.dbtype
    """
}

process MMSEQS_CLUST {
    label 'mmseqs'
    tag   "clust"

    input:
    path db_files
    path result_files

    output:
    path "cluster_result*", emit: cluster_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mmseqs clust combined_db search_result cluster_result \\
        --cluster-mode ${params.mmseqs_cluster_mode} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch cluster_result cluster_result.index cluster_result.dbtype
    """
}

process MMSEQS_CREATETSV {
    label      'mmseqs'
    tag        "createtsv"
    publishDir "${params.outdir}/${params.project}/mmseqs2", mode: 'copy'

    input:
    path db_files
    path cluster_files

    output:
    path "${params.project}.mmseqs2_clusters.tsv", emit: tsv

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    mmseqs createtsv combined_db combined_db cluster_result \\
        ${params.project}.mmseqs2_clusters.tsv
    """

    stub:
    """
    touch ${params.project}.mmseqs2_clusters.tsv
    """
}

workflow CLUSTER_MMSEQS2 {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    MMSEQS_CREATEDB(proteins_dir)
    MMSEQS_SEARCH(MMSEQS_CREATEDB.out.db_files)
    MMSEQS_CLUST(MMSEQS_CREATEDB.out.db_files, MMSEQS_SEARCH.out.result_files)
    MMSEQS_CREATETSV(MMSEQS_CREATEDB.out.db_files, MMSEQS_CLUST.out.cluster_files)

    emit:
    tsv = MMSEQS_CREATETSV.out.tsv
}

// ─────────────────────────────────────────────────────────────────────────────
// CLUSTER_MCL — DIAMOND all-vs-all blastp + MCL graph clustering
// ─────────────────────────────────────────────────────────────────────────────

process DIAMOND_MAKEDB {
    label    'diamond'
    tag      "makedb"
    storeDir "${params.outdir}/${params.project}/mcl/db"

    input:
    val proteins_dir

    output:
    path "combined.faa",  emit: combined_faa
    path "combined.dmnd", emit: db

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond
    cat "${proteins_dir}"/*.faa > combined.faa
    diamond makedb --in combined.faa -d combined --threads ${task.cpus}
    """

    stub:
    """
    touch combined.faa combined.dmnd
    """
}

process DIAMOND_BLASTP {
    label    'diamond'
    tag      "blastp"
    storeDir "${params.outdir}/${params.project}/mcl/search"

    input:
    path combined_faa
    path db

    output:
    path "blastp.tsv", emit: blastp_tsv

    script:
    def sens_flag = params.diamond_more_sensitive ? '--more-sensitive' : '--sensitive'
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond
    diamond blastp \\
        -d combined \\
        -q ${combined_faa} \\
        ${sens_flag} \\
        --max-target-seqs 500 \\
        -e ${params.mcl_evalue} \\
        --outfmt 6 \\
        --threads ${task.cpus} \\
        -o blastp.tsv
    """

    stub:
    """
    printf 'A\tB\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\n' > blastp.tsv
    """
}

process MCL_PREPARE {
    label 'mcl'
    tag   "mcl_prepare"

    input:
    path blastp_tsv

    output:
    path "mcl_input.abc", emit: abc

    script:
    """
    awk -F'\\t' '\$1 != \$2 {print \$1, \$2, \$12}' ${blastp_tsv} > mcl_input.abc
    """

    stub:
    """
    printf 'A B 200\n' > mcl_input.abc
    """
}

process MCL_RUN {
    label      'mcl'
    tag        "mcl"
    publishDir "${params.outdir}/${params.project}/mcl", mode: 'copy'

    input:
    path abc

    output:
    path "${params.project}.mcl_clusters.tsv", emit: tsv

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mcl
    mcl ${abc} --abc \\
        -I ${params.mcl_inflation} \\
        -te ${task.cpus} \\
        -o ${params.project}.mcl_clusters.tsv
    """

    stub:
    """
    touch ${params.project}.mcl_clusters.tsv
    """
}

workflow CLUSTER_MCL {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    DIAMOND_MAKEDB(proteins_dir)
    DIAMOND_BLASTP(DIAMOND_MAKEDB.out.combined_faa, DIAMOND_MAKEDB.out.db)
    MCL_PREPARE(DIAMOND_BLASTP.out.blastp_tsv)
    MCL_RUN(MCL_PREPARE.out.abc)

    emit:
    tsv = MCL_RUN.out.tsv
}

// ─────────────────────────────────────────────────────────────────────────────
// CLUSTER_ORTHOFINDER — OrthoFinder orthogroup inference
// ─────────────────────────────────────────────────────────────────────────────

process ORTHOFINDER_RUN {
    label    'orthofinder'
    tag      "orthofinder"
    storeDir "${params.outdir}/${params.project}/orthofinder/run"

    input:
    val proteins_dir

    output:
    path "orthofinder_out", emit: out_dir

    script:
    def msa_flag = params.orthofinder_msa ? '-M msa' : ''
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load orthofinder
    orthofinder -f "${proteins_dir}" \\
        -t ${task.cpus} -a ${task.cpus} \\
        ${msa_flag} \\
        -o orthofinder_out
    """

    stub:
    """
    mkdir -p orthofinder_out/Orthogroups orthofinder_out/Comparative_Genomics_Statistics
    printf 'Orthogroup\n' > orthofinder_out/Orthogroups/Orthogroups.tsv
    printf 'Number of species\t1\n' > orthofinder_out/Comparative_Genomics_Statistics/Statistics_Overall.tsv
    """
}

process ORTHOFINDER_PARSE {
    label      'orthofinder'
    tag        "parse"
    publishDir "${params.outdir}/${params.project}/orthofinder", mode: 'copy'

    input:
    path orthofinder_out

    output:
    path "Orthogroups.tsv",                     emit: orthogroups
    path "Orthogroups_UnassignedGenes.tsv",      emit: unassigned, optional: true
    path "Statistics_Overall.tsv",               emit: stats,      optional: true

    script:
    // Handle both OrthoFinder v2 (Results_<date>/ subdir) and v3 (flat layout).
    """
    OF_DIR="${orthofinder_out}"
    if ls "\${OF_DIR}"/Results_* >/dev/null 2>&1; then
        OF_DIR=\$(ls -d "\${OF_DIR}"/Results_* | sort | tail -1)
    fi
    cp "\${OF_DIR}/Orthogroups/Orthogroups.tsv" .
    [ -f "\${OF_DIR}/Orthogroups/Orthogroups_UnassignedGenes.tsv" ] && \\
        cp "\${OF_DIR}/Orthogroups/Orthogroups_UnassignedGenes.tsv" . || true
    [ -f "\${OF_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" ] && \\
        cp "\${OF_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv" . || true
    """

    stub:
    """
    touch Orthogroups.tsv Statistics_Overall.tsv
    """
}

workflow CLUSTER_ORTHOFINDER {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    ORTHOFINDER_RUN(proteins_dir)
    ORTHOFINDER_PARSE(ORTHOFINDER_RUN.out.out_dir)

    emit:
    orthogroups = ORTHOFINDER_PARSE.out.orthogroups
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────

workflow {
    if (!params.project) {
        error "Required parameter --project not specified.\n" +
              "  Example: --project MyComparison"
    }
    if (!params.group && !params.taxon) {
        error "At least one of --group or --taxon must be specified to select organisms.\n" +
              "  Examples: --group group.csv\n" +
              "            --taxon CLASS:Dothideomycetes\n" +
              "            --taxon PHYLUM:Ascomycota,ORDER:Hypocreales"
    }

    PREPARE_COMPARATIVE()

    if (params.run_mmseqs2) {
        CLUSTER_MMSEQS2(PREPARE_COMPARATIVE.out.proteins_dir)
    }

    if (params.run_mcl) {
        CLUSTER_MCL(PREPARE_COMPARATIVE.out.proteins_dir)
    }

    if (params.run_orthofinder) {
        CLUSTER_ORTHOFINDER(PREPARE_COMPARATIVE.out.proteins_dir)
    }
}
