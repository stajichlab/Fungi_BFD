#!/usr/bin/env nextflow

/*
 * PHYling Phylogenomics Pipeline
 * Multi-locus phylogeny from protein or CDS input using PHYling align/filter/tree.
 *
 * Usage (from project root):
 *   sbatch nextflow/run_phyling.sh
 *   nextflow run nextflow/phyling.nf -c nextflow/nextflow.config -profile phyling -resume
 *
 * Restrict to a taxonomic group (output goes to phylogeny/{VALUE}/{seq_type}/{markerset}/):
 *   sbatch nextflow/run_phyling.sh --taxon PHYLUM:Ascomycota
 *   sbatch nextflow/run_phyling.sh --taxon CLASS:Sordariomycetes --seq_type cds
 *
 * Multiple markersets run as independent branches in parallel:
 *   sbatch nextflow/run_phyling.sh --markerset "fungi_odb12,ascomycota_odb12"
 *
 * All params are defined in conf/profile_phyling.config; override any with --param value.
 */

// ════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ════════════════════════════════════════════════════════════════════════════

// Resolve each markerset: use the pre-extracted shared directory when available,
// otherwise extract the versioned tarball into the project's local markerset cache.
// storeDir means this runs at most once per markerset across all pipeline runs.
process MARKERSET_PREPARE {
    tag   "${markerset}"
    label 'setup'
    storeDir "${params.phylo_outdir}/markersets"

    input:
        val markerset

    output:
        tuple val(markerset), path("${markerset}"), emit: ready

    script:
    """
    src="${params.markerset_db}/${markerset}"
    if [ -d "\$src" ]; then
        ln -sfn "\$src" "${markerset}"
        echo "[INFO] Using pre-extracted markerset: \$src"
    else
        TARBALL=\$(ls "${params.markerset_db}/${markerset}".*.tar.gz 2>/dev/null | head -1)
        if [ -z "\$TARBALL" ]; then
            echo "ERROR: markerset '${markerset}' not found in ${params.markerset_db}" >&2
            exit 1
        fi
        echo "[INFO] Extracting \$TARBALL ..."
        mkdir -p "${markerset}"
        tar -xzf "\$TARBALL" -C "${markerset}" --strip-components=1
    fi
    """

    stub:
    """
    mkdir -p "${markerset}/hmms"
    touch "${markerset}/hmms/dummy.hmm"
    touch "${markerset}/links_to_ODB10.txt"
    """
}

// Stage input FASTAs with clean taxon names (strip .proteins / .cds-transcripts suffix),
// then run phyling align against the markerset HMMs.
process PHYLING_ALIGN {
    tag   "${markerset}"
    label 'phyling_align'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(markerset_dir), path(fastas)

    output:
        tuple val(markerset), val(taxon_slug), path("align"), emit: align

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    # Create clean-named symlinks so tree taxon labels = {Species_Strain}
    mkdir -p staged
    for f in *.fa *.fa.gz; do
        [ -e "\$f" ] || continue
        base="\$(basename "\$f" .gz)"
        base="\$(basename "\$base" .fa)"
        base="\${base%.proteins}"
        base="\${base%.cds-transcripts}"
        ext="\$(echo "\$f" | grep -oE '\\.fa(\\.gz)?$')"
        ln -sfn "\$(readlink -f "\$f")" "staged/\${base}\${ext}"
    done

    module load phyling
    phyling align \\
        -I staged \\
        -m ${markerset_dir} \\
        -o align \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p align
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > align/BUSCOmarker001.fa
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > align/BUSCOmarker002.fa
    """
}

// Select top-N markers by treeness/RCV score using phyling filter.
process PHYLING_FILTER {
    tag   "${markerset}"
    label 'phyling_filter'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(align_dir)

    output:
        tuple val(markerset), val(taxon_slug), path("filter"), emit: filtered

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    module load phyling
    phyling filter \\
        -I ${align_dir} \\
        -n ${params.top_n} \\
        -o filter \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p filter
    printf '>Taxon_A\\nACGTACGT\\n>Taxon_B\\nACGTACGT\\n' > filter/BUSCOmarker001.fa
    """
}

// Build a concatenated partitioned ML tree from the filtered marker MSAs.
// Uses phyling tree --concat --partition with the chosen method (ft/iqtree/raxml).
process PHYLING_TREE {
    tag   "${markerset}"
    label 'phyling_tree'
    publishDir { "${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/${markerset}" }, mode: 'copy'

    input:
        tuple val(markerset), val(taxon_slug), path(filter_dir)

    output:
        path("tree"), emit: tree

    script:
    def seqtype = params.seq_type == 'cds' ? 'dna' : 'pep'
    """
    module load phyling
    phyling tree \\
        -I ${filter_dir} \\
        -M ${params.tree_method} \\
        --concat \\
        --partition \\
        -o tree \\
        --seqtype ${seqtype} \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p tree
    echo '(Taxon_A:0.1,Taxon_B:0.1,(Taxon_C:0.05,Taxon_D:0.05):0.1);' > tree/concat.treefile
    touch tree/concat.partition
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {
    // ── Taxonomy filter ────────────────────────────────────────────────────────
    // Parse --taxon RANK:VALUE (e.g. --taxon PHYLUM:Ascomycota).
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { row -> true }
    }

    // ── Output slug ────────────────────────────────────────────────────────────
    // --taxon PHYLUM:Ascomycota → "Ascomycota";  no filter → "all"
    def taxon_slug = params.taxon
        ? (params.taxon as String).split(':', 2)[1].replaceAll(/\s+/, '_')
        : 'all'

    log.info "Output directory: ${params.phylo_outdir}/${taxon_slug}/${params.seq_type}/"

    // ── Input FASTA channel ────────────────────────────────────────────────────
    // Reads phylo.csv (default), applies taxon filter, collects all matching FASTAs.
    // phyling align needs the full set at once, so we collect() before passing it in.
    def input_dir = params.seq_type == 'cds' ? params.cds_dir : params.pep_dir
    def suffix    = params.seq_type == 'cds' ? '.cds-transcripts.fa' : '.proteins.fa'

    def fasta_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = (row.STRAIN?.trim() ?: '').split(';')[0].trim().replace("'", '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = [species, strain].findAll { it }.join('_').replaceAll(/[\s\/\#]+/, '_')
            def fasta    = file("${input_dir}/${basename}${suffix}", glob: false)
            if (!fasta.exists()) {
                log.warn "PHYling: skipping ${basename} (${locustag}): ${basename}${suffix} not found"
                return null
            }
            return fasta
        }
        .filter { it != null }
        .take(params.n_test > 0 ? params.n_test as int : -1)
        .collect()   // phyling align operates on a whole directory; collect all files first

    // ── Markerset channel ──────────────────────────────────────────────────────
    def markerset_ch = Channel.fromList(
        params.markerset.tokenize(',').collect { it.trim() }
    )

    // ── Pipeline ───────────────────────────────────────────────────────────────
    // Each markerset is prepared independently (storeDir caches extractions).
    MARKERSET_PREPARE(markerset_ch)

    // Broadcast the collected FASTA set to every markerset branch.
    // combine() with a single-element channel (from collect) creates one tuple per markerset.
    def align_input = MARKERSET_PREPARE.out.ready
        .map { markerset, markerset_dir -> tuple(markerset, taxon_slug, markerset_dir) }
        .combine(fasta_ch)

    PHYLING_ALIGN(align_input)
    PHYLING_FILTER(PHYLING_ALIGN.out.align)
    PHYLING_TREE(PHYLING_FILTER.out.filtered)
}
