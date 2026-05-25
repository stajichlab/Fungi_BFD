#!/usr/bin/env nextflow

/*
 * Sequence Statistics Pipeline
 * Computes per-species amino acid frequencies, codon frequencies,
 * gene structure tables, and chromosome info for all samples in samples.csv,
 * then consolidates each into a DuckDB-loadable <table>.csv.gz in tables/.
 *
 * Mirrors the logic of 08_sequence_stats.sh but runs each species in parallel
 * as a separate SLURM job.
 *
 * Usage (from project root):
 *   nextflow run nextflow/genome_seqstats.nf -c nextflow/nextflow.config -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/genome_seqstats.nf -c nextflow/nextflow.config -stub-run --n_test 2
 */

// ════════════════════════════════════════════════════════════════════════════
// SUBWORKFLOWS
// ════════════════════════════════════════════════════════════════════════════

// ── Amino acid frequency ──────────────────────────────────────────────────

workflow AA_FREQ {
    take: ch
    main:
        RUN_AA_FREQ(ch)
        MERGE_AA_FREQ(RUN_AA_FREQ.out.csv.collect())
    emit:
        merged = MERGE_AA_FREQ.out.csv
}

process RUN_AA_FREQ {
    tag   "${locustag}"
    label 'seqstats'
    publishDir "${params.stats_outdir}/aa_freq", mode: 'copy'

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.aa_freq.csv"), emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_AA_freq.py \\
        ${proteins} \\
        -o ${locustag}.aa_freq.csv
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' > ${locustag}.aa_freq.csv
    """
}

process MERGE_AA_FREQ {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(csvs)

    output:
        path("aa_freq.csv.gz"), emit: csv

    script:
    """
    first=\$(ls *.aa_freq.csv | head -1)
    head -1 "\$first" > aa_freq.csv
    for f in *.aa_freq.csv; do tail -n +2 "\$f" >> aa_freq.csv; done
    pigz aa_freq.csv
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > aa_freq.csv.gz
    """
}

// ── Codon frequency ───────────────────────────────────────────────────────

workflow CODON_FREQ {
    take: ch
    main:
        RUN_CODON_FREQ(ch)
        MERGE_CODON_FREQ(RUN_CODON_FREQ.out.csv.collect())
    emit:
        merged = MERGE_CODON_FREQ.out.csv
}

process RUN_CODON_FREQ {
    tag   "${locustag}"
    label 'seqstats'
    publishDir "${params.stats_outdir}/codon_freq", mode: 'copy'

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(cds)

    output:
        path("${locustag}.codon_freq.csv"), emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_codon_freq.py \\
        ${cds} \\
        -o ${locustag}.codon_freq.csv
    """

    stub:
    """
    printf 'species_prefix,codon,frequency\\n' > ${locustag}.codon_freq.csv
    """
}

process MERGE_CODON_FREQ {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(csvs)

    output:
        path("codon_freq.csv.gz"), emit: csv

    script:
    """
    first=\$(ls *.codon_freq.csv | head -1)
    head -1 "\$first" > codon_freq.csv
    for f in *.codon_freq.csv; do tail -n +2 "\$f" >> codon_freq.csv; done
    pigz codon_freq.csv
    """

    stub:
    """
    printf 'species_prefix,codon,frequency\\n' | gzip > codon_freq.csv.gz
    """
}

// ── Gene structure tables (build_genestats_bigquery.py) ──────────────────
// Each GFF3 + matching genome FASTA produces 7 per-species CSVs.
// build_genestats_bigquery.py is called with a single GFF3 and a local dna/
// directory containing a symlink named <basename>.scaffolds.fa → the staged
// genome file, matching the filename convention the script expects.

workflow GENE_STATS {
    take: ch   // tuple(locustag, basename, species, strain, gff3, genome)
    main:
        RUN_GENE_STATS(ch)
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_info.collect(),        'gene_info')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_transcripts.collect(), 'gene_transcripts')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_exons.collect(),       'gene_exons')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_CDS.collect(),         'gene_CDS')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_introns.collect(),     'gene_introns')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_trnas.collect(),       'gene_trnas')
        MERGE_GENE_STATS(RUN_GENE_STATS.out.gene_proteins.collect(),    'gene_proteins')
}

process RUN_GENE_STATS {
    tag   "${locustag}"
    label 'seqstats'
    publishDir "${params.stats_outdir}/gene_stats/${locustag}", mode: 'copy'

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(gff3), path(genome)

    output:
        path("${locustag}.gene_info.csv"),        emit: gene_info
        path("${locustag}.gene_transcripts.csv"), emit: gene_transcripts
        path("${locustag}.gene_exons.csv"),       emit: gene_exons
        path("${locustag}.gene_CDS.csv"),         emit: gene_CDS
        path("${locustag}.gene_introns.csv"),     emit: gene_introns
        path("${locustag}.gene_trnas.csv"),       emit: gene_trnas
        path("${locustag}.gene_proteins.csv"),    emit: gene_proteins

    script:
    """
    module load biopython
    mkdir -p dna tables
    ln -s "\$(readlink -f ${genome})" "dna/${basename}.scaffolds.fa"
    python3 ${params.scripts}/build_genestats_bigquery.py \\
        ${gff3} \\
        -d dna \\
        --outdir tables
    for tbl in gene_info gene_transcripts gene_exons gene_CDS gene_introns gene_trnas gene_proteins; do
        mv "tables/\${tbl}.csv" "${locustag}.\${tbl}.csv"
    done
    """

    stub:
    """
    for tbl in gene_info gene_transcripts gene_exons gene_CDS gene_introns gene_trnas gene_proteins; do
        printf 'gene_id\\n' > ${locustag}.\${tbl}.csv
    done
    """
}

process MERGE_GENE_STATS {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(csvs)
        val(table_name)

    output:
        path("${table_name}.csv.gz"), emit: csv

    script:
    """
    first=\$(ls *${table_name}.csv | head -1)
    head -1 "\$first" > ${table_name}.csv
    for f in *${table_name}.csv; do tail -n +2 "\$f" >> ${table_name}.csv; done
    pigz ${table_name}.csv
    """

    stub:
    """
    printf 'gene_id\\n' | gzip > ${table_name}.csv.gz
    """
}

// ── Chromosome / scaffold stats ───────────────────────────────────────────
// collect_chrom_info.py uses --run_with <line_index> to select a single row
// from samples.csv.  We fan out one job per row and merge the results.

workflow CHROM_INFO {
    take: ch   // tuple(line_idx, locustag)
    main:
        RUN_CHROM_INFO(ch)
        MERGE_CHROM_INFO(RUN_CHROM_INFO.out.csv.collect())
    emit:
        merged = MERGE_CHROM_INFO.out.csv
}

process RUN_CHROM_INFO {
    tag   "${locustag}"
    label 'seqstats'
    publishDir "${params.stats_outdir}/chrom_info", mode: 'copy'

    input:
        tuple val(line_idx), val(locustag)

    output:
        path("${locustag}.chrom_info.csv"), emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/collect_chrom_info.py \\
        --samples   ${params.samples} \\
        -d          ${params.genome_dir} \\
        --run_with  ${line_idx} \\
        --run_set   1 \\
        -o          ${locustag}.chrom_info.csv
    """

    stub:
    """
    printf 'LOCUSTAG,chrom_name,length,GC_percent,GC_count,left_50,right_50,lower_masked,N_masked\\n' > ${locustag}.chrom_info.csv
    """
}

process MERGE_CHROM_INFO {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(csvs)

    output:
        path("chrom_info.csv.gz"), emit: csv

    script:
    """
    first=\$(ls *.chrom_info.csv | head -1)
    head -1 "\$first" > chrom_info.csv
    for f in *.chrom_info.csv; do tail -n +2 "\$f" >> chrom_info.csv; done
    pigz chrom_info.csv
    """

    stub:
    """
    printf 'LOCUSTAG,chrom_name,length,GC_percent,GC_count,left_50,right_50,lower_masked,N_masked\\n' | gzip > chrom_info.csv.gz
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {
    // ── shared sample channel ──────────────────────────────────────────────
    // Builds tuple(locustag, basename, species, strain, <file>) per sample,
    // skipping rows whose input file is absent.

    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = (row.STRAIN?.trim() ?: '').split(';')[0].trim().replace("'", '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = [species, strain].findAll { it }.join('_').replaceAll(/[\s\/\#]+/, '_')
            [locustag, basename, species, strain]
        }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // protein FASTA channel for AA freq
    if (params.run_aa_freq) {
        def prot_ch = rows_ch.map { locustag, basename, species, strain ->
            def prot = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
            if (!prot.exists()) {
                log.warn "AA_FREQ: skipping ${basename} (${locustag}): protein file not found"
                return null
            }
            tuple(locustag, basename, species, strain, prot)
        }.filter { it != null }
        AA_FREQ(prot_ch)
    }

    // CDS FASTA channel for codon freq
    if (params.run_codon_freq) {
        def cds_ch = rows_ch.map { locustag, basename, species, strain ->
            def cds = file("${params.cds_dir}/${basename}.cds-transcripts.fa", glob: false)
            if (!cds.exists()) {
                log.warn "CODON_FREQ: skipping ${basename} (${locustag}): CDS file not found"
                return null
            }
            tuple(locustag, basename, species, strain, cds)
        }.filter { it != null }
        CODON_FREQ(cds_ch)
    }

    // GFF3 + genome channel for gene structure tables
    if (params.run_gene_stats) {
        def gff_ch = rows_ch.map { locustag, basename, species, strain ->
            def gff    = file("${params.gff_dir}/${basename}.gff3",           glob: false)
            def genome = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
            if (!gff.exists()) {
                log.warn "GENE_STATS: skipping ${basename} (${locustag}): GFF3 not found"
                return null
            }
            if (!genome.exists()) {
                log.warn "GENE_STATS: skipping ${basename} (${locustag}): genome FASTA not found"
                return null
            }
            tuple(locustag, basename, species, strain, gff, genome)
        }.filter { it != null }
        GENE_STATS(gff_ch)
    }

    // line-indexed channel for chrom info (one job per samples.csv row)
    if (params.run_chrom_info) {
        def chrom_ch = Channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .withIndex()   // emits [row, 0-based index]
            .map { row, idx ->
                def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
                tuple(idx, locustag)
            }
            .take(params.n_test > 0 ? params.n_test as int : -1)
        CHROM_INFO(chrom_ch)
    }
}
