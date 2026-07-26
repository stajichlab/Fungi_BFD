//
// BFD_GENOME_STATS — per-genome sequence and annotation statistics.
//
// Every process here is storeDir-cached per genome (BUSCO additionally per
// lineage, so switching to a finer-grained lineage produces a distinct file
// without invalidating the coarser one). clearIfStale() deletes any cached
// output older than its input first, so a re-annotated genome recomputes.
//
// The two frequency steps are batched rather than one-job-per-genome to amortise
// SLURM startup: planFreq() decides once, before anything runs, which genomes
// still need computing, and emits the cached ones separately so BFD_MERGE can
// build a merge list that does not race the batch jobs (see plan item 13).
//

include { BATCH_AA_FREQ    } from '../../modules/BFD/BATCH_AA_FREQ/main.nf'
include { BATCH_CODON_FREQ } from '../../modules/BFD/BATCH_CODON_FREQ/main.nf'
include { CALC_INTERGENIC  } from '../../modules/BFD/CALC_INTERGENIC/main.nf'
include { CALC_GENE_STATS  } from '../../modules/BFD/CALC_GENE_STATS/main.nf'
include { CALC_ASM_STATS   } from '../../modules/BFD/CALC_ASM_STATS/main.nf'
include { BUSCO_GENOME     } from '../../modules/BFD/BUSCO_GENOME/main.nf'
include { BUSCO_PEP        } from '../../modules/BFD/BUSCO_PEP/main.nf'

include { clearIfStale } from '../../modules/common/utils.nf'

// Split a per-genome frequency channel into genomes whose output already exists
// and genomes that still need computing.
//
// The exists() test happens exactly once, here, at channel-construction time and
// therefore before any BATCH_* task has run — so both the batch contents and the
// merge's file list are deterministic. Previously the merge rebuilt its file list
// by re-testing exists() on expected output paths while BATCH_* was still
// publishing, and publishDir copies asynchronously *after* a task completes, so
// MERGE_AA_FREQ / MERGE_CODON_FREQ fired or not depending on who won the race.
def planFreq(ch, String subdir, String suffix) {
    ch.map { locustag, basename, input ->
        def out = file("${params.genome_stats_outdir}/${subdir}/${basename}.${suffix}")
        clearIfStale(input, [out])
        tuple(locustag, basename, input, out, out.exists())
    }
    .branch {
        cached: it[4]
        todo  : !it[4]
    }
}

// Collate the not-yet-computed genomes into freq_batch_size jobs.
def batchFreq(todo_ch) {
    todo_ch
        .map { locustag, basename, input, _out, _cached -> tuple(locustag, basename, input) }
        .buffer(size: params.freq_batch_size as int, remainder: true)
        .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }, batch.collect { it[2] }) }
}

workflow BFD_GENOME_STATS {
    take:
    aa_freq_ch       // tuple(locustag, basename, proteins)
    codon_freq_ch    // tuple(locustag, basename, cds)
    intergenic_ch    // tuple(locustag, basename, gff)
    gene_stats_ch    // tuple(locustag, basename, gff, genome)
    asm_stats_ch     // tuple(asmid, basename, genome)
    busco_genome_ch  // tuple(locustag, basename, lineage, genome)
    busco_pep_ch     // tuple(locustag, basename, lineage, proteins)

    main:
    def stats = "${params.genome_stats_outdir}"

    def aa_plan    = null
    def codon_plan = null

    if (params.run_aa_freq.toBoolean()) {
        aa_plan = planFreq(aa_freq_ch, 'aa_freq', 'aa_freq.csv.gz')
        BATCH_AA_FREQ(batchFreq(aa_plan.todo))
    }
    if (params.run_codon_freq.toBoolean()) {
        codon_plan = planFreq(codon_freq_ch, 'codon_freq', 'codon_freq.csv.gz')
        BATCH_CODON_FREQ(batchFreq(codon_plan.todo))
    }

    if (params.run_intergenic.toBoolean())
        CALC_INTERGENIC(intergenic_ch.map { locustag, basename, gff ->
            clearIfStale(gff, [file("${stats}/intergenic_stats/${basename}.gene_intergenic_distances.csv.gz")])
            tuple(locustag, basename, gff)
        })

    if (params.run_gene_stats.toBoolean())
        CALC_GENE_STATS(gene_stats_ch.map { locustag, basename, gff, dna ->
            def outs = ['gene_info', 'gene_exons', 'gene_CDS', 'gene_introns',
                        'gene_transcripts', 'gene_trnas', 'gene_proteins'].collect {
                file("${stats}/gene_stats/${basename}.${it}.csv.gz")
            }
            // GFF3 (re-annotation) and genome FASTA are both inputs; either being
            // newer than the cached outputs should force a re-run.
            clearIfStale(gff, outs)
            clearIfStale(dna, outs)
            tuple(locustag, basename, gff, dna)
        })

    if (params.run_asm_stats.toBoolean())
        CALC_ASM_STATS(asm_stats_ch.map { asmid, basename, genome ->
            clearIfStale(genome, [file("${stats}/asm_stats/${asmid}.stats.txt")])
            tuple(asmid, basename, genome)
        })

    if (params.run_busco_genome.toBoolean())
        BUSCO_GENOME(busco_genome_ch.map { locustag, basename, lineage, genome ->
            clearIfStale(genome, [file("${stats}/BUSCO_genome/${basename}.BUSCO_summary.${lineage}.txt")])
            tuple(locustag, basename, lineage, genome)
        })

    if (params.run_busco_pep.toBoolean())
        BUSCO_PEP(busco_pep_ch.map { locustag, basename, lineage, prot ->
            clearIfStale(prot, [file("${stats}/BUSCO_protein/${basename}.BUSCO_summary.${lineage}.txt")])
            tuple(locustag, basename, lineage, prot)
        })

    emit:
    // Emitted for BFD_MERGE. The *_plan channels carry the genomes that were
    // already cached; the BATCH_* outputs carry the ones this run computed.
    aa_cached      = aa_plan    ? aa_plan.cached.map    { it[3] } : Channel.empty()
    codon_cached   = codon_plan ? codon_plan.cached.map { it[3] } : Channel.empty()
    aa_csv         = params.run_aa_freq.toBoolean()      ? BATCH_AA_FREQ.out.csv    : Channel.empty()
    codon_csv      = params.run_codon_freq.toBoolean()   ? BATCH_CODON_FREQ.out.csv : Channel.empty()
    intergenic_csv = params.run_intergenic.toBoolean()   ? CALC_INTERGENIC.out.csv  : Channel.empty()
    asm_stats      = params.run_asm_stats.toBoolean()    ? CALC_ASM_STATS.out.stats : Channel.empty()
    gene_stats     = params.run_gene_stats.toBoolean()
        ? CALC_GENE_STATS.out.gene_info
              .mix(CALC_GENE_STATS.out.gene_exons)
              .mix(CALC_GENE_STATS.out.gene_CDS)
              .mix(CALC_GENE_STATS.out.gene_introns)
              .mix(CALC_GENE_STATS.out.gene_transcripts)
              .mix(CALC_GENE_STATS.out.gene_trnas)
              .mix(CALC_GENE_STATS.out.gene_proteins)
        : Channel.empty()
}
