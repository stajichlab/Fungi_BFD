//
// FUNANNOTATE_GENOME_PREP — taxon DB setup, genome cleaning, repeat masking.
//
// Emits the genome each sample should be predicted from: the tantan soft-masked
// genome by default, or the clean unmasked genome under --run_repeatmasker false.
//
// Under --only_clean the emitted channel is empty, so the caller runs nothing
// downstream (the legacy workflow expressed this by wrapping everything after
// cleaning in `if (!params.only_clean)`).
//

include { SETUP_TAXONDB         } from '../../modules/funannotate/setup/SETUP_TAXONDB/main.nf'
include { GENOME_CLEAN          } from '../../modules/funannotate/genome/GENOME_CLEAN/main.nf'
include { GENOME_CLEAN_BATCH    } from '../../modules/funannotate/genome/GENOME_CLEAN_BATCH/main.nf'
include { MASKREPEAT_TANTAN_RUN } from '../../modules/funannotate/genome/MASKREPEAT_TANTAN_RUN/main.nf'

include { genomeFile } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_GENOME_PREP {
    take:
    jobs   // tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, gz, taxonid)

    main:
    // Ensure taxondb is populated before any GENOME_CLEAN task starts.
    // SETUP_TAXONDB uses storeDir so it runs at most once across all pipeline runs.
    SETUP_TAXONDB()
    def taxondb_ch = SETUP_TAXONDB.out.ready.map { params.taxondb }

    // Only clean genomes whose cleaned .fa does not already exist. This keeps batches
    // from being padded with finished genomes — a batch that is entirely cleaned is never
    // scheduled, so it never pays the ~30-min /dev/shm staging cost. (GENOME_CLEAN_BATCH
    // also re-checks per genome at runtime, which handles partial completion on retry.)
    def jobs_to_clean = jobs.filter { tup ->
        !genomeFile("${launchDir}/input_clean_genomes/${tup[1]}.fa").exists()
    }

    // Genome cleaning. The FCS-GX DB staging into /dev/shm costs ~30 min per task, so by
    // default we batch genomes (clean_batch_size, default 1000) into single SLURM jobs that
    // stage the DB once and then clean every genome in the batch. Set clean_batch_size = 0
    // to fall back to the original one-SLURM-job-per-genome GENOME_CLEAN process.
    // clean_done_ch gates downstream on cleaning finishing; ifEmpty([]) ensures it still
    // emits (so downstream runs) when every genome was already clean and nothing was scheduled.
    def clean_done_ch
    int clean_batch_size = params.clean_batch_size as int
    if (clean_batch_size > 0) {
        // Wrap each collated batch (a List of per-genome tuples) in a single-element list
        // so .combine() appends taxondb as the 2nd tuple element instead of spreading the
        // batch's rows into the tuple (which broke GENOME_CLEAN_BATCH's `tuple val(items),
        // val(taxondb)` declaration).
        def clean_batches = jobs_to_clean.collate(clean_batch_size).map { batch -> [ batch ] }
        GENOME_CLEAN_BATCH(clean_batches.combine(taxondb_ch), file(params.samples))
        clean_done_ch = GENOME_CLEAN_BATCH.out.manifest.collect().ifEmpty([])
        GENOME_CLEAN_BATCH.out.suppress
            .collectFile(name: 'TO_ADD_TO_SUPRESS.csv', storeDir: launchDir)
    } else {
        GENOME_CLEAN(jobs_to_clean.combine(taxondb_ch))
        clean_done_ch = GENOME_CLEAN.out.genome.map { it[8] }.collect().ifEmpty([])
        GENOME_CLEAN.out.suppress
            .collectFile(name: 'TO_ADD_TO_SUPRESS.csv', storeDir: launchDir)
    }

    // Re-attach the cleaned genome to its full per-sample metadata. The cleaned .fa
    // always lands at input_clean_genomes/<asmid>.fa, so we rebuild from the jobs
    // channel and gate on clean_done_ch (combine waits until all cleaning is done).
    // genome_fa is emitted as an absolute-path string so downstream val(genome_fa)
    // processes reference the file directly without Nextflow re-staging it.
    def clean_genome_ch = jobs
        .map { out, asmid, species, strain, locustag, busco, hlen, ttable, _gz, taxonid ->
            tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                  genomeFile("${launchDir}/input_clean_genomes/${asmid}.fa"), taxonid)
        }
        .combine(clean_done_ch)
        .map { it[0..9] }
        .filter { tup ->
            if (!tup[8].exists()) {
                log.warn "No cleaned genome for ${tup[0]} (asmid=${tup[1]}) — skipping downstream"
                return false
            }
            return true
        }
        .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
            tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                  genome_fa.toAbsolutePath().toString(), taxonid)
        }

    // ── Repeat masking ────────────────────────────────────────────────────────
    // predict_genome_ch carries the genome path to use for prediction — either
    // the tantan soft-masked genome (default) or the clean unmasked genome
    // (--run_repeatmasker false).
    //
    // Declared before the branch: a `def` inside an if-block is not visible to
    // the emit: section.
    def predict_genome_ch

    if (params.only_clean.toBoolean()) {
        // --only_clean: stop after cleaning, run nothing downstream.
        predict_genome_ch = Channel.empty()
    }
    else if (params.run_repeatmasker.toBoolean()) {
        MASKREPEAT_TANTAN_RUN(clean_genome_ch)
        predict_genome_ch = MASKREPEAT_TANTAN_RUN.out.masked
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, masked_fa, taxonid ->
                tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                    masked_fa.toAbsolutePath().toString(), taxonid)
            }
    }
    else {
        // --run_repeatmasker false: use masked genome if a prior run produced it, else unmasked.
        predict_genome_ch = clean_genome_ch
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                def masked = genomeFile("${launchDir}/input_clean_genomes/${asmid}.masked.fasta")
                def use_fa = masked.exists() ? masked.toString() : genome_fa
                if (params.debug.toBoolean()) {
                    log.info "[DEBUG] ${asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                }
                tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, use_fa, taxonid)
            }
    }

    emit:
    predict_genome = predict_genome_ch
}
