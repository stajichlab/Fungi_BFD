include { hashBucketForType } from '../../common/utils.nf'

// A single task processes a freq_batch_size-sized batch of genomes (amortises
// SLURM startup, see BFD_GENOME_STATS.nf), each hashing to a potentially
// different bucket -- so a single storeDir directive can't route them, and
// the script itself creates each genome's bucket subdirectory and writes
// directly into it. storeDir preserves the relative subdirectory structure
// captured by the "*/*.aa_freq.csv.gz" glob (verified via a standalone probe,
// same as publishDir), so the bucketed layout survives the copy into
// storeDir's root. Key is LOCUSTAG, matching planFreq()'s expected path
// (BFD_GENOME_STATS.nf) -- must stay in sync or every genome misclassifies as
// "todo" regardless of what's already on disk.
process BATCH_AA_FREQ {
    label      'genestats'
    // storeDir, not publishDir: publishDir copies asynchronously AFTER the task
    // completes, so `out.csv` emitted paths that did not exist yet and the
    // downstream manifest sampled exists()/size() on a file still being copied --
    // MERGE_AA_FREQ then fired or not at random. storeDir moves outputs into
    // the same directory and emits the stored paths, so they are guaranteed present.
    storeDir "${params.genome_stats_outdir}/aa_freq"

    input:
    tuple val(metas), path(prots)

    output:
    path "*/*.aa_freq.csv.gz", emit: csv

    script:
    def locustagList = (metas instanceof List ? metas : [metas]).collect { m -> m.locustag }
    def protList = prots     instanceof List ? prots     : [prots]
    def cmds = [locustagList, protList].transpose().collect { locustag, prot ->
        def bucket = hashBucketForType('aa_freq', locustag)
        "mkdir -p ${bucket} && python3 ${params.scripts}/calculate_AA_freq.py ${prot.name} -o ${bucket}/${locustag}.aa_freq.csv.gz"
    }.join('\n')
    """
    module load biopython
    ${cmds}
    """

    stub:
    def locustagList = (metas instanceof List ? metas : [metas]).collect { m -> m.locustag }
    """
    ${ locustagList.collect { locustag ->
        def bucket = hashBucketForType('aa_freq', locustag)
        "mkdir -p ${bucket} && printf 'species_prefix,amino_acid,frequency\\n' | gzip > ${bucket}/${locustag}.aa_freq.csv.gz"
    }.join('\n') }
    """
}
