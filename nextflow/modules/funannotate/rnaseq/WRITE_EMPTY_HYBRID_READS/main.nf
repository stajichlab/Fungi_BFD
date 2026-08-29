// Write zero-byte paired FASTQ placeholder files for hybrid-cross species, which
// never fetch their own SRA reads (see nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md
// -- they train against composite parent-transcript evidence via
// BUILD_HYBRID_COMPOSITE_TRINITY instead).
//
// This is a DELIBERATELY SEPARATE process/storeDir from WRITE_EMPTY_READS, not an
// alias of it. WRITE_EMPTY_READS writes to storeDir "${launchDir}/rnaseq_reads"
// using the exact filenames SRA_FETCH also writes for a species that DOES have real
// reads (${species_tag}_norm_R1/R2/SE.fastq.gz). Reusing that same storeDir+filename
// pair here would mean: for any hybrid species_tag that had already accumulated real
// reads under rnaseq_reads/ from BEFORE it was classified as hybrid (or from a run
// predating this feature), Nextflow's storeDir sees the file already exists and
// silently keeps serving those stale REAL reads instead of writing an empty one --
// confirmed 2026-08-28 against Saccharomyces_x_bayanus_NBRC1948: its real
// rnaseq_reads/Saccharomyces_x_bayanus_norm_R1.fastq.gz (dated 2026-08-06, ~147K
// read pairs, predating the composite-evidence feature) got silently reused,
// FUNANNOTATE_TRAIN then took the "PASA+PE" branch (composite trinity_fa AND real
// reads both present) instead of the intended readless "PASA only" branch, and
// funannotate train.py's own internal Trinity-GG re-assembly from those reads
// produced 0 transcripts and failed the whole task -- defeating the entire point of
// composite evidence (never touch this strain's/species' own reads). A dedicated
// storeDir under rnaseq_reads/hybrid_empty/ can never collide with SRA_FETCH's
// output paths, regardless of what any prior run (pre-dating this feature, or with
// it toggled differently) left behind.
process WRITE_EMPTY_HYBRID_READS {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads/hybrid_empty"

    cpus   1
    memory '1 GB'
    time   '5m'

    input:
    val(species_tag)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads

    script:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[INFO] Hybrid-cross species ${species_tag}: no own reads used by design (composite parent-transcript evidence instead) -- empty read placeholders"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    """
}
