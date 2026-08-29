// Build one composite transcript-evidence FASTA per hybrid-cross species_tag by
// concatenating its PARENT species' own (already-built, non-hybrid) Trinity-GG
// assemblies. See nextflow/docs/HYBRID_SPECIES_RNASEQ_SKIP_PLAN.md for the full
// rationale: interspecific hybrid genomes are mosaic (variable subgenome copy
// number/introgression per isolate), so sharing one hybrid strain's own admixed
// Trinity assembly across differently-recombined siblings assumes a similarity
// the biology doesn't support. PASA aligns transcripts locus-by-locus against the
// target genome, so a composite of the real parent transcriptomes lets each
// hybrid subgenome block align to whichever parent it actually descends from.
//
// Called once per hybrid species_tag (never per strain) -- FUNANNOTATE_RNASEQ.nf
// fans the resulting composite out to every strain of that hybrid cross the same
// way RNASEQ_PREPARE's representative-built Trinity fans out to ordinary species'
// siblings. storeDir-cached like RNASEQ_PREPARE -- inherits the same
// representative-identity invalidation gap noted in
// DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md Option 2 (no sentinel yet); if the
// resolved parent set for a hybrid ever changes, the cached composite needs
// manual invalidation until that sentinel exists.
process BUILD_HYBRID_COMPOSITE_TRINITY {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    cpus   1
    memory '4 GB'
    time   '30m'

    input:
    tuple val(species_tag), path(parent_fastas)

    output:
    tuple val(species_tag),
          path("${species_tag}.composite-parents.trinity-GG.fasta"), emit: shared

    script:
    """
    cat ${parent_fastas} > ${species_tag}.composite-parents.trinity-GG.fasta
    NPARENTS=\$(echo ${parent_fastas} | wc -w)
    NTRANSCRIPTS=\$(grep -c '^>' ${species_tag}.composite-parents.trinity-GG.fasta || true)
    echo "[INFO] Built composite parent-transcript evidence for ${species_tag} from \$NPARENTS parent Trinity file(s), \$NTRANSCRIPTS transcripts total"
    """

    stub:
    """
    echo ">stub_composite_${species_tag}" > ${species_tag}.composite-parents.trinity-GG.fasta
    """
}
