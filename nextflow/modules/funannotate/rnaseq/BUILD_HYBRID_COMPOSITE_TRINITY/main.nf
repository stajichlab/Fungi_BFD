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
    path("${species_tag}.composite-parents.parent_prefix_map.tsv"),  emit: prefix_map

    script:
    """
    # Trinity's genome-guided sequence IDs (Trinity_GG_<N>_c<X>_g<Y>_i<Z>) are only
    # unique WITHIN one assembly run -- every independent per-species Trinity run
    # restarts numbering from 1, so a naive \`cat\` of two parents' Trinity-GG FASTAs
    # collides immediately (confirmed 2026-08-28: PASA's upload_transcript_data.dbi
    # died on "Duplicate entry 'Trinity_GG_1_c0_g1_i1' for key 'cdna_acc_idx'" the
    # first time two parents both had that ID). Prefix every header with an 8-char
    # hash of its source parent's species_tag (derived from the staged filename,
    # which RNASEQ_PREPARE/the genus-wide fallback always name
    # \${species_tag}.trinity-GG.fasta) so IDs stay globally unique across the whole
    # composite regardless of parent count, without bloating headers with full
    # (sometimes very long, e.g. 4-way-cross) species names. The hash->species_tag
    # mapping is written alongside for traceability -- an md5sum, not a locustag,
    # since the specific representative strain that built each parent's own
    # Trinity-GG isn't threaded through this far (repr_ch picks it per-species but
    # only its Trinity path reaches here); revisit if per-strain provenance turns
    # out to matter more than per-species.
    : > ${species_tag}.composite-parents.trinity-GG.fasta
    : > ${species_tag}.composite-parents.parent_prefix_map.tsv
    for f in ${parent_fastas}; do
        parent_tag=\$(basename "\$f" .trinity-GG.fasta)
        prefix=\$(echo -n "\$parent_tag" | md5sum | cut -c1-8)
        echo -e "\${prefix}\\t\${parent_tag}" >> ${species_tag}.composite-parents.parent_prefix_map.tsv
        sed "s/^>/>\${prefix}__/" "\$f" >> ${species_tag}.composite-parents.trinity-GG.fasta
    done
    NPARENTS=\$(echo ${parent_fastas} | wc -w)
    NTRANSCRIPTS=\$(grep -c '^>' ${species_tag}.composite-parents.trinity-GG.fasta || true)
    NDUP=\$(grep '^>' ${species_tag}.composite-parents.trinity-GG.fasta | awk '{print \$1}' | sort | uniq -d | wc -l)
    if [ "\$NDUP" -gt 0 ]; then
        echo "[WARN] ${species_tag}: \$NDUP duplicate sequence ID(s) remain after prefixing -- either two parent tags hashed to the same 8-char prefix (vanishingly unlikely at this scale) or a parent Trinity file itself has internal duplicate IDs" >&2
    fi
    echo "[INFO] Built composite parent-transcript evidence for ${species_tag} from \$NPARENTS parent Trinity file(s), \$NTRANSCRIPTS transcripts total"
    """

    stub:
    """
    echo ">stub_composite_${species_tag}" > ${species_tag}.composite-parents.trinity-GG.fasta
    : > ${species_tag}.composite-parents.parent_prefix_map.tsv
    """
}
