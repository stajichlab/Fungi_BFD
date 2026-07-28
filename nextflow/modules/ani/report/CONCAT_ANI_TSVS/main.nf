// Concatenate multiple per-group ANI TSVs into one merged TSV and write a
// companion asmid-manifest so callers can detect dataset changes.
//
// Each input TSV:  query<TAB>ref<TAB>ANI[<TAB>AF_ref<TAB>AF_query]
// Output:          same 3-column format (query, ref, ANI) with header row.
// Companion manifest (all_pairs_merged.asmid_manifest.txt):
//   one asmid per line — the set of genomes that were in the input TSVs.
process CONCAT_ANI_TSVS {
    tag   "CONCAT_ANI_TSVS"
    label 'report'
    publishDir "${params.outdir}/${params.ani_method}/${params.compare}", mode: 'copy', overwrite: true

    input:
        path manifest  // file listing .ani.tsv paths, one per line
        val asmids     // list of ASMIDs that were used to produce the input TSVs

    output:
        path("all_pairs_merged.tsv"),           emit: out
        path("all_pairs_merged.asmid_manifest.txt"), emit: asmid_manifest

    script:
    // Dedupe/sort in Groovy (asmids is already available at template time) so the
    // shell side never has to round-trip asmid strings through word-splitting or
    // regex, and never needs a literal `$` inside this GString script block.
    def uniqAsmids = (asmids ?: [])
        .findAll { it }
        .collect { it.toString().trim() }
        .findAll { it }
        .unique()
        .sort()
    def asmidManifestCmd = uniqAsmids
        ? "cat > all_pairs_merged.asmid_manifest.txt << 'ASMID_EOF'\n${uniqAsmids.join('\n')}\nASMID_EOF"
        : "touch all_pairs_merged.asmid_manifest.txt"
    """
    # Write header, then all data rows from every file.
    printf 'query\\tref\\tANI\\n' > all_pairs_merged.tsv
    xargs -a "${manifest}" -I{} sh -c '[ -s "{}" ] && cat "{}"' >> all_pairs_merged.tsv

    # Write the asmid manifest so callers can detect when the dataset grew.
    ${asmidManifestCmd}
    """

    stub:
    """
    printf 'q\\tr\\t99.0\\n' > all_pairs_merged.tsv
    touch all_pairs_merged.asmid_manifest.txt
    """
}
