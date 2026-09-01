// Query NCBI SRA for available paired-end RNA-seq accessions per species.
// Lightweight: runs the esearch/efetch query only — no downloading.
// Records up to 5 candidates (ranked Illumina-first, then most-recent
// ReleaseDate, then spot count desc) in a per-species CSV.
// storeDir caches results so re-runs skip the network query.
// To invalidate the cache for a species, delete rnaseq_reads/sra_query/<species_tag>.sra_query.csv
//
// The runinfo fetch is capped to the first 250 UIDs (efetch -start 1 -stop
// 250 -- esearch itself has no -retmax flag) and given a 300s timeout: an
// uncapped query for a heavily-sequenced species (900+ matching SRA records)
// can take longer than a short timeout to fetch runinfo for, and would
// otherwise fail outright. esearch/sra also has no date -sort option, so
// "most recent first" is applied after the fetch, over whatever subset of
// UIDs esearch's default (not necessarily chronological) ordering returned
// in that first-250 window. See SRA_QUERY_BATCH/main.nf for the fuller
// rationale (both modules were fixed together, 2026-08-30).
//
// A best-effort keyword filter on LibraryName/SampleName excludes likely
// host-associated/co-infection or predation-model samples (mouse, macrophage,
// in vivo, blood, tissue, Acanthamoeba, Galleria, zebrafish, C. elegans,
// etc.) so Trinity isn't assembled against a mix of fungal + host/predator
// reads -- e.g. Acanthamoeba castellanii co-culture is a standard
// Cryptococcus virulence assay and would otherwise slip through untagged.
// This is a heuristic, not a guarantee -- SRA runinfo has no dedicated "pure
// culture vs host-associated" field; treat it as a first pass, not a
// substitute for a curated accession list on heavily-studied species.
//
// Second pass (added 2026-08-30): a BioSample docsum lookup against the
// final ~15-candidate short-list catches cases the LibraryName/SampleName
// filter can't -- confirmed for Trichophyton_rubrum, whose top-ranked
// SRR12096668/SRR12096667 had blank LibraryName and a bare GEO accession as
// SampleName (nothing to keyword-match), yet are actually THP-1
// human-macrophage + T. rubrum co-culture reads (BioSample SAMN15377118,
// title "CoTHP1_repIII", cell_type "human THP-1 macrophages and T.rubrum
// (co-culture)") -- 0.01% Hisat2 alignment rate, 0/300 exact-kmer hits
// against the genome. Only the short-list gets the extra lookup (~1.2s each)
// to stay well inside the timeout budget; fails open (keeps the candidate)
// on a lookup timeout, since that's more likely transient than a real signal.
process SRA_QUERY {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads/sra_query"

    cpus   1
    memory '4 GB'
    time   '30m'

    input:
    tuple val(species_tag), val(taxonid)

    output:
    tuple val(species_tag), path("${species_tag}.sra_query.csv"), emit: query_result

    script:
    """
    set -euo pipefail
    module load ncbi_edirect

    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${species_tag}.sra_query.csv

    timeout 300 bash -c "esearch -db sra \\
        -query 'txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])' | \\
        efetch -format runinfo -start 1 -stop 250" > _runinfo.tmp

    # ── BioSample-level host-association check (second pass) ────────────────
    # See SRA_QUERY_BATCH/main.nf header comment for the full rationale
    # (Trichophyton_rubrum/THP-1-macrophage co-culture case, 2026-08-30) --
    # LibraryName/SampleName in runinfo can both be blank/uninformative while
    # the real "co-culture"/"macrophage" signal only lives in the BioSample's
    # Title/Attributes. Only checked against the final ~15-candidate
    # short-list, not all 250 UIDs, to stay well inside the timeout budget.
    HOST_KEYWORD_RE='(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch|thp-?1)'
    biosample_is_host_associated() {
        local biosample="\$1" doc
        [ -z "\$biosample" ] && return 1
        doc=\$(timeout 15 bash -c "esearch -db biosample -query '\${biosample}[Accession]' | efetch -format docsum" 2>/dev/null)
        [ -z "\$doc" ] && return 1
        printf '%s' "\$doc" | tr '[:upper:]' '[:lower:]' | grep -qE "\$HOST_KEYWORD_RE"
    }

    # col 1=Run, col 2=ReleaseDate, col 4=spots, col 12=LibraryName, col 13=LibraryStrategy,
    # col 16=LibraryLayout, col 19=Platform, col 26=BioSample, col 30=SampleName.
    # Prepend a platform rank (0=Illumina, 1=BGI/other) so the short-list prefers Illumina,
    # then most-recent ReleaseDate first, then spot count desc as a final tiebreaker.
    # Best-effort keyword filter on LibraryName/SampleName excludes likely
    # host-associated/co-infection samples -- see module header comment. Widened to 15
    # candidates so the BioSample second-pass check below can reject a few and backfill.
    KEPT=0
    while IFS=',' read -r rank reldate acc spots platform biosample; do
        [ "\$KEPT" -ge 5 ] && break
        if biosample_is_host_associated "\$biosample"; then
            echo "[INFO] ${species_tag}: excluding \$acc (BioSample \$biosample looks host-associated on second-pass check)"
            continue
        fi
        printf '%s,%s,%s,%s,%s,PAIRED\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform" >> ${species_tag}.sra_query.csv
        KEPT=\$((KEPT + 1))
    done < <(
        awk -F',' '
            NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
                meta = " " tolower(\$12 " " \$30) " "
                if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch)/) next
                rank = (\$19 ~ /[Ii]llumina/) ? 0 : 1
                printf "%d,%s,%s,%s,%s,%s\\n", rank, \$2, \$1, \$4, \$19, \$26
            }' _runinfo.tmp | \\
            sort -t',' -k1,1n -k2,2r -k4,4rn | \\
            head -n 15
    )

    rm -f _runinfo.tmp
    NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
    echo "[INFO] Found \$NHITS paired-end SRA accessions for ${species_tag} (taxonid=${taxonid})"

    # SE fallback: if no PE hits found and enable_single_end is true, query SINGLE layout
    if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
        timeout 300 bash -c "esearch -db sra \\
            -query 'txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]' | \\
            efetch -format runinfo -start 1 -stop 250" > _runinfo_se.tmp
        SE_KEPT=0
        while IFS=',' read -r reldate acc spots platform biosample; do
            [ "\$SE_KEPT" -ge "${params.max_rnaseq_se_runs}" ] && break
            if biosample_is_host_associated "\$biosample"; then
                echo "[INFO] ${species_tag}: excluding SE \$acc (BioSample \$biosample looks host-associated on second-pass check)"
                continue
            fi
            printf '%s,%s,%s,%s,%s,SINGLE\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform" >> ${species_tag}.sra_query.csv
            SE_KEPT=\$((SE_KEPT + 1))
        done < <(
            awk -F',' '
                NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
                    meta = " " tolower(\$12 " " \$30) " "
                    if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch)/) next
                    printf "%s,%s,%s,%s,%s\\n", \$2, \$1, \$4, \$19, \$26
                }' _runinfo_se.tmp | \\
                sort -t',' -k1,1r -k3,3rn | \\
                head -n 15
        )
        rm -f _runinfo_se.tmp
        NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
        echo "[INFO] SE fallback: found \$NHITS single-end accessions for ${species_tag}"
    fi
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${species_tag}.sra_query.csv
    printf '%s,%s,SRR000001,1000000,ILLUMINA,PAIRED\n' "${species_tag}" "${taxonid}" >> ${species_tag}.sra_query.csv
    echo "[STUB] SRA_QUERY for ${species_tag}"
    """
}
