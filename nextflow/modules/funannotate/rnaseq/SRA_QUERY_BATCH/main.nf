// Batched SRA query: handles params.sra_query_batch_size species per SLURM job.
// maxForks 4 caps concurrent jobs to avoid overwhelming NCBI.
// Per-species esearch/efetch is retried up to 3 times inline with exponential
// backoff before writing an empty CSV.  Existing per-species CSVs already
// present in rnaseq_reads/sra_query/ are reused without re-querying.
//
// The runinfo fetch is capped to the first 250 UIDs (efetch -start 1 -stop
// 250 -- esearch itself has no -retmax flag) and given a 300s timeout per
// attempt. Without the cap, a heavily-sequenced species (e.g.
// Cryptococcus neoformans, Fusarium graminearum, Neurospora crassa -- each
// with 900+ matching SRA records) makes the runinfo fetch take well over the
// timeout, so all 3 attempts fail and an empty CSV gets written and cached
// forever via publishDir -- silently starving that species of RNA-seq even
// though real data exists (confirmed 2026-08-30: SRP609911 alone has 30+
// usable Cryptococcus neoformans PE runs that a capped query now finds).
//
// Candidates are ranked Illumina-over-BGI first, then most-recent ReleaseDate
// first, then spot count as a final tiebreaker -- so within whatever the
// 250-UID window returns, the newest runs win. esearch/sra has no -sort
// option for date order (checked 2026-08-30: only gene/geoprofiles/pubmed/
// sequences/snp dbs list sort choices), so this date-ranking happens after
// the fetch, over whatever subset of UIDs esearch's default ordering (not
// necessarily chronological) put in the first 250.
//
// A best-effort keyword filter on LibraryName/SampleName excludes likely
// host-associated/co-infection or predation-model samples (mouse, macrophage,
// in vivo, blood, tissue, Acanthamoeba, Galleria, zebrafish, C. elegans,
// etc.) so Trinity isn't assembled against a mix of fungal + host/predator
// reads -- e.g. Acanthamoeba castellanii co-culture is a standard
// Cryptococcus virulence assay and would otherwise slip through untagged.
// SRA runinfo has no dedicated "pure culture vs host-associated" field, so
// this is a heuristic, not a guarantee -- treat it as a first pass, not a
// substitute for a curated accession list on heavily-studied species where
// the automated top-5 pick matters (spot-check the resulting sra_query.csv;
// LibraryName/SampleName aren't retained in the CSV, so cross-reference
// accessions against SRA runinfo/BioSample directly if verifying).
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
process SRA_QUERY_BATCH {
    tag "${species_tags[0]}_+${species_tags.size() - 1}_more"

    publishDir "${launchDir}/rnaseq_reads/sra_query", mode: 'copy', overwrite: false

    maxForks 4
    cpus   1
    memory '4 GB'
    time   '4h'

    input:
    tuple val(species_tags), val(taxonids)

    output:
    path("*.sra_query.csv"), emit: query_results

    script:
    def cache_dir  = "${launchDir}/rnaseq_reads/sra_query"
    def batch_args = [species_tags, taxonids].transpose()
                         .collect { st, tid -> "${st}\\t${tid}" }
                         .join('\\n')
    """
    set -uo pipefail
    module load ncbi_edirect

    printf '${batch_args}\\n' > batch_input.tsv

    query_species() {
        local stag="\$1" tid="\$2" attempt

        for attempt in 1 2 3; do
            rm -f "_runinfo_\${stag}.tmp"
            if timeout 300 bash -c \\
                    "esearch -db sra -query 'txid\${tid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])' | efetch -format runinfo -start 1 -stop 250" \\
                    < /dev/null > "_runinfo_\${stag}.tmp"; then
                return 0
            fi
            echo "[WARN] Attempt \${attempt}/3 failed or timed out for \${stag}"
            [ "\${attempt}" -lt 3 ] && sleep \$((attempt * 30))
        done
        rm -f "_runinfo_\${stag}.tmp"
        return 1
    }

    # ── BioSample-level host-association check (second pass) ────────────────
    # The LibraryName/SampleName keyword filter above misses cases where a
    # run's runinfo LibraryName is blank and SampleName is just a bare GEO
    # accession (e.g. GSM4640578) -- the actual "co-culture"/"macrophage"
    # signal only lives in the BioSample's Title/Attributes, which runinfo
    # never exposes. Confirmed 2026-08-30: Trichophyton_rubrum's top-ranked
    # SRR12096668/SRR12096667 passed the LibraryName/SampleName filter clean
    # (both fields uninformative) but are actually THP-1 human-macrophage +
    # T. rubrum co-culture (BioSample SAMN15377118, title "CoTHP1_repIII",
    # cell_type "human THP-1 macrophages and T.rubrum (co-culture)") --
    # 0.01% Hisat2 alignment rate, 0/300 exact-kmer hits against the genome.
    # Only run against the final ranked candidate short-list (~15, not all
    # 250 UIDs) -- ~1.2s per BioSample docsum lookup, so this stays well
    # inside the per-species timeout budget. Fails OPEN (keeps the candidate)
    # on a lookup timeout/empty response, since that's far more likely to be
    # a transient NCBI hiccup than a real signal -- same spirit as the rest
    # of this heuristic: first pass, not a guarantee (see module header).
    HOST_KEYWORD_RE='(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch|thp-?1)'
    biosample_is_host_associated() {
        local biosample="\$1" doc
        [ -z "\$biosample" ] && return 1
        doc=\$(timeout 15 bash -c "esearch -db biosample -query '\${biosample}[Accession]' | efetch -format docsum" 2>/dev/null)
        [ -z "\$doc" ] && return 1
        printf '%s' "\$doc" | tr '[:upper:]' '[:lower:]' | grep -qE "\$HOST_KEYWORD_RE"
    }

    while IFS=\$(printf '\\t') read -r species_tag taxonid; do
        cached="${cache_dir}/\${species_tag}.sra_query.csv"
        if [ -s "\$cached" ]; then
            cp "\$cached" "\${species_tag}.sra_query.csv"
            echo "[INFO] Reusing cached result for \${species_tag}"
            continue
        fi

        if query_species "\${species_tag}" "\${taxonid}"; then
            printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
            # col 1=Run, col 2=ReleaseDate, col 4=spots, col 12=LibraryName, col 13=LibraryStrategy,
            # col 16=LibraryLayout, col 19=Platform, col 26=BioSample, col 30=SampleName.
            # Prepend a platform rank (0=Illumina, 1=BGI/other) so the short-list prefers
            # Illumina, then most-recent ReleaseDate first, then spot count desc as a final
            # tiebreaker. Best-effort keyword filter on LibraryName/SampleName excludes likely
            # host-associated/co-infection samples -- see module header comment. Widened to
            # 15 candidates (not 5) so the BioSample second-pass check below still has room
            # to reject a few and backfill from the rest.
            KEPT=0
            while IFS=',' read -r rank reldate acc spots platform biosample; do
                [ "\$KEPT" -ge 5 ] && break
                if biosample_is_host_associated "\$biosample"; then
                    echo "[INFO] \${species_tag}: excluding \$acc (BioSample \$biosample looks host-associated on second-pass check)"
                    continue
                fi
                printf '%s,%s,%s,%s,%s,PAIRED\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform" >> "\${species_tag}.sra_query.csv"
                KEPT=\$((KEPT + 1))
            done < <(
                awk -F',' '
                    NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
                        meta = " " tolower(\$12 " " \$30) " "
                        if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch)/) next
                        rank = (\$19 ~ /[Ii]llumina/) ? 0 : 1
                        printf "%d,%s,%s,%s,%s,%s\\n", rank, \$2, \$1, \$4, \$19, \$26
                    }' "_runinfo_\${species_tag}.tmp" | \\
                    sort -t',' -k1,1n -k2,2r -k4,4rn | \\
                    head -n 15
            )
            rm -f "_runinfo_\${species_tag}.tmp"
            NHITS=\$(awk 'END{print NR-1}' "\${species_tag}.sra_query.csv")
            echo "[INFO] Found \$NHITS paired-end accessions for \${species_tag} (taxonid=\${taxonid})"
            # SE fallback: if no PE hits and enable_single_end, query SINGLE layout
            if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
                rm -f "_runinfo_se_\${species_tag}.tmp"
                if timeout 300 bash -c \\
                        "esearch -db sra -query 'txid\${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]' | efetch -format runinfo -start 1 -stop 250" \\
                        < /dev/null > "_runinfo_se_\${species_tag}.tmp"; then
                    SE_KEPT=0
                    while IFS=',' read -r reldate acc spots platform biosample; do
                        [ "\$SE_KEPT" -ge "${params.max_rnaseq_se_runs}" ] && break
                        if biosample_is_host_associated "\$biosample"; then
                            echo "[INFO] \${species_tag}: excluding SE \$acc (BioSample \$biosample looks host-associated on second-pass check)"
                            continue
                        fi
                        printf '%s,%s,%s,%s,%s,SINGLE\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform" >> "\${species_tag}.sra_query.csv"
                        SE_KEPT=\$((SE_KEPT + 1))
                    done < <(
                        awk -F',' '
                            NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {
                                meta = " " tolower(\$12 " " \$30) " "
                                if (meta ~ /(mouse|murine| mice | rat |rabbit|macrophage|phagocyt|in.?vivo|infect|co.?infec|co.?cultur|amoeba|acanthamoeba|galleria|zebrafish|elegans| host |blood|serum|plasma|csf|cerebrospinal|lung|brain|spleen|kidney| liver |tissue|biopsy|patient|clinical|autopsy|necropsy|bronch)/) next
                                printf "%s,%s,%s,%s,%s\\n", \$2, \$1, \$4, \$19, \$26
                            }' "_runinfo_se_\${species_tag}.tmp" | \\
                            sort -t',' -k1,1r -k3,3rn | \\
                            head -n 15
                    )
                fi
                rm -f "_runinfo_se_\${species_tag}.tmp"
                NHITS=\$(awk 'END{print NR-1}' "\${species_tag}.sra_query.csv")
                echo "[INFO] SE fallback: \$NHITS single-end accessions for \${species_tag}"
            fi
        else
            printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
            echo "[WARN] All 3 attempts failed for \${species_tag}; writing empty CSV"
        fi
    done < batch_input.tsv
    """

    stub:
    def stub_args = [species_tags, taxonids].transpose()
                        .collect { st, tid -> "${st}\\t${tid}" }
                        .join('\\n')
    """
    printf '${stub_args}\\n' > batch_input.tsv
    while IFS=\$(printf '\\t') read -r species_tag taxonid; do
        printf 'species_tag,taxonid,sra_accession,spots,platform,layout\\n' > "\${species_tag}.sra_query.csv"
        printf '%s,%s,SRR000001,1000000,ILLUMINA,PAIRED\\n' "\${species_tag}" "\${taxonid}" >> "\${species_tag}.sra_query.csv"
    done < batch_input.tsv
    echo "[STUB] SRA_QUERY_BATCH (${species_tags.size()} species)"
    """
}
