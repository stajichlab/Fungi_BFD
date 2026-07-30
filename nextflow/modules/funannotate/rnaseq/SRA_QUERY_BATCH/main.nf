// Batched SRA query: handles params.sra_query_batch_size species per SLURM job.
// maxForks 4 caps concurrent jobs to avoid overwhelming NCBI.
// Per-species esearch/efetch is retried up to 3 times inline with exponential
// backoff before writing an empty CSV.  Existing per-species CSVs already
// present in rnaseq_reads/sra_query/ are reused without re-querying.
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
            if timeout 120 bash -c \\
                    "esearch -db sra -query 'txid\${tid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])' | efetch -format runinfo" \\
                    < /dev/null > "_runinfo_\${stag}.tmp"; then
                return 0
            fi
            echo "[WARN] Attempt \${attempt}/3 failed or timed out for \${stag}"
            [ "\${attempt}" -lt 3 ] && sleep \$((attempt * 30))
        done
        rm -f "_runinfo_\${stag}.tmp"
        return 1
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
            # col 1=Run, col 4=spots, col 13=LibraryStrategy, col 16=LibraryLayout, col 19=Platform
            # Prepend a platform rank (0=Illumina, 1=BGI/other) so the top 5 prefer Illumina,
            # then by spot count desc; BGI/other only fill remaining slots when Illumina runs run out.
            awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {rank=(\$19~/[Ii]llumina/)?0:1; printf "%d,%s,%s,%s\\n", rank, \$1, \$4, \$19}' "_runinfo_\${species_tag}.tmp" | \\
                sort -t',' -k1,1n -k3,3rn | \\
                head -n 5 | \\
                while IFS=',' read -r rank acc spots platform; do
                    printf '%s,%s,%s,%s,%s,PAIRED\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform"
                done >> "\${species_tag}.sra_query.csv"
            rm -f "_runinfo_\${species_tag}.tmp"
            NHITS=\$(awk 'END{print NR-1}' "\${species_tag}.sra_query.csv")
            echo "[INFO] Found \$NHITS paired-end accessions for \${species_tag} (taxonid=\${taxonid})"
            # SE fallback: if no PE hits and enable_single_end, query SINGLE layout
            if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
                rm -f "_runinfo_se_\${species_tag}.tmp"
                if timeout 120 bash -c \\
                        "esearch -db sra -query 'txid\${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]' | efetch -format runinfo" \\
                        < /dev/null > "_runinfo_se_\${species_tag}.tmp"; then
                    awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {printf "%s,%s,%s\\n", \$1, \$4, \$19}' "_runinfo_se_\${species_tag}.tmp" | \\
                        sort -t',' -k2 -rn | \\
                        head -n ${params.max_rnaseq_se_runs} | \\
                        while IFS=',' read -r acc spots platform; do
                            printf '%s,%s,%s,%s,%s,SINGLE\\n' "\${species_tag}" "\${taxonid}" "\$acc" "\$spots" "\$platform"
                        done >> "\${species_tag}.sra_query.csv"
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
