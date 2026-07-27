// Query NCBI SRA for available paired-end RNA-seq accessions per species.
// Lightweight: runs the esearch/efetch query only — no downloading.
// Records up to 5 candidates (sorted by spot count desc) in a per-species CSV.
// storeDir caches results so re-runs skip the network query.
// To invalidate the cache for a species, delete rnaseq_reads/sra_query/<species_tag>.sra_query.csv
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

    esearch -db sra \\
        -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])" | \\
        efetch -format runinfo > _runinfo.tmp

    # col 1=Run, col 4=spots, col 13=LibraryStrategy, col 16=LibraryLayout, col 19=Platform
    # Prepend a platform rank (0=Illumina, 1=BGI/other) so the top 5 prefer Illumina,
    # then by spot count desc; BGI/other only fill remaining slots when Illumina runs run out.
    awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {rank=(\$19~/[Ii]llumina/)?0:1; printf "%d,%s,%s,%s\\n", rank, \$1, \$4, \$19}' _runinfo.tmp | \\
        sort -t',' -k1,1n -k3,3rn | \\
        head -n 5 | \\
        while IFS=',' read -r rank acc spots platform; do
            printf '%s,%s,%s,%s,%s,PAIRED\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform"
        done >> ${species_tag}.sra_query.csv

    rm -f _runinfo.tmp
    NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
    echo "[INFO] Found \$NHITS paired-end SRA accessions for ${species_tag} (taxonid=${taxonid})"

    # SE fallback: if no PE hits found and enable_single_end is true, query SINGLE layout
    if [ "${params.enable_single_end}" = "true" ] && [ "\$NHITS" -eq 0 ]; then
        esearch -db sra \\
            -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND SINGLE[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND Illumina[Platform]" | \\
            efetch -format runinfo > _runinfo_se.tmp
        awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="SINGLE" && \$1~/^[SDE]RR/ && \$4+0>=250000 {printf "%s,%s,%s\\n", \$1, \$4, \$19}' _runinfo_se.tmp | \\
            sort -t',' -k2 -rn | \\
            head -n ${params.max_rnaseq_se_runs} | \\
            while IFS=',' read -r acc spots platform; do
                printf '%s,%s,%s,%s,%s,SINGLE\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots" "\$platform"
            done >> ${species_tag}.sra_query.csv
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
