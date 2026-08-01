// Download and normalize up to params.max_rnaseq_runs SRA accessions for species
// that have RNA-seq data. Accessions are read from the pre-queried per-species CSV
// produced by SRA_QUERY, so no NCBI network call is made here.
// Only invoked for species with data rows in their SRA_QUERY CSV; WRITE_EMPTY_READS
// handles the no-data case at the channel level.
process SRA_FETCH {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   32
    memory '96 GB'
    time   '2h'

    input:
    tuple val(species_tag), path(sra_query_csv)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads
    path("${species_tag}.se_candidates.csv"), optional: true, emit: se_candidates
    path("${species_tag}.blacklist_candidates.csv"), optional: true, emit: blacklist_candidates

    script:
    """
    module load sratoolkit
    module load parallel-fastq-dump
    module load parallel
    module load fastp
    module load BBTools
    module load workspace/scratch

    # Output files must always exist (storeDir requirement).
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz

    # Read pre-queried accessions from SRA_QUERY CSV (up to max_rnaseq_runs).
    # Consult the project override file for per-accession actions:
    #   skip           — exclude entirely
    #   rename_headers — parallel-fastq-dump as normal, then seqkit renumbers headers
    #                    with sequential integers so R1/R2 names match exactly.
    #                    Fixes BGI/MGISEQ read-name mismatch from block-parallel splits.
    # Accessions absent from the override file use the default parallel-fastq-dump path.
    BLACKLIST="${launchDir}/rnaseq_blacklist.csv"
    RAW_ACCESSIONS=\$(awk -F',' 'NR>1 {print \$3}' ${sra_query_csv} | head -n ${params.max_rnaseq_runs})
    TAXONID=\$(awk -F',' 'NR==2 {print \$2; exit}' ${sra_query_csv})

    # Helper: record an accession that yielded a non-empty _1 but no _2 (single-end data
    # mislabeled PAIRED in SRA). Writes a row in blacklist format so it can be pasted
    # straight into rnaseq_blacklist.csv (acc,species_tag,taxonid,SE_trinity); a trailing
    # spots column (SRA-reported, col 4 of the query CSV; empty if unknown) rides along as
    # info and is ignored by the blacklist parser. Collected to rnaseq_se_candidates.csv.
    flag_se_candidate() {
        local acc="\$1" spots
        spots=\$(awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$4; exit}' ${sra_query_csv})
        echo "[SE_CANDIDATE] \$acc ${species_tag} spots=\${spots:-NA} — add to rnaseq_blacklist.csv as SE_trinity"
        echo "\$acc,${species_tag},\$TAXONID,SE_trinity,\${spots}" >> ${species_tag}.se_candidates.csv
    }

    # Helper: record an accession whose download failed outright (parallel-fastq-dump and the
    # EBI FTP fallback both produced nothing usable). Writes a row in rnaseq_blacklist.csv
    # column order (sra_accession,species_tag,taxonid,action) with action=skip so it can be
    # pasted straight into the blacklist after review; the sra_query cache uses a different
    # column order, hence the re-ordering here. Trailing spots column is info only and is
    # ignored by the blacklist parser. Collected to rnaseq_blacklist_candidates.csv.
    flag_blacklist_candidate() {
        local acc="\$1" spots
        spots=\$(awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$4; exit}' ${sra_query_csv})
        echo "[BLACKLIST_CANDIDATE] \$acc ${species_tag} spots=\${spots:-NA} — download failed; add to rnaseq_blacklist.csv as skip"
        echo "\$acc,${species_tag},\$TAXONID,skip,\${spots}" >> ${species_tag}.blacklist_candidates.csv
    }

    # Helper: look up the explicit override action for an accession from the blacklist
    # (col 4 = action: skip | rename_headers).  Returns empty string if not listed.
    acc_action() {
        local acc="\$1"
        [ -f "\$BLACKLIST" ] && awk -F',' -v a="\$acc" 'NR>1 && \$1==a {print \$4; exit}' "\$BLACKLIST" || true
    }

    # Helper: look up sequencing platform from the per-species sra_query CSV (col 5).
    # Returns the raw SRA platform string, e.g. BGISEQ, ILLUMINA.
    acc_platform() {
        local acc="\$1"
        awk -F',' -v a="\$acc" 'NR>1 && \$3==a {print \$5; exit}' ${sra_query_csv}
    }

    # Resolve the effective download strategy for an accession:
    #   explicit blacklist action takes precedence;
    #   BGISEQ platform auto-maps to rename_headers;
    #   everything else is empty (default parallel-fastq-dump path).
    acc_strategy() {
        local acc="\$1" action platform
        action=\$(acc_action "\$acc")
        if [ -n "\$action" ]; then
            echo "\$action"
        else
            platform=\$(acc_platform "\$acc")
            case "\$platform" in
                BGISEQ|BGI*|MGI*) echo "rename_headers" ;;
                *)                 echo "" ;;
            esac
        fi
    }

    # Compute the EBI FTP directory URL for an SRA accession.
    # Path formula: vol1/fastq/<first-6-chars>/<subdir>/<acc>
    # subdir is omitted for ≤9-char accessions; for longer ones it zero-pads
    # chars 10-onward of the accession string into a 3-char field.
    ebi_ftp_dir() {
        local acc="\$1" len prefix base
        len="\${#acc}"
        prefix="\${acc:0:6}"
        base="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/\${prefix}"
        if   [ "\$len" -le 9  ]; then printf '%s/%s'      "\$base" "\$acc"
        elif [ "\$len" -eq 10 ]; then printf '%s/00%s/%s'  "\$base" "\${acc:9:1}" "\$acc"
        elif [ "\$len" -eq 11 ]; then printf '%s/0%s/%s'   "\$base" "\${acc:9:2}" "\$acc"
        else                          printf '%s/%s/%s'    "\$base" "\${acc:9:3}" "\$acc"
        fi
    }

    # Partition accessions into skip / active.
    # SE_trinity entries are treated as skip here: they route to SRA_FETCH_SE instead.
    ACCESSIONS=""
    for ACC in \$RAW_ACCESSIONS; do
        STRATEGY=\$(acc_strategy "\$ACC")
        ACTION=\$(acc_action "\$ACC")
        if [ "\$STRATEGY" = "skip" ] || [ "\$ACTION" = "SE_trinity" ]; then
            echo "[INFO] Skipping accession \$ACC for ${species_tag} (action: \${ACTION:-\$STRATEGY})"
        else
            ACCESSIONS="\$ACCESSIONS \$ACC"
        fi
    done
    ACCESSIONS=\$(echo "\$ACCESSIONS" | xargs)   # trim whitespace

    if [ -z "\$ACCESSIONS" ]; then
        echo "[INFO] No paired-end RNA-seq runs found for ${species_tag} (no accessions in query CSV or all skipped)"
    else
        echo "[INFO] SRA accessions for ${species_tag}: \$ACCESSIONS"
        TMPDIR=\${SCRATCH:-/tmp}
        mkdir -p reads

        # Download and concatenate in accession order so R1/R2 stay matched.
        for ACC in \$ACCESSIONS; do
            STRATEGY=\$(acc_strategy "\$ACC")
            echo "[INFO] Downloading \$ACC (strategy: \${STRATEGY:-default}) ..."
	    maxspot=""
	    if [ ${params.max_rnaseq_spot} -gt "0" ]; then
	   	maxspot="-X ${params.max_rnaseq_spot}"
	    fi

            if [ "\$STRATEGY" = "rename_headers" ]; then
                # BGI/MGISEQ: parallel-fastq-dump block-splitting desynchronises spot IDs
                # between R1 and R2.  Download identically to the default path, then pipe
                # each file through seqkit to replace every header with a sequential integer
                # before handing off — guarantees R1 read N and R2 read N have matching names.
                module load seqkit
                parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \
                    --outdir reads/ --split-files --gzip --tmpdir \$TMPDIR || {
                    echo "[WARN] Download failed for \$ACC (\$maxspot), skipping"
                    flag_blacklist_candidate "\$ACC"
                    continue
                }
                if [ -f reads/\${ACC}_1.fastq.gz ] && [ -f reads/\${ACC}_2.fastq.gz ]; then
                    seqkit replace -j ${task.cpus} -p '.+' -r '{nr}' reads/\${ACC}_1.fastq.gz \
                        | ${params.fastq_hdr_script} --read 1 /dev/stdin \
                            --max-reads ${params.max_rnaseq_reads} \
                        | pigz -c >> \$TMPDIR/${species_tag}_R1.fastq.gz
                    seqkit replace -j ${task.cpus} -p '.+' -r '{nr}' reads/\${ACC}_2.fastq.gz \
                        | ${params.fastq_hdr_script} --read 2 /dev/stdin \
                            --max-reads ${params.max_rnaseq_reads} \
                        | pigz -c >> \$TMPDIR/${species_tag}_R2.fastq.gz
                    rm reads/\${ACC}_[12].fastq.gz
                elif [ -s reads/\${ACC}_1.fastq.gz ]; then
                    flag_se_candidate "\$ACC"
                    rm -f reads/\${ACC}_1.fastq.gz
                else
                    echo "[WARN] Missing pair for \$ACC after download, skipping"
                    flag_blacklist_candidate "\$ACC"
                fi
            else
                parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \
                    --outdir reads/ --split-files \$maxspot --gzip --tmpdir \$TMPDIR || {
                    echo "[WARN] Download failed for \$ACC (\$maxspot), skipping"
                    flag_blacklist_candidate "\$ACC"
                    continue
                }
                if [ -f reads/\${ACC}_1.fastq.gz ] && [ -f reads/\${ACC}_2.fastq.gz ]; then
                    parallel -j 2 ${params.fastq_hdr_script} --read {} reads/\${ACC}_{}.fastq.gz \
		    --max-reads ${params.max_rnaseq_reads} \
		    \\| pigz -c \\>\\> \$TMPDIR/${species_tag}_R{}.fastq.gz  ::: 1 2
                    rm reads/\${ACC}_[12].fastq.gz
                elif [ -s reads/\${ACC}_1.fastq.gz ]; then
                    flag_se_candidate "\$ACC"
                    rm -f reads/\${ACC}_1.fastq.gz
                else
                    echo "[WARN] Missing pair for \$ACC after download, skipping"
                    flag_blacklist_candidate "\$ACC"
                fi
            fi
        done
        rm -rf reads
        ENFORCE="${projectDir}/bin/enforce_seqpair_readlen"
        [[ -x "\$ENFORCE" ]] || { echo "[ERROR] enforce_seqpair_readlen not found or not executable at \$ENFORCE"; exit 1; }
        if ! "\$ENFORCE" in=\$TMPDIR/${species_tag}_R1.fastq.gz \
                in2=\$TMPDIR/${species_tag}_R2.fastq.gz \
                out=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz \
                out2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz minlen=75; then
            echo "[WARN] enforce_seqpair_readlen failed for ${species_tag} (likely parallel-fastq-dump read-name mismatch); retrying via EBI FTP..."
            module load aria2
            rm -f \$TMPDIR/${species_tag}_R1.fastq.gz \$TMPDIR/${species_tag}_R2.fastq.gz \
                  \$TMPDIR/${species_tag}_trunc_R1.fastq.gz \$TMPDIR/${species_tag}_trunc_R2.fastq.gz
            mkdir -p reads_ebi
            for ACC in \$ACCESSIONS; do
                EBI_DIR=\$(ebi_ftp_dir "\$ACC")
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \
                    "\${EBI_DIR}/\${ACC}_1.fastq.gz" -d reads_ebi/ -o "\${ACC}_1.fastq.gz" || true
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \
                    "\${EBI_DIR}/\${ACC}_2.fastq.gz" -d reads_ebi/ -o "\${ACC}_2.fastq.gz" || true
                if [ ! -s reads_ebi/\${ACC}_2.fastq.gz ]; then
                    if [ -s reads_ebi/\${ACC}_1.fastq.gz ]; then
                        echo "[WARN] \$ACC: no paired-end R2 at EBI but R1 present (single-end data)"
                        flag_se_candidate "\$ACC"
                    else
                        echo "[WARN] \$ACC: no paired-end R2 at EBI (single-end or accession absent); skipping"
                        flag_blacklist_candidate "\$ACC"
                    fi
                    rm -f reads_ebi/\${ACC}_1.fastq.gz reads_ebi/\${ACC}_2.fastq.gz
                    continue
                fi
                if [ ! -s reads_ebi/\${ACC}_1.fastq.gz ]; then
                    echo "[WARN] \$ACC: R1 missing at EBI; skipping"
                    flag_blacklist_candidate "\$ACC"
                    rm -f reads_ebi/\${ACC}_2.fastq.gz
                    continue
                fi
                parallel -j 2 ${params.fastq_hdr_script} --read {} reads_ebi/\${ACC}_{}.fastq.gz \
                    --max-reads ${params.max_rnaseq_reads} \
                    \\| pigz -c \\>\\> \$TMPDIR/${species_tag}_R{}.fastq.gz  ::: 1 2
                rm -f reads_ebi/\${ACC}_[12].fastq.gz
            done
            rm -rf reads_ebi
            "\$ENFORCE" in=\$TMPDIR/${species_tag}_R1.fastq.gz \
                in2=\$TMPDIR/${species_tag}_R2.fastq.gz \
                out=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz \
                out2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz minlen=75 \
                || { echo "[ERROR] enforce_seqpair_readlen failed for ${species_tag} even after EBI FTP fallback"; exit 1; }
        fi
        bbnorm.sh in=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz in2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz \
            out1=\$TMPDIR/${species_tag}_norm_R1.fastq.gz \
            out2=\$TMPDIR/${species_tag}_norm_R2.fastq.gz target=30 ecc=t

        fastp   --in1 \$TMPDIR/${species_tag}_norm_R1.fastq.gz \
                --in2 \$TMPDIR/${species_tag}_norm_R2.fastq.gz \
                --out1 ${species_tag}_norm_R1.fastq.gz --out2 ${species_tag}_norm_R2.fastq.gz \
                --thread ${task.cpus} --detect_adapter_for_pe \
                --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5 \
                --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5 \
                --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
                --length_required 25

        rm \$TMPDIR/${species_tag}_*

        NPAIRS=\$(zcat ${species_tag}_norm_R1.fastq.gz 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0)
        echo "[INFO] Combined \$NPAIRS normalized read pairs for ${species_tag}"

        mkdir -p "${launchDir}/rnaseq_reads"
    fi
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[STUB] SRA_FETCH for ${species_tag}"
    """
}
