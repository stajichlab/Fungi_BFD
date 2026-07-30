// Download and normalize single-end RNA-seq for species with no paired-end data,
// or where rnaseq_blacklist.csv lists accessions with action=SE_trinity (SRA metadata
// says PAIRED but the data is actually single-end).
// Outputs zero-byte PE stubs so all read channels carry the same 4-tuple shape.
// storeDir note: if a species was previously fetched via SRA_FETCH (PE), re-routing
// it to SRA_FETCH_SE requires deleting rnaseq_reads/${species_tag}_norm_* files first.
process SRA_FETCH_SE {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   32
    memory '96 GB'
    time   '2h'

    input:
    tuple val(species_tag), path(sra_query_csv)

    output:
    tuple val(species_tag),
          path("${species_tag}_norm_R1.fastq.gz"),
          path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads

    script:
    """
    module load sratoolkit
    module load parallel-fastq-dump
    module load fastp
    module load BBTools
    module load workspace/scratch

    # Zero-byte PE stubs (this process produces SE output only).
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz

    BLACKLIST="${launchDir}/rnaseq_blacklist.csv"

    acc_action() {
        local acc="\$1"
        [ -f "\$BLACKLIST" ] && awk -F',' -v a="\$acc" 'NR>1 && \$1==a {print \$4; exit}' "\$BLACKLIST" || true
    }

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

    # Collect SE accessions from the query CSV:
    #   SE_trinity entries: col 6 (layout) may be PAIRED in SRA but blacklist says SE_trinity.
    #                       Download with --split-3 and take _1 (real SE mate) if present.
    #   SINGLE layout entries: col 6 == SINGLE; --split-3 gives ACC.fastq.gz (no _1/_2).
    SE_ACCESSIONS=""
    while IFS=',' read -r stag tid acc spots platform layout rest; do
        [ "\$stag" = "species_tag" ] && continue
        ACTION=\$(acc_action "\$acc")
        if [ "\$ACTION" = "SE_trinity" ] || [ "\${layout:-PAIRED}" = "SINGLE" ]; then
            SE_ACCESSIONS="\$SE_ACCESSIONS \$acc"
        fi
    done < ${sra_query_csv}
    SE_ACCESSIONS=\$(echo "\$SE_ACCESSIONS" | tr ' ' '\\n' | grep -v '^\$' | head -n ${params.max_rnaseq_se_runs} | tr '\\n' ' ' | xargs)

    if [ -z "\$SE_ACCESSIONS" ]; then
        echo "[INFO] No SE accessions found for ${species_tag}"
        exit 0
    fi

    echo "[INFO] SE accessions for ${species_tag}: \$SE_ACCESSIONS"
    TMPDIR=\${SCRATCH:-/tmp}
    mkdir -p reads

    for ACC in \$SE_ACCESSIONS; do
        ACTION=\$(acc_action "\$ACC")
        echo "[INFO] Downloading \$ACC (SE mode, action: \${ACTION:-default}) ..."

        # Download with --split-3 (NOT --split-files): --split-files forces paired
        # splitting and fails on genuine single-end runs, which have no second mate to
        # split out. --split-3 routes reads by what each spot actually contains:
        #   genuine SINGLE layout        -> all reads land in ACC.fastq.gz (no _1/_2)
        #   SE_trinity, 2 reads per spot -> ACC_1.fastq.gz + ACC_2.fastq.gz (take _1)
        #   SE_trinity, 1 read per spot  -> unpaired reads land in ACC.fastq.gz
        # The SE_FILE detection below picks _1 first (SE_trinity mate), else ACC.fastq.gz.
        parallel-fastq-dump --sra-id "\$ACC" --threads ${task.cpus} \\
            --outdir reads/ --split-3 --gzip --tmpdir "\$TMPDIR" || {
            echo "[WARN] pfd failed for \$ACC; trying EBI FTP..."
            module load aria2
            EBI_DIR=\$(ebi_ftp_dir "\$ACC")
            # SE_trinity (mislabeled PAIRED): EBI has ACC_1.fastq.gz
            # Genuine SINGLE: EBI has ACC.fastq.gz
            aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \\
                "\${EBI_DIR}/\${ACC}_1.fastq.gz" -d reads/ -o "\${ACC}_1.fastq.gz" 2>/dev/null || true
            if [ ! -s "reads/\${ACC}_1.fastq.gz" ]; then
                aria2c --max-connection-per-server=4 --min-split-size=1M --max-tries=3 --retry-wait=5 \\
                    "\${EBI_DIR}/\${ACC}.fastq.gz" -d reads/ -o "\${ACC}.fastq.gz" 2>/dev/null || true
            fi
        }

        SE_FILE=""
        if   [ -s "reads/\${ACC}_1.fastq.gz" ]; then SE_FILE="reads/\${ACC}_1.fastq.gz"
        elif [ -s "reads/\${ACC}.fastq.gz"   ]; then SE_FILE="reads/\${ACC}.fastq.gz"
        fi

        if [ -n "\$SE_FILE" ]; then
            ${params.fastq_hdr_script} --read 1 "\$SE_FILE" \\
                --max-reads ${params.max_rnaseq_reads} \\
                | pigz -c >> "\$TMPDIR/${species_tag}_SE.fastq.gz"
            rm -f "reads/\${ACC}_1.fastq.gz" "reads/\${ACC}_2.fastq.gz" "reads/\${ACC}.fastq.gz"
        else
            echo "[WARN] No SE file found for \$ACC after download; skipping"
            rm -f "reads/\${ACC}"*.fastq.gz
        fi
    done
    rm -rf reads

    if [ ! -s "\$TMPDIR/${species_tag}_SE.fastq.gz" ]; then
        echo "[WARN] No SE reads downloaded for ${species_tag}; leaving empty SE placeholder"
        exit 0
    fi

    bbnorm.sh in="\$TMPDIR/${species_tag}_SE.fastq.gz" \\
        out="\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz" target=30 ecc=f

    fastp --in1 "\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz" \\
          --out1 ${species_tag}_norm_SE.fastq.gz \\
          --thread ${task.cpus} \\
          --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5 \\
          --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5 \\
          --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \\
          --length_required 25

    rm -f "\$TMPDIR/${species_tag}_SE.fastq.gz" "\$TMPDIR/${species_tag}_norm_SE_bbn.fastq.gz"

    NREADS=\$(zcat ${species_tag}_norm_SE.fastq.gz 2>/dev/null | awk 'NR%4==1' | wc -l || echo 0)
    echo "[INFO] \$NREADS normalized SE reads for ${species_tag}"

    mkdir -p "${launchDir}/rnaseq_reads"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[STUB] SRA_FETCH_SE for ${species_tag}"
    """
}
