#!/usr/bin/env nextflow

/*
 * SOURCE: ../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf
 * Last synced: 2026-05-23
 * Changes vs source: removed nextflow.enable.dsl=2; params block moved to
 *                    conf/profile_funannotate.config.
 *
 * Usage (from project root):
 *   sbatch nextflow/run_funannotate.sh
 *   nextflow run nextflow/funannotate.nf -c nextflow/nextflow.config -profile funannotate -resume
 */

// Metadata tuple order used throughout:
//   val(out), val(asmid), val(species), val(strain), val(locustag),
//   val(busco_lineage), val(header_length), val(transl_table)
// GENOME_CLEAN receives: ..., path(genome_gz), val(taxonid), val(taxondb)
//   → emits: ..., path(genome_fa), val(taxonid)   [storeDir moves .fa; workflow maps to abs string]
//   → writes <asmid>.fa to input_clean_genomes/ (storeDir; skip check targets this file)
//   → purge/FCS intermediates written as side effects to input_clean_genomes/clean/
// MASKREPEAT_TANTAN_RUN receives: ..., val(genome_fa), val(taxonid)
//   → emits: ..., path(masked_fa), val(taxonid)   [storeDir caches input_clean_genomes/<asmid>.masked.fasta]
//   [skipped unless --run_repeatmasker; masked_fa falls back to unmasked .fa if .masked.fasta absent]
// SRA_FETCH receives: val(species_tag), val(taxonid)   [only when --run_sra_fetch; one per species]
//   → emits: val(species_tag), path(norm_R1.fastq.gz), path(norm_R2.fastq.gz)
//   → storeDir caches normalized reads at rnaseq_reads/<species_tag>_norm_{R1,R2}.fastq.gz
//   → empty files (0 bytes) written when no RNA-seq found; downstream checks size to skip
//   → SRA_FETCH handles: download → fastp trim → bbnorm normalization internally
// --stop_after_sra_fetch: when true, pipeline halts after SRA_FETCH (skips RNASEQ_PREPARE,
//   FUNANNOTATE_TRAIN, FUNANNOTATE_PREDICT and all downstream steps).
// RNASEQ_PREPARE receives: ..., val(genome_fa), path(norm_r1), path(norm_r2)   [representative only]
//   → emits: val(species_tag), path(trinity-GG.fasta)   [storeDir caches in rnaseq_data/]
//   → normalized reads stay in rnaseq_reads/ and are NOT re-emitted from RNASEQ_PREPARE
// FUNANNOTATE_TRAIN receives: ..., val(genome_fa), path(norm_r1), path(norm_r2), path(trinity_fa)
//   → norm reads come directly from SRA_FETCH; trinity_fa from RNASEQ_PREPARE
//   → emits: ..., val(genome_fa)
// FUNANNOTATE_PREDICT receives: ..., val(genome_fa)   [from TRAIN or directly after masking/clean]

// Download and extract NCBI taxdump once; storeDir caches it at params.taxondb so
// subsequent runs skip this entirely.
process SETUP_TAXONDB {
    storeDir params.taxondb

    cpus   1
    memory '4 GB'
    time   '1h'

    output:
    path "names.dmp",    emit: ready
    path "nodes.dmp"
    path "merged.dmp"
    path "delnodes.dmp"
    path "division.dmp"
    path "gencode.dmp"
    path "citations.dmp"

    script:
    """
    set -euo pipefail
    wget --no-verbose https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar zxf taxdump.tar.gz
    rm taxdump.tar.gz
    """

    stub:
    """
    for f in names.dmp nodes.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp citations.dmp; do
        touch \$f
    done
    """
}

process GENOME_CLEAN {
    tag "$asmid"

    // container '/rhome/jstajich/projects/AAFTF/AAFTF_v0.6.1-signed.sif'

    // Nextflow skips this task when input_clean_genomes/<asmid>.fa already exists.
    storeDir "${launchDir}/input_clean_genomes"

    cpus   16
    memory '450 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_gz), val(taxonid), val(taxondb)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.fa"), val(taxonid), emit: genome

    script:
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
    source ${launchDir}/scripts/setup_fcs_shm.sh
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    TAXONKIT_DB=${taxondb}
    module load taxonkit
    phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
    	phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
	# weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
    fi
    module unload taxonkit
    echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    module load AAFTF
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
        -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
        -o \$SCRATCH/${asmid}.purge.fasta \
        -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
    mkdir -p ${launchDir}/input_clean_genomes/clean
    cat \$SCRATCH/${asmid}.purge.fasta | \
        ${params.clean_script} --len ${params.min_contig_len} > ${asmid}.fa
    echo "[INFO] Clean genome written: ${asmid}.fa (\$(du -sh ${asmid}.fa | cut -f1))"
    pigz \$SCRATCH/${asmid}.purge.fasta
    pigz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv
    mv \$SCRATCH/${asmid}.purge.fasta.gz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ${launchDir}/input_clean_genomes/clean/
    """

    stub:
    """
    echo ">stub_${asmid}" > ${asmid}.fa
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}

// Soft-mask each assembly using funannotate mask with tantan.
// storeDir caches the masked FASTA alongside the clean genome.
process MASKREPEAT_TANTAN_RUN {
    tag "$asmid"

    storeDir "${launchDir}/input_clean_genomes"

    cpus   8
    memory '16 GB'
    time   '2h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(taxonid)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.masked.fasta"), val(taxonid), emit: masked

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate
    funannotate mask -i ${genome_fa} -o ${asmid}.masked.fasta -m tantan --cpus ${task.cpus}
    """

    stub:
    """
    echo ">stub_${asmid}_masked" > ${asmid}.masked.fasta
    """
}

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

    printf 'species_tag,taxonid,sra_accession,spots\n' > ${species_tag}.sra_query.csv

    esearch -db sra \\
        -query "txid${taxonid}[Organism:noexp] AND RNA-Seq[Strategy] AND PAIRED[Layout] AND 00000000075[ReadLength] : 00000000300[ReadLength] AND (BGISEQ[Platform] OR Illumina[Platform])" | \\
        efetch -format runinfo > _runinfo.tmp

    awk -F',' 'NR>1 && \$13=="RNA-Seq" && \$16=="PAIRED" && \$1~/^[SDE]RR/ && \$4+0>=250000 {printf "%s,%s\\n", \$1, \$4}' _runinfo.tmp | \\
        sort -t',' -k2 -rn | \\
        head -n 5 | \\
        while IFS=',' read -r acc spots; do
            printf '%s,%s,%s,%s\\n' "${species_tag}" "${taxonid}" "\$acc" "\$spots"
        done >> ${species_tag}.sra_query.csv

    rm -f _runinfo.tmp
    NHITS=\$(awk 'END{print NR-1}' ${species_tag}.sra_query.csv)
    echo "[INFO] Found \$NHITS SRA accessions for ${species_tag} (taxonid=${taxonid})"
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots\n' > ${species_tag}.sra_query.csv
    printf '%s,%s,SRR000001,1000000\n' "${species_tag}" "${taxonid}" >> ${species_tag}.sra_query.csv
    echo "[STUB] SRA_QUERY for ${species_tag}"
    """
}

// Merge all per-species SRA query CSVs into a single named manifest.
// Output: {stem}.rnaseq_sra.csv written alongside the input samples file.
// Columns: species_tag, taxonid, sra_accession, spots
process COLLECT_SRA_QUERY {
    publishDir { file(params.samples).parent.toAbsolutePath().toString() }, mode: 'copy'

    cpus   1
    memory '1 GB'
    time   '10m'

    input:
    path(query_csvs)
    val(stem)

    output:
    path("${stem}.rnaseq_sra.csv"), emit: manifest

    script:
    """
    printf 'species_tag,taxonid,sra_accession,spots\n' > ${stem}.rnaseq_sra.csv
    for f in ${query_csvs}; do
        tail -n +2 "\$f" >> ${stem}.rnaseq_sra.csv
    done
    NSPECIES=\$(awk -F',' 'NR>1{print \$1}' ${stem}.rnaseq_sra.csv | sort -u | wc -l)
    NACCESSIONS=\$(awk 'NR>1' ${stem}.rnaseq_sra.csv | wc -l)
    echo "[INFO] ${stem}.rnaseq_sra.csv: \$NACCESSIONS accessions across \$NSPECIES species with RNA-seq data"
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots\n' > ${stem}.rnaseq_sra.csv
    """
}

// Write zero-byte paired FASTQ placeholder files for species with no SRA data.
// Called only for species whose SRA_QUERY CSV has no data rows, avoiding a
// SLURM job allocation for what would be an immediate empty-file write.
process WRITE_EMPTY_READS {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   1
    memory '1 GB'
    time   '5m'

    input:
    val(species_tag)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"), emit: reads

    script:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    echo "[INFO] No SRA data for ${species_tag}; created empty read placeholders"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    """
}

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
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"), emit: reads

    script:
    """
    module load sratoolkit
    module load parallel-fastq-dump
    module load fastp
    module load BBTools
    module load workspace/scratch

    # Output files must always exist (storeDir requirement).
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz

    # Read pre-queried accessions from SRA_QUERY CSV (up to max_rnaseq_runs).
    ACCESSIONS=\$(awk -F',' 'NR>1 {print \$3}' ${sra_query_csv} | head -n ${params.max_rnaseq_runs} | tr '\n' ' ')
    TAXONID=\$(awk -F',' 'NR==2 {print \$2; exit}' ${sra_query_csv})

    if [ -z "\$(echo \$ACCESSIONS | tr -d ' ')" ]; then
        echo "[INFO] No paired-end RNA-seq runs found for ${species_tag} (no accessions in query CSV)"
    else
        echo "[INFO] SRA accessions for ${species_tag}: \$ACCESSIONS"
        TMPDIR=\${SCRATCH:-/tmp}
        mkdir -p reads

        # Download and concatenate in accession order so R1/R2 stay matched.
        for ACC in \$ACCESSIONS; do
            echo "[INFO] Downloading \$ACC ..."
            parallel-fastq-dump --sra-id \$ACC --threads ${task.cpus} \
                --outdir reads/ --split-files --gzip --tmpdir \$TMPDIR || {
                echo "[WARN] Download failed for \$ACC, skipping"
                continue
            }
            if [ -f reads/\${ACC}_1.fastq.gz ] && [ -f reads/\${ACC}_2.fastq.gz ]; then
                # if we could run these in parallel?
                parallel -j 2 ${params.fastq_hdr_script} --read {} reads/\${ACC}_{}.fastq.gz \
		\\| head -n ${params.max_rnaseq_reads} \\|
                    \\| pigz -c \\>\\> \$TMPDIR/${species_tag}_R{}.fastq.gz  ::: 1 2
                rm reads/\${ACC}_[12].fastq.gz
            else
                echo "[WARN] Missing pair for \$ACC after download, skipping"
            fi
        done
        rm -rf reads
	${launchDir}/scripts/enforce_seqpair_readlen in=\$TMPDIR/${species_tag}_R1.fastq.gz \
	in2=\$TMPDIR/${species_tag}_R2.fastq.gz out=\$TMPDIR/${species_tag}_trunc_R1.fastq.gz \
	out2=\$TMPDIR/${species_tag}_trunc_R2.fastq.gz minlen=75
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

        # Append provenance manifest.
        mkdir -p "${launchDir}/rnaseq_reads"
        MANIFEST="${launchDir}/rnaseq_reads/rnaseq_manifest.tsv"
        if [ ! -f "\$MANIFEST" ]; then
            printf "species_tag\ttaxonid\taccessions\ttimestamp\n" > "\$MANIFEST"
        fi
        printf "%s\t%s\t%s\t%s\n" \
            "${species_tag}" "\$TAXONID" \
            "\$(echo \$ACCESSIONS | tr '[:space:]' ',' | sed 's/,\$//')" \
            "\$(date -Iseconds)" >> "\$MANIFEST"
    fi
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    echo "[STUB] SRA_FETCH for ${species_tag}"
    """
}

// Run funannotate train on the representative (first) assembly of each species, then
// archive the Trinity-GG transcripts (normalized reads are in rnaseq_reads)
// reads into rnaseq_data/ so all other strains can skip those expensive steps.
// storeDir skips this process entirely if all five output files already exist.
process RNASEQ_PREPARE {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(species_tag), val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(r1), path(r2)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared

    script:
    """
    # ── Empty-reads sentinel: no RNA-seq found by SRA_FETCH ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${species_tag}; writing empty shared markers"
        touch ${species_tag}.trinity-GG.fasta
        exit 0
    fi

    # ── If representative was already trained, just extract shared files ──────
    TRAIN_GFF3="${params.target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; extracting shared files to rnaseq_data"
        TRAINDIR="${params.target}/${out}/training"
        TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
        if [ -n "\$TRINITY_FA" ]; then
            cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
        else
            touch ${species_tag}.trinity-GG.fasta
        fi
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate
    module load fastp

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ── Run full funannotate train on the representative genome ───────────────
    # Use SCRATCH for the funannotate output dir so Trinity/HISAT2/normalize
    # intermediates land on fast local storage and don't consume project quota.
    echo "[INFO] RNASEQ_PREPARE: running funannotate train for representative ${out} (species: ${species_tag})"

    funannotate train -i ${genome_fa} -o \$SCRATCH/${out} \\
        --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
        --species "${species}" --strain "${strain}" \\
        --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
        --header_length ${header_length} \\
        --jaccard_clip --no-progress --min_coverage 4 \\
        --max_intronlen ${params.max_intronlen} \\
        --stop_after_trinity --no_trimmomatic

    # ── Copy shared outputs to rnaseq_data/ ──────────────────────────────────
    TRAINDIR="\$SCRATCH/${out}/training"
    TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
    if [ -n "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
    else
        echo "[WARN] No trinity.fasta found under \$TRAINDIR for ${out}"
        touch ${species_tag}.trinity-GG.fasta
    fi

    # ── Clean up scratch output dir (all intermediates were temporary) ────────
    rm -rf "\$SCRATCH/${out}"
    echo "[INFO] RNASEQ_PREPARE complete for ${species_tag}"
    """

    stub:
    """
    echo ">stub_trinity_${species_tag}" > ${species_tag}.trinity-GG.fasta
    mkdir -p ${params.target}/${out}/training
    touch ${params.target}/${out}/training/funannotate_train.pasa.gff3
    """
}

// For non-representative strains: funannotate train --trinity <shared_fasta> runs only
// PASA (skips Trimmomatic, normalization, HISAT2, and Trinity-GG assembly).
// Falls back to a full train when no shared Trinity is available (e.g. species with
// a single strain or when run_sra_fetch is false).
process FUNANNOTATE_TRAIN {
    tag "$out"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), path(r1), path(r2), path(trinity_fa)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    script:
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if training output already present ───────────────────────────────
    TRAIN_GFF3="${params.target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; skipping"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # ──  may be unnecessary if overridden by -B option later? ──
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/mysql_db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    if [ -s "${trinity_fa}" ]; then
        echo "[INFO] Running funannotate train (PASA only) for ${out} using shared Trinity from rnaseq_data"
        funannotate train -i ${genome_fa} -o ${params.target}/${out} \\
            --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (no shared Trinity) for ${out} using pre-normalized reads"
        funannotate train -i ${genome_fa} -o ${params.target}/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi

    # ── Remove large intermediates not needed for predict or update ─────────────
    # Keeps: *.bam, *.bai, *.pasa.gff3, *.stringtie.gtf, *.transcripts.gff3
    TRAINDIR="${params.target}/${out}/training"
    echo "[INFO] Removing large training intermediates in \$TRAINDIR"
    rm -rf "\$TRAINDIR/hisat2"
    rm -rf "\$TRAINDIR/trinity_gg"
    echo "[INFO] Training cleanup complete for ${out}"
    echo "mysql is ${params.pasa_mysql}"
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
    mkdir -p ${params.target}/${out}/training
    touch ${params.target}/${out}/training/funannotate_train.pasa.gff3
    """
}

process FUNANNOTATE_PREDICT {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '32h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa)

    output:
    path("${out}"), emit: dir
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table), emit: metadata

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    if [ "${params.debug}" = "true" ]; then
        echo "[DEBUG] out          = ${out}"
        echo "[DEBUG] asmid        = ${asmid}"
        echo "[DEBUG] species      = ${species}"
        echo "[DEBUG] strain       = ${strain}"
        echo "[DEBUG] locustag     = ${locustag}"
        echo "[DEBUG] busco        = ${busco_lineage}"
        echo "[DEBUG] transl_table = ${transl_table}"
        echo "[DEBUG] proteins     = ${params.proteins}"
        echo "[DEBUG] genome_fa    = ${genome_fa}"
        echo "[DEBUG] TMPDIR       = \$TMPDIR"
        echo "[DEBUG] pwd          = \$(pwd)"
    fi

    # Link training data into work dir so funannotate predict finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.target}/${out}/training" ]; then
        ln -sfn "${params.target}/${out}/training" "${out}/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    funannotate predict --name ${locustag} -i ${genome_fa} --strain "${strain}" \\
        -o ${out} -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 glimmerhmm:0 \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table}

    EXPECTED_GBK="${out}/predict_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate predict did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    if [ -d "${out}/predict_misc/ab_initio_parameters" ]; then
        mv ${out}/predict_misc/ab_initio_parameters ${out}
        mv ${out}/predict_misc/trnascan.no-overlaps.gff3 ${out}
        rm -rf ${out}/predict_misc
        mkdir -p ${out}/predict_misc
        mv ${out}/ab_initio_parameters ${out}/trnascan.no-overlaps.gff3 ${out}/predict_misc
    fi
    find ${out}/predict_results/ -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    # Remove the training symlink so publishDir does not overwrite the real training dir.
    rm -f "${out}/training"
    sync
    """

    stub:
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || { echo "ERROR: genome not found at ${genome_fa}" >&2; exit 1; }
    mkdir -p ${out}/predict_results ${out}/predict_misc
    touch ${out}/predict_results/${out}.gbk ${out}/predict_results/${out}.proteins.fa
    """
}

process ANTISMASH_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/antismash_local/**")

    script:
    def gbk = "${params.target}/${out}/predict_results/${out}.gbk"
    """
    if [ ! -f "${gbk}" ]; then
        echo "ERROR: predict GBK not found: ${gbk}" >&2
        exit 1
    fi
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load antismash
    mkdir -p ${out}/antismash_local
    antismash --taxon ${params.antismash_taxon} \\
        --output-dir ${out}/antismash_local \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        ${gbk}
    pigz ${out}/antismash_local/*.json
    """

    stub:
    """
    mkdir -p ${out}/antismash_local
    touch ${out}/antismash_local/${out}.json.gz
    touch ${out}/antismash_local/index.html
    """
}

// IPRSCAN5
process INTERPROSCAN_RUN {
    tag "$out"

    cpus   8
    memory '32 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/iprscan.xml")

    script:
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    mkdir -p ${out}/annotate_misc
    module load interproscan
    interproscan.sh -i ${proteins} -f XML -o ${out}/annotate_misc/iprscan.xml \\
        -dp -goterms -pa -t p -cpu ${task.cpus}
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}

process SIGNALP_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '12h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/signalp.results.txt")

    script:
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    module load signalp/6-gpu
    TMPDIR=\${SCRATCH:-/tmp}
    signalp6 -od \$TMPDIR/${out}_signalp \\
        -org euk --mode fast -format txt \\
        -fasta ${proteins} \\
        --write_procs ${task.cpus} -bs 16
    mkdir -p ${out}/annotate_misc
    cp \$TMPDIR/${out}_signalp/prediction_results.txt ${out}/annotate_misc/signalp.results.txt
    rm -rf \$TMPDIR/${out}_signalp
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    """
}

process FUNANNOTATE_ANNOTATE {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '48h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}.annotate.done"), emit: marker

    script:
    def antiSm    = file("${params.target}/${out}/antismash_local/${out}.gbk")
    def antiSmArg = antiSm.exists() ? "--antismash ${antiSm}" : ""
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    funannotate annotate -i ${params.target}/${out} -o ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        --header_length ${header_length} \\
        ${antiSmArg} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR

    EXPECTED_GBK="${params.target}/${out}/annotate_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate annotate did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    touch ${out}.annotate.done
    """

    stub:
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${params.target}/${out}/annotate_results ${params.target}/${out}/annotate_misc
    touch ${params.target}/${out}/annotate_results/${out}.gbk
    touch ${out}.annotate.done
    """
}

workflow {
    def suppressSet = (params.suppress && file(params.suppress).exists())
        ? file(params.suppress).readLines()
              .collect { it.trim().split(',')[0].trim() }
              .findAll { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }

    // ── Taxonomy filter ───────────────────────────────────────────────────────
    // Parse --taxon RANK:VALUE (e.g. --taxon PHYLUM:Ascomycota).
    // taxonFilter is a closure applied after splitCsv on the raw row map.
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { row -> true }
    }

    // ── Prediction pipeline ───────────────────────────────────────────────────
    def jobs = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species       = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
            def strain        = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
            strain = strain.replaceAll(/;.*$/, '').trim()
            def out           = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_').replaceAll(/[\[\]\*\?\{\}]/, '_')
            def asmid         = row.ASMID?.trim()
            def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def busco         = row.BUSCO_LINEAGE?.trim()
            def header_length = 24
            def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
            def taxonid       = row.NCBI_TAXONID?.trim()
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid -> out && asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, _tid ->
            if (suppressSet.contains(asmid)) {
                log.info "Suppressing ${out} (asmid=${asmid})"
                return false
            }
            return true
        }
        .map { out, asmid, species, strain, locustag, busco, header_length, transl_table, taxonid ->
            def gz = file("${params.source}/${asmid}/${asmid}_genomic.fna.gz")
            tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)
        }
        .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt, gz, _tid ->
            if (!gz.exists()) {
                log.warn "Missing genome for ${out} (asmid=${asmid}): ${gz}"
                return false
            }
            if (params.debug) {
                log.info "Queuing ${out}: genome=${gz} (${gz.size()} bytes)"
            }
            return true
        }

    if (params.debug) {
        jobs.view { t -> "[CHANNEL] Submitting: out=${t[0]}, asmid=${t[1]}, transl_table=${t[7]}, gz=${t[8]}" }
    }

    // Ensure taxondb is populated before any GENOME_CLEAN task starts.
    // SETUP_TAXONDB uses storeDir so it runs at most once across all pipeline runs.
    SETUP_TAXONDB()
    def taxondb_ch = SETUP_TAXONDB.out.ready.map { params.taxondb }
    GENOME_CLEAN(jobs.combine(taxondb_ch))

    if (!params.only_clean) {
        // Convert path output to absolute-path string so downstream val(genome_fa) processes
        // can reference the file directly without Nextflow re-staging it per-process.
        def clean_genome_ch = GENOME_CLEAN.out.genome
            .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                      genome_fa.toAbsolutePath().toString(), taxonid)
            }

        // ── Repeat masking ────────────────────────────────────────────────────────
        // predict_genome_ch carries the genome path to use for prediction — either
        // the tantan soft-masked genome (default) or the clean unmasked genome
        // (--run_repeatmasker false).
        def predict_genome_ch
        if (params.run_repeatmasker) {
            MASKREPEAT_TANTAN_RUN(clean_genome_ch)
            predict_genome_ch = MASKREPEAT_TANTAN_RUN.out.masked
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, masked_fa, taxonid ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable,
                        masked_fa.toAbsolutePath().toString(), taxonid)
                }
        } else {
            // --run_repeatmasker false: use masked genome if a prior run produced it, else unmasked.
            predict_genome_ch = clean_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    def masked = file("${launchDir}/input_clean_genomes/${asmid}.masked.fasta")
                    def use_fa = masked.exists() ? masked.toString() : genome_fa
                    if (params.debug) {
                        log.info "[DEBUG] ${asmid}: genome_fa=${use_fa} (masked=${masked.exists()})"
                    }
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, use_fa, taxonid)
                }
        }

        // FUNANNOTATE_PREDICT input tuple drops taxonid (not needed after masking/clean).
        // When SRA is enabled: SRA_FETCH fetches reads once per species; RNASEQ_PREPARE runs
        // funannotate train on the representative assembly and archives Trinity-GG, trimmed, and
        // normalized reads to rnaseq_data/; all other strains run FUNANNOTATE_TRAIN --trinity.
        def predict_input_ch
        def reads_ch = Channel.empty()
        if (params.run_sra_fetch) {
            // Build per-species input: group assemblies, keep first taxonid per species.
            def sra_input = predict_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    def species_tag = species.replaceAll(/\s+/, '_')
                    tuple(species_tag, taxonid)
                }
                .groupTuple(by: 0)
                .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

            // Step 1: Lightweight per-species SRA query (cacheable via storeDir)
            SRA_QUERY(sra_input)

            // Step 2: Collect all per-species results into {stem}.rnaseq_sra.csv
            def stem = file(params.samples).baseName
            COLLECT_SRA_QUERY(
                SRA_QUERY.out.query_result.map { _stag, csv -> csv }.collect(),
                stem
            )

            if (!params.stop_after_sra_query) {
            // Step 3: Branch — species with hits go to SRA_FETCH; species without hits
            // get zero-byte placeholder files written by WRITE_EMPTY_READS without a
            // SLURM job for download.
            def branched_sra = SRA_QUERY.out.query_result
                .branch {
                    has_data: it[1].readLines().size() > 1
                    no_data:  true
                }

            SRA_FETCH(branched_sra.has_data)
            WRITE_EMPTY_READS(branched_sra.no_data.map { stag, _csv -> stag })
            reads_ch = SRA_FETCH.out.reads.mix(WRITE_EMPTY_READS.out.reads)

            if (!params.stop_after_sra_fetch) {
            // Build per-assembly channel keyed by species_tag with SRA reads joined.
            def assembly_with_reads = predict_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
                    def species_tag = species.replaceAll(/\s+/, '_')
                    tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
                }
                .combine(reads_ch, by: 0)

            // RNASEQ_PREPARE: run funannotate train --stop_after_trinity once per species on
            // the representative (first) assembly, then cache the Trinity-GG FASTA in rnaseq_data/
            // so all other strains share it. Normalized reads stay in rnaseq_reads/ (SRA_FETCH storeDir).
            // pasa.gff3 is NOT produced here (--stop_after_trinity stops before PASA);
            // it is produced by FUNANNOTATE_TRAIN for every strain including the representative.
            def repr_ch = assembly_with_reads
                .groupTuple(by: 0)
                .map { species_tag, outs, asmids, species_list, strains, locustags,
                       buscos, hlens, ttables, genomes, r1s, r2s ->
                    tuple(species_tag, outs[0], asmids[0], species_list[0], strains[0],
                          locustags[0], buscos[0], hlens[0], ttables[0], genomes[0], r1s[0], r2s[0])
                }
            RNASEQ_PREPARE(repr_ch)

            // Join shared Trinity from rnaseq_data back to every assembly for FUNANNOTATE_TRAIN.
            // Normalized reads (r1/r2) come from SRA_FETCH via assembly_with_reads; they are NOT
            // re-emitted by RNASEQ_PREPARE (they live in rnaseq_reads/ via storeDir).
            def train_input = assembly_with_reads
                .combine(RNASEQ_PREPARE.out.shared, by: 0)
                .map { species_tag, out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, trinity_fa ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa, r1, r2, trinity_fa)
                }

            // Branch on r1 (index 9) and trinity_fa (index 11) file sizes.
            // Assemblies with no RNA-seq bypass FUNANNOTATE_TRAIN entirely and go straight to predict.
            def branched = train_input.branch {
                has_rnaseq: it[9].size() > 0 || it[11].size() > 0
                no_rnaseq:  true
            }
            def predict_no_rnaseq = branched.no_rnaseq
                .map { out, asmid, sp, st, lt, bl, hl, tt, genome_fa, _r1, _r2, _tf ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
                }

            // Skip TRAIN at the channel level when pasa.gff3 already exists and is non-empty.
            // pasa.gff3 is produced by FUNANNOTATE_TRAIN (not RNASEQ_PREPARE) for every strain,
            // including the representative. Size check guards against zero-byte incomplete files.
            def train_todo = branched.has_rnaseq.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _gfa, _r1, _r2, _tf ->
                def gff3 = file("${params.target}/${out}/training/funannotate_train.pasa.gff3")
                !gff3.exists() || gff3.size() == 0
            }
            def train_done = branched.has_rnaseq
                .filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _gfa, _r1, _r2, _tf ->
                    def gff3 = file("${params.target}/${out}/training/funannotate_train.pasa.gff3")
                    gff3.exists() && gff3.size() > 0
                }
                .map { out, asmid, sp, st, lt, bl, hl, tt, genome_fa, _r1, _r2, _tf ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
                }
            FUNANNOTATE_TRAIN(train_todo)
            predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
            } // end if (!params.stop_after_sra_fetch)
            } // end if (!params.stop_after_sra_query)
        } else {
            predict_input_ch = predict_genome_ch
                .map { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, _taxonid ->
                    tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)
                }
        }

        if ((!params.stop_after_sra_fetch && !params.stop_after_sra_query) || !params.run_sra_fetch) {
        def predict_ch = predict_input_ch
            .filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt, _gfa ->
                def f = file("${params.target}/${out}/predict_results/${out}.gbk")
                !f.exists() || f.size() == 0
            }
        FUNANNOTATE_PREDICT(predict_ch)

        // ── Post-predict steps and annotation ────────────────────────────────────
        // postpredict: all samples with a completed predict_results/*.gbk, whether
        // produced in this run or a prior one. This is the source for all optional
        // pre-annotate steps and for FUNANNOTATE_ANNOTATE itself.
        def postpredict = channel.fromPath(params.samples)
            .splitCsv(header: true)
            .filter(taxonFilter)
            .map { row ->
                def species       = row.SPECIES?.trim()?.replaceAll(/['"]/, '')
                def strain        = row.STRAIN?.trim()?.replaceAll(/['"]/, '')
                strain = strain.replaceAll(/;.*$/, '').trim()
                def out           = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_').replaceAll(/[\[\]\*\?\{\}]/, '_')
                def asmid         = row.ASMID?.trim()
                def locustag      = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
                def busco         = row.BUSCO_LINEAGE?.trim()
                def header_length = 24
                def transl_table  = row.TRANSL_TABLE?.trim() ?: '1'
                tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table)
            }
            .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt -> out && asmid }
            .take((params.n_test as int) > 0 ? params.n_test as int : -1)
            .filter { out, asmid, _sp, _st, _lt, _bl, _hl, _tt -> !suppressSet.contains(asmid) }
            .filter { out, _asmid, _sp, _st, _lt, _bl, _hl, _tt ->
                def f = file("${params.target}/${out}/predict_results/${out}.gbk")
                f.exists() && f.size() > 0
            }

        // annotate_ready_ch threads through optional pre-annotate steps. Each optional
        // step splits the channel into "needs to run" vs "already done", processes the
        // former, then mixes the freshly-completed items back. FUNANNOTATE_ANNOTATE only
        // fires once all requested optional steps are complete for a given sample.
        // Joining ANTISMASH/INTERPRO/SIGNALP output back through postpredict reconstructs
        // the metadata tuple while encoding the dependency edge in the channel DAG.
        def annotate_ready_ch = postpredict

        if (params.run_antismash) {
            def as_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                def asDir = file("${params.target}/${out}/antismash_local")
                !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
            }
            def as_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                def asDir = file("${params.target}/${out}/antismash_local")
                asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') }
            }
            ANTISMASH_RUN(as_todo)
            def as_completed = ANTISMASH_RUN.out
                .map { out, _files -> tuple(out, 'done') }
                .join(postpredict)
                .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
            annotate_ready_ch = as_completed.mix(as_done)
        }

        if (params.run_interpro) {
            def ipr_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                !file("${params.target}/${out}/annotate_misc/iprscan.xml").exists()
            }
            def ipr_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                file("${params.target}/${out}/annotate_misc/iprscan.xml").exists()
            }
            INTERPROSCAN_RUN(ipr_todo)
            def ipr_completed = INTERPROSCAN_RUN.out
                .map { out, _xml -> tuple(out, 'done') }
                .join(postpredict)
                .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
            annotate_ready_ch = ipr_completed.mix(ipr_done)
        }

        if (params.run_signalp) {
            def sp_todo = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                !file("${params.target}/${out}/annotate_misc/signalp.results.txt").exists()
            }
            def sp_done = annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                file("${params.target}/${out}/annotate_misc/signalp.results.txt").exists()
            }
            SIGNALP_RUN(sp_todo)
            def sp_completed = SIGNALP_RUN.out
                .map { out, _txt -> tuple(out, 'done') }
                .join(postpredict)
                .map { out, _flag, asmid, sp, st, lt, bl, hl, tt -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
            annotate_ready_ch = sp_completed.mix(sp_done)
        }

        if (params.run_update) {
            if (params.run_sra_fetch) {
                // UPDATE runs from predict results in parallel with antismash/interpro/signalp.
                // Reads are joined from SRA_FETCH (storeDir-cached, so prior-run reads are reused).
                // The join on upd_signal gates annotate_ready_ch so ANNOTATE waits for UPDATE.
                def upd_input = postpredict
                    .map { out, asmid, species, strain, locustag, busco, hlen, ttable ->
                        def species_tag = species.replaceAll(/\s+/, '_')
                        tuple(species_tag, out, asmid, species, strain, locustag, busco, hlen, ttable)
                    }
                    .combine(reads_ch, by: 0)
                    .map { _st, out, asmid, species, strain, locustag, busco, hlen, ttable, r1, r2 ->
                        tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, r1, r2)
                    }
                def upd_todo = upd_input.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _r1, _r2 ->
                    !file("${params.target}/${out}/update_results/${out}.gbk").exists()
                }
                def upd_done_signal = upd_input
                    .filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _r1, _r2 ->
                        file("${params.target}/${out}/update_results/${out}.gbk").exists()
                    }
                    .map { out, _a, _sp, _st, _lt, _bl, _hl, _tt, _r1, _r2 -> tuple(out, 'upd') }
                FUNANNOTATE_UPDATE(upd_todo)
                def upd_signal = FUNANNOTATE_UPDATE.out
                    .map { out, _a, _sp, _st, _lt, _bl, _hl, _tt -> tuple(out, 'upd') }
                    .mix(upd_done_signal)
                annotate_ready_ch = annotate_ready_ch
                    .join(upd_signal)
                    .map { out, asmid, sp, st, lt, bl, hl, tt, _flag -> tuple(out, asmid, sp, st, lt, bl, hl, tt) }
            } else {
                log.warn "run_update=true but run_sra_fetch=false; funannotate update skipped (no reads available)"
            }
        }

        if (params.run_annotate) {
            FUNANNOTATE_ANNOTATE(annotate_ready_ch.filter { out, _a, _sp, _st, _lt, _bl, _hl, _tt ->
                def f = file("${params.target}/${out}/annotate_results/${out}.gbk")
                !f.exists() || f.size() == 0
            })
        }
        } // end if (!params.stop_after_sra_fetch || !params.run_sra_fetch)
    }
}

process FUNANNOTATE_UPDATE {
    tag "$out"

    cpus   16
    memory '96 GB'
    time   '48h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(r1), path(r2)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    script:
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no reads (empty marker file from SRA_FETCH) ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate update"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load funannotate

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Link training data into work dir so funannotate update finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.target}/${out}/training" ]; then
        ln -sfn "${params.target}/${out}/training" "${out}/training"
    fi

    # r1/r2 are pre-normalized reads from SRA_FETCH (fastp-trimmed + bbnorm-normalized).
    # funannotate update will still run its internal alignment step against these.
    echo "[INFO] Running funannotate update for ${out}"
    funannotate update -i ${params.target}/${out} \\
        --left ${r1} --right ${r2} \\
        --cpus ${task.cpus} \\
        \$pasa_db_arg
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    EXPECTED="${params.target}/${out}/update_results/${out}.gbk"
    if [ ! -f "\$EXPECTED" ]; then
        echo "ERROR: funannotate update did not produce expected GBK: \$EXPECTED" >&2
        exit 1
    fi
    """

    stub:
    """
    echo "[STUB] FUNANNOTATE_UPDATE stub for ${out} (r1=${r1}, r2=${r2})"
    mkdir -p ${params.target}/${out}/update_results
    touch ${params.target}/${out}/update_results/${out}.tbl
    touch ${params.target}/${out}/update_results/${out}.gbk
    touch ${params.target}/${out}/update_results/${out}.gff3
    """
}
