#!/usr/bin/env nextflow
// DSL2 is the default from Nextflow 22+ — no explicit declaration needed

/*
 * Fungi 5K Functional Annotation Pipeline
 * Runs Pfam, CAZy, MEROPS, SignalP, TMHMM, TargetP, IDP (AIUPred), WoLF PSORT,
 * and predGPI on all species in samples.csv, then consolidates each tool's results
 * into a DuckDB-loadable <tool>.csv.gz in tables/.
 *
 * Usage (from project root):
 *   sbatch nextflow/run_functional.sh
 *   nextflow run nextflow/BFD.nf -c nextflow/nextflow.config -profile BFD -resume
 *
 * Dry-run / testing:
 *   nextflow run nextflow/BFD.nf -c nextflow/nextflow.config -profile BFD -stub-run --n_test 2
 *
 * Key params (all defined in nextflow.config / conf/profile_BFD.config):
 *   --merge_all   (default true)  Merge ALL result files in the output subdirs,
 *                                 not just those produced in this run. Use false
 *                                 to restrict to the current run only.
 *   --skip_merge  (default false) Skip all MERGE_* steps entirely.
 */

// All params (samples, pep_dir, outdir, tables, scripts, run_*, merge_all,
// skip_merge, n_test) are defined in nextflow.config — do not redeclare defaults here.

// ════════════════════════════════════════════════════════════════════════════
// RUN PROCESSES  (storeDir → skip automatically if all outputs already exist)
// ════════════════════════════════════════════════════════════════════════════

// ── Pfam ─────────────────────────────────────────────────────────────────────
process RUN_PFAM {
    tag        "${locustag}"
    label      'pfam'
    storeDir   "${params.outdir}/pfam_hmmscan"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.pfam.gz"),    emit: domtbl
        path("${basename}.tblout.gz"),  emit: tblout

    script:
    def mpi_launch = params.pfam_tasks > 1 ? "srun -N ${params.pfam_nodes} -n ${params.pfam_tasks}" : ""
    def mpi_flag   = params.pfam_tasks > 1 ? "--mpi" : ""
    """
    # PFAM_DB and hmmer module loaded by beforeScript; version recorded in trace
    if [ ! -z "${mpi_flag}" ]; then
        module load hmmer/3.4-mpi
    else
        module load hmmer/3.4
    fi
    module load db-pfam
    ${mpi_launch} hmmsearch ${mpi_flag} --cut_ga --noali --cpu ${task.cpus} \\
        --domtbl    ${basename}.pfam \\
        --tblout    ${basename}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${basename}.pfam ${basename}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${basename}.pfam.gz
    printf '' | gzip     > ${basename}.tblout.gz
    """
}

// ── CAZy ─────────────────────────────────────────────────────────────────────
process RUN_CAZY {
    tag        "${locustag}"
    label      'cazy'
    storeDir   { "${params.outdir}/cazy/${basename}" }

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.overview.tsv.gz"),   emit: overview
        path("${basename}.cazymes.tsv.gz"),    emit: cazymes
        path("${basename}.substrates.tsv.gz"), emit: substrates

    script:
    """
    module load dbcanlight
    # the version of dbCAN should be recorded in the metadata for reproducibility
    mkdir -p ${basename}
    dbcanlight search -i ${proteins} -m cazyme -o ${basename} -t ${task.cpus}
    dbcanlight search -i ${proteins} -m sub    -o ${basename} -t ${task.cpus}
    dbcanlight conclude ${basename}
    pigz -f ${basename}/cazymes.tsv ${basename}/substrates.tsv ${basename}/overview.tsv
    mv ${basename}/overview.tsv.gz   ${basename}.overview.tsv.gz
    mv ${basename}/cazymes.tsv.gz    ${basename}.cazymes.tsv.gz
    mv ${basename}/substrates.tsv.gz ${basename}.substrates.tsv.gz
    """

    stub:
    """
    printf 'Gene_ID\\tEC\\tcazyme_fam\\tsub_fam\\tdiamond_fam\\tSubstrate\\t#ofTools\\n' | gzip > ${basename}.overview.tsv.gz
    printf 'HMM_Profile\\tProfile_Length\\tGene_ID\\tGene_Length\\tEvalue\\tProfile_Start\\tProfile_End\\tGene_Start\\tGene_End\\tCoverage\\n' | gzip > ${basename}.cazymes.tsv.gz
    printf 'Gene_ID\\tSubstrate\\n' | gzip > ${basename}.substrates.tsv.gz
    """
}

// ── MEROPS ────────────────────────────────────────────────────────────────────
process RUN_MEROPS {
    tag        "${locustag}"
    label      'merops'
    storeDir   "${params.outdir}/merops"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.blasttab.gz"), emit: blasttab

    script:
    """
    module load ncbi-blast
    module load db-merops
    # the version of MEROPS should be recorded in the metadata for reproducibility
    blastp -query ${proteins} \\
        -db \$MEROPS_DB/merops_scan.lib \\
        -out ${basename}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -outfmt 6 \\
        -use_sw_tback
    pigz ${basename}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${basename}.blasttab.gz
    """
}

// ── SignalP ───────────────────────────────────────────────────────────────────
process RUN_SIGNALP {
    tag        "${locustag}"
    label      'signalp'
    storeDir   "${params.outdir}/signalp"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.signalp.gff3.gz"),         emit: gff3
        path("${basename}.signalp.results.txt.gz"),  emit: results

    script:
    """
    module load signalp/6-gpu
    OUTD=\$(mktemp -d)
    signalp6 -od \$OUTD -org euk --mode fast -format txt \\
        -fasta ${proteins} --write_procs ${task.cpus} -bs 100
    pigz -c \$OUTD/output.gff3             > ${basename}.signalp.gff3.gz
    pigz -c \$OUTD/prediction_results.txt  > ${basename}.signalp.results.txt.gz
    rm -rf \$OUTD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${basename}.signalp.gff3.gz
    printf '# SignalP-6.0\\n'   | gzip > ${basename}.signalp.results.txt.gz
    """
}

// ── TMHMM ─────────────────────────────────────────────────────────────────────
process RUN_TMHMM {
    tag        "${locustag}"
    label      'tmhmm'
    storeDir   "${params.outdir}/tmhmm"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.tmhmm_short.tsv.gz"),   emit: short_tsv
        path("${basename}.tmhmm_results.tsv.gz"), emit: full_tsv

    script:
    """
    module load tmhmm
    tmhmm --noplot         < ${proteins} > ${basename}.tmhmm_results.tsv
    tmhmm --short --noplot < ${proteins} > ${basename}.tmhmm_short.tsv
    pigz ${basename}.tmhmm_results.tsv ${basename}.tmhmm_short.tsv
    """

    stub:
    """
    printf '# TMHMM\\n' | gzip > ${basename}.tmhmm_results.tsv.gz
    printf '# TMHMM\\n' | gzip > ${basename}.tmhmm_short.tsv.gz
    """
}

// ── TargetP ───────────────────────────────────────────────────────────────────
process RUN_TARGETP {
    tag        "${locustag}"
    label      'targetp'
    storeDir   "${params.outdir}/targetP"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}_summary.targetp2.gz"), emit: summary

    script:
    """
    TMPD=\$(mktemp -d)
    module load targetp
    targetp -batch 50 -tmp \$TMPD -format short \\
        -fasta ${proteins} -org non-pl -prefix ${basename}
    pigz -f ${basename}_summary.targetp2
    rm -rf \$TMPD
    """

    stub:
    """
    printf '# TargetP-2.0\\n' | gzip > ${basename}_summary.targetp2.gz
    """
}

// ── IDP (AIUPred) ─────────────────────────────────────────────────────────────
process RUN_IDP {
    tag        "${locustag}"
    label      'idp'
    storeDir   "${params.outdir}/aiupred"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.aiupred.txt.gz"),     emit: raw
        path("${basename}.idp.csv.gz"),         emit: idp_csv
        path("${basename}.idp_summary.csv.gz"), emit: idp_summary_csv

    script:
    """
    module load aiupred
    aiupred.py -i ${proteins} -o ${basename}.aiupred.txt
    pigz ${basename}.aiupred.txt
    python3 ${params.scripts}/gather_AIUPred.py ${basename}.aiupred.txt.gz \\
        --outfile      ${basename}.idp.csv \\
        --outfilesum   ${basename}.idp_summary.csv
    pigz ${basename}.idp.csv ${basename}.idp_summary.csv
    """

    stub:
    """
    printf '' | gzip > ${basename}.aiupred.txt.gz
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > ${basename}.idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > ${basename}.idp_summary.csv.gz
    """
}

// ── WoLF PSORT ────────────────────────────────────────────────────────────────
process RUN_WOLFPSORT {
    tag        "${locustag}"
    label      'wolfpsort'
    storeDir   "${params.outdir}/wolfpsort"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.wolfpsort.results.txt.gz"), emit: results

    script:
    """
    module load wolfpsort
    cat ${proteins} | runWolfPsortSummary fungi > ${basename}.wolfpsort.results.txt
    pigz ${basename}.wolfpsort.results.txt
    """

    stub:
    """
    printf '# WoLF PSORT\\n' | gzip > ${basename}.wolfpsort.results.txt.gz
    """
}

// ── predGPI ───────────────────────────────────────────────────────────────────
process RUN_PREDGPI {
    tag        "${locustag}"
    label      'predgpi'
    storeDir   "${params.outdir}/predgpi"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.predgpi.gff3.gz"), emit: gff3

    script:
    """
    module load predgpi
    predgpi.py -f ${proteins} -m gff3 -o ${basename}.predgpi.gff3
    pigz ${basename}.predgpi.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${basename}.predgpi.gff3.gz
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MERGE PROCESSES  (publishDir; inputs are staged file lists)
// ════════════════════════════════════════════════════════════════════════════

process MERGE_PFAM {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(domtbls)

    output:
        path("pfam.csv.gz"), emit: csv

    script:
    """
    python3 ${params.scripts}/pfamtbl_to_long.py \\
        --outfile pfam.csv \\
        ${domtbls}
    pigz pfam.csv
    """

    stub:
    """
    printf 'protein_id,hmm_id,hmm_acc,hmm_len,full_seq_e_value,full_seq_score,full_seq_bias,domain_num,domain_num_of,domain_c_evalue,domain_i_evalue,domain_score,domain_bias,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to\\n' | gzip > pfam.csv.gz
    """
}

process MERGE_CAZY {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(overviews)
        path(cazymes)

    output:
        path("cazy.overview.csv.gz"),    emit: overview_csv
        path("cazy.cazymes_hmm.csv.gz"), emit: hmm_csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_cazy.py \\
        --overviews ${overviews} \\
        --cazymes   ${cazymes} \\
        --out-overview cazy.overview.csv \\
        --out-hmm      cazy.cazymes_hmm.csv
    pigz cazy.overview.csv cazy.cazymes_hmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,EC,cazyme_fam,sub_fam,diamond_fam,substrate,toolcount\\n' | gzip > cazy.overview.csv.gz
    printf 'species_prefix,HMM_id,profile_length,protein_id,protein_length,evalue,q_start,q_end,s_start,s_end,coverage\\n' | gzip > cazy.cazymes_hmm.csv.gz
    """
}

process MERGE_MEROPS {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("merops.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_merops.py -o merops.csv ${blasttabs}
    pigz merops.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,merops_id,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' | gzip > merops.csv.gz
    """
}

process MERGE_SIGNALP {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(gff3s)

    output:
        path("signalp.signal_peptide.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_signalp.py -o signalp.signal_peptide.csv ${gff3s}
    pigz signalp.signal_peptide.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,peptide_start,peptide_end,probability\\n' | gzip > signalp.signal_peptide.csv.gz
    """
}

process MERGE_TMHMM {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(tsvs)

    output:
        path("tmhmm.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_tmhmm.py -o tmhmm.csv ${tsvs}
    pigz tmhmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,len,ExpAA,First60,PredHel,Topology\\n' | gzip > tmhmm.csv.gz
    """
}

process MERGE_TARGETP {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(summaries)

    output:
        path("targetP.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_targetp.py -o targetP.csv ${summaries}
    pigz targetP.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,prediction,probability,cleavage_position_start,cleavage_position_end,cleavage_probability,motif\\n' | gzip > targetP.csv.gz
    """
}

process MERGE_IDP {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(idp_files)
        path(idp_sums)

    output:
        path("idp.csv.gz"),         emit: idp
        path("idp_summary.csv.gz"), emit: summary

    script:
    """
    # Concatenate per-species tables keeping one header line each
    first_idp=\$(ls *.idp.csv.gz | head -1)
    zcat \$first_idp | head -1 > idp.csv
    for f in *.idp.csv.gz; do zcat "\$f" | tail -n +2 >> idp.csv; done
    pigz idp.csv

    first_sum=\$(ls *.idp_summary.csv.gz | head -1)
    zcat \$first_sum | head -1 > idp_summary.csv
    for f in *.idp_summary.csv.gz; do zcat "\$f" | tail -n +2 >> idp_summary.csv; done
    pigz idp_summary.csv
    """

    stub:
    """
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > idp_summary.csv.gz
    """
}

process MERGE_WOLFPSORT {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(results)

    output:
        path("wolfpsort.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_wolfpsort.py -o wolfpsort.csv ${results}
    pigz wolfpsort.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,localization,score\\n' | gzip > wolfpsort.csv.gz
    """
}

process MERGE_PREDGPI {
    label      'merge'
    publishDir "${params.tables}", mode: 'copy'

    input:
        path(gff3s)

    output:
        path("predgpi.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_predgpi.py -o predgpi.csv ${gff3s}
    pigz predgpi.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,source,feature,start,end,score,strand,phase,attributes\\n' | gzip > predgpi.csv.gz
    """
}

// ════════════════════════════════════════════════════════════════════════════
// PER-GENOME STATISTICS PROCESSES  (storeDir-cached per LOCUSTAG)
// ════════════════════════════════════════════════════════════════════════════

// ── AA Frequency ──────────────────────────────────────────────────────────────
process CALC_AA_FREQ {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/aa_freq"

    input:
    tuple val(locustag), path(proteins_faa)

    output:
    path "${locustag}.aa_freq.csv.gz", emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_AA_freq.py \\
        ${proteins_faa} -o ${locustag}.aa_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > ${locustag}.aa_freq.csv.gz
    """
}

process MERGE_AA_FREQ {
    label      'merge'
    publishDir path: {
        params.taxon
            ? "${params.tables}/${params.taxon.split(':',2)[1].replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
            : "${params.tables}/All_Taxa"
    }, mode: 'copy'

    input:
    path 'inputs/*'

    output:
    path "aa_freq.csv.gz", emit: csv

    script:
    """
    first=1
    for f in inputs/*.aa_freq.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > aa_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > aa_freq.csv.gz
    """
}

// ── Codon Frequency ───────────────────────────────────────────────────────────
process CALC_CODON_FREQ {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/codon_freq"

    input:
    tuple val(locustag), path(cds_faa)

    output:
    path "${locustag}.codon_freq.csv.gz", emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_codon_freq.py \\
        ${cds_faa} -o ${locustag}.codon_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,codon,frequency\\n' | gzip > ${locustag}.codon_freq.csv.gz
    """
}

process MERGE_CODON_FREQ {
    label      'merge'
    publishDir path: {
        params.taxon
            ? "${params.tables}/${params.taxon.split(':',2)[1].replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
            : "${params.tables}/All_Taxa"
    }, mode: 'copy'

    input:
    path 'inputs/*'

    output:
    path "codon_freq.csv.gz", emit: csv

    script:
    """
    first=1
    for f in inputs/*.codon_freq.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > codon_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,codon,frequency\\n' | gzip > codon_freq.csv.gz
    """
}

// ── Intergenic Distances ──────────────────────────────────────────────────────
process CALC_INTERGENIC {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/intergenic_stats"

    input:
    tuple val(locustag), path(gff_file)

    output:
    path "${locustag}.gene_intergenic_distances.csv.gz", emit: csv

    script:
    """
    module load biopython
    python3 ${params.scripts}/calculate_intergenic.py \\
        ${gff_file} -o .
    pigz gene_pairwise_distances.csv
    mv gene_pairwise_distances.csv.gz ${locustag}.gene_intergenic_distances.csv.gz
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > ${locustag}.gene_intergenic_distances.csv.gz
    """
}

process MERGE_INTERGENIC {
    label      'merge'
    publishDir path: {
        params.taxon
            ? "${params.tables}/${params.taxon.split(':',2)[1].replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
            : "${params.tables}/All_Taxa"
    }, mode: 'copy'

    input:
    path 'inputs/*'

    output:
    path "gene_intergenic_distances.csv.gz", emit: csv

    script:
    """
    first=1
    for f in inputs/*.gene_intergenic_distances.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > gene_intergenic_distances.csv.gz
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > gene_intergenic_distances.csv.gz
    """
}

// ── Gene Statistics ───────────────────────────────────────────────────────────
process CALC_GENE_STATS {
    label    'genestats'
    tag      locustag
    storeDir "${params.genome_stats_outdir}/gene_stats"

    input:
    tuple val(locustag), path(gff_file), path(dna_file)

    output:
    path "${locustag}.gene_info.csv.gz",        emit: gene_info
    path "${locustag}.gene_exons.csv.gz",       emit: gene_exons
    path "${locustag}.gene_CDS.csv.gz",         emit: gene_CDS
    path "${locustag}.gene_introns.csv.gz",     emit: gene_introns
    path "${locustag}.gene_transcripts.csv.gz", emit: gene_transcripts
    path "${locustag}.gene_trnas.csv.gz",       emit: gene_trnas
    path "${locustag}.gene_proteins.csv.gz",    emit: gene_proteins

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load biopython
    python3 ${params.scripts}/build_genestats_table.py \\
        ${gff_file} \\
        -d . \\
        -o .
    pigz gene_info.csv gene_exons.csv gene_CDS.csv gene_introns.csv \\
         gene_transcripts.csv gene_trnas.csv gene_proteins.csv
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        mv \${f}.csv.gz ${locustag}.\${f}.csv.gz
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > ${locustag}.\${f}.csv.gz
    done
    """
}

process MERGE_GENE_STATS {
    label      'merge'
    publishDir path: {
        params.taxon
            ? "${params.tables}/${params.taxon.split(':',2)[1].replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
            : "${params.tables}/All_Taxa"
    }, mode: 'copy'

    input:
    path 'inputs/*'

    output:
    path "gene_info.csv.gz",        emit: gene_info
    path "gene_exons.csv.gz",       emit: gene_exons
    path "gene_CDS.csv.gz",         emit: gene_CDS
    path "gene_introns.csv.gz",     emit: gene_introns
    path "gene_transcripts.csv.gz", emit: gene_transcripts
    path "gene_trnas.csv.gz",       emit: gene_trnas
    path "gene_proteins.csv.gz",    emit: gene_proteins

    script:
    """
    for type in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        first=1
        for f in inputs/*.\${type}.csv.gz; do
            if [ "\$first" = "1" ]; then zcat "\$f"; first=0
            else zcat "\$f" | tail -n +2; fi
        done | gzip > \${type}.csv.gz
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > \${f}.csv.gz
    done
    """
}

// ════════════════════════════════════════════════════════════════════════════
// INPUT SETUP (symlinks from genome_annotation → input/)
// ════════════════════════════════════════════════════════════════════════════

workflow SETUP_INPUT {
    take: ch   // tuple(locustag, basename, species, strain)
    main:
        // Write all locustag+basename pairs to a TSV file for the bash loop
        rows_file = ch
            .map { locustag, basename, species, strain -> "${locustag}\t${basename}" }
            .collectFile(name: 'setup_rows.tsv', newLine: true)
        SETUP_SYMLINKS(rows_file)
        // Re-emit original rows gated by process completion
        done_ch = SETUP_SYMLINKS.out.done
            .combine(ch)
            .map { _flag, locustag, basename, species, strain ->
                tuple(locustag, basename, species, strain)
            }
    emit:
        done = done_ch
}

process SETUP_SYMLINKS {
    label 'setup'

    input:
        path(rows_file)

    output:
        val(true), emit: done

    script:
    """
    mkdir -p "${params.pep_dir}" "${params.cds_dir}" "${params.gff_dir}" \\
             "${params.genome_dir}" "${params.trna_dir}"

    make_link() {
        local target=\$1 linkname=\$2
        if [ ! -e "\$target" ]; then
            echo "[WARN] source not found, skipping: \$target" >&2
            return 0
        fi
        if [[ ! -L "\$linkname" || ! -e "\$linkname" ]]; then
            ln -sfn "\$target" "\$linkname"
            echo "[INFO] linked \$linkname -> \$target"
        else
            echo "[INFO] symlink already valid, skipping: \$linkname"
        fi
    }

    while IFS=\$'\\t' read -r locustag basename; do
        src="${params.genome_annotation}/\${basename}/predict_results"
        misc="${params.genome_annotation}/\${basename}/predict_misc"

        if [ ! -d "\$src" ]; then
            echo "[WARN] predict_results not found for \${basename}: \$src" >&2
            continue
        fi

        make_link "\$src/\${basename}.proteins.fa"        "${params.pep_dir}/\${basename}.proteins.fa"
        make_link "\$src/\${basename}.cds-transcripts.fa" "${params.cds_dir}/\${basename}.cds-transcripts.fa"
        make_link "\$src/\${basename}.gff3"               "${params.gff_dir}/\${basename}.gff3"
        make_link "\$src/\${basename}.scaffolds.fa"       "${params.genome_dir}/\${basename}.scaffolds.fa"

        if [ -f "\$misc/trnascan.no-overlaps.gff3" ]; then
            make_link "\$misc/trnascan.no-overlaps.gff3" "${params.trna_dir}/\${basename}.trna.gff3"
        else
            echo "[INFO] no trnascan GFF3 for \${basename}, skipping trna symlink"
        fi
    done < ${rows_file}
    """

    stub:
    """
    echo "[STUB] SETUP_SYMLINKS: \$(wc -l < ${rows_file}) species"
    """
}

// ════════════════════════════════════════════════════════════════════════════
// HELPERS
// ════════════════════════════════════════════════════════════════════════════

// Build a gated glob channel.
//   sync_ch — emits one value when the RUN step is done (or Channel.of(true))
//   glob    — shell-style glob relative to params.outdir
// Returns a channel of matching Path objects, or empty if none found.
def gatedGlob(sync_ch, String glob) {
    sync_ch
        .flatMap { files("${params.outdir}/${glob}") }
        .filter  { it.size() > 0 }
        .collect()
        .filter  { !it.isEmpty() }
}

// Like gatedGlob but roots the glob in params.genome_stats_outdir instead.
// Used for merge_all=true when no --taxon filter is active.
def gatedGlobStats(sync_ch, String glob) {
    sync_ch
        .flatMap { files("${params.genome_stats_outdir}/${glob}") }
        .filter  { it.size() > 0 }
        .collect()
        .filter  { !it.isEmpty() }
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {
    // ── Taxonomy filter ────────────────────────────────────────────────────────
    // Parse --taxon RANK:VALUE (e.g. --taxon PHYLUM:Ascomycota).
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

    // ── Base sample channel ────────────────────────────────────────────────────
    // Emits tuple(locustag, basename, species, strain) per row.
    // STRAIN: first ';'-delimited token; single quotes stripped.
    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = (row.STRAIN?.trim() ?: '').split(';')[0].trim().replace("'", '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = [species, strain].findAll { it }.join('_').replaceAll(/[\s\/\#]+/, '_')
            tuple(locustag, basename, species, strain)
        }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // ── Input setup: symlink predict_results files into input/ subdirs ─────────
    def ready_ch
    if (params.run_setup.toBoolean()) {
        SETUP_INPUT(rows_ch)
        ready_ch = SETUP_INPUT.out.done
    } else {
        ready_ch = rows_ch
    }

    // ── Protein channel ────────────────────────────────────────────────────────
    proteins_ch = ready_ch
        .map { locustag, basename, species, strain ->
            def prot = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
            if (!prot.exists()) {
                log.warn "Skipping ${basename} (${locustag}): protein file not found"
                return null
            }
            return tuple(locustag, basename, species, strain, prot)
        }
        .filter { it != null }

    // ── Per-genome statistics input channels ───────────────────────────────────
    aa_freq_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
        f.exists() ? tuple(locustag, f) : null
    }.filter { it != null }

    codon_freq_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.cds_dir}/${basename}.cds-transcripts.fa", glob: false)
        f.exists() ? tuple(locustag, f) : null
    }.filter { it != null }

    intergenic_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.gff_dir}/${basename}.gff3", glob: false)
        f.exists() ? tuple(locustag, f) : null
    }.filter { it != null }

    gene_stats_ch = ready_ch.map { locustag, basename, species, strain ->
        def gff = file("${params.gff_dir}/${basename}.gff3", glob: false)
        def dna = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
        (gff.exists() && dna.exists()) ? tuple(locustag, gff, dna) : null
    }.filter { it != null }

    // ── Per-species RUN steps ──────────────────────────────────────────────────
    // storeDir means Nextflow skips a species automatically if all its output
    // files already exist in the store directory — no -resume needed.
    if (params.run_pfam.toBoolean())      RUN_PFAM(proteins_ch)
    if (params.run_cazy.toBoolean())      RUN_CAZY(proteins_ch)
    if (params.run_merops.toBoolean())    RUN_MEROPS(proteins_ch)
    if (params.run_signalp.toBoolean())   RUN_SIGNALP(proteins_ch)
    if (params.run_tmhmm.toBoolean())     RUN_TMHMM(proteins_ch)
    if (params.run_targetp.toBoolean())   RUN_TARGETP(proteins_ch)
    if (params.run_idp.toBoolean())       RUN_IDP(proteins_ch)
    if (params.run_wolfpsort.toBoolean()) RUN_WOLFPSORT(proteins_ch)
    if (params.run_predgpi.toBoolean())   RUN_PREDGPI(proteins_ch)

    // ── MERGE steps ────────────────────────────────────────────────────────────
    //
    // merge_all=true  (default): collect ALL result files from the output dirs,
    //                            including those from previous runs.  A sync
    //                            barrier (collect on the current run's outputs)
    //                            ensures newly generated files are on disk before
    //                            globbing.  If a tool was not run this session the
    //                            barrier is an immediate Channel.of(true).
    //                            MERGE is skipped when the glob finds no files.
    //
    // merge_all=false:           merge only the files produced in this run.
    //                            MERGE is skipped when the tool was not run.
    //
    // skip_merge=true:           skip all MERGE steps unconditionally.

    if (!params.skip_merge.toBoolean()) {

        if (params.merge_all.toBoolean()) {

            if (params.run_pfam.toBoolean()) {
                def sync = RUN_PFAM.out.domtbl.collect()
                MERGE_PFAM( gatedGlob(sync, "pfam_hmmscan/*.pfam.gz") )
            } else {
                MERGE_PFAM( gatedGlob(Channel.of(true), "pfam_hmmscan/*.pfam.gz") )
            }

            if (params.run_cazy.toBoolean()) {
                def ov_sync = RUN_CAZY.out.overview.collect()
                def ca_sync = RUN_CAZY.out.cazymes.collect()
                MERGE_CAZY(
                    gatedGlob(ov_sync, "cazy/*/*.overview.tsv.gz"),
                    gatedGlob(ca_sync, "cazy/*/*.cazymes.tsv.gz")
                )
            } else {
                MERGE_CAZY(
                    gatedGlob(Channel.of(true), "cazy/*/*.overview.tsv.gz"),
                    gatedGlob(Channel.of(true), "cazy/*/*.cazymes.tsv.gz")
                )
            }

            if (params.run_merops.toBoolean()) {
                def sync = RUN_MEROPS.out.blasttab.collect()
                MERGE_MEROPS( gatedGlob(sync, "merops/*.blasttab.gz") )
            } else {
                MERGE_MEROPS( gatedGlob(Channel.of(true), "merops/*.blasttab.gz") )
            }

            if (params.run_signalp.toBoolean()) {
                def sync = RUN_SIGNALP.out.gff3.collect()
                MERGE_SIGNALP( gatedGlob(sync, "signalp/*.signalp.gff3.gz") )
            } else {
                MERGE_SIGNALP( gatedGlob(Channel.of(true), "signalp/*.signalp.gff3.gz") )
            }

            if (params.run_tmhmm.toBoolean()) {
                def sync = RUN_TMHMM.out.short_tsv.collect()
                MERGE_TMHMM( gatedGlob(sync, "tmhmm/*.tmhmm_short.tsv.gz") )
            } else {
                MERGE_TMHMM( gatedGlob(Channel.of(true), "tmhmm/*.tmhmm_short.tsv.gz") )
            }

            if (params.run_targetp.toBoolean()) {
                def sync = RUN_TARGETP.out.summary.collect()
                MERGE_TARGETP( gatedGlob(sync, "targetP/*_summary.targetp2.gz") )
            } else {
                MERGE_TARGETP( gatedGlob(Channel.of(true), "targetP/*_summary.targetp2.gz") )
            }

            if (params.run_idp.toBoolean()) {
                def idp_sync = RUN_IDP.out.idp_csv.collect()
                def sum_sync = RUN_IDP.out.idp_summary_csv.collect()
                MERGE_IDP(
                    gatedGlob(idp_sync, "aiupred/*.idp.csv.gz"),
                    gatedGlob(sum_sync, "aiupred/*.idp_summary.csv.gz")
                )
            } else {
                MERGE_IDP(
                    gatedGlob(Channel.of(true), "aiupred/*.idp.csv.gz"),
                    gatedGlob(Channel.of(true), "aiupred/*.idp_summary.csv.gz")
                )
            }

            if (params.run_wolfpsort.toBoolean()) {
                def sync = RUN_WOLFPSORT.out.results.collect()
                MERGE_WOLFPSORT( gatedGlob(sync, "wolfpsort/*.wolfpsort.results.txt.gz") )
            } else {
                MERGE_WOLFPSORT( gatedGlob(Channel.of(true), "wolfpsort/*.wolfpsort.results.txt.gz") )
            }

            if (params.run_predgpi.toBoolean()) {
                def sync = RUN_PREDGPI.out.gff3.collect()
                MERGE_PREDGPI( gatedGlob(sync, "predgpi/*.predgpi.gff3.gz") )
            } else {
                MERGE_PREDGPI( gatedGlob(Channel.of(true), "predgpi/*.predgpi.gff3.gz") )
            }

        } else {
            // merge_all=false: merge only the outputs produced in this run.
            if (params.run_pfam.toBoolean())
                MERGE_PFAM(RUN_PFAM.out.domtbl.collect())
            if (params.run_cazy.toBoolean())
                MERGE_CAZY(RUN_CAZY.out.overview.collect(), RUN_CAZY.out.cazymes.collect())
            if (params.run_merops.toBoolean())
                MERGE_MEROPS(RUN_MEROPS.out.blasttab.collect())
            if (params.run_signalp.toBoolean())
                MERGE_SIGNALP(RUN_SIGNALP.out.gff3.collect())
            if (params.run_tmhmm.toBoolean())
                MERGE_TMHMM(RUN_TMHMM.out.short_tsv.collect())
            if (params.run_targetp.toBoolean())
                MERGE_TARGETP(RUN_TARGETP.out.summary.collect())
            if (params.run_idp.toBoolean())
                MERGE_IDP(RUN_IDP.out.idp_csv.collect(), RUN_IDP.out.idp_summary_csv.collect())
            if (params.run_wolfpsort.toBoolean())
                MERGE_WOLFPSORT(RUN_WOLFPSORT.out.results.collect())
            if (params.run_predgpi.toBoolean())
                MERGE_PREDGPI(RUN_PREDGPI.out.gff3.collect())
        }
    }

    // ── Per-genome statistics + MERGE ──────────────────────────────────────────
    // When merge_all=true and no --taxon is active, glob ALL files from storeDir
    // so that species not in this run's samples.csv are still included.
    // When --taxon is set, always use only current-run outputs (already filtered).
    def use_glob = params.merge_all.toBoolean() && !params.taxon

    if (params.run_aa_freq.toBoolean())    CALC_AA_FREQ(aa_freq_ch)
    if (params.run_codon_freq.toBoolean()) CALC_CODON_FREQ(codon_freq_ch)
    if (params.run_intergenic.toBoolean()) CALC_INTERGENIC(intergenic_ch)
    if (params.run_gene_stats.toBoolean()) CALC_GENE_STATS(gene_stats_ch)

    if (!params.skip_merge.toBoolean()) {
        if (use_glob) {
            if (params.run_aa_freq.toBoolean()) {
                MERGE_AA_FREQ(gatedGlobStats(CALC_AA_FREQ.out.csv.collect(), "aa_freq/*.aa_freq.csv.gz"))
            } else {
                MERGE_AA_FREQ(gatedGlobStats(Channel.of(true), "aa_freq/*.aa_freq.csv.gz"))
            }
            if (params.run_codon_freq.toBoolean()) {
                MERGE_CODON_FREQ(gatedGlobStats(CALC_CODON_FREQ.out.csv.collect(), "codon_freq/*.codon_freq.csv.gz"))
            } else {
                MERGE_CODON_FREQ(gatedGlobStats(Channel.of(true), "codon_freq/*.codon_freq.csv.gz"))
            }
            if (params.run_intergenic.toBoolean()) {
                MERGE_INTERGENIC(gatedGlobStats(CALC_INTERGENIC.out.csv.collect(), "intergenic_stats/*.gene_intergenic_distances.csv.gz"))
            } else {
                MERGE_INTERGENIC(gatedGlobStats(Channel.of(true), "intergenic_stats/*.gene_intergenic_distances.csv.gz"))
            }
            if (params.run_gene_stats.toBoolean()) {
                MERGE_GENE_STATS(gatedGlobStats(
                    CALC_GENE_STATS.out.gene_info
                        .mix(CALC_GENE_STATS.out.gene_exons)
                        .mix(CALC_GENE_STATS.out.gene_CDS)
                        .mix(CALC_GENE_STATS.out.gene_introns)
                        .mix(CALC_GENE_STATS.out.gene_transcripts)
                        .mix(CALC_GENE_STATS.out.gene_trnas)
                        .mix(CALC_GENE_STATS.out.gene_proteins)
                        .collect(),
                    "gene_stats/*.csv.gz"
                ))
            } else {
                MERGE_GENE_STATS(gatedGlobStats(Channel.of(true), "gene_stats/*.csv.gz"))
            }
        } else {
            // current-run outputs only (merge_all=false, or --taxon active)
            if (params.run_aa_freq.toBoolean())    MERGE_AA_FREQ(CALC_AA_FREQ.out.csv.collect())
            if (params.run_codon_freq.toBoolean()) MERGE_CODON_FREQ(CALC_CODON_FREQ.out.csv.collect())
            if (params.run_intergenic.toBoolean()) MERGE_INTERGENIC(CALC_INTERGENIC.out.csv.collect())
            if (params.run_gene_stats.toBoolean()) {
                MERGE_GENE_STATS(
                    CALC_GENE_STATS.out.gene_info
                        .mix(CALC_GENE_STATS.out.gene_exons)
                        .mix(CALC_GENE_STATS.out.gene_CDS)
                        .mix(CALC_GENE_STATS.out.gene_introns)
                        .mix(CALC_GENE_STATS.out.gene_transcripts)
                        .mix(CALC_GENE_STATS.out.gene_trnas)
                        .mix(CALC_GENE_STATS.out.gene_proteins)
                        .collect()
                )
            }
        }
    }
}
