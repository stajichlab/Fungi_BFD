#!/usr/bin/env nextflow
// DSL2 is the default from Nextflow 22+ — no explicit declaration needed

/*
 * Fungi 5K Functional Annotation Pipeline
 * Runs Pfam, CAZy, MEROPS, SignalP, TMHMM, TargetP, IDP (AIUPred), WoLF PSORT,
 * and predGPI on all species in samples.csv, then consolidates each tool's results
 * into a DuckDB-loadable <tool>.csv.gz in tables/.
 *
 * Usage (from project root):
 *   sbatch pipeline/nextflow/run_functional.sh
 *   nextflow run pipeline/nextflow/genome_functional.nf -c pipeline/nextflow/nextflow.config -resume
 *
 * Dry-run / testing:
 *   nextflow run pipeline/nextflow/genome_functional.nf -stub-run --n_test 2
 */

// All params (samples, pep_dir, outdir, tables, scripts, run_*, n_test) are
// defined in nextflow.config — do not redeclare defaults here.

// ════════════════════════════════════════════════════════════════════════════
// SUBWORKFLOWS
// ════════════════════════════════════════════════════════════════════════════

// ── Pfam ─────────────────────────────────────────────────────────────────────
workflow PFAM {
    take: ch
    main:
        RUN_PFAM(ch)
        MERGE_PFAM(RUN_PFAM.out.domtbl.collect())
    emit:
        merged = MERGE_PFAM.out.csv
}

process RUN_PFAM {
    tag        "${locustag}"
    label      'pfam'
    publishDir "${params.outdir}/pfam_hmmscan", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.pfam.gz"),    emit: domtbl
        path("${locustag}.tblout.gz"),  emit: tblout

    script:
    """
    hmmscan --cut_ga --cpu ${task.cpus} \\
        --domtblout ${locustag}.pfam \\
        --tblout    ${locustag}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${locustag}.pfam ${locustag}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${locustag}.pfam.gz
    printf '' | gzip     > ${locustag}.tblout.gz
    """
}

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

// ── CAZy ─────────────────────────────────────────────────────────────────────
workflow CAZY {
    take: ch
    main:
        RUN_CAZY(ch)
        MERGE_CAZY(
            RUN_CAZY.out.overview.collect(),
            RUN_CAZY.out.cazymes.collect()
        )
    emit:
        merged_overview = MERGE_CAZY.out.overview_csv
        merged_hmm      = MERGE_CAZY.out.hmm_csv
}

process RUN_CAZY {
    tag        "${locustag}"
    label      'cazy'
    publishDir { "${params.outdir}/cazy/${locustag}" }, mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.overview.tsv.gz"),   emit: overview
        path("${locustag}.cazymes.tsv.gz"),    emit: cazymes
        path("${locustag}.substrates.tsv.gz"), emit: substrates

    script:
    """
    mkdir -p ${locustag}
    dbcanlight search -i ${proteins} -m cazyme -o ${locustag} -t ${task.cpus}
    dbcanlight search -i ${proteins} -m sub    -o ${locustag} -t ${task.cpus}
    dbcanlight conclude ${locustag}
    pigz -f ${locustag}/cazymes.tsv ${locustag}/substrates.tsv ${locustag}/overview.tsv
    mv ${locustag}/overview.tsv.gz   ${locustag}.overview.tsv.gz
    mv ${locustag}/cazymes.tsv.gz    ${locustag}.cazymes.tsv.gz
    mv ${locustag}/substrates.tsv.gz ${locustag}.substrates.tsv.gz
    """

    stub:
    """
    printf 'Gene_ID\\tEC\\tcazyme_fam\\tsub_fam\\tdiamond_fam\\tSubstrate\\t#ofTools\\n' | gzip > ${locustag}.overview.tsv.gz
    printf 'HMM_Profile\\tProfile_Length\\tGene_ID\\tGene_Length\\tEvalue\\tProfile_Start\\tProfile_End\\tGene_Start\\tGene_End\\tCoverage\\n' | gzip > ${locustag}.cazymes.tsv.gz
    printf 'Gene_ID\\tSubstrate\\n' | gzip > ${locustag}.substrates.tsv.gz
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

// ── MEROPS ────────────────────────────────────────────────────────────────────
workflow MEROPS {
    take: ch
    main:
        RUN_MEROPS(ch)
        MERGE_MEROPS(RUN_MEROPS.out.blasttab.collect())
    emit:
        merged = MERGE_MEROPS.out.csv
}

process RUN_MEROPS {
    tag        "${locustag}"
    label      'merops'
    publishDir "${params.outdir}/merops", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.blasttab.gz"), emit: blasttab

    script:
    """
    blastp -query ${proteins} \\
        -db \$MEROPS_DB/merops_scan.lib \\
        -out ${locustag}.blasttab \\
        -num_threads ${task.cpus} \\
        -seg yes -soft_masking true \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -outfmt 6 \\
        -use_sw_tback
    pigz ${locustag}.blasttab
    """

    stub:
    """
    printf '' | gzip > ${locustag}.blasttab.gz
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
    merge_merops.py -o merops.csv ${blasttabs}
    pigz merops.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,merops_id,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' | gzip > merops.csv.gz
    """
}

// ── SignalP ───────────────────────────────────────────────────────────────────
workflow SIGNALP {
    take: ch
    main:
        RUN_SIGNALP(ch)
        MERGE_SIGNALP(RUN_SIGNALP.out.gff3.collect())
    emit:
        merged = MERGE_SIGNALP.out.csv
}

process RUN_SIGNALP {
    tag        "${locustag}"
    label      'signalp'
    publishDir "${params.outdir}/signalp", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.signalp.gff3.gz"),         emit: gff3
        path("${locustag}.signalp.results.txt.gz"),  emit: results

    script:
    """
    OUTD=\$(mktemp -d)
    signalp6 -od \$OUTD -org euk --mode fast -format txt \\
        -fasta ${proteins} --write_procs ${task.cpus} -bs 100
    pigz -c \$OUTD/output.gff3             > ${locustag}.signalp.gff3.gz
    pigz -c \$OUTD/prediction_results.txt  > ${locustag}.signalp.results.txt.gz
    rm -rf \$OUTD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${locustag}.signalp.gff3.gz
    printf '# SignalP-6.0\\n'   | gzip > ${locustag}.signalp.results.txt.gz
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
    merge_signalp.py -o signalp.signal_peptide.csv ${gff3s}
    pigz signalp.signal_peptide.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,peptide_start,peptide_end,probability\\n' | gzip > signalp.signal_peptide.csv.gz
    """
}

// ── TMHMM ─────────────────────────────────────────────────────────────────────
workflow TMHMM {
    take: ch
    main:
        RUN_TMHMM(ch)
        MERGE_TMHMM(RUN_TMHMM.out.short_tsv.collect())
    emit:
        merged = MERGE_TMHMM.out.csv
}

process RUN_TMHMM {
    tag        "${locustag}"
    label      'tmhmm'
    publishDir "${params.outdir}/tmhmm", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.tmhmm_short.tsv.gz"),   emit: short_tsv
        path("${locustag}.tmhmm_results.tsv.gz"), emit: full_tsv

    script:
    """
    tmhmm --noplot         < ${proteins} > ${locustag}.tmhmm_results.tsv
    tmhmm --short --noplot < ${proteins} > ${locustag}.tmhmm_short.tsv
    pigz ${locustag}.tmhmm_results.tsv ${locustag}.tmhmm_short.tsv
    """

    stub:
    """
    printf '# TMHMM\\n' | gzip > ${locustag}.tmhmm_results.tsv.gz
    printf '# TMHMM\\n' | gzip > ${locustag}.tmhmm_short.tsv.gz
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
    merge_tmhmm.py -o tmhmm.csv ${tsvs}
    pigz tmhmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,len,ExpAA,First60,PredHel,Topology\\n' | gzip > tmhmm.csv.gz
    """
}

// ── TargetP ───────────────────────────────────────────────────────────────────
workflow TARGETP {
    take: ch
    main:
        RUN_TARGETP(ch)
        MERGE_TARGETP(RUN_TARGETP.out.summary.collect())
    emit:
        merged = MERGE_TARGETP.out.csv
}

process RUN_TARGETP {
    tag        "${locustag}"
    label      'targetp'
    publishDir "${params.outdir}/targetP", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}_summary.targetp2.gz"), emit: summary

    script:
    """
    TMPD=\$(mktemp -d)
    targetp -batch 50 -tmp \$TMPD -format short \\
        -fasta ${proteins} -org non-pl -prefix ${locustag}
    pigz -f ${locustag}_summary.targetp2
    rm -rf \$TMPD
    """

    stub:
    """
    printf '# TargetP-2.0\\n' | gzip > ${locustag}_summary.targetp2.gz
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
    merge_targetp.py -o targetP.csv ${summaries}
    pigz targetP.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,prediction,probability,cleavage_position_start,cleavage_position_end,cleavage_probability,motif\\n' | gzip > targetP.csv.gz
    """
}

// ── IDP (AIUPred) ─────────────────────────────────────────────────────────────
workflow IDP {
    take: ch
    main:
        RUN_IDP(ch)
        MERGE_IDP(
            RUN_IDP.out.idp_csv.collect(),
            RUN_IDP.out.idp_summary_csv.collect()
        )
    emit:
        merged_idp     = MERGE_IDP.out.idp
        merged_summary = MERGE_IDP.out.summary
}

process RUN_IDP {
    tag        "${locustag}"
    label      'idp'
    publishDir "${params.outdir}/aiupred", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.aiupred.txt.gz"),     emit: raw
        path("${locustag}.idp.csv.gz"),         emit: idp_csv
        path("${locustag}.idp_summary.csv.gz"), emit: idp_summary_csv

    script:
    """
    aiupred.py -i ${proteins} -o ${locustag}.aiupred.txt
    pigz ${locustag}.aiupred.txt
    python3 ${params.scripts}/gather_AIUPred.py ${locustag}.aiupred.txt.gz \\
        --outfile      ${locustag}.idp.csv \\
        --outfilesum   ${locustag}.idp_summary.csv
    pigz ${locustag}.idp.csv ${locustag}.idp_summary.csv
    """

    stub:
    """
    printf '' | gzip > ${locustag}.aiupred.txt.gz
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > ${locustag}.idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > ${locustag}.idp_summary.csv.gz
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

// ── WoLF PSORT ────────────────────────────────────────────────────────────────
workflow WOLFPSORT {
    take: ch
    main:
        RUN_WOLFPSORT(ch)
        MERGE_WOLFPSORT(RUN_WOLFPSORT.out.results.collect())
    emit:
        merged = MERGE_WOLFPSORT.out.csv
}

process RUN_WOLFPSORT {
    tag        "${locustag}"
    label      'wolfpsort'
    publishDir "${params.outdir}/wolfpsort", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.wolfpsort.results.txt.gz"), emit: results

    script:
    """
    cat ${proteins} | runWolfPsortSummary fungi > ${locustag}.wolfpsort.results.txt
    pigz ${locustag}.wolfpsort.results.txt
    """

    stub:
    """
    printf '# WoLF PSORT\\n' | gzip > ${locustag}.wolfpsort.results.txt.gz
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
    merge_wolfpsort.py -o wolfpsort.csv ${results}
    pigz wolfpsort.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,localization,score\\n' | gzip > wolfpsort.csv.gz
    """
}

// ── predGPI ───────────────────────────────────────────────────────────────────
workflow PREDGPI {
    take: ch
    main:
        RUN_PREDGPI(ch)
        MERGE_PREDGPI(RUN_PREDGPI.out.gff3.collect())
    emit:
        merged = MERGE_PREDGPI.out.csv
}

process RUN_PREDGPI {
    tag        "${locustag}"
    label      'predgpi'
    publishDir "${params.outdir}/predgpi", mode: 'copy'

    input:
        tuple val(locustag), val(label), val(species), val(strain), path(proteins)

    output:
        path("${locustag}.predgpi.gff3.gz"), emit: gff3

    script:
    """
    predgpi.py -f ${proteins} -m gff3 -o ${locustag}.predgpi.gff3
    pigz ${locustag}.predgpi.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${locustag}.predgpi.gff3.gz
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
    merge_predgpi.py -o predgpi.csv ${gff3s}
    pigz predgpi.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,source,feature,start,end,score,strand,phase,attributes\\n' | gzip > predgpi.csv.gz
    """
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {
    // ── input channel ──────────────────────────────────────────────────────────
    // Filename convention (matching 1KFG): {SPECIES}_{STRAIN}.proteins.fa
    //   STRAIN: first ';'-delimited token; single quotes and extra whitespace stripped
    proteins_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = (row.STRAIN?.trim() ?: '').split(';')[0].trim().replace("'", '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def label    = [species, strain].findAll { it }.join('_').replaceAll(/\s+/, '_')
            def prot     = file("${params.pep_dir}/${label}.proteins.fa", glob: false)
            if (!prot.exists()) {
                log.warn "Skipping ${label} (${locustag}): protein file not found"
                return null
            }
            return tuple(locustag, label, species, strain, prot)
        }
        .filter { it != null }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    if (params.run_pfam)      PFAM(proteins_ch)
    if (params.run_cazy)      CAZY(proteins_ch)
    if (params.run_merops)    MEROPS(proteins_ch)
    if (params.run_signalp)   SIGNALP(proteins_ch)
    if (params.run_tmhmm)     TMHMM(proteins_ch)
    if (params.run_targetp)   TARGETP(proteins_ch)
    if (params.run_idp)       IDP(proteins_ch)
    if (params.run_wolfpsort) WOLFPSORT(proteins_ch)
    if (params.run_predgpi)   PREDGPI(proteins_ch)
}
