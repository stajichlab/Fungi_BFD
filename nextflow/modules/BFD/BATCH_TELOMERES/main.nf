include { hashBucketForType; globAnyOf } from '../../common/utils.nf'

// Telomere detection for a batch of genomes in one task. Per genome it runs
//   1. find_telomeres.py          known-motif search (both orientations)
//   2. denovo_telomere_repeats.py motif-free terminal tandem-repeat search
// and writes three files into the genome's hash bucket:
//   <asmid>.telomeres.tsv.gz              motif tracts (schema unchanged)
//   <asmid>.telomere_denovo.tsv           de novo unit summary + status row
//   <asmid>.telomere_denovo_ends.tsv.gz   one row per end with a tandem tract
//
// Batched like BATCH_AA_FREQ: a task runs ~telomere_batch_size genomes to
// amortise SLURM start-up. The batch's ASMIDs are named in the output patterns
// so storeDir does not treat the task as cached because of another genome's
// files (see the 2026-10-30 note in BATCH_AA_FREQ). Output goes to a new
// directory (telomeres_v2) so results from the old FIND_TELOMERES are not reused.
process BATCH_TELOMERES {
    label      'telomeres'
    storeDir   "${params.genome_stats_outdir}/telomeres_v2"

    input:
    tuple val(metas), path(genomes)

    output:
    path "*/${globAnyOf(metas.collect { m -> m.asmid })}.telomeres.tsv.gz",            emit: tsv
    path "*/${globAnyOf(metas.collect { m -> m.asmid })}.telomere_denovo.tsv",          emit: denovo
    path "*/${globAnyOf(metas.collect { m -> m.asmid })}.telomere_denovo_ends.tsv.gz",  emit: denovo_ends

    script:
    def metaList   = metas   instanceof List ? metas   : [metas]
    def genomeList = genomes instanceof List ? genomes : [genomes]
    def patterns_arg = params.telomere_patterns ? "-p ${params.telomere_patterns.join(' ')}" : ''
    def fuzzy_arg    = params.telomere_fuzzy.toBoolean() ? '--fuzzy' : ''
    def both_ends    = params.telomere_both_ends.toBoolean() ? '--both-ends' : ''
    def allow_internal = params.telomere_allow_internal.toBoolean() ? '--allow-internal' : ''
    def cmds = [metaList, genomeList].transpose().collect { meta, genome ->
        def asmid  = meta.asmid
        def bucket = hashBucketForType('telomeres_v2', asmid as String)
        """
        mkdir -p ${bucket}
        python3 ${projectDir}/bin/find_telomeres.py \\
            ${patterns_arg} \\
            -n ${params.telomere_min_repeats} \\
            -l ${params.telomere_min_length} \\
            -w ${params.telomere_flank_window} \\
            -s ${params.telomere_search_window} \\
            ${fuzzy_arg} \\
            --max-mismatch ${params.telomere_max_mismatch} \\
            --max-indel ${params.telomere_max_indel} \\
            --terminal-tolerance ${params.telomere_terminal_tolerance} \\
            ${both_ends} \\
            ${allow_internal} \\
            ${genome.name} \\
            -o ${bucket}/${asmid}.telomeres.tsv
        gzip -n -f ${bucket}/${asmid}.telomeres.tsv
        python3 ${projectDir}/bin/denovo_telomere_repeats.py ${genome.name} \\
            --asmid ${asmid} \\
            --window ${params.telomere_denovo_window} \\
            --tolerance ${params.telomere_denovo_tolerance} \\
            --purity ${params.telomere_denovo_purity} \\
            --min-length ${params.telomere_denovo_min_length} \\
            --background ${params.telomere_denovo_background} \\
            --min-ends ${params.telomere_denovo_min_ends} \\
            ${params.telomere_denovo_args} \\
            --summary-out ${bucket}/${asmid}.telomere_denovo.tsv \\
            --ends-out ${bucket}/${asmid}.telomere_denovo_ends.tsv.gz
        """.stripIndent()
    }.join('\n')
    """
    module load biopython
    ${cmds}
    """

    stub:
    def metaList = metas instanceof List ? metas : [metas]
    """
    ${ metaList.collect { meta ->
        def asmid  = meta.asmid
        def bucket = hashBucketForType('telomeres_v2', asmid as String)
        """
        mkdir -p ${bucket}
        printf 'scaffold\\tend\\tstrand\\tmonomer\\trepeat_count\\ttract_length\\tstart\\tend_coord\\tterminal\\tdistance_to_end\\ttract_seq\\tflank_seq\\nscaffold_1\\t5prime\\t+\\tTTAGGG\\t12\\t72\\t0\\t72\\tTrue\\t0\\tttagggttagggttaggg\\tACGT\\n' | gzip -n > ${bucket}/${asmid}.telomeres.tsv.gz
        printf 'ASMID\\tunit\\tperiod\\tn_ends\\tn_ends_examined\\tfrac_ends\\tmedian_tract_len\\tbg_windows\\tbg_hits\\tenrichment\\tstrand_biased_ends\\tstatus\\n${asmid}\\t.\\t0\\t0\\t2\\t0\\t0\\t0\\t0\\tNA\\t0\\tno_terminal_repeat\\n' > ${bucket}/${asmid}.telomere_denovo.tsv
        printf 'ASMID\\tscaffold\\tscaffold_len\\tend\\tunit_canonical\\tunit_observed\\tperiod\\tstart\\tend_coord\\ttract_length\\tpurity\\tca_frac_100\\tmethod\\n' | gzip -n > ${bucket}/${asmid}.telomere_denovo_ends.tsv.gz
        """.stripIndent()
    }.join('\n') }
    """
}
