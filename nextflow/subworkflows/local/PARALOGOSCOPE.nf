include { WGD_DMD } from '../../modules/comparative/wgd/WGD_DMD/main.nf'
include { WGD_DMD_BUNDLE } from '../../modules/comparative/wgd/WGD_DMD_BUNDLE/main.nf'
include { WGD_KSD } from '../../modules/comparative/wgd/WGD_KSD/main.nf'
include { WGD_PAIR_KS } from '../../modules/comparative/wgd/WGD_PAIR_KS/main.nf'
include { WGD_SYN } from '../../modules/comparative/wgd/WGD_SYN/main.nf'
include { MERGE_WGD_KSD } from '../../modules/comparative/wgd/MERGE_WGD_KSD/main.nf'

// Expected CPU-seconds of the pair stage for one genome, from its wgd dmd families file.
// Constants are measured: about 9 ms per pair for the first pass (alignment + identity) and about 0.15 s per
// codeml run (Boletaceae: 1.88 M pairs, 6 k codeml runs, 17.7 k CPU-s at K = 5; Glomus: 119 codeml runs/s on 32 CPUs).
// With node sampling the codeml runs are at most 5 per node; without it, about 65% of the pairs.
def estimateCost(fam) {
    long pairs = 0
    long nodes = 0
    boolean first = true
    fam.eachLine { line ->
        if (first) { first = false; return }
        def i = line.indexOf('\t')
        if (i < 0) return
        long n = line.substring(i + 1).count(',') + 1
        pairs += n * (n - 1) / 2
        nodes += n - 1
    }
    def ns = (params.wgd_speedups.toString().toBoolean() ? (params.wgd_node_sample as int) : 0)
    double codeml = ns > 0 ? Math.min(pairs * 0.65, ns * (double) nodes) : pairs * 0.65
    return 0.009 * pairs + 0.15 * codeml + 60.0
}

// First-fit-decreasing: genomes (tuple(meta, families, cds)) into bundles of at most cap expected CPU-seconds.
// A genome above the cap gets its own bundle.
def packBundles(List items, double cap) {
    def costed = items.collect { it -> [it, estimateCost(it[1])] }.sort { a, b -> b[1] <=> a[1] }
    def bins = []
    costed.each { c ->
        def bin = bins.find { b -> b[0] + c[1] <= cap }
        if (bin) { bin[0] += c[1]; bin[1] << c[0] } else { bins << [c[1], [c[0]]] }
    }
    return bins.collect { it[1] }
}

workflow PARALOGOSCOPE {
    take:
    species_ch   // channel: tuple(meta, cds) per species
    gff_ch       // channel: tuple(meta, gff) per species (may be empty)

    main:
    // wgd dmd: per genome (wgd_dmd_bundle = 0) or bundled (wgd_dmd_bundle genomes per task).
    // Both give tuple(meta, cds, families) per genome.
    def dmd_fams
    if ((params.wgd_dmd_bundle as int) > 0) {
        WGD_DMD_BUNDLE(
            species_ch
                .buffer(size: (params.wgd_dmd_bundle as int), remainder: true)
                .map { lst -> tuple(lst.collect { it[0] }, lst.collect { it[1] }) }
        )
        // one tuple per genome again; genomes whose families file is missing (failed) are dropped
        dmd_fams = WGD_DMD_BUNDLE.out.families
            .flatMap { metas, fams ->
                def fam_by = (fams instanceof List ? fams : [fams]).collectEntries { f ->
                    [(f.name.replaceAll(/\.cds-transcripts\.fa\.tsv$/, '')): f] }
                metas.findAll { m -> fam_by.containsKey(m.id) }
                     .collect { m -> tuple(m, fam_by[m.id]) }
            }
            .join(species_ch)
            .map { meta, families, cds -> tuple(meta, cds, families) }
    } else {
        WGD_DMD(species_ch)
        dmd_fams = WGD_DMD.out.families
    }

    // WGD_DMD echoes (meta, cds, families); strip cds for the next step.
    def fam_only = dmd_fams
        .map { meta, cds, families -> tuple(meta, families) }

    def fam_cds = fam_only.join(species_ch)
        .map { meta, families, cds -> tuple(meta, families, cds) }

    // params.wgd_stage = 'dmd': stop after wgd dmd (families cached in outdir/wgd_dmd) -- used to build
    // the families for a whole sample sheet first and count pairs per genome before any dating.
    if (!(params.wgd_stage in ['all', 'dmd'])) {
        error "--wgd_stage must be all or dmd (got '${params.wgd_stage}')"
    }
    def ks_ch = Channel.empty()
    def ks_table_ch = Channel.empty()
    if (params.wgd_stage == 'all') {
        // params.wgd_ks_method: 'pair' (WGD_PAIR_KS, pair-level codeml; default), 'ksd_pairwise'
        // (wgd ksd --pairwise) or 'ksd' (wgd ksd default, family-level gap stripping).
        if (!(params.wgd_ks_method in ['pair', 'ksd_pairwise', 'ksd'])) {
            error "--wgd_ks_method must be pair, ksd_pairwise or ksd (got '${params.wgd_ks_method}')"
        }
        if (params.wgd_ks_method == 'pair') {
            // bundles of about 1 h: params.wgd_bundle_mode = load packs genomes by expected CPU-seconds
            // (params.wgd_bundle_cpu_seconds per bundle); count puts params.wgd_pair_bundle genomes in each
            def grouped = (params.wgd_bundle_mode == 'load')
                ? fam_cds.toList().flatMap { items -> packBundles(items, params.wgd_bundle_cpu_seconds as double) }
                : fam_cds.buffer(size: (params.wgd_pair_bundle as int), remainder: true)
            def bundled = grouped
                .map { lst -> tuple(lst.collect { it[0] }, lst.collect { it[1] }, lst.collect { it[2] }) }
            WGD_PAIR_KS(bundled)
            ks_ch = WGD_PAIR_KS.out.ks.flatten()
        } else {
            WGD_KSD(fam_cds)
            ks_ch = WGD_KSD.out.ks
        }

        if (params.run_wgd_syn.toBoolean()) {
            WGD_SYN(
                fam_only.join(gff_ch)
                    .map { meta, families, gff -> tuple(meta, families, gff) }
            )
        }

        // Final merge: after all of this invocation's ksds finish (or immediately
        // when all are cached), glob every ks.tsv published so far off disk into
        // tables/wgd.ks.parquet. outdir must be absolute (task workdirs don't share
        // the launch dir's relative path).
        def ks_outdir = params.outdir.startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"
        MERGE_WGD_KSD(ks_outdir, ks_ch.collect().map { it }, species_ch.count())
        ks_table_ch = MERGE_WGD_KSD.out.ks_table

        MERGE_WGD_KSD.out.check.subscribe { f ->
            def t = f.text.trim()
            if (t != 'OK') log.error "PARALOGOSCOPE merge check found problems:\n${t}"
        }

    }

    emit:
    families = dmd_fams
    ks       = ks_ch
    anchors  = (params.run_wgd_syn.toBoolean()) ? WGD_SYN.out.anchors : Channel.empty()
    ks_table = ks_table_ch
}
