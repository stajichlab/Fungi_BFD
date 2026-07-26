//
// BFD Functional Annotation Pipeline — main workflow definition
//
// All process logic lives in modules/BFD/ and modules/common/.
// This file orchestrates the channel wiring and conditional execution.
//

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// ── Channel helper functions (need DSL scope — cannot live in lib/ Groovy class) ──

def gatedGlobIn(sync_ch, String baseDir, String glob) {
    sync_ch
        .flatMap { files("${baseDir}/${glob}") }
        .filter  { it.size() > 0 }
        .collect()
        .filter  { !it.isEmpty() }
}

def gatedGlob(sync_ch, String glob) {
    gatedGlobIn(sync_ch, params.outdir, glob)
}

def gatedGlobStats(sync_ch, String glob) {
    gatedGlobIn(sync_ch, params.genome_stats_outdir, glob)
}

def toManifest(ch, String name) {
    ch.flatten()
      .map { f -> "${f.toString()}\t${f.lastModified()}\t${f.size()}" }
      .collectFile(name: name, newLine: true, sort: true)
}

// ── BFD modules ─────────────────────────────────────────────────────────────
include { RUN_PFAM       } from '../modules/BFD/PFAM/main.nf'
include { RUN_CAZY       } from '../modules/BFD/CAZY/main.nf'
include { RUN_MEROPS     } from '../modules/BFD/MEROPS/main.nf'
include { RUN_SIGNALP    } from '../modules/BFD/SIGNALP/main.nf'
include { RUN_TMHMM      } from '../modules/BFD/TMHMM/main.nf'
include { RUN_TARGETP    } from '../modules/BFD/TARGETP/main.nf'
include { RUN_IDP        } from '../modules/BFD/IDP/main.nf'
include { RUN_WOLFPSORT  } from '../modules/BFD/WOLFPSORT/main.nf'
include { RUN_PREDGPI    } from '../modules/BFD/PREDGPI/main.nf'

// ── Per-genome statistics modules ───────────────────────────────────────────
include { BATCH_AA_FREQ    } from '../modules/BFD/BATCH_AA_FREQ/main.nf'
include { BATCH_CODON_FREQ } from '../modules/BFD/BATCH_CODON_FREQ/main.nf'
include { CALC_INTERGENIC  } from '../modules/BFD/CALC_INTERGENIC/main.nf'
include { CALC_GENE_STATS  } from '../modules/BFD/CALC_GENE_STATS/main.nf'
include { CALC_ASM_STATS   } from '../modules/BFD/CALC_ASM_STATS/main.nf'
include { BUSCO_GENOME     } from '../modules/BFD/BUSCO_GENOME/main.nf'
include { BUSCO_PEP        } from '../modules/BFD/BUSCO_PEP/main.nf'

// ── Merge modules ───────────────────────────────────────────────────────────
include { MERGE_PFAM       } from '../modules/BFD/MERGE_PFAM/main.nf'
include { MERGE_CAZY       } from '../modules/BFD/MERGE_CAZY/main.nf'
include { MERGE_MEROPS     } from '../modules/BFD/MERGE_MEROPS/main.nf'
include { MERGE_SIGNALP    } from '../modules/BFD/MERGE_SIGNALP/main.nf'
include { MERGE_TMHMM      } from '../modules/BFD/MERGE_TMHMM/main.nf'
include { MERGE_TARGETP    } from '../modules/BFD/MERGE_TARGETP/main.nf'
include { MERGE_IDP        } from '../modules/BFD/MERGE_IDP/main.nf'
include { MERGE_WOLFPSORT  } from '../modules/BFD/MERGE_WOLFPSORT/main.nf'
include { MERGE_PREDGPI    } from '../modules/BFD/MERGE_PREDGPI/main.nf'
include { MERGE_AA_FREQ    } from '../modules/BFD/MERGE_AA_FREQ/main.nf'
include { MERGE_CODON_FREQ } from '../modules/BFD/MERGE_CODON_FREQ/main.nf'
include { MERGE_INTERGENIC } from '../modules/BFD/MERGE_INTERGENIC/main.nf'
include { MERGE_GENE_STATS } from '../modules/BFD/MERGE_GENE_STATS/main.nf'
include { MERGE_ASM_STATS  } from '../modules/BFD/MERGE_ASM_STATS/main.nf'
include { MERGE_SAMPLES    } from '../modules/BFD/MERGE_SAMPLES/main.nf'

// ── Subworkflows ────────────────────────────────────────────────────────────
include { INPUT_SETUP } from '../subworkflows/local/INPUT_SETUP.nf'

// ── Shared helper functions ─────────────────────────────────────────────────
include { clearIfStale } from '../modules/common/utils.nf'

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow BFD {

    // ── --help: print parameter documentation and exit ──────────────────────────
    if (params.help) {
        log.info paramsHelp(command: "nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile BFD")
        return
    }
    // ── Validate params and samples.csv structure (fail fast) ───────────────────
    validateParameters()
    log.info paramsSummaryLog(workflow)
    samplesheetToList(params.samples, "${projectDir}/assets/schema_input.json")

    // ── Taxonomy filter ────────────────────────────────────────────────────────
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
    def rows_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def species  = row.SPECIES?.trim() ?: ''
            def strain   = SampleUtils.cleanStrain(row.STRAIN?.trim() ?: '')
            def locustag = row.LOCUSTAG?.replaceAll(/[\r\n]/, '')?.trim()
            def basename = SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            tuple(locustag, basename, species, strain)
        }
        .take((params.n_test as int) > 0 ? (params.n_test as int) : -1)

    // ── Input setup: symlink predict_results files into input/ subdirs ─────────
    def ready_ch
    if (params.run_setup.toBoolean()) {
        INPUT_SETUP(rows_ch)
        ready_ch = INPUT_SETUP.out.done
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
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    codon_freq_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.cds_dir}/${basename}.cds-transcripts.fa", glob: false)
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    intergenic_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.gff_dir}/${basename}.gff3", glob: false)
        f.exists() ? tuple(locustag, basename, f) : null
    }.filter { it != null }

    gene_stats_ch = ready_ch.map { locustag, basename, species, strain ->
        def gff = file("${params.gff_dir}/${basename}.gff3", glob: false)
        def dna = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
        (gff.exists() && dna.exists()) ? tuple(locustag, basename, gff, dna) : null
    }.filter { it != null }

    asmid_by_basename_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .map { row ->
            def basename = SampleUtils.makeSampleTag(row.SPECIES?.trim() ?: '', row.STRAIN?.trim() ?: '')
            tuple(basename, row.ASMID?.trim() ?: '')
        }
        .filter { basename, asmid -> asmid }

    asm_stats_ch = ready_ch
        .map { locustag, basename, species, strain ->
            def f = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
            f.exists() ? tuple(basename, f) : null
        }
        .filter { it != null }
        .join(asmid_by_basename_ch)
        .map { basename, genome, asmid -> tuple(asmid, basename, genome) }

    // ── BUSCO input channels ───────────────────────────────────────────────────
    busco_genome_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.genome_dir}/${basename}.scaffolds.fa", glob: false)
        f.exists() ? tuple(locustag, basename, params.busco_lineage, f) : null
    }.filter { it != null }

    busco_pep_ch = ready_ch.map { locustag, basename, species, strain ->
        def f = file("${params.pep_dir}/${basename}.proteins.fa", glob: false)
        f.exists() ? tuple(locustag, basename, params.busco_lineage, f) : null
    }.filter { it != null }

    // ── Per-species RUN steps ──────────────────────────────────────────────────
    if (params.run_pfam.toBoolean())
        RUN_PFAM(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/pfam_hmmscan/${basename}.domtblout.gz"),
                file("${params.outdir}/pfam_hmmscan/${basename}.tblout.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_cazy.toBoolean())
        RUN_CAZY(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/cazy/${basename}/${basename}.overview.tsv.gz"),
                file("${params.outdir}/cazy/${basename}/${basename}.cazymes.tsv.gz"),
                file("${params.outdir}/cazy/${basename}/${basename}.substrates.tsv.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_merops.toBoolean())
        RUN_MEROPS(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/merops/${basename}.blasttab.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_signalp.toBoolean())
        RUN_SIGNALP(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/signalp/${basename}.signalp.gff3.gz"),
                file("${params.outdir}/signalp/${basename}.signalp.results.txt.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_tmhmm.toBoolean())
        RUN_TMHMM(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/tmhmm/${basename}.tmhmm_short.tsv.gz"),
                file("${params.outdir}/tmhmm/${basename}.tmhmm_results.tsv.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_targetp.toBoolean())
        RUN_TARGETP(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/targetP/${basename}_summary.targetp2.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_idp.toBoolean())
        RUN_IDP(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/aiupred/${basename}.aiupred.txt.gz"),
                file("${params.outdir}/aiupred/${basename}.idp.csv.gz"),
                file("${params.outdir}/aiupred/${basename}.idp_summary.csv.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_wolfpsort.toBoolean())
        RUN_WOLFPSORT(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/wolfpsort/${basename}.wolfpsort.results.txt.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })
    if (params.run_predgpi.toBoolean())
        RUN_PREDGPI(proteins_ch.map { locustag, basename, sp, st, prot ->
            clearIfStale(prot, [
                file("${params.outdir}/predgpi/${basename}.predgpi.gff3.gz")
            ])
            tuple(locustag, basename, sp, st, prot)
        })

    // ── MERGE steps ────────────────────────────────────────────────────────────
    if (!params.skip_merge.toBoolean()) {

        if (params.merge_all.toBoolean() && !params.taxon) {

            if (params.run_pfam.toBoolean()) {
                def sync = RUN_PFAM.out.domtbl.collect()
                MERGE_PFAM( gatedGlob(sync, "pfam_hmmscan/*.domtblout.gz") )
            } else {
                MERGE_PFAM( gatedGlob(Channel.of(true), "pfam_hmmscan/*.domtblout.gz") )
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
    def use_glob = params.merge_all.toBoolean() && !params.taxon

    if (params.run_aa_freq.toBoolean()) {
        def aa_batch_ch = aa_freq_ch
            .map { locustag, basename, prot ->
                clearIfStale(prot, [file("${params.genome_stats_outdir}/aa_freq/${basename}.aa_freq.csv.gz")])
                file("${params.genome_stats_outdir}/aa_freq/${basename}.aa_freq.csv.gz").exists()
                    ? null : tuple(locustag, basename, prot)
            }
            .filter { it != null }
            .buffer(size: params.freq_batch_size as int, remainder: true)
            .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }, batch.collect { it[2] }) }
        BATCH_AA_FREQ(aa_batch_ch)
    }
    if (params.run_codon_freq.toBoolean()) {
        def codon_batch_ch = codon_freq_ch
            .map { locustag, basename, prot ->
                clearIfStale(prot, [file("${params.genome_stats_outdir}/codon_freq/${basename}.codon_freq.csv.gz")])
                file("${params.genome_stats_outdir}/codon_freq/${basename}.codon_freq.csv.gz").exists()
                    ? null : tuple(locustag, basename, prot)
            }
            .filter { it != null }
            .buffer(size: params.freq_batch_size as int, remainder: true)
            .map { batch -> tuple(batch.collect { it[0] }, batch.collect { it[1] }, batch.collect { it[2] }) }
        BATCH_CODON_FREQ(codon_batch_ch)
    }
    if (params.run_intergenic.toBoolean())
        CALC_INTERGENIC(intergenic_ch.map { locustag, basename, gff ->
            clearIfStale(gff, [
                file("${params.genome_stats_outdir}/intergenic_stats/${basename}.gene_intergenic_distances.csv.gz")
            ])
            tuple(locustag, basename, gff)
        })
    if (params.run_gene_stats.toBoolean())
        CALC_GENE_STATS(gene_stats_ch.map { locustag, basename, gff, dna ->
            def outs = ['gene_info', 'gene_exons', 'gene_CDS', 'gene_introns',
                        'gene_transcripts', 'gene_trnas', 'gene_proteins'].collect {
                file("${params.genome_stats_outdir}/gene_stats/${basename}.${it}.csv.gz")
            }
            clearIfStale(gff, outs)
            clearIfStale(dna, outs)
            tuple(locustag, basename, gff, dna)
        })

    if (params.run_asm_stats.toBoolean())
        CALC_ASM_STATS(asm_stats_ch.map { asmid, basename, genome ->
            clearIfStale(genome, [
                file("${params.genome_stats_outdir}/asm_stats/${asmid}.stats.txt")
            ])
            tuple(asmid, basename, genome)
        })

    if (params.run_busco_genome.toBoolean())
        BUSCO_GENOME(busco_genome_ch.map { locustag, basename, lineage, genome ->
            clearIfStale(genome, [
                file("${params.genome_stats_outdir}/BUSCO_genome/${basename}.BUSCO_summary.${lineage}.txt")
            ])
            tuple(locustag, basename, lineage, genome)
        })

    if (params.run_busco_pep.toBoolean())
        BUSCO_PEP(busco_pep_ch.map { locustag, basename, lineage, prot ->
            clearIfStale(prot, [
                file("${params.genome_stats_outdir}/BUSCO_protein/${basename}.BUSCO_summary.${lineage}.txt")
            ])
            tuple(locustag, basename, lineage, prot)
        })

    if (!params.skip_merge.toBoolean()) {
        def matched_keys_file = ready_ch
            .map { locustag, basename, species, strain -> (locustag ?: '').trim() }
            .filter { it }
            .collectFile(name: 'matched_locustags.txt', newLine: true)
        MERGE_SAMPLES(file(params.samples), matched_keys_file)

        if (use_glob) {
            if (params.run_aa_freq.toBoolean()) {
                MERGE_AA_FREQ(toManifest(gatedGlobStats(BATCH_AA_FREQ.out.csv.flatten().collect().ifEmpty([]), "aa_freq/*.aa_freq.csv.gz"), 'aa_freq.manifest.txt'))
            } else {
                MERGE_AA_FREQ(toManifest(gatedGlobStats(Channel.of(true), "aa_freq/*.aa_freq.csv.gz"), 'aa_freq.manifest.txt'))
            }
            if (params.run_codon_freq.toBoolean()) {
                MERGE_CODON_FREQ(toManifest(gatedGlobStats(BATCH_CODON_FREQ.out.csv.flatten().collect().ifEmpty([]), "codon_freq/*.codon_freq.csv.gz"), 'codon_freq.manifest.txt'))
            } else {
                MERGE_CODON_FREQ(toManifest(gatedGlobStats(Channel.of(true), "codon_freq/*.codon_freq.csv.gz"), 'codon_freq.manifest.txt'))
            }
            if (params.run_intergenic.toBoolean()) {
                MERGE_INTERGENIC(toManifest(gatedGlobStats(CALC_INTERGENIC.out.csv.collect(), "intergenic_stats/*.gene_intergenic_distances.csv.gz"), 'intergenic.manifest.txt'))
            } else {
                MERGE_INTERGENIC(toManifest(gatedGlobStats(Channel.of(true), "intergenic_stats/*.gene_intergenic_distances.csv.gz"), 'intergenic.manifest.txt'))
            }
            if (params.run_gene_stats.toBoolean()) {
                MERGE_GENE_STATS(toManifest(gatedGlobStats(
                    CALC_GENE_STATS.out.gene_info
                        .mix(CALC_GENE_STATS.out.gene_exons)
                        .mix(CALC_GENE_STATS.out.gene_CDS)
                        .mix(CALC_GENE_STATS.out.gene_introns)
                        .mix(CALC_GENE_STATS.out.gene_transcripts)
                        .mix(CALC_GENE_STATS.out.gene_trnas)
                        .mix(CALC_GENE_STATS.out.gene_proteins)
                        .collect(),
                    "gene_stats/*.csv.gz"
                ), 'gene_stats.manifest.txt'))
            } else {
                MERGE_GENE_STATS(toManifest(gatedGlobStats(Channel.of(true), "gene_stats/*.csv.gz"), 'gene_stats.manifest.txt'))
            }
            if (params.run_asm_stats.toBoolean()) {
                MERGE_ASM_STATS(toManifest(gatedGlobStats(CALC_ASM_STATS.out.stats.collect(), "asm_stats/*.stats.txt"), 'asm_stats.manifest.txt'))
            } else {
                MERGE_ASM_STATS(toManifest(gatedGlobStats(Channel.of(true), "asm_stats/*.stats.txt"), 'asm_stats.manifest.txt'))
            }
        } else {
            if (params.run_aa_freq.toBoolean()) {
                def aa_sync  = BATCH_AA_FREQ.out.csv.flatten().collect().map { true }.ifEmpty(true)
                def aa_paths = aa_freq_ch.map { locustag, basename, _ignored ->
                    file("${params.genome_stats_outdir}/aa_freq/${basename}.aa_freq.csv.gz")
                }
                MERGE_AA_FREQ(toManifest(
                    aa_sync.combine(aa_paths)
                           .map { _s, f -> f }
                           .filter { it.exists() },
                    'aa_freq.manifest.txt'
                ))
            }
            if (params.run_codon_freq.toBoolean()) {
                def codon_sync  = BATCH_CODON_FREQ.out.csv.flatten().collect().map { true }.ifEmpty(true)
                def codon_paths = codon_freq_ch.map { locustag, basename, _ignored ->
                    file("${params.genome_stats_outdir}/codon_freq/${basename}.codon_freq.csv.gz")
                }
                MERGE_CODON_FREQ(toManifest(
                    codon_sync.combine(codon_paths)
                              .map { _s, f -> f }
                              .filter { it.exists() },
                    'codon_freq.manifest.txt'
                ))
            }
            if (params.run_intergenic.toBoolean()) MERGE_INTERGENIC(toManifest(CALC_INTERGENIC.out.csv.collect(), 'intergenic.manifest.txt'))
            if (params.run_gene_stats.toBoolean()) {
                MERGE_GENE_STATS(toManifest(
                    CALC_GENE_STATS.out.gene_info
                        .mix(CALC_GENE_STATS.out.gene_exons)
                        .mix(CALC_GENE_STATS.out.gene_CDS)
                        .mix(CALC_GENE_STATS.out.gene_introns)
                        .mix(CALC_GENE_STATS.out.gene_transcripts)
                        .mix(CALC_GENE_STATS.out.gene_trnas)
                        .mix(CALC_GENE_STATS.out.gene_proteins)
                        .collect(),
                    'gene_stats.manifest.txt'
                ))
            }
            if (params.run_asm_stats.toBoolean()) MERGE_ASM_STATS(toManifest(CALC_ASM_STATS.out.stats.collect(), 'asm_stats.manifest.txt'))
        }
    }
}
