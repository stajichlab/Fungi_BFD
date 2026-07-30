//
// ANI_COMPARE_METHOD — sketch + compare a group of genomes with the selected
// --ani_method, emitting one normalized ANI table per group.
//
// All four methods converge on the same output contract:
//     tuple( group_name, path(<group>.ani.tsv) )   cols: query <TAB> reference <TAB> ANI
//
// Isolating the dispatch here keeps workflows/compare_ANI.nf to channel
// construction and reporting (REFACTOR_NEXTFLOW_PLAN.md §2.5).
//

include { MASH_SKETCH      } from '../../modules/ani/sketch/MASH_SKETCH/main.nf'
include { SOURMASH_SKETCH  } from '../../modules/ani/sketch/SOURMASH_SKETCH/main.nf'
include { SKANI_COMPARE    } from '../../modules/ani/compare/SKANI_COMPARE/main.nf'
include { MASH_COMPARE     } from '../../modules/ani/compare/MASH_COMPARE/main.nf'
include { SOURMASH_COMPARE } from '../../modules/ani/compare/SOURMASH_COMPARE/main.nf'
include { FASTANI_COMPARE  } from '../../modules/ani/compare/FASTANI_COMPARE/main.nf'
include { MASH_PREFILTER   } from '../../modules/ani/compare/MASH_PREFILTER/main.nf'
include { MASH_COMPONENTS  } from '../../modules/ani/compare/MASH_COMPONENTS/main.nf'
include { MERGE_ANI        } from '../../modules/ani/report/MERGE_ANI/main.nf'

workflow ANI_COMPARE_METHOD {
    take:
    genome_ch      // tuple(group_name, [genome files])
    method         // validated, lower-case: skani | mash | sourmash | fastani

    main:
    // Flattened (group, genome) for the per-genome sketchers.
    def genome_flat = genome_ch.flatMap { group_name, genomes ->
        genomes.collect { g -> tuple(group_name, g) }
    }

    def ani_tsv_ch

    if (method == 'skani') {
        // skani 0.3.x: sketch all genomes in one pass, then skani triangle.
        // Chunking was for the old per-genome .sketch API; skani triangle handles
        // large groups efficiently with a single sketch directory.
        // Pass n_genomes as an explicit val so the process resource directive can
        // use it reliably (channel inputs aren't accessible in closures at submit time).
        def skani_input = genome_ch.map { gn, genomes -> tuple(gn, genomes.size(), genomes) }
        ani_tsv_ch = SKANI_COMPARE(skani_input)

    } else if (method == 'mash') {
        ani_tsv_ch = MASH_SKETCH(genome_flat)
            .groupTuple()
            | MASH_COMPARE

    } else if (method == 'sourmash') {
        ani_tsv_ch = SOURMASH_SKETCH(genome_flat)
            .groupTuple()
            | SOURMASH_COMPARE

    } else { // fastani
        def part_ch

        if (params.fastani_prefilter as boolean) {
            // (b) mash pre-cluster, then fastANI all-vs-all within each component.
            def grouped_msk = MASH_SKETCH(genome_flat).groupTuple()
            def comps_ch    = MASH_PREFILTER(grouped_msk) | MASH_COMPONENTS

            // Explode components.tsv (comp_id <TAB> genome_filename) into rows
            // keyed by (group, filename) so we can re-attach the genome file path.
            def comp_rows = comps_ch
                .splitCsv(elem: 1, sep: '\t')
                .map { gn, row -> tuple(tuple(gn, row[1]), tuple(gn, row[0])) }

            def genome_lut = genome_ch.flatMap { gn, genomes ->
                genomes.collect { g -> tuple(tuple(gn, g.name), g) }
            }

            def per_comp = comp_rows
                .combine(genome_lut, by: 0)                       // (key, (gn,comp), gpath)
                .map { _key, gncomp, gpath -> tuple(gncomp[0], gncomp[1], gpath) }
                .groupTuple(by: [0, 1])                           // (gn, comp, [gpaths])
                .filter { _gn, _comp, gs -> gs.size() >= params.min_group_size as int }
                .map { gn, comp, gs -> tuple(gn, gs, gs, "c${comp}", gs.size()) }

            part_ch = FASTANI_COMPARE(per_comp)

        } else {
            def batchSize = params.ani_batch_size as int

            if (batchSize > 0) {
                // Upper-triangle batch pairs; carry total group size for CPU scaling.
                def batch_pairs_ch = genome_ch.flatMap { group_name, genomes ->
                    def N      = genomes.size()
                    def nBatch = (N + batchSize - 1).intdiv(batchSize)
                    def batches = (0..<nBatch).collect { i ->
                        def end = (i + 1) * batchSize < N ? (i + 1) * batchSize : N
                        genomes.subList(i * batchSize, end)
                    }
                    (0..<nBatch).collectMany { i ->
                        (i..<nBatch).collect { j ->
                            tuple(group_name, batches[i], batches[j], "b${i}_b${j}", N)
                        }
                    }
                }
                part_ch = FASTANI_COMPARE(batch_pairs_ch)

            } else {
                // Single job per group; query == ref == all genomes.
                def single_ch = genome_ch.map { group_name, genomes ->
                    tuple(group_name, genomes, genomes, "full", genomes.size())
                }
                part_ch = FASTANI_COMPARE(single_ch)
            }
        }

        // Merge per-group partials (1 for single job, many for batch/prefilter).
        ani_tsv_ch = part_ch
            .groupTuple()
            .map { gn, tsv_list -> tuple(gn, tsv_list.flatten()) }
            | MERGE_ANI
    }

    emit:
    ani_tsv = ani_tsv_ch
}
