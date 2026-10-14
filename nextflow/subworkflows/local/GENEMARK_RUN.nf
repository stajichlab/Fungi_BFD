// DRAFT. Drop-in replacement for the GENEMARK_RUN process: same take (one 11-field row per
// genome), same emits (gtf, other_gff, mod), now a subworkflow. Callers keep calling
// GENEMARK_RUN(...) and GENEMARK_RUN.out.{gtf,other_gff,mod}; only the include path changes.
//
//   GENEMARK_TRAIN ──┬─ monolithic / skip route ─────────────────────────────► gtf, other_gff, mod
//                    └─ chunked route (params.genemark_chunked, >= N contigs)
//                         chunks ► GENEMARK_PREDICT_CHUNK (one per chunk, --predict_with the model)
//                                    ► groupTuple(groupKey(out, n_chunks)) ► GENEMARK_MERGE ► gtf, other_gff, mod
//
// With params.genemark_chunked=false (default) GENEMARK_TRAIN never takes the chunked route and
// this behaves exactly like the old single process.

include { GENEMARK_TRAIN }         from '../../modules/funannotate/predict/GENEMARK_TRAIN/main.nf'
include { GENEMARK_PREDICT_CHUNK } from '../../modules/funannotate/predict/GENEMARK_PREDICT_CHUNK/main.nf'
include { GENEMARK_MERGE }         from '../../modules/funannotate/predict/GENEMARK_MERGE/main.nf'

workflow GENEMARK_RUN {
    take:
    rows   // tuple(out, asmid, species, strain, genome_fa, transl_table, mode, training_bam, force_independent, shared_mod, is_microsporidia)

    main:
    GENEMARK_TRAIN(rows)

    // per-genome facts needed later: species, genome path, translation table
    meta = rows.map { out, _asmid, species, _strain, genome_fa, transl_table, _mode, _bam, _fi, _sm, _micro ->
        tuple(out, species, genome_fa, transl_table)
    }

    // model per genome in the chunked route: fresh (just trained) or the stored shared one
    model = GENEMARK_TRAIN.out.train_mod
        .map { out, _species, mod -> tuple(out, mod, true) }
        .mix(GENEMARK_TRAIN.out.reuse_mod.map { out, mod -> tuple(out, mod, false) })

    // scatter: one row per chunk file, in chunk order
    chunk_rows = GENEMARK_TRAIN.out.chunks
        .flatMap { out, files ->
            // sort on the integer in chunk_<NNNN>.fa.gz, not on the name, so ordering never depends on zero padding
            def fl = (files instanceof List ? files : [files]).sort { (it.name =~ /chunk_(\d+)\.fa\.gz/)[0][1] as int }
            fl.withIndex().collect { f, i -> tuple(out, String.format('%04d', i), fl.size(), f) }
        }
        .combine(model, by: 0)
        .combine(meta, by: 0)
        .map { out, idx, n, chunk, mod, _fresh, _species, _genome, tt -> tuple(out, idx, n, chunk, mod, tt) }

    GENEMARK_PREDICT_CHUNK(chunk_rows)

    // gather: wait for all n chunks of a genome (groupKey carries the expected count)
    gathered = GENEMARK_PREDICT_CHUNK.out.gtf
        .map { out, _idx, n, gtf -> tuple(groupKey(out, n as int), gtf) }
        .groupTuple()
        .map { key, gtfs -> tuple(key.getGroupTarget(), key.getGroupSize(), gtfs) }

    // failOnDuplicate: two rows with the same `out` would cross-pair chunks and models. (Not failOnMismatch: one
    // permanently failed chunk must only drop its own genome, not stop the run.)
    merge_in = gathered
        .join(model, failOnDuplicate: true)
        .join(GENEMARK_TRAIN.out.other_pending, failOnDuplicate: true)
        .join(meta, failOnDuplicate: true)
        .map { out, n, gtfs, mod, fresh, pending, species, genome, _tt ->
            tuple(out, n, species, genome, fresh, gtfs, mod, pending)
        }

    GENEMARK_MERGE(merge_in)

    emit:
    gtf       = GENEMARK_TRAIN.out.gtf.mix(GENEMARK_MERGE.out.gtf)
    other_gff = GENEMARK_TRAIN.out.other_gff.mix(GENEMARK_MERGE.out.other_gff)
    mod       = GENEMARK_TRAIN.out.mod.mix(GENEMARK_MERGE.out.mod)
}
