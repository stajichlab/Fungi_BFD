# Chunked GeneMark (train once, scatter predict, gather)

`GENEMARK_RUN` is now a subworkflow with the same name, take and emits (`gtf`, `other_gff`, `mod`) as the old
process. With `params.genemark_chunked = false` (the default) every genome takes the old single-task route
(`GENEMARK_TRAIN`, the old process body renamed) and produces the old files.

## Route
- Genomes with at least `genemark_chunk_min_contigs` (25,000) contigs, in ES or fast-reuse mode, stop after
  training (or immediately, for a stored model) and write `chunks/chunk_NNNN.fa.gz` (about
  `genemark_chunk_contigs` = 20,000 contigs each, cut by contig count, genome order).
- `GENEMARK_PREDICT_CHUNK` runs `gmes_petap.pl --predict_with <model>` per chunk (8 CPUs, node-local `$SCRATCH`).
  Attempt 1 uses the 2 h `short` queue; attempt 2 `preempt` 6 h; attempt 3 the default queue 12 h; then `ignore`.
  `maxForks = 32`.
- `GENEMARK_MERGE` merges in chunk order, renumbers gene and transcript IDs, validates against the whole genome,
  and publishes to `genemark_store/<out>/` (the only place the chunked route publishes; a freshly trained model is
  published only after the merge).
- ET rows with a training BAM stay on the single-task route (chunked ET is untested).

## Why it is safe
Same model, chunked vs monolithic `--predict_with`: byte-identical GTF on Racocetra (43,967 contigs), a 110-contig
genome with a 9.47 Mb contig, a 1,759-contig genome, and Jimgerdemannia_flammicorona_AD002 (26,180 contigs, compared with
the GTF a production run stored on 2026-10-03). Massospora_platypediae (241,057 contigs): all 233,526 whole-contig
pieces that the killed monolithic job had finished match the merged GTF. ES training is NOT deterministic (two
identical runs give different `gmhmm.mod`), so the stored `.mod` is the reference.

Production shell is `/bin/bash -l` (no `set -e`): the chunk and merge scripts set `set -euo pipefail`; the train
script keeps its old options and checks every exit 0 with `route_check`.

## After a run
- `genemark_chunked.tsv` in `params.target`: one row per chunk (contigs and bases) for every genome that took the route.
- `bin/genemark_chunked_missing.py --chunked-tsv <target>/genemark_chunked.tsv --store <genemark_store>` lists
  chunked genomes with no GTF in the store (a permanently failed chunk or merge is `ignore`d and leaves no GTF).
- `bin/list_failed_genomes.py` reports failed chunks by genome with a `chunk` column.

## Known limits
- Chunks are cut by contig count only. Two chunks of equal contig count but 3.5x different bases ran 30 and 52 min
  (Jimgerdemannia); chunk sizes are logged for review.
- Not tested: genetic-code tables 6/26 through the chunked route, ET, a chunk on an abu_dhabi node.
- `maxForks` applies per process instance, so the `_SIB` alias can add its own 32.
