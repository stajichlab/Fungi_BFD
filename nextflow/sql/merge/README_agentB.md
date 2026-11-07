# DuckDB merge templates: CAZY, PREDGPI, TARGETP, INFERNAL, TF (2026-11-06)

Runner: `bin/duckdb_merge.py --sql sql/merge/<tpl> --files <list> --out <parquet>`.
Each template replaces one `merge_*.py` step plus its `read_csv_auto(sample_size=-1)` step.
All write `(FORMAT PARQUET, COMPRESSION ZSTD)`.

| Module | Template | Table | Extra input |
|---|---|---|---|
| MERGE_CAZY | cazy_overview.sql | cazy.overview | none (list of `*.overview.tsv.gz`) |
| MERGE_CAZY | cazy_cazymes_hmm.sql | cazy.cazymes_hmm | none (list of `*.cazymes.tsv.gz`) |
| MERGE_PREDGPI | predgpi.sql | predgpi | none |
| MERGE_TARGETP | targetp.sql | targetP | none |
| MERGE_INFERNAL | infernal.sql | infernal_rfam | none |
| MERGE_TF | tf.sql | tf_inventory | `--param tf_domains=<path>/assets/fungal_tf_pfam_domains.csv` (raw path, the template adds the quotes) |

MERGE_CAZY needs two runner calls, one per table, each with its own file list.

## Verification

Method: run the Python script as the module does, load its CSV with
`read_csv_auto(sample_size=-1)` (reference), run the template on the same files,
then compare row count, column names and types, and `EXCEPT ALL` in both directions.
Two rounds: 30 files (login node), then up to 1,000 files (sbatch, epyc, 4 CPUs, `--threads 4`).
Files were spread across the folder. Pfam/predGPI samples include one 9-character locustag (F3DC5624A).

Results (1,000-file round; the 30-file round also gave 0 / 0 for every table):

| Table | Files | Rows | ref EXCEPT new | new EXCEPT ref | Types equal | Python + CSV->Parquet (s) | Template (s) |
|---|---|---|---|---|---|---|---|
| cazy.overview | 1,000 | 400,731 | 0 | 0 | yes | 19.5 + 0.9 | 1.2 |
| cazy.cazymes_hmm | 1,000 | 394,538 | 0 | 0 | yes | 17.9 + 1.0 | 1.0 |
| predgpi | 1,000 | 9,924,221 | 0 | 0 | yes | 64.5 + 15.5 | 10.9 |
| targetP | 125 (all) | 178,535 | 0 | 0 | yes | 6.1 + 0.6 | 1.5 |
| infernal_rfam | 15 (all) | 3,925 | 0 | 0 | yes | 0.4 + 0.1 | 0.4 |
| tf_inventory | 1,000 | 483,388 | 0 | 0 | yes | 89.7 + 1.3 | 18.0 |

Template time includes DuckDB start-up and the runner's row count.
The CAZY Python times are for one table each (the other list was given one file).
Column types equal the current tables in `Fungi_BFD_runs/tables/` (tf_inventory has no table yet; compared to the Python output only).

## Caveats

1. Row order. With `--threads` > 1, DuckDB 1.1.3 can write whole input files in a
   different order than the list. Rows of one file stay together and keep their
   file order (checked: 30/30 files contiguous, same internal order). With `--threads 1`
   the order equals the Python order exactly (0 mismatches for all six tables).
   This is a property of multi-file `read_csv`, so it applies to pfam.sql and the others too.
   Also note: the Python scripts sort by the staged basename. The runner keeps the list order,
   so the list must be sorted the same way to get the same file order.
2. INFERNAL output is wrong today, and the template copies it on purpose.
   `merge_infernal.py` uses `--fmt 1` column indexes, but the files are cmscan `--fmt 2`
   (leading `idx` column and a `clan name` column). So: target_name = idx (BIGINT),
   target_acc = target name, query_name and contig = Rfam accession, mdl_from = clan,
   mdl_to = 'cm', seq_from/seq_to = model from/to, strand = seq from, score = GC fraction,
   evalue = bias. Correct fmt-2 fields (1-based after whitespace split): 2 target name,
   3 accession, 4 query (contig), 8-9 mdl from/to, 10-11 seq from/to, 12 strand,
   17 score, 18 E-value. A fix changes the table schema; it is not done here.
   Also, species_prefix = basename up to the first '.', so names with '.' are cut
   (Absidia_glauca_CBS_101.48... -> Absidia_glauca_CBS_101) and two genomes can share one
   prefix (Acarospora_aff._strigata and Acarospora_aff._strigata_T1882 -> Acarospora_aff).
3. TF domain list: the PF11955 (WOPR) row in `assets/fungal_tf_pfam_domains.csv` has an
   unquoted comma. Python keeps tf_family = "WOPR DNA-binding domain (fungal"; the template
   does the same (map checked against the Python dict: 25 entries, same values).
   PF11955 had no hits in the 1,000 test genomes, so this row was not exercised by data.
   The TF template uses a regex pre-filter on the listed accessions; without it the template
   was slower than Python (9.2 s vs 4.8 s on 30 files).
4. Edge cases not seen in the data: Python `csv.reader` (CAZY, predGPI) treats '"' as a quote
   character; the templates do not. No '"' was found in the 30-file samples of these inputs.
   Empty fields become NULL, as they do after the CSV round trip.
