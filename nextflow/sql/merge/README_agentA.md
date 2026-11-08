# DuckDB merge templates: MEROPS, SIGNALP, TMHMM, WOLFPSORT

Tested 2026-11-06 with DuckDB CLI v1.1.3 and `bin/duckdb_merge.py`.

## Method

- Reference: the module's `merge_<tool>.py` writes a CSV. DuckDB loads it with
  `read_csv_auto(sample_size=-1)` into Parquet. This is the current module path.
- New: `duckdb_merge.py --sql sql/merge/<tool>.sql` on the same files, in `LC_ALL=C` sorted
  order (the Python scripts sort their inputs).
- Checks: row count, column names and types, `EXCEPT ALL` in both directions,
  and `EXCEPT ALL` with a row number added (row order).
- Inputs: files named `<8-hex LOCUSTAG>.<suffix>` only.
  Sample A = 30 files (`awk 'NR%700==1'`). Sample B = 301 files (`awk 'NR%78==1'`).
  Sample B ran in an sbatch job (epyc, 4 CPU). Sample A ran on the interactive node.
- Empty inputs: 6 MEROPS files are empty, and 2 SignalP files contain only the
  GFF header. Both were tested mixed with normal files. Results were identical.

## Results

| Tool | Template | Table | Rows A / B | EXCEPT ALL (ref-new, new-ref, order) | Types equal | Python path A / B (s) | DuckDB A / B (s) |
|---|---|---|---|---|---|---|---|
| MEROPS | merops.sql | merops.parquet | 60,091 / 543,227 | 0, 0, 0 | yes | 1.3 / 16.2 | 0.5 / 1.6 |
| SIGNALP | signalp.sql | signalp.signal_peptide.parquet | 23,922 / 239,376 | 0, 0, 0 | yes | 0.5 / 13.5 | 0.3 / 0.9 |
| TMHMM | tmhmm.sql | tmhmm.parquet | 61,796 / 567,222 | 0, 0, 0 | yes | 3.7 / 24.1 | 0.9 / 3.4 |
| WOLFPSORT | wolfpsort.sql | wolfpsort.parquet | 331,328 / 2,988,120 | 0, 0, 0 | yes | 2.8 / 21.3 | 1.6 / 7.2 |

"Python path" = merge_<tool>.py + CSV load to Parquet. "DuckDB" = duckdb_merge.py wall time.
All output is identical, including row order. Each module writes one table, so there is one template per tool.

## Parsing rules kept

- species_prefix = protein_id up to the first `_`.
- MEROPS: all 12 blast columns, every non-empty line.
- SIGNALP: skip `#` lines. protein_id = GFF column 1 up to the first space. Columns 4, 5, 6.
- TMHMM: value after `=` for each field. Drop proteins with PredHel = 0.
- WOLFPSORT: first (best) localization and its score only. Skip `#` lines.
- Each line is read whole (`delim = chr(1)`, no quote, no escape) and split in SQL.

## Caveat: row order needs threads = 1

DuckDB 1.1.3 `read_csv` over a list of files does not keep the file order when
threads > 1. Rows of one file stay together and in order, but the files are
interleaved. On 30 MEROPS files, 26,339 of 60,091 rows were at a different
position with 4 threads. Each template therefore starts with `SET threads = 1;`.
This overrides `--threads` from duckdb_merge.py.

Cost on 301 WoLF PSORT files (interactive node): 4 threads, unordered: 10.1 s, 199 MB.
1 thread, ordered: 18.8 s, 58 MB. An explicit `ORDER BY file_idx, row_number() OVER ()`
gave the right order but took 18.5 s and 1.4 GB, so it was not used.

The same reader is used in pfam.sql and swissprot_hits.sql. I did not test them.
They probably also lose file order when threads > 1.

## Other caveats

- The Python scripts sort input paths. duckdb_merge.py keeps the given order. The file list
  must be sorted to get the same order as before.
- Content differences are only possible for malformed lines (for example a line with only
  spaces, or a TMHMM line with a different field count). No such lines were found in the
  tested files. I did not scan all 23,407 files for them.
