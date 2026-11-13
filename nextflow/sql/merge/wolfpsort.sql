-- MERGE_WOLFPSORT: per-genome WoLF PSORT output (*.wolfpsort.results.txt.gz) -> wolfpsort.parquet
-- Same rows, columns and types as bin/merge_wolfpsort.py + read_csv_auto(sample_size=-1):
-- line format "<protein_id> <loc> <score>, <loc> <score>, ..."; only the first
-- (best) localization is kept. '#' lines and blank lines are skipped.
-- Row order: DuckDB 1.1.3 read_csv over a file list with threads > 1 does NOT keep
-- the file order (rows of a file stay together and in order, but files are
-- interleaved; measured 2026-11-06). threads = 1 keeps the given file order and
-- was still faster than the Python merge (see README_agentA.md).
SET threads = 1;
COPY (
  WITH l AS (
    SELECT trim(line) AS line
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = chr(1),
                  quote = '', escape = '', auto_detect = false)
    WHERE NOT starts_with(line, '#') AND len(trim(line)) > 0
  ), b AS (
    SELECT split_part(line, ' ', 1) AS pid,
           split_part(substr(line, strpos(line, ' ') + 1), ', ', 1) AS best
    FROM l
  )
  SELECT split_part(pid, '_', 1) AS species_prefix, pid AS protein_id,
         split_part(best, ' ', 1) AS localization, split_part(best, ' ', 2)::DOUBLE AS score
  FROM b
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
