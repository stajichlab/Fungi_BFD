-- MERGE_TMHMM: per-genome TMHMM short output (*.tmhmm_short.tsv.gz) -> tmhmm.parquet
-- Same rows, columns and types as bin/merge_tmhmm.py + read_csv_auto(sample_size=-1):
-- fields are KEY=VALUE after the protein id; only the value is kept.
-- Proteins with PredHel=0 are dropped (as in merge_tmhmm.py). '#' and blank lines skipped.
-- Row order: DuckDB 1.1.3 read_csv over a file list with threads > 1 does NOT keep
-- the file order (rows of a file stay together and in order, but files are
-- interleaved; measured 2026-11-06). threads = 1 keeps the given file order and
-- was still faster than the Python merge (see README_agentA.md).
SET threads = 1;
COPY (
  WITH l AS (
    SELECT string_split(line, chr(9)) AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = chr(1),
                  quote = '', escape = '', auto_detect = false)
    WHERE len(line) > 0 AND NOT starts_with(line, '#')
  ), v AS (
    SELECT f[1] AS pid, split_part(f[2], '=', 2) AS len_s, split_part(f[3], '=', 2) AS expaa_s,
           split_part(f[4], '=', 2) AS first60_s, split_part(f[5], '=', 2) AS predhel_s,
           split_part(f[6], '=', 2) AS topology
    FROM l
  )
  SELECT split_part(pid, '_', 1) AS species_prefix, pid AS protein_id,
         len_s::BIGINT AS len, expaa_s::DOUBLE AS ExpAA, first60_s::DOUBLE AS First60,
         predhel_s::BIGINT AS PredHel, topology AS Topology
  FROM v
  WHERE predhel_s <> '0'
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
