-- MERGE_MEROPS: per-genome MEROPS blastp -outfmt 6 tables (*.blasttab.gz) -> merops.parquet
-- Same rows, columns and types as bin/merge_merops.py + read_csv_auto(sample_size=-1):
-- every non-empty line is one row; species_prefix = protein_id up to the first '_'.
-- Each line is read whole (delim chr(1) never occurs) and split on TAB, so blank lines and empty files are safe.
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
    WHERE len(line) > 0
  )
  SELECT split_part(f[1], '_', 1) AS species_prefix, f[1] AS protein_id, f[2] AS merops_id,
         f[3]::DOUBLE AS percent_identity, f[4]::BIGINT AS aln_length,
         f[5]::BIGINT AS mismatches, f[6]::BIGINT AS gap_openings,
         f[7]::BIGINT AS q_start, f[8]::BIGINT AS q_end,
         f[9]::BIGINT AS s_start, f[10]::BIGINT AS s_end,
         f[11]::DOUBLE AS evalue, f[12]::DOUBLE AS bitscore
  FROM l
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
