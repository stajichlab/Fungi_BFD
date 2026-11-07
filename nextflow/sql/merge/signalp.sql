-- MERGE_SIGNALP: per-genome SignalP 6 GFF3 (*.signalp.gff3.gz) -> signalp.signal_peptide.parquet
-- Same rows, columns and types as bin/merge_signalp.py + read_csv_auto(sample_size=-1):
-- lines starting with '#' and blank lines are skipped; protein_id = GFF column 1
-- up to the first space; species_prefix = protein_id up to the first '_';
-- peptide_start/end/probability = GFF columns 4/5/6.
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
  )
  SELECT split_part(split_part(f[1], ' ', 1), '_', 1) AS species_prefix,
         split_part(f[1], ' ', 1) AS protein_id,
         f[4]::BIGINT AS peptide_start, f[5]::BIGINT AS peptide_end,
         f[6]::DOUBLE AS probability
  FROM l
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
