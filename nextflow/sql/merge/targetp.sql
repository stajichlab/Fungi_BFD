-- MERGE_TARGETP: <LOCUSTAG>_summary.targetp2.gz files -> targetP.parquet
-- Same rows, columns and types as bin/merge_targetp.py followed by
-- read_csv_auto(sample_size=-1):
--   * '#' lines and lines with fewer than 3 TAB fields are skipped;
--   * noTP rows are dropped;
--   * probability = SP column for SP, mTP column for mTP, else NULL;
--   * the CS Position text is parsed with the same regex; when it does not
--     match, start/end = 0, probability = 0.0 and motif = NULL.
-- Python strips leading/trailing whitespace from each line before splitting,
-- so a trailing empty CS field is absent, not ''.
COPY (
  WITH l AS (
    SELECT string_split(regexp_replace(line, '^\s+|\s+$', '', 'g'), E'\t') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = E'\x01',
                  quote = '', escape = '', auto_detect = false)
    WHERE line IS NOT NULL AND line NOT LIKE '#%'
  ), p AS (
    SELECT f[1] AS pid, f[2] AS prediction,
           CASE f[2] WHEN 'SP' THEN coalesce(f[4], '') WHEN 'mTP' THEN coalesce(f[5], '') ELSE '' END AS prob,
           coalesce(f[6], '') AS cs
    FROM l
    WHERE len(f) >= 3 AND f[2] <> 'noTP'
  ), m AS (
    SELECT *, regexp_matches(cs, '^CS pos:\s+(\d+)-(\d+)\.\s+(\S+)\.\s+Pr:\s+(\S+)') AS hit,
           regexp_extract(cs, '^CS pos:\s+(\d+)-(\d+)\.\s+(\S+)\.\s+Pr:\s+(\S+)',
                          ['s', 'e', 'motif', 'pr']) AS g
    FROM p
  )
  SELECT NULLIF(split_part(pid, '_', 1), '') AS species_prefix,
         NULLIF(pid, '') AS protein_id,
         NULLIF(prediction, '') AS prediction,
         NULLIF(prob, '')::DOUBLE AS probability,
         (CASE WHEN hit THEN g.s ELSE '0' END)::BIGINT AS cleavage_position_start,
         (CASE WHEN hit THEN g.e ELSE '0' END)::BIGINT AS cleavage_position_end,
         (CASE WHEN hit THEN g.pr ELSE '0.0' END)::DOUBLE AS cleavage_probability,
         NULLIF(CASE WHEN hit THEN g.motif ELSE '' END, '') AS motif
  FROM m
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
