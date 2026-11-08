-- MERGE_PREDGPI: <LOCUSTAG>.predgpi.gff3.gz files -> predgpi.parquet
-- Same rows, columns and types as bin/merge_predgpi.py followed by
-- read_csv_auto(sample_size=-1). Lines starting with '#' and empty lines are
-- skipped. The 9 GFF3 columns are kept; species_prefix is the part of the
-- protein ID before the first '_'.
COPY (
  WITH l AS (
    SELECT string_split(line, E'\t') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = E'\x01',
                  quote = '', escape = '', auto_detect = false)
    WHERE line IS NOT NULL AND line <> '' AND line NOT LIKE '#%'
  )
  SELECT NULLIF(split_part(f[1], '_', 1), '') AS species_prefix,
         NULLIF(f[1], '') AS protein_id,
         NULLIF(f[2], '') AS source,
         NULLIF(f[3], '') AS feature,
         NULLIF(f[4], '')::BIGINT AS "start",
         NULLIF(f[5], '')::BIGINT AS "end",
         NULLIF(f[6], '')::DOUBLE AS score,
         NULLIF(f[7], '') AS strand,
         NULLIF(f[8], '') AS phase,
         NULLIF(f[9], '') AS attributes
  FROM l
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
