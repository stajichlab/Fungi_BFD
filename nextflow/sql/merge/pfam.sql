-- MERGE_PFAM: hmmscan/hmmsearch --domtblout files -> pfam.parquet
-- Same columns and values as bin/pfamtbl_to_long.py (one row per domain hit).
-- domtblout columns are separated by runs of spaces and the last field (target
-- description) contains spaces, so each line is read whole and split on
-- whitespace; only fields 1-21 are used. Lines starting with '#' are comments.
COPY (
  WITH l AS (
    SELECT regexp_split_to_array(trim(line), '\s+') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = '\t',
                  quote = '', escape = '', auto_detect = false)
    WHERE line NOT LIKE '#%' AND len(trim(line)) > 0
  )
  SELECT f[1] AS protein_id, f[4] AS pfam_id, f[5] AS pfam_acc, f[6]::BIGINT AS pfam_len,
         f[7]::DOUBLE AS full_seq_e_value, f[8]::DOUBLE AS full_seq_score, f[9]::DOUBLE AS full_seq_bias,
         f[10]::BIGINT AS domain_num, f[11]::BIGINT AS domain_num_of,
         f[12]::DOUBLE AS domain_c_evalue, f[13]::DOUBLE AS domain_i_evalue,
         f[14]::DOUBLE AS domain_score, f[15]::DOUBLE AS domain_bias,
         f[16]::BIGINT AS hmm_from, f[17]::BIGINT AS hmm_to, f[18]::BIGINT AS ali_from,
         f[19]::BIGINT AS ali_to, f[20]::BIGINT AS env_from, f[21]::BIGINT AS env_to
  FROM l
  WHERE len(f) >= 2
  -- sorted by protein: lets queries by protein or genome skip row groups. On
  -- 1,018 genomes, 4 threads + this sort took 128 s (Python path: 175 s).
  ORDER BY protein_id
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
