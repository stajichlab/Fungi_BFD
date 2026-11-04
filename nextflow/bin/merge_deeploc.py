#!/usr/bin/env python3
"""Merge per-species DeepLoc 2.1 result CSVs into one Parquet table.

Input: *.deeploc_<model>.csv.gz files from RUN_DEEPLOC. Each is the
results_<timestamp>.csv that `deeploc2` writes, gzipped and renamed. Header:

  Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Nucleus,
  Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,
  Lysosome/Vacuole,Golgi apparatus,Peroxisome,Peripheral,Transmembrane,
  Lipid anchor,Soluble

Multi-value fields (Localizations, Signals, Membrane types) are '|'-joined.

DeepLoc writes one row per protein, so the merged table has one row per
protein in BFD (about 200 M rows at full scale). A Python CSV pass would write
a large intermediate file, so this script only builds one DuckDB COPY
statement and runs the duckdb CLI on it. Use --sql-only to print the SQL.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

# (column name in the DeepLoc CSV, column name in the Parquet table)
PROB_COLUMNS = [
    ("Cytoplasm", "p_cytoplasm"),
    ("Nucleus", "p_nucleus"),
    ("Extracellular", "p_extracellular"),
    ("Cell membrane", "p_cell_membrane"),
    ("Mitochondrion", "p_mitochondrion"),
    ("Plastid", "p_plastid"),
    ("Endoplasmic reticulum", "p_endoplasmic_reticulum"),
    ("Lysosome/Vacuole", "p_lysosome_vacuole"),
    ("Golgi apparatus", "p_golgi"),
    ("Peroxisome", "p_peroxisome"),
    ("Peripheral", "p_peripheral"),
    ("Transmembrane", "p_transmembrane"),
    ("Lipid anchor", "p_lipid_anchor"),
    ("Soluble", "p_soluble"),
]
TEXT_COLUMNS = [
    ("Protein_ID", "protein_id"),
    ("Localizations", "localizations"),
    ("Signals", "signals"),
    ("Membrane types", "membrane_types"),
]
EXPECTED_HEADER = [src for src, _ in TEXT_COLUMNS + PROB_COLUMNS]
MODEL_RE = re.compile(r"\.deeploc_(fast|accurate)\.csv\.gz$")


def sql_str(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def build_sql(files: list[str], outfile: str, memory_limit: str = "", threads: int = 0) -> str:
    raw_cols = ", ".join(
        f"{sql_str(dst)}: 'VARCHAR'" for _, dst in TEXT_COLUMNS
    ) + ", " + ", ".join(f"{sql_str(dst)}: 'DOUBLE'" for _, dst in PROB_COLUMNS)
    file_list = ",\n    ".join(sql_str(f) for f in files)
    probs = ",\n    ".join(f"round({dst}, 4)::FLOAT AS {dst}" for _, dst in PROB_COLUMNS)
    settings = "SET temp_directory = 'duckdb_tmp';\n"
    if memory_limit:
        settings += f"SET memory_limit = {sql_str(memory_limit)};\n"
    if threads:
        settings += f"SET threads = {int(threads)};\n"
    return settings + f"""COPY (
SELECT
    split_part(protein_id, '_', 1)               AS species_prefix,
    protein_id,
    regexp_extract(filename, '\\.deeploc_([a-z]+)\\.csv\\.gz$', 1) AS model,
    localizations,
    signals,
    membrane_types,
    coalesce(contains(signals, 'Signal peptide'), FALSE)                AS has_signal_peptide,
    coalesce(contains(signals, 'Mitochondrial transit peptide'), FALSE) AS has_mito_transit_peptide,
    {probs}
FROM read_csv([
    {file_list}
  ],
  header = true, delim = ',', quote = '"', filename = true,
  columns = {{{raw_cols}}})
ORDER BY species_prefix, protein_id
) TO {sql_str(outfile)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""


def check_inputs(files: list[str]) -> None:
    """Fail early on a wrong file name or a changed DeepLoc header."""
    import gzip
    models = set()
    for f in files:
        m = MODEL_RE.search(f)
        if not m:
            sys.exit(f"ERROR: {f}: name does not match *.deeploc_<fast|accurate>.csv.gz")
        models.add(m.group(1))
        with gzip.open(f, "rt") as fh:
            header = fh.readline().rstrip("\n").split(",")
        if header != EXPECTED_HEADER:
            sys.exit(f"ERROR: {f}: unexpected header {header}")
    if len(models) > 1:
        sys.exit(f"ERROR: inputs mix DeepLoc models {sorted(models)}; merge one model at a time")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="+", help="*.deeploc_<model>.csv.gz files")
    ap.add_argument("-o", "--outfile", default="deeploc.parquet")
    ap.add_argument("--duckdb", default="duckdb", help="duckdb CLI to run")
    ap.add_argument("--memory-limit", default="",
                    help="DuckDB memory_limit, e.g. 20GB (set it below the job limit)")
    ap.add_argument("--threads", type=int, default=0, help="DuckDB threads (0 = default)")
    ap.add_argument("--sql-only", action="store_true", help="print the SQL and exit")
    args = ap.parse_args()

    files = sorted(str(Path(f)) for f in args.inputs)
    check_inputs(files)
    sql = build_sql(files, args.outfile, args.memory_limit, args.threads)
    if args.sql_only:
        print(sql)
        return
    # The SQL can be ~1 MB for 18k inputs; pass it on stdin, not argv.
    subprocess.run([args.duckdb], input=sql, text=True, check=True)
    print(f"Written: {args.outfile} from {len(files)} files", file=sys.stderr)


if __name__ == "__main__":
    main()
