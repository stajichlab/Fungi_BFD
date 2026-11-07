#!/usr/bin/env python3
"""schema_doc.py -- write a Markdown schema document (tables, columns, row
counts, join keys, Mermaid ER diagram) from a built BFD DuckDB file.

Usage:
    schema_doc.py --db db/BFD.duckdb --out docs/BFD_SCHEMA.md

Reads only information_schema and duckdb_tables() (estimated_size is the row
count DuckDB stores for a base table; views show no count). GitHub renders the
```mermaid block as a diagram.
"""
import argparse, csv, datetime, io, os, subprocess, sys

# Join keys. species_prefix in the per-genome tables holds the LOCUSTAG value.
KEYS = ["LOCUSTAG", "species_prefix", "ASMID", "gene_id", "transcript_id", "protein_id"]
# (child column, parent table, parent column, label)
LINKS = [
    ("LOCUSTAG", "species", "LOCUSTAG", "genome"),
    ("species_prefix", "species", "LOCUSTAG", "genome"),
    ("gene_id", "gene_info", "gene_id", "gene"),
    ("transcript_id", "gene_transcripts", "transcript_id", "transcript"),
    ("protein_id", "gene_proteins", "protein_id", "protein"),
]


def q(db, sql):
    r = subprocess.run(["duckdb", "-readonly", db, "-csv", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode != 0 or "Error" in r.stderr:
        sys.exit(f"ERROR: duckdb: {r.stderr.strip()[:1000]}")
    return list(csv.DictReader(io.StringIO(r.stdout)))


def mm_type(t):
    return t.split("(")[0].replace(" ", "_").replace("[]", "_list")


def mm_col(c):
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cols = q(a.db, "SELECT c.table_name, t.table_type, c.column_name, c.data_type "
                   "FROM information_schema.columns c JOIN information_schema.tables t "
                   "USING (table_schema, table_name) WHERE c.table_schema = 'main' "
                   "AND c.table_name NOT LIKE '\\_%' ESCAPE '\\' "
                   "ORDER BY c.table_name, c.ordinal_position")
    rows = {r["table_name"]: int(r["estimated_size"])
            for r in q(a.db, "SELECT table_name, estimated_size FROM duckdb_tables()")}
    tables, kind = {}, {}
    for r in cols:
        tables.setdefault(r["table_name"], []).append((r["column_name"], r["data_type"]))
        kind[r["table_name"]] = "view" if r["table_type"] == "VIEW" else "table"

    def group(t):
        if t in ("species", "asm_stats", "busco_genome") or t.startswith("telomere") \
                or t.endswith("_frequency") or t.startswith("v_telomere"):
            return "Genome level"
        if t.startswith("gene_"):
            return "Gene models"
        if t.startswith("v_"):
            return "Summary views"
        if t.startswith("wgd"):
            return "Paralog Ks"
        return "Protein function"

    order = ["Genome level", "Gene models", "Protein function", "Paralog Ks", "Summary views"]
    stamp = datetime.date.today().isoformat()
    out = [f"# BFD DuckDB schema",
           "",
           f"Generated {stamp} by `nextflow/bin/schema_doc.py` from `{os.path.basename(a.db)}`. "
           f"Do not edit by hand; rerun the script after a rebuild.",
           "",
           f"{sum(1 for t in kind if kind[t] == 'table')} tables and "
           f"{sum(1 for t in kind if kind[t] == 'view')} views. Row counts are DuckDB "
           "`estimated_size` for base tables. Views have no stored count.",
           "",
           "## Join keys",
           "",
           "| Key | Meaning | Parent |",
           "|---|---|---|",
           "| `LOCUSTAG` | genome key | `species` |",
           "| `species_prefix` | the LOCUSTAG of the genome, in per-genome result tables | `species.LOCUSTAG` |",
           "| `gene_id` | gene model | `gene_info` |",
           "| `transcript_id` | mRNA / tRNA transcript | `gene_transcripts` |",
           "| `protein_id` | predicted protein | `gene_proteins` |",
           "",
           "Protein and CDS sequences are not stored in the database. They are in the "
           "bgzip + `samtools faidx` files under `seqdb/` (see `seqdb/genomes.tsv`).",
           "",
           "## Diagram",
           "",
           "Each box shows only the key columns. The full column lists are below.",
           "",
           "```mermaid",
           "erDiagram"]
    for t in sorted(tables):
        if kind[t] == "view" and t.startswith("v_"):
            continue
        keys = [(c, ty) for c, ty in tables[t] if c in KEYS]
        out.append(f"    {t} {{")
        for c, ty in keys:
            out.append(f"        {mm_type(ty)} {mm_col(c)}")
        out.append("    }")
    seen = set()
    for t in sorted(tables):
        if t.startswith("v_"):
            continue
        names = {c for c, _ in tables[t]}
        for child, parent, pcol, label in LINKS:
            if child in names and parent in tables and parent != t and (t, parent) not in seen:
                # a table links to its nearest parent only (protein > transcript > gene > genome)
                if label == "genome" and names & {"protein_id", "transcript_id", "gene_id"} \
                        and t not in ("gene_info", "gene_transcripts", "gene_proteins"):
                    continue
                if label == "gene" and names & {"protein_id", "transcript_id"} and t != "gene_transcripts" \
                        and t != "gene_proteins":
                    continue
                seen.add((t, parent))
                out.append(f"    {parent} ||--o{{ {t} : {label}")
    out += ["```", ""]

    for g in order:
        ts = [t for t in sorted(tables) if group(t) == g]
        if not ts:
            continue
        out += [f"## {g}", ""]
        for t in ts:
            n = rows.get(t)
            out.append(f"### `{t}` ({kind[t]}" + (f", {n:,} rows)" if n is not None else ")"))
            out += ["", "| Column | Type |", "|---|---|"]
            out += [f"| `{c}` | {ty} |" for c, ty in tables[t]]
            out.append("")
    open(a.out, "w").write("\n".join(out))
    print(f"{a.out}: {len(tables)} tables/views")


if __name__ == "__main__":
    main()
