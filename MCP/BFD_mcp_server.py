#!/usr/bin/env python3
"""
MCP server for the BFD protein/genome functional database.

Exposes read-only tools for agentic investigation of protein composition,
domain architecture, and genome characteristics stored in DuckDB.

Usage:
    python BFD_mcp_server.py
    PROTEIN_DB_PATH=/path/to/BFD.duckdb python BFD_mcp_server.py

Configure in Claude Desktop (~/.claude/claude_desktop_config.json) or
claude-code (.mcp.json) as an stdio server.
"""

import json
import os
import re
import sys
from contextlib import asynccontextmanager
from typing import Any, Optional

import duckdb
from mcp.server.fastmcp import Context, FastMCP
from pydantic import BaseModel, ConfigDict, Field, field_validator

# ── Constants ─────────────────────────────────────────────────────────────────

DB_PATH = os.environ.get("PROTEIN_DB_PATH", "db/BFD.duckdb")
DEFAULT_LIMIT = 1000
MAX_LIMIT = 5000
# Keywords whose presence indicates a mutating statement.
_WRITE_PATTERN = re.compile(
    r"\b(INSERT|UPDATE|DELETE|CREATE|DROP|ALTER|TRUNCATE|REPLACE|COPY|EXPORT"
    r"|PRAGMA|ATTACH|DETACH|VACUUM|CHECKPOINT|LOAD)\b",
    re.IGNORECASE,
)

# ── Lifespan: open a single read-only DuckDB connection ───────────────────────


@asynccontextmanager
async def _lifespan(server: FastMCP):
    if not os.path.exists(DB_PATH):
        print(
            f"ERROR: DuckDB not found at '{DB_PATH}'. "
            "Set PROTEIN_DB_PATH or run scripts/build_duckdb.sh first.",
            file=sys.stderr,
        )
        yield {"db": None}
        return

    conn = duckdb.connect(DB_PATH, read_only=True)
    yield {"db": conn}
    conn.close()


mcp = FastMCP("BFD_db_mcp", lifespan=_lifespan)


# ── Helpers ───────────────────────────────────────────────────────────────────


def _get_conn(ctx) -> duckdb.DuckDBPyConnection:
    conn = ctx.request_context.lifespan_context.get("db")
    if conn is None:
        raise RuntimeError(
            f"Database not available. Check that '{DB_PATH}' exists "
            "and run scripts/build_duckdb.sh if needed."
        )
    return conn


def _reject_write(sql: str) -> Optional[str]:
    """Return an error string if the SQL contains mutating keywords, else None."""
    m = _WRITE_PATTERN.search(sql)
    if m:
        return f"Error: '{m.group()}' is not allowed — this server is read-only."
    return None


def _rows_to_json(cursor: duckdb.DuckDBPyConnection, limit: int) -> dict[str, Any]:
    cols = [d[0] for d in cursor.description]
    rows = cursor.fetchmany(limit)
    truncated = cursor.fetchone() is not None  # peek for more rows
    return {
        "columns": cols,
        "rows": [dict(zip(cols, row)) for row in rows],
        "row_count": len(rows),
        "truncated": truncated,
        "note": f"Results capped at {limit}. Use LIMIT/OFFSET in SQL for pagination."
        if truncated
        else None,
    }


# ── Input models ──────────────────────────────────────────────────────────────


class SqlInput(BaseModel):
    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    sql: str = Field(
        ...,
        description="Read-only SQL SELECT statement. "
        "Do NOT include LIMIT in your SQL — use the `limit` parameter instead.",
        min_length=6,
    )
    limit: int = Field(
        DEFAULT_LIMIT,
        description=f"Maximum rows to return (1–{MAX_LIMIT}). Default {DEFAULT_LIMIT}.",
        ge=1,
        le=MAX_LIMIT,
    )

    @field_validator("sql")
    @classmethod
    def strip_trailing_semicolon(cls, v: str) -> str:
        return v.rstrip(";").strip()


class TableNameInput(BaseModel):
    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    table_name: str = Field(
        ...,
        description="Exact table or view name (from db_list_tables).",
        min_length=1,
        max_length=100,
        pattern=r"^[A-Za-z_][A-Za-z0-9_]*$",
    )


class SampleInput(TableNameInput):
    n: int = Field(5, description="Number of rows to preview.", ge=1, le=50)


class CorrelationInput(BaseModel):
    model_config = ConfigDict(str_strip_whitespace=True, extra="forbid")

    x_column: str = Field(
        ...,
        description="Numeric column for X axis. Must exist in the FROM table or join.",
        min_length=1,
        max_length=100,
    )
    y_column: str = Field(
        ...,
        description="Numeric column for Y axis. Must exist in the FROM table or join.",
        min_length=1,
        max_length=100,
    )
    from_sql: str = Field(
        ...,
        description="SQL fragment after FROM that produces the x_column and y_column. "
        "Example: 'v_species_summary' or "
        "'v_species_summary WHERE PHYLUM = ''Ascomycota'''",
    )
    group_by: Optional[str] = Field(
        None,
        description="Optional column to stratify by (e.g. 'PHYLUM'). "
        "Returns one Pearson r per group.",
        max_length=100,
    )


# ── Tools ─────────────────────────────────────────────────────────────────────


@mcp.tool(
    name="db_execute_sql",
    annotations={
        "title": "Execute read-only SQL",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_execute_sql(params: SqlInput, ctx: Context) -> str:
    """Execute a read-only SQL SELECT against the protein/genome DuckDB.

    This is the primary exploration tool. Use it to query any table or view,
    compute aggregates, join across tables, or filter by taxonomy.

    The database contains 22,412 fungal genomes with protein annotations.
    Key join keys:
      - species.LOCUSTAG = asm_stats.LOCUSTAG = gene_proteins.species_prefix
        = pfam.species_prefix = cazy.species_prefix = idp_summary.species_prefix ...
      - gene_proteins.protein_id = pfam.protein_id = signalp.protein_id ...

    Pre-built views for common cross-table joins:
      - v_species_summary  — one row per species: genome stats + feature counts
      - v_protein_annotation — one row per protein: all annotation flags

    Args:
        params.sql   (str): SELECT statement (no trailing semicolon, no LIMIT clause).
        params.limit (int): Max rows returned (default 1000, max 5000).

    Returns:
        JSON with keys: columns, rows, row_count, truncated, note.

    Examples:
        - "How many species per phylum?"
          sql="SELECT PHYLUM, count(*) AS n FROM species GROUP BY PHYLUM ORDER BY n DESC"
        - "Top 10 species by CAZyme count"
          sql="SELECT SPECIES, cazy_proteins FROM v_species_summary ORDER BY cazy_proteins DESC"
          limit=10
        - "IDP fraction vs genome size for Ascomycota"
          sql="SELECT TOTAL_LENGTH, idp_high_disorder*1.0/gene_count AS idp_rate
               FROM v_species_summary WHERE PHYLUM='Ascomycota'"
    """
    err = _reject_write(params.sql)
    if err:
        return err

    wrapped = f"SELECT * FROM ({params.sql}) __q LIMIT {params.limit + 1}"
    try:
        conn = _get_conn(ctx)
        cur = conn.execute(wrapped)
        result = _rows_to_json(cur, params.limit)
        return json.dumps(result, default=str, indent=2)
    except duckdb.Error as e:
        return f"Error: DuckDB query failed — {e}\nSQL attempted:\n{params.sql}"
    except RuntimeError as e:
        return f"Error: {e}"


@mcp.tool(
    name="db_list_tables",
    annotations={
        "title": "List all tables and views",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_list_tables(ctx: Context) -> str:
    """List every table and view in the protein database with estimated row counts.

    Returns a JSON array. Each entry has:
      - name        (str)  : table or view name
      - type        (str)  : "table" or "view"
      - row_count   (int)  : estimated rows (NULL for views — query to count)
      - description (str)  : brief annotation

    Call db_describe_table(name) for column details, or db_schema_overview()
    for the full relationship map.

    Returns:
        JSON array of {name, type, row_count, description}.
    """
    TABLE_DOCS: dict[str, str] = {
        "species": "One row per genome: taxonomy (PHYLUM→SPECIES), LOCUSTAG, ASMID, BUSCO_LINEAGE",
        "asm_stats": "Assembly statistics per genome: TOTAL_LENGTH, GC_PERCENT, N50_bp, contig_count, masked_pct",
        "gene_info": "Gene coordinates: chrom, start, end, strand, gene_type",
        "gene_proteins": "Protein catalogue: protein_id, transcript_id, length, peptide, md5checksum",
        "gene_transcripts": "Transcript coordinates and attributes",
        "gene_exons": "Per-exon coordinates",
        "gene_CDS": "Per-CDS coordinates and GC content",
        "gene_introns": "Intron coordinates and attributes",
        "gene_trna": "tRNA gene annotations",
        "gene_intergenic_distances": "Pairwise intergenic distances between adjacent genes (38M rows)",
        "pfam": "Pfam domain hits: protein_id, pfam_id, pfam_acc, e-values, alignment coords",
        "cazy": "CAZyme HMM hits: HMM_id, protein_id, evalue, coverage",
        "cazy_overview": "CAZyme summary per protein: family, EC number, substrate",
        "merops": "MEROPS peptidase hits: merops_id, percent_identity, evalue",
        "signalp": "Signal peptide predictions: protein_id, cleavage site",
        "targetp": "TargetP subcellular targeting: prediction (SP/mTP/CS), probability",
        "wolfpsort": "WoLF PSORT localisation: localization, score",
        "predgpi": "GPI-anchor predictions: protein_id, chain coordinates",
        "tmhmm": "TMHMM transmembrane helices: PredHel count, topology",
        "idp_summary": "Intrinsic disorder per protein: IDP_fraction, IDP_residues",
        "idp": "Per-region disorder coordinates: IDP_start, IDP_end",
        "aa_frequency": "Amino-acid usage frequency per species",
        "codon_frequency": "Codon usage frequency per species",
        "v_species_summary": "VIEW — per-species genome stats + functional feature counts (recommended starting point)",
        "v_protein_annotation": "VIEW — per-protein annotation flags: domains, secretome, disorder, localisation",
    }

    try:
        conn = _get_conn(ctx)

        tables = conn.execute(
            "SELECT table_name, estimated_size FROM duckdb_tables() ORDER BY table_name"
        ).fetchall()
        views = conn.execute(
            "SELECT view_name FROM duckdb_views() ORDER BY view_name"
        ).fetchall()

        result = []
        for name, est in tables:
            result.append(
                {
                    "name": name,
                    "type": "table",
                    "row_count": est,
                    "description": TABLE_DOCS.get(name, ""),
                }
            )
        for (name,) in views:
            result.append(
                {
                    "name": name,
                    "type": "view",
                    "row_count": None,
                    "description": TABLE_DOCS.get(name, ""),
                }
            )

        return json.dumps(result, indent=2)
    except (duckdb.Error, RuntimeError) as e:
        return f"Error: {e}"


@mcp.tool(
    name="db_describe_table",
    annotations={
        "title": "Describe table columns",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_describe_table(params: TableNameInput, ctx: Context) -> str:
    """Return column names, types, and null counts for a table or view.

    Also returns a sample of distinct values for low-cardinality VARCHAR columns
    to help formulate WHERE clauses (e.g. PHYLUM values, prediction categories).

    Args:
        params.table_name (str): Exact table or view name from db_list_tables.

    Returns:
        JSON with keys:
          - table_name  (str)
          - columns     (list of {name, type, null_count, sample_values?})
          - row_count   (int, for tables; null for views)
    """
    try:
        conn = _get_conn(ctx)

        # Validate name exists to give a clean error
        names = {
            r[0]
            for r in conn.execute(
                "SELECT table_name FROM duckdb_tables() "
                "UNION SELECT view_name FROM duckdb_views()"
            ).fetchall()
        }
        if params.table_name not in names:
            available = sorted(names)
            return (
                f"Error: '{params.table_name}' not found. "
                f"Available: {available}"
            )

        desc = conn.execute(f"DESCRIBE {params.table_name}").fetchall()
        col_info = []
        for row in desc:
            col_name, col_type = row[0], row[1]
            entry: dict[str, Any] = {"name": col_name, "type": col_type}

            # For low-cardinality text columns, fetch distinct values
            if col_type in ("VARCHAR", "TEXT") and col_name.upper() in (
                "PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY",
                "BUSCO_LINEAGE", "PREDICTION", "LOCALIZATION", "GENE_TYPE",
                "TARGETP_CLASS",
            ):
                # Always quote the column name to handle reserved words like ORDER
                quoted = f'"{col_name}"'
                try:
                    vals = conn.execute(
                        f"SELECT DISTINCT {quoted} FROM {params.table_name} "
                        f"WHERE {quoted} IS NOT NULL ORDER BY 1 LIMIT 30"
                    ).fetchall()
                    entry["sample_values"] = [v[0] for v in vals]
                except duckdb.Error:
                    pass

            col_info.append(entry)

        # Row count for tables only
        row_count = None
        try:
            rc = conn.execute(
                f"SELECT estimated_size FROM duckdb_tables() "
                f"WHERE table_name = '{params.table_name}'"
            ).fetchone()
            if rc:
                row_count = rc[0]
        except duckdb.Error:
            pass

        return json.dumps(
            {"table_name": params.table_name, "row_count": row_count, "columns": col_info},
            indent=2,
        )
    except (duckdb.Error, RuntimeError) as e:
        return f"Error: {e}"


@mcp.tool(
    name="db_sample_rows",
    annotations={
        "title": "Preview table rows",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_sample_rows(params: SampleInput, ctx: Context) -> str:
    """Return the first N rows of a table or view.

    Useful for inspecting data format and column contents before writing
    more complex SQL. For filtered previews, use db_execute_sql instead.

    Args:
        params.table_name (str): Exact table or view name.
        params.n          (int): Rows to return (default 5, max 50).

    Returns:
        JSON with keys: columns, rows, row_count.
    """
    try:
        conn = _get_conn(ctx)
        cur = conn.execute(f"SELECT * FROM {params.table_name} LIMIT {params.n}")
        cols = [d[0] for d in cur.description]
        rows = cur.fetchall()
        return json.dumps(
            {
                "columns": cols,
                "rows": [dict(zip(cols, r)) for r in rows],
                "row_count": len(rows),
            },
            default=str,
            indent=2,
        )
    except duckdb.Error as e:
        return f"Error: {e}"
    except RuntimeError as e:
        return f"Error: {e}"


@mcp.tool(
    name="db_pearson_correlation",
    annotations={
        "title": "Pearson correlation between two columns",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_pearson_correlation(params: CorrelationInput, ctx: Context) -> str:
    """Compute Pearson r between two numeric columns, optionally stratified.

    Uses DuckDB's built-in corr() aggregate. Returns r, n, and p-value
    approximation (t-test on r√(n-2)/√(1-r²), two-tailed).

    Args:
        params.x_column (str): Numeric X column name.
        params.y_column (str): Numeric Y column name.
        params.from_sql (str): SQL fragment after FROM providing both columns.
                               Can include WHERE clauses.
                               Example: "v_species_summary WHERE PHYLUM='Ascomycota'"
        params.group_by (str, optional): Column to stratify by.

    Returns:
        JSON with {correlation_r, n, approx_p_value} or list thereof when group_by set.

    Examples:
        x_column="GC_PERCENT", y_column="pfam_proteins", from_sql="v_species_summary"
        → Pearson r between GC% and Pfam domain count across all species

        x_column="TOTAL_LENGTH", y_column="idp_high_disorder",
        from_sql="v_species_summary", group_by="PHYLUM"
        → r per phylum between genome size and high-disorder protein count
    """
    err = _reject_write(params.from_sql)
    if err:
        return err

    try:
        conn = _get_conn(ctx)

        if params.group_by:
            sql = (
                f"SELECT {params.group_by} AS group_val, "
                f"  corr({params.x_column}, {params.y_column}) AS r, "
                f"  count(*) AS n "
                f"FROM {params.from_sql} "
                f"WHERE {params.x_column} IS NOT NULL "
                f"  AND {params.y_column} IS NOT NULL "
                f"GROUP BY {params.group_by} "
                f"ORDER BY abs(r) DESC NULLS LAST "
                f"LIMIT 100"
            )
        else:
            sql = (
                f"SELECT "
                f"  corr({params.x_column}, {params.y_column}) AS r, "
                f"  count(*) AS n "
                f"FROM {params.from_sql} "
                f"WHERE {params.x_column} IS NOT NULL "
                f"  AND {params.y_column} IS NOT NULL"
            )

        rows = conn.execute(sql).fetchall()

        def _p_approx(r: Optional[float], n: int) -> Optional[float]:
            """Two-tailed p-value approximation via t-distribution (df=n-2)."""
            if r is None or n < 3:
                return None
            import math
            t = r * math.sqrt(n - 2) / math.sqrt(max(1e-12, 1 - r * r))
            # Abramowitz & Stegun approximation for two-tailed p
            x = abs(t)
            df = n - 2
            # Use regularized incomplete beta function approximation
            try:
                from scipy.stats import t as tdist  # type: ignore
                return float(2 * tdist.sf(x, df))
            except ImportError:
                # Rough approximation without scipy
                z = x / math.sqrt(1 + x * x / df)
                p = math.exp(-0.717 * z - 0.416 * z * z)
                return round(min(1.0, 2 * p), 6)

        if params.group_by:
            result = []
            for row in rows:
                group_val, r, n = row
                result.append(
                    {
                        params.group_by: group_val,
                        "r": round(r, 6) if r is not None else None,
                        "n": n,
                        "approx_p_value": _p_approx(r, n),
                    }
                )
            return json.dumps(
                {
                    "x": params.x_column,
                    "y": params.y_column,
                    "stratified_by": params.group_by,
                    "results": result,
                },
                indent=2,
            )
        else:
            r_val, n_val = rows[0]
            return json.dumps(
                {
                    "x": params.x_column,
                    "y": params.y_column,
                    "r": round(r_val, 6) if r_val is not None else None,
                    "n": n_val,
                    "approx_p_value": _p_approx(r_val, n_val),
                    "from": params.from_sql,
                },
                indent=2,
            )
    except duckdb.Error as e:
        return f"Error: DuckDB query failed — {e}"
    except RuntimeError as e:
        return f"Error: {e}"


@mcp.tool(
    name="db_schema_overview",
    annotations={
        "title": "Database schema and join key overview",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def db_schema_overview(ctx: Context) -> str:
    """Return a concise schema overview including join keys and table relationships.

    Use this first to orient yourself before writing SQL. Returns:
      - Join key map (how tables connect to each other)
      - All table and view names grouped by topic
      - Notes on large tables that need filtering

    Returns:
        Markdown-formatted schema reference string.
    """
    overview = """# BFD Protein/Genome Database Schema

## Scale
- 22,412 fungal genomes across Ascomycota, Basidiomycota, and other phyla
- ~1.2M proteins, ~1.2M Pfam domain hits, ~38M intergenic distance records

## Primary Join Keys

### Genome-level joins (species ↔ assembly ↔ composition)
```
species.LOCUSTAG = asm_stats.LOCUSTAG
species.ASMID    = asm_stats.ASMID
species.LOCUSTAG = aa_frequency.species_prefix
species.LOCUSTAG = codon_frequency.species_prefix
```

### Species → Protein (most functional tables)
```
species.LOCUSTAG  = gene_proteins.species_prefix
                  = pfam.species_prefix
                  = cazy.species_prefix / cazy_overview.species_prefix
                  = merops.species_prefix
                  = signalp.species_prefix
                  = targetp.species_prefix
                  = wolfpsort.species_prefix
                  = predgpi.species_prefix
                  = tmhmm.species_prefix
                  = idp_summary.species_prefix / idp.species_prefix
                  = gene_intergenic_distances.species_prefix
```

### Protein → annotation (protein_id key)
```
gene_proteins.protein_id = pfam.protein_id
                         = cazy_overview.protein_id
                         = merops.protein_id
                         = signalp.protein_id
                         = targetp.protein_id
                         = wolfpsort.protein_id
                         = predgpi.protein_id
                         = tmhmm.protein_id
                         = idp_summary.protein_id
```

### Gene structure
```
gene_info.gene_id      = gene_proteins.gene_id = gene_transcripts.gene_id
gene_transcripts.transcript_id = gene_exons.transcript_id
                               = gene_introns.transcript_id
                               = gene_CDS.transcript_id
```

## Pre-built Views (use these for cross-table analysis)

| View | Description |
|------|-------------|
| v_species_summary | One row per species: genome stats + pfam/cazy/merops/signalp/tmhmm/idp counts |
| v_protein_annotation | One row per protein: all annotation flags joined |

## Tables by Topic

### Taxonomy & Assembly
- `species` — LOCUSTAG, ASMID, PHYLUM, CLASS, "ORDER", FAMILY, GENUS, SPECIES
- `asm_stats` — TOTAL_LENGTH, GC_PERCENT, N50_bp, contig_count, masked_pct

### Gene Structure
- `gene_info`, `gene_proteins`, `gene_transcripts`, `gene_exons`, `gene_CDS`, `gene_introns`, `gene_trna`
- `gene_intergenic_distances` ⚠️ 38M rows — always filter by species_prefix

### Protein Domains
- `pfam` — 1.2M hits; pfam_id, pfam_acc, e-values, alignment coords
- `cazy` — CAZyme HMM hits; HMM_id (e.g. "GH16_1", "GT25")
- `cazy_overview` — per-protein CAZyme summary with EC and substrate
- `merops` — peptidase hits; merops_id

### Secretome / Targeting
- `signalp` — classical secretion signal peptides
- `targetp` — prediction: SP (secretory), mTP (mitochondrial), CS (chloroplast)
- `wolfpsort` — predicted localisation: cyto, nucl, plas, E.R., golgi, mito, pero, vacu, extr
- `predgpi` — GPI-anchor chain predictions
- `tmhmm` — transmembrane helix count (PredHel)

### Intrinsic Disorder
- `idp_summary` — IDP_fraction, IDP_residues per protein
- `idp` — per-region disorder coordinates

### Composition
- `aa_frequency` — per-species amino-acid usage (amino_acid, frequency)
- `codon_frequency` — per-species codon usage (codon, frequency)

## Taxonomy Hierarchy
PHYLUM → SUBPHYLUM → CLASS → ORDER → FAMILY → GENUS → SPECIES

Common PHYLUMs: Ascomycota, Basidiomycota, Mucoromycota, Chytridiomycota
Use db_execute_sql to query `SELECT DISTINCT PHYLUM FROM species` for the full list.
"""
    return overview


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run()
