#!/usr/bin/env python3
"""Drop the statements in a SQL file whose targets are not present, so an
incomplete tables/ directory degrades instead of failing the build.

build_BFD_duckDB.sh is written to tolerate a missing table: mkview prints
"SKIP:" and moves on, and generate_bfd_catalog_sql.py emits "-- SKIP:" for the
same reason -- DuckDB binds a view eagerly and errors at CREATE time when the
Parquet behind it is absent. The three SQL files moved in from BFD_dashboard
carry no such guards, because the dashboard only ever ran them against a
complete production tables/ directory. Applied unguarded they make roughly ten
optional tables mandatory:

  * COMMENT ON VIEW <name> raises "Catalog Error: View with name <name> does
    not exist!" -- the duckdb CLI runs the remaining statements and then exits
    1, so under `set -e` the builder dies with the file already half applied.
  * CREATE OR REPLACE VIEW ... read_parquet('<missing>') fails the same way.

A fresh release directory built with params_merge_only.yaml has no
infernal_rfam.parquet (MERGE_INFERNAL is in skip_merges) and no
wgd.ks.nodes.parquet, so without this filter BUILD_DUCKDB dies at the last
step, after hours of merges.

Statements are kept unless we can prove the target is absent: this filter is a
guard against known-optional data, not a validator. Anything it does not
recognise is passed through so DuckDB still reports real errors.

Usage:
    filter_sql_for_db.py --sql FILE --existing-file NAMES [--tables-dir DIR]

NAMES is one catalog name per line (views and tables already in the database),
produced by the caller with the same duckdb binary that will apply the output,
so this script needs no duckdb module and cannot hit a storage-version
mismatch. Kept statements go to stdout; skips are reported on stderr.
"""
import argparse
import glob
import os
import re
import sys


def split_statements(sql):
    """Split on semicolons that are not inside a string literal or a comment.

    Handles '...' with '' escaping, "..." identifiers, -- line comments and
    /* */ block comments. Returns (statement_text, had_terminator) pairs with
    original text preserved, so comments and formatting survive.
    """
    out, buf = [], []
    i, n = 0, len(sql)
    while i < n:
        c = sql[i]
        if c == "'" or c == '"':
            q = c
            buf.append(c)
            i += 1
            while i < n:
                if sql[i] == q:
                    if i + 1 < n and sql[i + 1] == q:   # '' / "" escape
                        buf.append(sql[i:i + 2])
                        i += 2
                        continue
                    buf.append(q)
                    i += 1
                    break
                buf.append(sql[i])
                i += 1
            continue
        if c == "-" and sql.startswith("--", i):
            j = sql.find("\n", i)
            j = n if j == -1 else j
            buf.append(sql[i:j])
            i = j
            continue
        if c == "/" and sql.startswith("/*", i):
            j = sql.find("*/", i + 2)
            j = n if j == -1 else j + 2
            buf.append(sql[i:j])
            i = j
            continue
        if c == ";":
            out.append(("".join(buf) + ";", True))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    tail = "".join(buf)
    if tail:
        # Keep even a blank tail: it is the file's trailing newline, and a
        # filter run that changes nothing must reproduce its input exactly.
        out.append((tail, False))
    return out


def strip_comments(stmt):
    """Statement text with comments blanked, for pattern matching only."""
    out = []
    i, n = 0, len(stmt)
    while i < n:
        if stmt.startswith("--", i):
            j = stmt.find("\n", i)
            i = n if j == -1 else j
            continue
        if stmt.startswith("/*", i):
            j = stmt.find("*/", i + 2)
            i = n if j == -1 else j + 2
            continue
        c = stmt[i]
        if c == "'":
            out.append("''")                      # collapse literals
            i += 1
            while i < n:
                if stmt[i] == "'":
                    if i + 1 < n and stmt[i + 1] == "'":
                        i += 2
                        continue
                    i += 1
                    break
                i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def strip_comments_keep_literals(stmt):
    """Comments removed, string literals intact -- for pulling paths out.

    strip_comments() deliberately collapses literals to '' so a ';' or a
    keyword inside quoted text cannot be mistaken for syntax. Paths live
    inside those literals, so path extraction has to use this instead.
    """
    out = []
    i, n = 0, len(stmt)
    while i < n:
        if stmt.startswith("--", i):
            j = stmt.find("\n", i)
            i = n if j == -1 else j
            continue
        if stmt.startswith("/*", i):
            j = stmt.find("*/", i + 2)
            i = n if j == -1 else j + 2
            continue
        c = stmt[i]
        if c == "'":
            j = i + 1
            while j < n:
                if stmt[j] == "'":
                    if j + 1 < n and stmt[j + 1] == "'":
                        j += 2
                        continue
                    j += 1
                    break
                j += 1
            out.append(stmt[i:j])
            i = j
            continue
        out.append(c)
        i += 1
    return "".join(out)


_COMMENT_ON = re.compile(r"^\s*COMMENT\s+ON\s+(VIEW|TABLE|COLUMN)\s+([A-Za-z_][\w]*)", re.I)
_CREATE = re.compile(r"^\s*CREATE\s+(?:OR\s+REPLACE\s+)?(VIEW|TABLE)\s+([A-Za-z_][\w]*)", re.I)
# a bare identifier source: FROM foo / JOIN foo, not FROM read_parquet(...) or a subquery
_FROM_IDENT = re.compile(r"\b(?:FROM|JOIN)\s+([A-Za-z_][\w]*)\s*(?:AS\s+\w+\s*)?(?=[,);\s]|$)", re.I)
_READ_PATH = re.compile(r"read_(?:parquet|csv|csv_auto|json|json_auto)\(\s*'([^']+)'", re.I)

SQL_KEYWORDS = {"select", "values", "lateral", "unnest", "generate_series", "range"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sql", required=True)
    ap.add_argument("--existing-file", required=True,
                    help="one catalog name per line, from the caller's duckdb")
    ap.add_argument("--tables-dir", default=None,
                    help="only used to report paths relative to it")
    args = ap.parse_args()

    with open(args.existing_file, encoding="utf-8") as fh:
        existing = {ln.strip() for ln in fh if ln.strip()}

    with open(args.sql, encoding="utf-8") as fh:
        text = fh.read()

    kept, skipped = [], []
    # A statement may create a name that a later statement in the same file
    # depends on, so track what this file will add as we go.
    will_exist = set(existing)

    for raw, _term in split_statements(text):
        bare = strip_comments(raw)
        paths_src = strip_comments_keep_literals(raw)
        if not bare.strip():
            kept.append(raw)                       # comment-only chunk
            continue

        reason = None

        m = _COMMENT_ON.match(bare)
        if m and m.group(2) not in will_exist:
            reason = f"{m.group(1).lower()} {m.group(2)} not in the database"

        if reason is None:
            for path in _READ_PATH.findall(paths_src):
                if not os.path.isabs(path):
                    continue
                if glob.glob(path):
                    continue
                reason = f"no file matches {path}"
                break

        if reason is None:
            mc = _CREATE.match(bare)
            if mc and not _READ_PATH.search(paths_src):
                for src in _FROM_IDENT.findall(bare):
                    if src.lower() in SQL_KEYWORDS:
                        continue
                    if src in will_exist:
                        continue
                    reason = f"source {src} not in the database"
                    break

        if reason is None:
            mc = _CREATE.match(bare)
            if mc:
                will_exist.add(mc.group(2))
            kept.append(raw)
        else:
            target = (m.group(2) if m else (mc.group(2) if (mc := _CREATE.match(bare)) else "?"))
            skipped.append((target, reason))

    sys.stdout.write("".join(kept))
    for target, reason in skipped:
        print(f"SKIP: {os.path.basename(args.sql)}: {target} ({reason})", file=sys.stderr)
    print(f"{os.path.basename(args.sql)}: {len(kept)} statements kept, "
          f"{len(skipped)} skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
