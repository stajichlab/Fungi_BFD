#!/usr/bin/env python3
"""
build_runtime_metadata_table.py — build a durable, appendable table of
per-task funannotate pipeline runtimes joined with genome size / repeat
content, for plotting runtime vs. assembly characteristics.

Data sources (durable, on-disk; no dependency on .nextflow/cache retention
or sacct accounting windows):
  1. logs/nextflow/funannotate_trace*.txt — Nextflow -with-trace files,
     already written on every run of this pipeline. Each row is one task
     attempt (hash), with realtime/cpu/rss. The same hash can reappear in
     later trace files as CACHED (a resumed run reusing that result); those
     are deduped down to a single observation per hash.
  2. results/genome_stats_by_name/<Genus>/<species>.asm_stats.stats.txt —
     PRIMARY genome-metadata source. This is a shared, curated per-species
     assembly-stats table (symlinked identically into every project root),
     keyed by species name rather than by Nextflow work dir, so it survives
     work dir cleanup. Covers ~75% of species with a FUNANNOTATE_PREDICT
     observation. Joined onto every process row for a given species
     (FUNANNOTATE_PREDICT, FUNANNOTATE_TRAIN, GENEMARK_RUN, RNASEQ_PREPARE,
     ...), not just PREDICT/TRAIN, since it's genome-level metadata.
  3. work/<subdir>/<hash>/.command.log — FALLBACK, only used for
     FUNANNOTATE_PREDICT/FUNANNOTATE_TRAIN hashes whose species has no
     asm_stats file (e.g. a species run through the pipeline but not yet
     added to the shared genome_stats DB). The funannotate stdout log
     contains a "Genome loaded: N scaffolds; M bp; P% repeats masked" line.
     Best-effort: work dirs get cleaned up over time, so this fallback finds
     little for old hashes once their work dir is gone.

A distinct hash = a distinct real execution attempt, so the same assembly
run through the same workflow step multiple times (e.g. a failed attempt
that got retried, or a repredict after retraining) contributes multiple
rows — this is intentional; the table is observations, not a per-assembly
summary.

Output: an Apache Parquet table (plus optional CSV) that can be re-run
against new trace files over time; existing rows are kept and new hashes
are appended (dedup key: hash, since a hash is intrinsically unique to one
task attempt regardless of which root/species it was; root+species+process
are carried along as columns for grouping/plotting).

Usage:
  build_runtime_metadata_table.py \
      --roots do_annotation do_annotation_asco do_annotation_Rhod \
      --outdir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis \
      --also-csv
"""

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

TRACE_GLOB_PATTERNS = [
    'logs/nextflow/funannotate_trace*.txt',
]
EXPECTED_HEADER = {'task_id', 'hash', 'name', 'status', 'exit', 'realtime',
                   '%cpu', 'rss', 'tag'}

TS_RE = re.compile(r'(\d{4})-(\d{2})-(\d{2})_(\d{2})_(\d{2})_(\d{2})')
REALTIME_RE = re.compile(r'([\d.]+)\s*(ms|s|m|h|d)')
PCT_RE = re.compile(r'([\d.]+)\s*%')
RSS_RE = re.compile(r'([\d.]+)\s*(B|KB|MB|GB|TB)')
SKIP_EXIT_THRESHOLD_S = 10.0
# Processes that wrap a funannotate-family tool with its own existing-output
# idempotency check: the tool checks the target folder before doing any real
# work and exits immediately (exit 0, COMPLETED at the Nextflow level, since
# it's a *fresh* task hash rather than a Nextflow-level CACHED reuse) if the
# expected output is already there. Confirmed 2026-09-11 (project owner) for
# FUNANNOTATE_PREDICT after the runtime-scaling analysis flagged a cluster of
# ~1.4-2.2s "COMPLETED" PREDICT attempts with large-genome metadata but no
# real work done. Applied here to the other funannotate-family processes too
# since they wrap tools with the same output-already-exists convention.
SKIP_EXIT_PROCESSES = {'FUNANNOTATE_PREDICT', 'FUNANNOTATE_PREDICT_SIB',
                        'FUNANNOTATE_TRAIN', 'GENEMARK_RUN', 'GENEMARK_RUN_SIB'}

GENOME_STATS_RE = re.compile(
    r'Genome loaded:\s*([\d,]+)\s*scaffolds;\s*([\d,]+)\s*bp;\s*([\d.]+)%\s*repeats masked')

# Only these processes' .command.log carry the funannotate "Genome loaded"
# banner; the command.log fallback lookup is skipped entirely for everything
# else. The asm_stats.stats.txt join (primary source) applies to ALL
# processes with a species tag, since it's genome-level metadata.
COMMAND_LOG_FALLBACK_PROCESSES = {'FUNANNOTATE_PREDICT', 'FUNANNOTATE_TRAIN'}

ASM_STATS_FIELD_RE = {
    'num_scaffolds': re.compile(r'CONTIG COUNT\s*=\s*([\d,]+)'),
    'genome_size_bp': re.compile(r'TOTAL LENGTH\s*=\s*([\d,]+)'),
    'gc_pct': re.compile(r'GC%\s*=\s*([\d.]+)'),
    'n50': re.compile(r'^\s*N50\s*=\s*([\d,]+)', re.MULTILINE),
    'repeat_masked_pct': re.compile(r'PERCENT MASKED\s*=\s*([\d.]+)'),
}


def parse_realtime_seconds(s):
    if not s or s == '-':
        return None
    matches = REALTIME_RE.findall(s)
    if not matches:
        return None
    units = {'ms': 1 / 1000, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}
    return sum(float(val) * units[unit] for val, unit in matches)


def parse_pct(s):
    if not s or s == '-':
        return None
    m = PCT_RE.search(s)
    return float(m.group(1)) if m else None


def parse_rss_mb(s):
    if not s or s == '-':
        return None
    m = RSS_RE.search(s)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2)
    factor = {'B': 1 / (1024 * 1024), 'KB': 1 / 1024, 'MB': 1, 'GB': 1024, 'TB': 1024 * 1024}
    return val * factor[unit]


def trace_datetime_from_name(path):
    m = TS_RE.search(str(path))
    if not m:
        return None
    try:
        return datetime(*(int(g) for g in m.groups()))
    except ValueError:
        return None


def process_basename(name_full):
    """FUNANNOTATE:FUNANNOTATE_PREDICTION:FUNANNOTATE_PREDICT (Sp) -> FUNANNOTATE_PREDICT"""
    name_no_tag = re.sub(r'\s*\(.*$', '', name_full).strip()
    return name_no_tag.split(':')[-1].strip()


def species_from_name(name_full, tag):
    if tag and tag != '-':
        return tag
    m = re.search(r'\((.*)\)\s*$', name_full)
    return m.group(1) if m else None


def collect_trace_files(roots):
    files = []
    for root in roots:
        root_path = Path(root)
        if not root_path.is_dir():
            print(f"Warning: root not found, skipping: {root}", file=sys.stderr)
            continue
        for pattern in TRACE_GLOB_PATTERNS:
            for p in root_path.glob(pattern):
                if p.is_file():
                    files.append((root_path, p))
    return sorted(set(files), key=lambda rp: rp[1])


def load_trace_rows(trace_files):
    """Long-form rows, one per (root, hash) task attempt, deduped by first
    chronological sighting of that hash within a root."""
    rows_by_key = {}
    files_parsed = 0
    for root_path, trace_path in trace_files:
        trace_dt = trace_datetime_from_name(trace_path)
        try:
            with open(trace_path) as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if reader.fieldnames is None or not EXPECTED_HEADER.issubset(set(reader.fieldnames)):
                    print(f"Warning: unexpected header, skipping: {trace_path}", file=sys.stderr)
                    continue
                for row in reader:
                    h = row.get('hash', '')
                    if not h:
                        continue
                    key = (str(root_path), h)
                    if key in rows_by_key:
                        continue  # already have this task attempt from an earlier trace file
                    name_full = row.get('name', '')
                    proc = process_basename(name_full)
                    if not proc:
                        continue
                    rt = parse_realtime_seconds(row.get('realtime', ''))
                    rows_by_key[key] = {
                        'root': str(root_path),
                        'hash': h,
                        'process': proc,
                        'species': species_from_name(name_full, row.get('tag')),
                        'status': row.get('status', ''),
                        'exit_code': row.get('exit', ''),
                        'realtime_s': rt,
                        'cpu_pct': parse_pct(row.get('%cpu', '')),
                        'rss_mb': parse_rss_mb(row.get('rss', '')),
                        'first_seen_trace': trace_path.name,
                        'first_seen_date': trace_dt.date().isoformat() if trace_dt else None,
                    }
            files_parsed += 1
        except Exception as e:
            print(f"Warning: failed to parse {trace_path}: {e}", file=sys.stderr)
    print(f"Parsed {files_parsed}/{len(trace_files)} trace files -> "
          f"{len(rows_by_key)} unique task attempts", file=sys.stderr)
    return list(rows_by_key.values())


def find_workdir(root, hash_):
    """work/<subdir>/<2-char>/<rest>* — <subdir> varies (funannotate/, ANI/, ...)."""
    if '/' not in hash_:
        return None
    prefix, rest = hash_.split('/', 1)
    work_root = Path(root) / 'work'
    if not work_root.is_dir():
        return None
    for subdir in work_root.iterdir():
        candidate_parent = subdir / prefix
        if not candidate_parent.is_dir():
            continue
        matches = sorted(candidate_parent.glob(f'{rest}*'))
        if matches:
            return matches[0]
    return None


def extract_genome_stats(workdir):
    if workdir is None:
        return None, None, None
    for fname in ('.command.log', '.command.out'):
        p = workdir / fname
        if not p.is_file():
            continue
        try:
            with open(p, errors='replace') as fh:
                for line in fh:
                    m = GENOME_STATS_RE.search(line)
                    if m:
                        scaffolds = int(m.group(1).replace(',', ''))
                        bp = int(m.group(2).replace(',', ''))
                        pct_repeat = float(m.group(3))
                        return bp, scaffolds, pct_repeat
        except Exception:
            continue
    return None, None, None


def find_asm_stats_file(genome_stats_root, species):
    if not species or '_' not in species:
        return None
    genus = species.split('_', 1)[0]
    p = Path(genome_stats_root) / genus / f'{species}.asm_stats.stats.txt'
    return p if p.is_file() else None


def find_gene_info_file(genome_stats_root, species):
    if not species or '_' not in species:
        return None
    genus = species.split('_', 1)[0]
    p = Path(genome_stats_root) / genus / f'{species}.gene_stats.gene_info.csv.gz'
    return p if p.is_file() else None


def count_genes(path):
    import gzip
    try:
        with gzip.open(path, 'rt') as fh:
            n = sum(1 for _ in fh) - 1  # header row
        return max(n, 0)
    except Exception:
        return None


def parse_asm_stats(path):
    try:
        text = path.read_text(errors='replace')
    except Exception:
        return {}
    out = {}
    for field, pat in ASM_STATS_FIELD_RE.items():
        m = pat.search(text)
        if not m:
            continue
        val = m.group(1).replace(',', '')
        out[field] = float(val) if field in ('gc_pct', 'repeat_masked_pct') else int(val)
    return out


def build_asm_stats_lookup(species_list, genome_stats_root):
    """One asm_stats + gene-count lookup per distinct species (not per row)."""
    lookup = {}
    n_found = 0
    n_genes_found = 0
    for sp in species_list:
        f = find_asm_stats_file(genome_stats_root, sp)
        stats = {} if f is None else parse_asm_stats(f)
        if f is not None:
            stats['_source_file'] = str(f)
            n_found += 1

        gf = find_gene_info_file(genome_stats_root, sp)
        if gf is not None:
            gene_count = count_genes(gf)
            if gene_count is not None:
                stats['gene_count'] = gene_count
                n_genes_found += 1

        lookup[sp] = stats or None
    print(f"asm_stats DB: matched {n_found}/{len(species_list)} distinct species "
          f"under {genome_stats_root}", file=sys.stderr)
    print(f"gene_info DB: matched {n_genes_found}/{len(species_list)} distinct species "
          f"(gene counts)", file=sys.stderr)
    return lookup


def add_genome_stats(rows, genome_stats_root):
    species_list = sorted({r['species'] for r in rows if r['species']})
    asm_lookup = build_asm_stats_lookup(species_list, genome_stats_root)

    n_asm_db = 0
    n_gene_count = 0
    n_command_log = 0
    n_command_log_attempted = 0
    for r in rows:
        stats = asm_lookup.get(r['species']) or {}
        r['genome_size_bp'] = stats.get('genome_size_bp')
        r['num_scaffolds'] = stats.get('num_scaffolds')
        r['repeat_masked_pct'] = stats.get('repeat_masked_pct')
        r['gc_pct'] = stats.get('gc_pct')
        r['n50'] = stats.get('n50')
        r['gene_count'] = stats.get('gene_count')
        r['genome_stats_source'] = 'asm_stats_db' if r['genome_size_bp'] is not None else None
        r['workdir'] = None
        if r['genome_size_bp'] is not None:
            n_asm_db += 1
        if r['gene_count'] is not None:
            n_gene_count += 1

        if r['genome_size_bp'] is not None or r['process'] not in COMMAND_LOG_FALLBACK_PROCESSES:
            continue
        n_command_log_attempted += 1
        wd = find_workdir(r['root'], r['hash'])
        bp, scaffolds, pct_repeat = extract_genome_stats(wd)
        r['workdir'] = str(wd) if wd else None
        if bp is not None:
            r['genome_size_bp'] = bp
            r['num_scaffolds'] = scaffolds
            r['repeat_masked_pct'] = pct_repeat
            r['genome_stats_source'] = 'command_log'
            n_command_log += 1

    print(f"Genome stats joined: {n_asm_db} rows from asm_stats DB, "
          f"{n_command_log}/{n_command_log_attempted} additional rows from "
          f"command.log fallback (PREDICT/TRAIN only, work dir still present); "
          f"{n_gene_count} rows with gene_count", file=sys.stderr)
    return rows


def merge_with_existing(new_df, parquet_path):
    if parquet_path.is_file():
        old_df = pd.read_parquet(parquet_path)
        combined = pd.concat([old_df, new_df], ignore_index=True)
        before = len(combined)
        combined = combined.drop_duplicates(subset=['root', 'hash'], keep='last')
        print(f"Merged with existing table: {len(old_df)} old + {len(new_df)} new -> "
              f"{len(combined)} rows ({before - len(combined)} duplicate hashes resolved)",
              file=sys.stderr)
        return combined
    return new_df


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--roots', nargs='+', required=True,
                     help='Run root dirs to scan for logs/nextflow/funannotate_trace*.txt')
    ap.add_argument('--outdir',
                     default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis',
                     help='Output directory for the parquet table')
    ap.add_argument('--genome-stats-root',
                     default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/genome_stats_by_name',
                     help='genome_stats_by_name dir (species -> asm_stats.stats.txt), '
                          'shared across all project roots')
    ap.add_argument('--name', default='funannotate_runtime_metadata',
                     help='Base filename (no extension) for the output table')
    ap.add_argument('--also-csv', action='store_true',
                     help='Also write a sibling .csv of the full table')
    ap.add_argument('--no-merge', action='store_true',
                     help='Overwrite the existing parquet instead of merging into it')
    args = ap.parse_args()

    trace_files = collect_trace_files(args.roots)
    if not trace_files:
        print("No trace files found under roots.", file=sys.stderr)
        return 1
    print(f"Found {len(trace_files)} trace files across {len(args.roots)} roots", file=sys.stderr)

    rows = load_trace_rows(trace_files)
    if not rows:
        print("No task attempts parsed.", file=sys.stderr)
        return 1
    rows = add_genome_stats(rows, args.genome_stats_root)

    df = pd.DataFrame(rows)
    df['likely_skip_exit'] = (
        df['process'].isin(SKIP_EXIT_PROCESSES)
        & (df['status'] == 'COMPLETED')
        & (df['realtime_s'] < SKIP_EXIT_THRESHOLD_S)
    )
    col_order = ['root', 'species', 'process', 'hash', 'status', 'exit_code',
                 'realtime_s', 'likely_skip_exit', 'cpu_pct', 'rss_mb', 'genome_size_bp',
                 'num_scaffolds', 'repeat_masked_pct', 'gc_pct', 'n50', 'gene_count',
                 'genome_stats_source', 'first_seen_trace', 'first_seen_date', 'workdir']
    df = df[col_order]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    parquet_path = outdir / f'{args.name}.parquet'

    if not args.no_merge:
        df = merge_with_existing(df, parquet_path)

    df.to_parquet(parquet_path, index=False, engine='pyarrow', compression='snappy')
    print(f"Wrote {len(df)} rows -> {parquet_path}", file=sys.stderr)

    if args.also_csv:
        csv_path = outdir / f'{args.name}.csv'
        df.to_csv(csv_path, index=False)
        print(f"Wrote {csv_path}", file=sys.stderr)

    print("\nObservations per process:", file=sys.stderr)
    print(df.groupby('process').size().sort_values(ascending=False).to_string(), file=sys.stderr)

    print("\nGenome-stats coverage per process (rows with genome_size_bp populated):",
          file=sys.stderr)
    coverage = (df.assign(has_stats=df['genome_size_bp'].notna())
                  .groupby('process')['has_stats']
                  .agg(['sum', 'count']))
    coverage['pct'] = (100 * coverage['sum'] / coverage['count']).round(1)
    coverage = coverage.sort_values('count', ascending=False)
    print(coverage.rename(columns={'sum': 'with_stats', 'count': 'n'}).to_string(),
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
