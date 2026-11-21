#!/usr/bin/env python3
"""Build samples.csv from the NCBI accession + taxonomy tables.

Reproducible, dependency-light replacement for the original ad-hoc script.

Pipeline (all upstream files are produced by ../../1KFG/2026/NCBI_fungi/Makefile):
    ncbi_accessions.csv          (ACCESSION,SPECIES,STRAIN,NCBI_TAXID,BIOPROJECT,...,ASM_NAME,ASM_FOLDER)
    ncbi_accessions_taxonomy.csv (ASM_ACCESSION,NCBI_TAXID,SPECIES_IN,STRAIN,<lineage>,SPECIES)
        |
        v   sanitize (scripts/sample_sanitize.py)  +  curation (data/curation/)
        v
    samples.csv

ASMID is the upstream ASM_FOLDER value: the filesystem-safe (no '#', spaces,
commas, parens) '{ACCESSION}_{ASM_NAME}' that names the on-disk genome folder
and its '<ASM_FOLDER>_genomic.fna.gz' file, and is also the taxonomy join key
(ncbi_accessions_taxonomy.csv ASM_ACCESSION == ASM_FOLDER).  This keeps ASMID
in lock-step with the path nextflow/funannotate.nf builds for input genomes.
Older accession CSVs without an ASM_FOLDER column fall back to the previous
behaviour: ASMID = clean_asmid(raw) joined to taxonomy on the raw base.

Curation is data, not manual post-edits:
    data/curation/exclude_asmids.txt  hard removals (one raw ASMID per line; '\\t reason'; '#' comments)
    data/curation/keep_dupes.csv      species+strain isolates intentionally kept as >1 assembly
    data/curation/overrides.csv       optional per-ASMID field corrections (ASMID,FIELD,VALUE)

Dedup: known duplicate assemblies are handled by the exclude list.  Any *new*
species+strain collision not covered by curation is resolved by a default
tie-breaker (prefer RefSeq GCF, else newest GCA accession) and logged, unless
the isolate is in keep_dupes.csv.

LOCUSTAG is the md5-derived tag of the *raw* '{ACCESSION}_{ASM_NAME}' (stable
across display cleaning); collisions get deterministic A/B/C suffixes after a
stable sort.
"""

import argparse
import csv
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sample_sanitize import clean_species, clean_strain, clean_asmid, backfill_strain

DEF_NCBI = "../../1KFG/2026/NCBI_fungi"

OUT_FIELDS = ['ASMID', 'SPECIES_IN', 'STRAIN', 'BIOPROJECT', 'NCBI_TAXONID',
              'BUSCO_LINEAGE', 'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS',
              'ORDER', 'FAMILY', 'GENUS', 'SPECIES', 'TRANSL_TABLE', 'LOCUSTAG']
LINEAGE_COLS = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS']

# --- translation tables ----------------------------------------------------
TRANSL_TABLE_12_FAMILIES = {'Debaryomycetaceae'}
TRANSL_TABLE_12_GENERA = {'Clavispora'}
TRANSL_TABLE_26_GENERA = {'Pachysolen'}


def get_transl_table(family, genus):
    if family in TRANSL_TABLE_12_FAMILIES or genus in TRANSL_TABLE_12_GENERA:
        return 12
    if genus in TRANSL_TABLE_26_GENERA:
        return 26
    return 1


def busco_lineage(phylum):
    if phylum in ('Ascomycota', 'Basidiomycota'):
        return 'dikarya'
    if phylum == 'Mucoromycota' or phylum == "Glomeromyocta" or phylum == "Calcarisporiellomycota":
        return 'fungi'
    if phylum == 'Microsporidia':
        return 'microsporidia'
    return 'fungi'


def get_locus(raw_asm_base):
    return 'F' + hashlib.md5(raw_asm_base.encode('utf-8')).hexdigest().upper()[-7:]


def accnum(asmid):
    """(numeric, version) of a GCA_/GCF_ accession for 'newest' comparison."""
    parts = asmid.split('_')
    try:
        base, ver = parts[1].split('.')
        return (int(base), int(ver))
    except Exception:
        return (0, 0)


def prefix(asmid):
    return asmid.split('_', 1)[0]


def log(msg):
    print(msg, file=sys.stderr)


def clean_strain_tag(raw):
    """Port of cleanStrain() in nextflow/modules/common/utils.nf, step for step.

    Keep the steps and their ORDER identical to the Groovy. The pipeline names
    each results folder with the Groovy version; this copy only predicts that
    name for the clash check. If the two disagree, the check passes and two
    genomes still land in one folder -- the failure it exists to prevent."""
    s = (raw or '').strip()
    s = re.sub(r'[\'"]', '', s)             # .replaceAll(/['"]/, '')
    s = s.split(';')[0].strip()             # .split(';')[0].trim()
    s = s.replace(':', ' ')                 # .replace(':', ' ')
    s = re.sub(r'^\s*\*+', '', s)           # leading '*' -> removed
    s = re.sub(r'\*+\s*$', '', s)           # trailing '*' -> removed
    s = re.sub(r'\s*\*+\s*', '-', s)        # '*' between words -> '-'
    s = re.sub(r'[&`$|<>()]+', ' ', s)      # shell metacharacters -> space
    return s.strip()


def make_sample_tag(species, strain):
    """Port of makeSampleTag() in nextflow/modules/common/utils.nf.

    This is the per-genome folder key `out` under Fungi_BFD_runs/genome_annotation.
    Two samples.csv rows with the same tag would share one results folder."""
    sp = re.sub(r'[\'"]', '', (species or '').strip())
    st = clean_strain_tag(strain)
    return re.sub(r'[\s/#\[\]?{}]+', '_', '_'.join(x for x in (sp, st) if x))


def load_out_suffix(path):
    """ASMID -> (SUFFIX, REASON). SUFFIX is appended to STRAIN as '<strain>_<suffix>'
    (or is the whole strain when it is empty). SUFFIX 'HOLD' keeps the row out of
    samples.csv until someone picks a real suffix."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline='') as f:
        lines = (l for l in f if l.strip() and not l.lstrip().startswith('#'))
        for row in csv.DictReader(lines):
            out[row['ASMID'].strip()] = ((row.get('SUFFIX') or '').strip(),
                                         (row.get('REASON') or '').strip())
    return out


# --- loaders ---------------------------------------------------------------
def load_exclude(path):
    excl = {}
    if not path or not os.path.exists(path):
        return excl
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split('\t')
            asmid = parts[0].strip()
            reason = parts[1].strip() if len(parts) > 1 else ''
            if asmid:
                excl[asmid] = reason
    return excl


def load_preferred(path):
    """Set of ASMIDs to force as the winner of their species+strain dedup group.

    Entries must be the ASMID exactly as it appears in samples.csv (the upstream
    ASM_FOLDER value); one per line, '#' comments, optional '\\t reason'.  Stored
    verbatim so it matches rec['ASMID']."""
    pref = set()
    if not path or not os.path.exists(path):
        return pref
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            asmid = line.split('\t')[0].strip()
            if asmid:
                pref.add(asmid)
    return pref


def load_keep_dupes(path):
    keep = set()
    if not path or not os.path.exists(path):
        return keep
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        for row in r:
            keep.add((clean_species(row.get('SPECIES', '')).lower(),
                      clean_strain(row.get('STRAIN', '')).lower()))
    return keep


def load_overrides(path):
    ov = {}
    if not path or not os.path.exists(path):
        return ov
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        for row in r:
            ov.setdefault(row['ASMID'], {})[row['FIELD']] = row['VALUE']
    return ov


def load_taxonomy(path):
    """asm_base -> dict with lineage cols + binomial SPECIES (deduped, first wins)."""
    tax = {}
    n_lines = n_conflict = 0
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        for row in r:
            n_lines += 1
            asm = row['ASM_ACCESSION']
            rec = {c: (row.get(c, '') or '').strip() for c in LINEAGE_COLS}
            rec['SPECIES'] = (row.get('SPECIES', '') or '').strip()
            if asm in tax:
                if tax[asm] != rec:
                    n_conflict += 1
                continue
            tax[asm] = rec
    log(f"[taxonomy] {n_lines} lines -> {len(tax)} unique assemblies "
        f"({n_lines - len(tax)} duplicate lines, {n_conflict} conflicting)")
    return tax


# --- main ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--accessions', default=f'{DEF_NCBI}/ncbi_accessions.csv')
    ap.add_argument('--taxonomy', default=f'{DEF_NCBI}/ncbi_accessions_taxonomy.csv')
    ap.add_argument('--exclude', default='data/curation/exclude_asmids.txt')
    ap.add_argument('--keep-dupes', default='data/curation/keep_dupes.csv')
    ap.add_argument('--prefer', default='data/curation/preferred_asmids.txt')
    ap.add_argument('--overrides', default='data/curation/overrides.csv')
    ap.add_argument('--out-suffix', default='data/curation/out_suffix.csv',
                    help='ASMID,SUFFIX,REASON: disambiguate haplotype/clashing genome folders')
    ap.add_argument('--clash-report', default=None,
                    help='write rows whose results-folder key (`out`) is shared (held back) to this CSV')
    ap.add_argument('--allow-out-clash', action='store_true',
                    help='keep rows that share an `out` key (default: hold them all back)')
    ap.add_argument('--allow-missing-haplotype', action='store_true',
                    help='proceed when the accessions CSV has no HAPLOTYPE column (old '
                         'snapshots). hap1/hap2 pairs of one isolate may then be collapsed '
                         'by the dedup rule.')
    ap.add_argument('--outfile', default='samples.csv')
    ap.add_argument('--no-dedup-rule', action='store_true',
                    help='do not auto-resolve uncurated species+strain collisions')
    args = ap.parse_args()

    for p in (args.accessions, args.taxonomy):
        if not os.path.exists(p):
            ap.error(f'input not found: {p}')

    exclude = load_exclude(args.exclude)
    keep_dupes = load_keep_dupes(args.keep_dupes)
    preferred = load_preferred(args.prefer)
    overrides = load_overrides(args.overrides)
    out_suffix = load_out_suffix(args.out_suffix)
    held_suffix = []
    tax = load_taxonomy(args.taxonomy)

    # ----- read accessions (master), sanitize, apply exclude ---------------
    records = []
    n_in = n_excluded = n_no_tax = 0
    seen_asmbase = set()
    with open(args.accessions, newline='') as f:
        r = csv.DictReader(f)
        need = {'ACCESSION', 'SPECIES', 'STRAIN', 'NCBI_TAXID', 'BIOPROJECT', 'ASM_NAME'}
        missing = need - set(r.fieldnames or [])
        if missing:
            ap.error(f'accessions file missing columns: {sorted(missing)}')
        # HAPLOTYPE carries the hap1/hap2 role and is part of the dedup key, so
        # the two haplotypes of one isolate stay two genomes. Without the column
        # every row looks like haplotype '', the pair collapses to a single
        # species+strain group and the tie-breaker drops one half -- silently,
        # because nothing else in the run distinguishes them. There is no safe
        # fallback (unlike ASM_FOLDER below, which has one), so refuse instead.
        # Old snapshots predate the column; extract_hand_edits.py regenerates
        # against those on purpose and passes --allow-missing-haplotype.
        if 'HAPLOTYPE' not in (r.fieldnames or []) and not args.allow_missing_haplotype:
            ap.error('accessions file has no HAPLOTYPE column: hap1/hap2 assemblies of one '
                     'isolate would be silently collapsed by the dedup rule. Re-export the '
                     'accessions table, or pass --allow-missing-haplotype to accept that.')
        # ASM_FOLDER is the upstream-sanitized (filesystem-safe; no '#', spaces,
        # commas, parens) '{ACCESSION}_{ASM_NAME}'.  It is the canonical on-disk
        # folder + genome-file prefix and is the join key in the taxonomy table
        # (ASM_ACCESSION == ASM_FOLDER), so we adopt it verbatim as ASMID.  Older
        # accession CSVs predate the column: fall back to clean_asmid(raw), which
        # reproduces the previous ASMID/join behaviour.
        has_asm_folder = 'ASM_FOLDER' in (r.fieldnames or [])
        for row in r:
            n_in += 1
            raw_asm_base = f"{row['ACCESSION']}_{row['ASM_NAME']}"
            asm_folder = (row.get('ASM_FOLDER') or '').strip() if has_asm_folder else ''
            asmid = asm_folder or clean_asmid(raw_asm_base)
            if raw_asm_base in seen_asmbase:
                continue  # defensive: duplicate accession line
            seen_asmbase.add(raw_asm_base)
            if raw_asm_base in exclude:
                n_excluded += 1
                continue

            species_in_raw = row['SPECIES']
            # New CSVs: taxonomy ASM_ACCESSION == ASM_FOLDER == asmid.
            # Old CSVs: taxonomy was keyed on the raw '{ACCESSION}_{ASM_NAME}'.
            t = tax.get(asmid if has_asm_folder else raw_asm_base)
            if t is None:
                n_no_tax += 1
                t = {c: '' for c in LINEAGE_COLS}
                t['SPECIES'] = ''
            binomial = clean_species(t['SPECIES'] or species_in_raw)
            species_in = clean_species(species_in_raw)
            strain = backfill_strain(clean_strain(row['STRAIN']), species_in, binomial)

            rec = {
                'ASMID': asmid,
                'RAW_ASM_BASE': raw_asm_base,
                'SPECIES_IN': species_in,
                'STRAIN': strain,
                'BIOPROJECT': row['BIOPROJECT'],
                'NCBI_TAXONID': row['NCBI_TAXID'],
                'BUSCO_LINEAGE': busco_lineage(t['PHYLUM']),
                'SPECIES': binomial,
                'TRANSL_TABLE': get_transl_table(t['FAMILY'], t['GENUS']),
            }
            for c in LINEAGE_COLS:
                rec[c] = t[c]
            # haplotype role from the accessions table (hap1/hap2/''); dedup key only
            rec['HAPLOTYPE'] = (row.get('HAPLOTYPE') or '').strip()
            # per-ASMID overrides
            for fld, val in overrides.get(rec['ASMID'], {}).items():
                rec[fld] = val
            # results-folder disambiguation (after overrides, before dedup)
            suf = out_suffix.get(rec['ASMID'])
            if suf:
                if suf[0].upper() == 'HOLD':
                    held_suffix.append((rec['ASMID'], suf[1]))
                    continue
                rec['STRAIN'] = f"{rec['STRAIN']}_{suf[0]}" if rec['STRAIN'] else suf[0]
            records.append(rec)

    log(f"[accessions] {n_in} rows -> {len(records)} kept "
        f"({n_excluded} excluded by curation, {n_no_tax} without taxonomy)")

    # ----- dedup tie-breaker for uncurated species+strain collisions -------
    if not args.no_dedup_rule:
        groups = {}
        for rec in records:
            # haplotype role is part of the key: hap1 and hap2 are separate genomes
            k = (rec['SPECIES'].lower(), rec['STRAIN'].lower(), rec['HAPLOTYPE'])
            groups.setdefault(k, []).append(rec)
        drop = set()
        n_collapsed = 0
        for k, recs in groups.items():
            if not k[1] or len(recs) == 1:
                continue            # no strain, or singleton -> nothing to resolve
            if k[:2] in keep_dupes:
                continue            # intentional multi-assembly isolate
            # winner selection:
            #   1. an explicitly preferred ASMID, else
            #   2. RefSeq GCF, else
            #   3. newest accession
            pref = [x for x in recs if x['ASMID'] in preferred]
            if pref:
                pool = pref
            else:
                gcf = [x for x in recs if prefix(x['ASMID']) == 'GCF']
                pool = gcf if gcf else recs
            winner = max(pool, key=lambda x: accnum(x['ASMID']))
            for x in recs:
                if x is not winner:
                    drop.add(id(x))
                    n_collapsed += 1
                    log(f"[dedup] {k[0]} | {k[1]}: drop {x['ASMID']} "
                        f"(kept {winner['ASMID']})")
        if drop:
            records = [x for x in records if id(x) not in drop]
        log(f"[dedup] auto-collapsed {n_collapsed} uncurated duplicate assemblies")

    # ----- results-folder (`out`) clash check -------------------------------
    # Two rows with the same `out` would share one genome_annotation folder, and
    # -resume would reuse the first genome's predictions for the second. Hold every
    # clashing row back (loudly) unless out_suffix.csv has separated them.
    by_out = {}
    for rec in records:
        by_out.setdefault(make_sample_tag(rec['SPECIES'], rec['STRAIN']), []).append(rec)
    clashes = {o: rs for o, rs in by_out.items() if len(rs) > 1}
    # keep_dupes isolates are intentionally >1 assembly. They still share a folder,
    # so warn, but do not hold them back (that would drop rows already in use).
    for o, rs in list(clashes.items()):
        if all((x['SPECIES'].lower(), x['STRAIN'].lower()) in keep_dupes for x in rs):
            log(f"[out-clash] WARNING {o}: keep_dupes isolate, {len(rs)} assemblies share one "
                f"folder key ({', '.join(x['ASMID'] for x in rs)}); give one an out_suffix.csv row")
            del clashes[o]
    clash_rows = []
    for o, rs in sorted(clashes.items()):
        for x in rs:
            clash_rows.append([x['ASMID'], o, x['HAPLOTYPE'], len(rs), x['SPECIES'], x['STRAIN'],
                               'kept' if args.allow_out_clash else 'held'])
            log(f"[out-clash] {o}: {x['ASMID']} ({len(rs)} rows share this folder key)")
    for asmid, reason in held_suffix:
        clash_rows.append([asmid, '', '', 0, '', '', 'held (out_suffix HOLD): ' + reason])
        log(f"[out-hold] {asmid}: held by out_suffix.csv ({reason})")
    if clashes and not args.allow_out_clash:
        held_ids = {id(r) for rs in clashes.values() for r in rs}
        records = [r for r in records if id(r) not in held_ids]
    log(f"[out-clash] {len(clashes)} shared folder keys; {len(held_suffix)} rows on HOLD")
    if args.clash_report:
        with open(args.clash_report, 'w', newline='') as f:
            w = csv.writer(f, lineterminator='\n')
            w.writerow(['ASMID', 'OUT', 'HAPLOTYPE', 'N_IN_GROUP', 'SPECIES', 'STRAIN', 'STATUS'])
            w.writerows(clash_rows)

    # ----- deterministic LOCUSTAG ------------------------------------------
    records.sort(key=lambda x: x['RAW_ASM_BASE'])
    seen = {}
    for rec in records:
        locus = get_locus(rec['RAW_ASM_BASE'])
        if locus in seen:
            suffixed = locus + chr(ord('A') + seen[locus] - 1)
            log(f"[locus] collision {locus} -> {suffixed} for {rec['ASMID']}")
            seen[locus] += 1
            locus = suffixed
        else:
            seen[locus] = 1
        rec['LOCUSTAG'] = locus

    # ----- sort by species for output ----------------------------------------
    records.sort(key=lambda x: (x['SPECIES'].lower(), x['STRAIN'].lower(), x['ASMID']))

    # ----- write -----------------------------------------------------------
    with open(args.outfile, 'wt', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(OUT_FIELDS)
        for rec in records:
            w.writerow([rec[c] for c in OUT_FIELDS])
    log(f"[write] {len(records)} assemblies -> {args.outfile}")


if __name__ == '__main__':
    main()
