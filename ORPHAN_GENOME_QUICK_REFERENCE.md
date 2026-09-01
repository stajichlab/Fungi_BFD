# Orphaned Genome Quick Reference

## Problem
When ASMIDs are removed from `samples.csv`, they may still exist on disk in various locations, creating "orphaned" files that:
- Waste disk space
- May cause confusion about what's actually being used
- Can slow down file system operations
- May conflict with newer versions

## Three-Layer Solution

### Layer 1: Tracking (After each regeneration)
**File: `data/curation/removed_asmids.csv`**

Documents every ASMID removal with:
- ASMID (what was removed)
- REASON (version_update, superseded, excluded, etc.)
- DATE_REMOVED (when it was removed)
- SPECIES & STRAIN (for identification)
- REPLACEMENT_ASMID (if updated to newer version)
- NOTES (any additional context)

**Example:**
```csv
GCA_036244875.1_ASM3624487v1,version_update,2026-08-03,Rhodotorula paludigena,,GCA_036244875.1_Rhopal_B15_v1.0,New assembly name format
```

### Layer 2: Detection (Find orphans)
**Script: `scripts/check_orphaned_genomes.py`**

Automatically finds genomes on disk not in samples.csv:

```bash
# Run after regenerating samples.csv
python3 scripts/check_orphaned_genomes.py

# Output shows:
# ✓ EXPECTED: GCA_036244875.1_ASM3624487v1 (in removed_asmids.csv)
# ✗ UNEXPECTED: GCA_900000000.1_Unknown (NOT in removed_asmids.csv - investigate!)
```

### Layer 3: Recovery (Git history)
**Using git to track changes:**

```bash
# See what was in samples.csv at previous commit
git log --oneline -- samples.csv | head -10

# Show samples.csv before recent changes
git show HEAD~1:samples.csv | grep GCA_036244875.1_ASM3624487v1

# Restore samples.csv to previous version
git checkout HEAD~1 -- samples.csv
```

## Quick Start Workflow

### When regenerating samples.csv:

```bash
# 1. Regenerate
python3 scripts/create_samples_file.py --outfile samples.csv

# 2. Find what changed (git shows removed vs added)
git diff samples.csv | grep "^<" | cut -d',' -f1 > /tmp/removed.txt
git diff samples.csv | grep "^>" | cut -d',' -f1 > /tmp/added.txt

# 3. Manually update removed_asmids.csv
# Edit: data/curation/removed_asmids.csv
# Add: new rows for each removed ASMID with reason

# 4. Verify no unexpected orphans
python3 scripts/check_orphaned_genomes.py

# 5. Commit changes
git add samples.csv data/curation/removed_asmids.csv
git commit -m "Regenerate samples.csv: +366 new, -20 old (8 updated, 12 removed)

See data/curation/removed_asmids.csv for removal details"
```

## Example Scenarios

### Scenario 1: Version Update (most common)
```
OLD: GCA_036986815.1_ASM3698681v1 (Saccharomyces cerevisiae)
NEW: GCA_036986815.1_AH465_1.0 (same species, improved assembly name)

REASON: version_update
REPLACEMENT_ASMID: GCA_036986815.1_AH465_1.0
NOTES: Improved assembly name
```
**Action:** Safe to remove old version, new one is in samples.csv

### Scenario 2: Superseded (old version still exists)
```
OLD: GCA_000827195.1_Laccaria_amethystina_LaAM-08-1_v1.0 (Laccaria amethystina)
NEW: GCA_000827195.2_ASM82719v2 (v2 available)

REASON: superseded
REPLACEMENT_ASMID: GCA_000827195.2_ASM82719v2
NOTES: v2 available
```
**Action:** v1 no longer in use, can archive/delete

### Scenario 3: Unexpected Orphan (investigate!)
```
ORPHANED: GCA_900000000.1_Unknown
STATUS: ✗ UNEXPECTED (not in removed_asmids.csv)

ACTION:
1. Is this a custom genome you added?
2. Is it a pipeline output that should be deleted?
3. Should it be added back to samples.csv?
4. Document decision in removed_asmids.csv or restore it
```

## Common Cleanup Tasks

### Find all orphaned genomes in data/
```bash
python3 scripts/check_orphaned_genomes.py

# Or manually:
tail -n +2 samples.csv | cut -d',' -f1 | sort > /tmp/current.txt
find data/genomes -type d -name 'GC*' | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | sort -u | \
  while read asmid; do
    grep -q "^$asmid" /tmp/current.txt || echo "ORPHANED: $asmid"
  done
```

### Archive removed genomes
```bash
# Save to external storage
tar czf archive/removed_genomes_2026-08-03.tar.gz \
  data/genomes/GCA_036244875.1_ASM3624487v1 \
  data/genomes/GCA_036245355.1_ASM3624535v1

# Then delete
rm -rf data/genomes/GCA_036244875.1_ASM3624487v1 \
       data/genomes/GCA_036245355.1_ASM3624535v1

# Commit
git add archive/
git commit -m "Archive removed genomes from 2026-08-03 regeneration"
```

### Recover a removed genome
```bash
# Find in git history
git log --all --full-history -- "data/genomes/GCA_036244875.1_*"

# Restore from specific commit
git checkout <commit-hash> -- data/genomes/GCA_036244875.1_ASM3624487v1/

# Or restore entire samples.csv from previous version
git show <commit-hash>:samples.csv > samples.csv.old
# Then manually edit samples.csv to re-add the ASMID
```

## Key Files

| File | Purpose | Update When |
|------|---------|-------------|
| `data/curation/removed_asmids.csv` | Master record of removals | After samples.csv regeneration |
| `scripts/check_orphaned_genomes.py` | Automated orphan detection | Monthly or after regeneration |
| `data/curation/orphan_audit.log` | Append-only log of checks | Auto-generated by script |
| `DATA_MANIFEST.md` | High-level summary of changes | After samples.csv regeneration |

## Integration with samples.csv Workflow

**BEFORE regenerating:**
1. Read this guide
2. Check codebase for any custom genome additions

**AFTER regenerating:**
1. Run `python3 scripts/check_orphaned_genomes.py`
2. Update `data/curation/removed_asmids.csv` with removed ASMIDs
3. Verify no unexpected orphans (all marked as "✓ EXPECTED")
4. Commit with detailed message about removals
5. Optional: Archive and delete old genomes to free space

**MONTHLY or QUARTERLY:**
1. Run orphan check: `python3 scripts/check_orphaned_genomes.py`
2. Review for unexpected orphans
3. Archive/delete if confirmed
4. Document in `DATA_MANIFEST.md`

## Troubleshooting

### Q: I found unexpected orphans. What do I do?
A: Check if they should be:
1. In `samples.csv` → add them back via `data/curation/preferred_asmids.txt`
2. Deleted → add to `data/curation/removed_asmids.csv` with reason
3. In a different location → move them and document
4. Archived → run archive script and remove

### Q: I accidentally removed an ASMID. How do I recover it?
A: Restore from git:
```bash
git show HEAD~1:samples.csv | grep GCA_xxxxx
# If it was there, restore it
git checkout HEAD~1 -- samples.csv
# Then re-apply other changes you wanted to keep
```

### Q: How do I find out why an ASMID was removed?
A: Check these in order:
1. `data/curation/removed_asmids.csv` - explicit reason
2. `git log -p -- samples.csv` - see commit that removed it
3. `data/curation/exclude_asmids.txt` - was it explicitly excluded?
4. `data/curation/keep_dupes.csv` - did it fail dedup logic?

### Q: Orphan script says "UNEXPECTED". Should I be worried?
A: Check:
1. Is it a pipeline output? (analysis/, misc/) → probably safe to delete
2. Is it data you added? → add to samples.csv or remove_asmids.csv
3. Is it old test data? → delete and update removed_asmids.csv
4. Unsure? → ask before deleting

## References

- Full strategy guide: `ORPHAN_GENOME_STRATEGIES.md`
- Samples regeneration report: `SAMPLES_CSV_REGENERATION_REPORT.txt`
- Git history: `git log -- samples.csv`
