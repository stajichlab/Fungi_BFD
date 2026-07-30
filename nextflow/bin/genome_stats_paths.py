#!/usr/bin/env python3
"""genome_stats_paths.py — hash-bucket path helpers for genome_stats/function outputs.

genome_stats/function outputs are stored one file per genome (or one file per
genome per data column), which reaches tens of thousands of files in a single
flat directory (see todo/genome_stats_storage_reorg.md, T-014) — slow to
list/back up on this cluster's NFS-backed storage. hash_bucket() fans a stable
key (ASMID or LOCUSTAG, never the derived Genus_species_strain tag — see the
plan for why) out into a hex subdirectory so no single directory holds more
than a few hundred files.

The Groovy equivalent is nextflow/modules/common/utils.nf::hashBucket()/
bucketWidthFor()/hashBucketForType(). Both must produce identical output for
identical (key, width) — a mismatch here is a silent path miss, not a crash.
Verified in sync via nextflow/tests/test_hash_bucket_parity.py.

Usage as a library:
    from genome_stats_paths import hash_bucket, bucket_width_for, hash_bucket_for_type, genome_stats_path

    hash_bucket_for_type("asm_stats", "GCA_000149445.2_ASM14944v2")
    genome_stats_path("results/genome_stats", "asm_stats", "GCA_000149445.2_ASM14944v2", ext="stats.txt")

Usage as a CLI (for the Groovy/Python parity test):
    python3 genome_stats_paths.py --key GCA_000149445.2_ASM14944v2 --width 2
"""

import argparse
import hashlib

# Per-type bucket width (todo/genome_stats_storage_reorg.md §A): most
# genome_stats types get 256 buckets (2 hex chars); gene_stats and
# intergenic_stats get 4096 buckets (3 hex chars) because they carry ~7x more
# files per genome than the rest. Must stay in sync with
# nextflow/modules/common/utils.nf::bucketWidthFor().
BUCKET_WIDTH = {
    "asm_stats": 2,
    "asm_reports": 2,
    "BUSCO_genome": 2,
    "BUSCO_protein": 2,
    "aa_freq": 2,
    "codon_freq": 2,
    "gene_stats": 3,
    "intergenic_stats": 3,
}

# Default width for any type not in BUCKET_WIDTH (e.g. a function/* subfolder
# not yet migrated). Must match the Groovy default in bucketWidthFor().
DEFAULT_WIDTH = 2


def hash_bucket(key: str, width: int) -> str:
    """Return the leading `width` hex chars of sha1(key)."""
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()
    return digest[:width]


def bucket_width_for(type_: str) -> int:
    """Resolve the bucket width for a genome_stats/function subfolder type."""
    return BUCKET_WIDTH.get(type_, DEFAULT_WIDTH)


def hash_bucket_for_type(type_: str, key: str) -> str:
    """Convenience wrapper: look up the right width for `type_` and hash `key`."""
    return hash_bucket(key, bucket_width_for(type_))


def genome_stats_path(root: str, type_: str, key: str, ext: str) -> str:
    """Build the full bucketed path: {root}/{type_}/{bucket}/{key}.{ext}."""
    bucket = hash_bucket_for_type(type_, key)
    return f"{root}/{type_}/{bucket}/{key}.{ext}"


def main():
    parser = argparse.ArgumentParser(
        description="Hash-bucket helper CLI, primarily for the Groovy/Python parity test."
    )
    parser.add_argument("--key", required=True, help="ASMID or LOCUSTAG key to bucket")
    parser.add_argument("--width", type=int, required=True, help="bucket width in hex chars")
    args = parser.parse_args()
    print(hash_bucket(args.key, args.width))


if __name__ == "__main__":
    main()
