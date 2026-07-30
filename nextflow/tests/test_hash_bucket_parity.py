#!/usr/bin/env python3
"""test_hash_bucket_parity.py — cross-check Groovy hashBucket() against Python hash_bucket().

Generates a fixed set of realistic ASMID/LOCUSTAG-like keys, runs them through
both the Groovy implementation (nextflow/modules/common/utils.nf::hashBucket(),
via the throwaway nextflow/tests/hash_bucket_probe.nf probe script) and the
Python implementation (nextflow/bin/genome_stats_paths.py::hash_bucket()), and
asserts the two agree bit-for-bit for every key and both bucket widths (2 and
3 hex chars).

A mismatch here is a silent path miss in production (a file written by
Nextflow at one bucket path, looked up by a Python consumer script at a
different one) -- not a crash -- so this is a required check before either
implementation ships, not an optional nicety.

Usage:
    python3 nextflow/tests/test_hash_bucket_parity.py
"""

import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent.parent
sys.path.insert(0, str(TESTS_DIR.parent / "bin"))

from genome_stats_paths import hash_bucket  # noqa: E402


def build_test_keys():
    """Realistic ASMID/LOCUSTAG-shaped keys, plus edge cases, deterministic (seeded)."""
    rng = random.Random(20260730)
    keys = [
        "GCA_000149445.2_ASM14944v2",
        "GCF_010015735.1_Aaoar1",
        "FF5840CF",
        "F2EE6837",
        "",  # edge case: empty key must still hash deterministically, not crash
        "a",  # edge case: single char
    ]
    for _ in range(500):
        acc = rng.choice(["GCA", "GCF"])
        digits = "".join(rng.choice("0123456789") for _ in range(9))
        ver = rng.randint(1, 3)
        suffix = "".join(rng.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789") for _ in range(8))
        keys.append(f"{acc}_{digits}.{ver}_{suffix}")
    return keys


def run_groovy_probe(keys, widths):
    if shutil.which("nextflow") is None:
        print("SKIP: nextflow not found on PATH — cannot run the Groovy side of the parity test", file=sys.stderr)
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        keys_file = Path(tmpdir) / "keys.tsv"
        with keys_file.open("w") as fh:
            for key in keys:
                for width in widths:
                    fh.write(f"{key}\t{width}\n")

        probe = TESTS_DIR / "hash_bucket_probe.nf"
        result = subprocess.run(
            ["nextflow", "run", str(probe), "--keys_file", str(keys_file)],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            print(result.stdout, file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"nextflow probe run failed (exit {result.returncode})")

        groovy_results = {}
        for line in result.stdout.splitlines():
            if not line.startswith("BUCKET_RESULT:"):
                continue
            payload = line[len("BUCKET_RESULT:"):]
            key, width, bucket = payload.split("\t")
            groovy_results[(key, int(width))] = bucket
        return groovy_results


def main():
    keys = build_test_keys()
    widths = [2, 3]

    groovy_results = run_groovy_probe(keys, widths)

    mismatches = []
    checked = 0
    for key in keys:
        for width in widths:
            py_bucket = hash_bucket(key, width)
            assert len(py_bucket) == width, f"Python hash_bucket returned wrong width for {key!r}"
            checked += 1
            if groovy_results is not None:
                groovy_bucket = groovy_results.get((key, width))
                if groovy_bucket != py_bucket:
                    mismatches.append((key, width, py_bucket, groovy_bucket))

    if mismatches:
        print(f"FAIL: {len(mismatches)} Groovy/Python hash_bucket mismatches:", file=sys.stderr)
        for key, width, py_bucket, groovy_bucket in mismatches[:20]:
            print(f"  key={key!r} width={width} python={py_bucket} groovy={groovy_bucket}", file=sys.stderr)
        sys.exit(1)

    mode = "Groovy+Python" if groovy_results is not None else "Python-only (nextflow unavailable)"
    print(f"OK: {checked} key/width combinations checked ({mode}), no mismatches")


if __name__ == "__main__":
    main()
