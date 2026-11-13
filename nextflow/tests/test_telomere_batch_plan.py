"""Test the telomere batch plan (planFreq + batchFreq in BFD_GENOME_STATS.nf).

Runs the real Groovy functions through nextflow/tests/telomere_plan_probe.nf on
fake genomes and checks that:
  - every genome that needs computing is in exactly one batch, none in two;
  - no batch holds more than telomere_batch_size genomes;
  - genomes with a current output in the right bucket are not batched and are
    emitted on the cached channel (the one BFD_MERGE adds to the merge list);
  - an output older than its genome is deleted and the genome is batched;
  - an output in the wrong bucket does not count as cached.

Needs `nextflow` on PATH (module load nextflow); skipped otherwise.
Run: python3 -m pytest nextflow/tests/test_telomere_batch_plan.py
"""
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent / "bin"))
from genome_stats_paths import hash_bucket_for_type  # noqa: E402

SUBDIR = "telomeres_v2"
BATCH_SIZE = 5


def store_path(root, asmid):
    return root / SUBDIR / hash_bucket_for_type(SUBDIR, asmid) / f"{asmid}.telomeres.tsv.gz"


def run_probe(tmp_path, genome_dir, stats_dir):
    result = subprocess.run(
        ["nextflow", "run", "-ansi-log", "false", str(TESTS_DIR / "telomere_plan_probe.nf"),
         "--genome_dir", str(genome_dir),
         "--genome_stats_outdir", str(stats_dir),
         "--telomere_batch_size", str(BATCH_SIZE)],
        cwd=tmp_path, capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    # Search, do not anchor: Nextflow can print a marker on the same line as
    # its own log text.
    cached, batches = {}, []
    for kind, payload in re.findall(r"(CACHED|BATCH):([^\n]*)", result.stdout):
        if kind == "CACHED":
            asmid, path = payload.split("\t")
            cached[asmid] = path
        else:
            batches.append(payload.split(","))
    return cached, batches


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")
def test_plan_assigns_each_genome_once_and_passes_cached(tmp_path):
    genome_dir = tmp_path / "dna"
    stats_dir = tmp_path / "genome_stats"
    genome_dir.mkdir()
    asmids = [f"GCA_{i:09d}.1_TEST{i}" for i in range(23)]
    old = time.time() - 3600
    for a in asmids:
        g = genome_dir / f"{a}.fa"
        g.write_text(f">{a}\nACGT\n")
        os.utime(g, (old, old))

    current = asmids[0:4]      # output newer than genome: cached
    stale = asmids[4:6]        # output older than genome: deleted, recomputed
    misplaced = asmids[6:7]    # output in the wrong bucket: not cached
    for a in current:
        p = store_path(stats_dir, a)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"x")
    for a in stale:
        p = store_path(stats_dir, a)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"x")
        os.utime(p, (old - 3600, old - 3600))
    for a in misplaced:
        right = hash_bucket_for_type(SUBDIR, a)
        wrong = "00" if right != "00" else "01"
        p = stats_dir / SUBDIR / wrong / f"{a}.telomeres.tsv.gz"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"x")

    cached, batches = run_probe(tmp_path, genome_dir, stats_dir)

    todo = set(asmids) - set(current)
    batched = [a for b in batches for a in b]
    assert len(batched) == len(set(batched)), "a genome is in more than one batch"
    assert set(batched) == todo, "batched genomes differ from the genomes to compute"
    assert all(1 <= len(b) <= BATCH_SIZE for b in batches)
    assert len(batches) == math.ceil(len(todo) / BATCH_SIZE)

    assert set(cached) == set(current)
    assert not set(cached) & set(batched)
    for a, path in cached.items():
        assert Path(path) == store_path(stats_dir, a)

    for a in stale:
        assert not store_path(stats_dir, a).exists(), "stale output was not deleted"
