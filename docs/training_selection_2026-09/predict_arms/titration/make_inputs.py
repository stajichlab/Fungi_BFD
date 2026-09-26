#!/usr/bin/env python
"""getBestModel (code_new4) on each genome's r1_r2 PASA output, then the pool of
complete models on the TRAIN chromosomes, which subsets are drawn from.
Writes inputs/<genome>.best.gff3, inputs/<genome>.pool.gff3, availability.tsv, tasks.tsv."""
import os
from titration_lib import A, X, T, GENOMES, NS, draws_for, lib, train, gzip, shutil

os.makedirs(f"{T}/inputs", exist_ok=True)
avail = {}
for g, asm in GENOMES.items():
    d = f"{X}/r1_r2/{g}.rust"
    full = f"{A}/{g}/genome.fa"
    best = f"{T}/inputs/{g}.best.gff3"
    train.getBestModel(f"{d}/pasa.step1.gff3", full, f"{d}/kallisto.tsv", best, pasa_alignment_overlap=30)
    keep = set(open(f"{A}/{g}/split.train_chroms.txt").read().split())
    G = {}
    G = lib.gff2dict(best, full, G)
    pool = {k: v for k, v in G.items() if v["contig"] in keep and lib.is_complete_model(v)}
    lib.dict2gff3(pool, f"{T}/inputs/{g}.pool.gff3")
    avail[g] = len(pool)
    print(g, "getBestModel", len(G), "complete on train chroms", len(pool), flush=True)
with open(f"{T}/availability.tsv", "w") as f:
    f.write("genome\tcomplete_on_train_chroms\n")
    for g, n in avail.items():
        f.write(f"{g}\t{n}\n")
with open(f"{T}/tasks.tsv", "w") as f:
    for g, n_av in avail.items():
        for n in NS:
            if n > n_av:
                continue
            for d in range(1, draws_for(n) + 1):
                f.write(f"{g}\t{n}\t{d}\n")
print("tasks:", sum(1 for _ in open(f"{T}/tasks.tsv")))
