#!/usr/bin/env python
"""Draw N complete models from inputs/<genome>.pool.gff3 with a fixed seed.
Usage: subset.py GENOME N DRAW OUT_GFF3"""
import random, sys
from titration_lib import A, T, lib, seed_for

g, n, draw, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
G = {}
G = lib.gff2dict(f"{T}/inputs/{g}.pool.gff3", f"{A}/{g}/genome.fa", G)
ids = sorted(G)
pick = random.Random(seed_for(g, n, draw)).sample(ids, n)
lib.dict2gff3({k: G[k] for k in pick}, out)
print(len(pick))
