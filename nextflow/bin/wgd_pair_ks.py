#!/usr/bin/env python3
"""Pair-level dS for paralog families (alternative to `wgd ksd`).

For every gene pair inside each wgd dmd family:
  1. global BLOSUM62 protein alignment of the two proteins (no family-wide MSA),
  2. back-translation to a codon alignment, protein identity, NG86 dS (cheap, no codeml),
  3. pre-filter: pairs below --min-identity (or --min-ng-codons aligned codons, or above
     --max-ngds when > 0) are written without a codeml estimate,
  4. pairs that pass: one pairwise codeml (runmode -2) on the gap-free codon alignment.

Output: TSV (gzip when the name ends in .gz) with the wgd ks.tsv columns that
merge_wgd_ks.py reads (pair, family, g1, g2, gene1, gene2, N, S, dN, dS, dN/dS, t,
alignmentlength, strippedalignmentlength, alignmentcoverage, alignmentidentity) plus
ng_dS, ng_codons, status, node. N+S = 3 x codons used by codeml.

node: wgd mix and build_wgd_ksd_summary.py average dS over (family, node). Here the node of a
pair is its last common ancestor in a neighbour-joining tree of the family built on protein
identity distance (1 - identity, every pair, no family alignment), midpoint-rooted. Pairs
without an identity use the family maximum distance. Validation on simulated families with a
known tree: eval/sim_controls/report.md (section "Node-method validation"). The node labels
are only moderately accurate (ARI 0.55-0.62 vs the true nodes); the node-mean dS distribution
is close to the truth (Wasserstein-1 0.045-0.057 vs 0.41-0.43 with no averaging).

Identical sequences: genes of one family with a byte-identical CDS are collapsed to one
representative (--no-collapse turns this off). Only representative pairs are aligned and run
through codeml; their dS applies to every gene pair between the two groups (copies1 x copies2 =
npairs). Each group with 2 or more copies gets one row with status "identical" (dS = 0, no
computation; npairs = n(n-1)/2). The mapping from representative to member genes is written to
<out>.groups.tsv.gz. Without --no-collapse the output has one row per sequence pair, not per
gene pair; the sum of npairs equals the number of gene pairs.

Node sampling (--node-sample K, K > 0): wgd mix style statistics average dS per (family, node),
so every pair of a node measures the same duplication event. Here every pair is aligned and
scored for identity (cheap), the NJ tree gives every pair its node, and pairs of each node
(among the pairs that pass the pre-filter, in a seeded random order, so a rerun gives the same
pairs) go to codeml in rounds until the node has K gated pairs (dS <= --gate-max-ds and at least
--gate-min-codons codons) or --node-max-draws pairs were run or no pair is left. Without the
rounds a node whose first K pairs all fall outside the gate would lose its value (Glomus: 2%
of the nodes and a 3% lower median; with up to 20 draws: 0.5% and 1%). The sampled pairs are written as rows with node, node_pairs (gene pairs in
the node) and node_sampled; all nodes are written to <out>.nodes.tsv.gz. Unsampled pairs have
no row. K = 0 (default) keeps the all-pairs output. Tested offline on Glomus: K = 5 reproduces
the node-mean distribution closely (eval/sampling_test/offline_k_glomus.txt).

Exit codes: 0 ok; 2 empty/invalid families file; 3 more than 1% of pairs have missing genes;
4 more than 5% of attempted codeml runs failed (timeout/no_out/no_estimate).

Evaluation: eval/sim_controls/report.md (section "Cap, re-cluster and pair-level
pre-filter tests"). Run it inside the wgd container (Biopython 1.76, codeml in PATH).
Standard genetic code by default; --table / --icode are the hook for other codes
(not validated).

Usage: wgd_pair_ks.py --cds cds.fa[.gz] --families fam.tsv --out out.ks.tsv.gz --threads 8
"""
import argparse, gzip, itertools, math, os, re, shutil, subprocess, sys, tempfile
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--cds", required=True)
ap.add_argument("--families", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--min-identity", type=float, default=0.4, help="pre-filter: protein identity over aligned columns")
ap.add_argument("--max-ngds", type=float, default=0.0, help="pre-filter: NG86 dS upper bound (0 = off)")
ap.add_argument("--min-ng-codons", type=int, default=30, help="pre-filter: aligned codons")
ap.add_argument("--table", type=int, default=1, help="NCBI genetic code for translation / NG86")
ap.add_argument("--icode", default="0", help="codeml icode")
ap.add_argument("--tmpdir", default=os.environ.get("TMPDIR", "/tmp"))
ap.add_argument("--node-sample", type=int, default=0, help="pairs per node sent to codeml (0 = all pairs)")
ap.add_argument("--gate-max-ds", type=float, default=2.5, help="node sampling: a pair counts as gated when dS <= this (0 = no limit)")
ap.add_argument("--gate-min-codons", type=int, default=30, help="node sampling: a pair counts as gated with at least this many codons used by codeml")
ap.add_argument("--node-max-draws", type=int, default=20, help="node sampling: at most this many pairs per node go to codeml")
ap.add_argument("--seed", type=int, default=1, help="seed for node sampling")
ap.add_argument("--no-collapse", action="store_true", help="do not collapse identical CDS")
ap.add_argument("--max-pairs", type=int, default=0, help="debug: stop after this many pairs")
A = ap.parse_args()

TABLE = CodonTable.unambiguous_dna_by_id[A.table]
FWD = dict(TABLE.forward_table)
STOPS = set(TABLE.stop_codons)
NT = "ACGT"
ALL_CODONS = ["".join(p) for p in itertools.product(NT, repeat=3)]


def syn_sites(codon):
    if codon in STOPS or codon not in FWD:
        return None
    aa = FWD[codon]
    s = 0.0
    for i in range(3):
        n_syn = 0
        for b in NT:
            if b == codon[i]:
                continue
            c2 = codon[:i] + b + codon[i + 1:]
            if c2 in STOPS:
                continue
            if FWD[c2] == aa:
                n_syn += 1
        s += n_syn / 3.0
    return s


SYN = {c: syn_sites(c) for c in ALL_CODONS}


def diff_counts(c1, c2):
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    paths = []
    for order in itertools.permutations(pos):
        cur = c1
        sd = nd = 0
        ok = True
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            if nxt in STOPS:
                ok = False
                break
            if FWD[cur] == FWD[nxt]:
                sd += 1
            else:
                nd += 1
            cur = nxt
        if ok:
            paths.append((sd, nd))
    if not paths:
        return float("nan"), float("nan")
    return sum(p[0] for p in paths) / len(paths), sum(p[1] for p in paths) / len(paths)


def jc(p):
    if p >= 0.75:
        return float("inf")
    return -0.75 * math.log(1 - 4 * p / 3)


def ng86_dS(cod1, cod2):
    """NG86 synonymous distance (Jukes-Cantor corrected) and number of usable codons."""
    S = Sd = 0.0
    L = 0
    for a, b in zip(cod1, cod2):
        if SYN.get(a) is None or SYN.get(b) is None:
            continue
        sd, nd = diff_counts(a, b)
        if math.isnan(sd):
            continue
        S += (SYN[a] + SYN[b]) / 2
        Sd += sd
        L += 1
    if S <= 0:
        return L, float("nan")
    return L, jc(Sd / S)


cds = {}
opener = gzip.open if A.cds.endswith(".gz") else open
with opener(A.cds, "rt") as fh:
    for rec in SeqIO.parse(fh, "fasta"):
        s = str(rec.seq).upper()
        cds[rec.id] = s[: len(s) - len(s) % 3]

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = "global"

CTL = """seqfile = aln.phy
outfile = out
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = {icode}
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = .4
cleandata = 1
"""
RX = re.compile(r"t=\s*([\d.]+)\s+S=\s*([\d.]+)\s+N=\s*([\d.]+)\s+dN/dS=\s*([\d.]+)\s+dN =\s*([\d.]+)\s+dS =\s*([\d.]+)")
HDR = ["pair", "family", "g1", "g2", "gene1", "gene2", "N", "S", "dN", "dS", "dN/dS", "t",
       "alignmentlength", "strippedalignmentlength", "alignmentcoverage", "alignmentidentity",
       "ng_dS", "ng_codons", "status", "node", "copies1", "copies2", "npairs", "node_pairs", "node_sampled"]


def fmt(x):
    if isinstance(x, float):
        if x != x:
            return "NaN"
        if x == float("inf"):
            return "inf"
        return f"{x:.6g}"
    return str(x)


def work(job):
    fam, g1, g2, n1, n2, kind = job
    pair = f"{g1}__{g2}"
    row = {"pair": pair, "family": fam, "g1": g1, "g2": g2, "gene1": g1, "gene2": g2,
           "copies1": float(n1), "copies2": float(n2), "node_pairs": float("nan"), "node_sampled": float("nan"),
           "npairs": float(n1 * (n1 - 1) // 2) if kind == "ident" else float(n1 * n2)}
    nan = float("nan")
    for k in ("N", "S", "dN", "dS", "dN/dS", "t", "alignmentlength", "strippedalignmentlength",
              "alignmentcoverage", "alignmentidentity", "ng_dS", "ng_codons"):
        row[k] = nan
    if kind == "ident":
        L = float(len(cds[g1]))
        row.update({"status": "identical", "dN": 0.0, "dS": 0.0, "t": 0.0, "alignmentlength": L,
                    "strippedalignmentlength": L, "alignmentcoverage": 1.0, "alignmentidentity": 1.0,
                    "ng_dS": 0.0, "ng_codons": L / 3.0})
        return row
    if g1 not in cds or g2 not in cds:
        row["status"] = "missing"
        return row
    s1, s2 = cds[g1], cds[g2]
    p1 = str(Seq(s1).translate(table=A.table))
    p2 = str(Seq(s2).translate(table=A.table))
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    p1 = p1.replace("*", "X")
    p2 = p2.replace("*", "X")
    if not p1 or not p2:
        row["status"] = "empty"
        return row
    aln = aligner.align(p1, p2)[0]
    c1, c2, ncols, nident = [], [], 0, 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            i, j = a0 + k, b0 + k
            ncols += 1
            if p1[i] == p2[j]:
                nident += 1
            x = s1[3 * i: 3 * i + 3]
            y = s2[3 * j: 3 * j + 3]
            if len(x) == 3 == len(y) and set(x) <= set(NT) and set(y) <= set(NT) and x not in STOPS and y not in STOPS:
                c1.append(x)
                c2.append(y)
    row["alignmentlength"] = float(3 * ncols)
    row["alignmentcoverage"] = ncols / max(len(p1), len(p2)) if ncols else nan
    ident = nident / ncols if ncols else nan
    row["alignmentidentity"] = ident
    L, ngds = ng86_dS(c1, c2)
    row["ng_codons"] = float(L)
    row["ng_dS"] = ngds
    if not (ident >= A.min_identity) or L < A.min_ng_codons or (A.max_ngds > 0 and not (ngds <= A.max_ngds)):
        row["status"] = "prefiltered"
        return row
    n = len(c1)
    if n < 10:
        row["status"] = "too_short"
        return row
    d = tempfile.mkdtemp(prefix="pk_", dir=A.tmpdir)
    try:
        with open(f"{d}/aln.phy", "w") as fh:
            fh.write(f" 2 {3 * n}\ns1  {''.join(c1)}\ns2  {''.join(c2)}\n")
        with open(f"{d}/c.ctl", "w") as fh:
            fh.write(CTL.format(icode=A.icode))
        try:
            subprocess.run(["codeml", "c.ctl"], cwd=d, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=60)
        except subprocess.TimeoutExpired:
            row["status"] = "timeout"
            return row
        try:
            txt = open(f"{d}/out").read()
        except FileNotFoundError:
            row["status"] = "no_out"
            return row
        m = RX.search(txt)
        if not m:
            row["status"] = "no_estimate"
            return row
        t, S, N, w, dN, dS = (float(x) for x in m.groups())
        row.update({"t": t, "S": S, "N": N, "dN/dS": w, "dN": dN, "dS": dS,
                    "strippedalignmentlength": float(3 * n), "status": "ok"})
        return row
    finally:
        shutil.rmtree(d, ignore_errors=True)


def nj_midpoint_nodes(genes, dist):
    """Neighbour-joining tree of `genes` from the distance matrix `dist` (numpy, n x n), rooted at
    the midpoint of the longest tip-to-tip path. Returns {(i, j): node_label} for i < j (tip indexes)."""
    import numpy as np
    n = len(genes)
    if n < 2:
        return {}
    if n == 2:
        return {(0, 1): "n0"}
    adj = {i: {} for i in range(n)}  # node -> {neighbour: length}
    D = {i: {j: float(dist[i, j]) for j in range(n) if j != i} for i in range(n)}
    active = list(range(n))
    nxt = n
    while len(active) > 3:
        m = len(active)
        r = {i: sum(D[i][j] for j in active if j != i) for i in active}
        best, bi, bj = None, None, None
        for a in range(m):
            i = active[a]
            for b in range(a + 1, m):
                j = active[b]
                q = (m - 2) * D[i][j] - r[i] - r[j]
                if best is None or q < best:
                    best, bi, bj = q, i, j
        dij = D[bi][bj]
        li = 0.5 * dij + (r[bi] - r[bj]) / (2.0 * (m - 2))
        lj = dij - li
        u = nxt
        nxt += 1
        adj[u] = {bi: li, bj: lj}
        adj[bi][u] = li
        adj[bj][u] = lj
        D[u] = {}
        for k in active:
            if k in (bi, bj):
                continue
            d = 0.5 * (D[bi][k] + D[bj][k] - dij)
            D[u][k] = d
            D[k][u] = d
        for k in (bi, bj):
            for x in list(D[k]):
                D[x].pop(k, None)
            del D[k]
        active = [k for k in active if k not in (bi, bj)] + [u]
    a, b, c = active
    la = 0.5 * (D[a][b] + D[a][c] - D[b][c])
    lb = 0.5 * (D[a][b] + D[b][c] - D[a][c])
    lc = 0.5 * (D[a][c] + D[b][c] - D[a][b])
    center = nxt
    adj[center] = {a: la, b: lb, c: lc}
    adj[a][center] = la
    adj[b][center] = lb
    adj[c][center] = lc

    def far(src):
        dd, prev, st = {src: 0.0}, {src: None}, [src]
        while st:
            x = st.pop()
            for y, w in adj[x].items():
                if y not in dd:
                    dd[y] = dd[x] + w
                    prev[y] = x
                    st.append(y)
        t = max((i for i in range(n) if i != src), key=lambda i: dd[i])  # farthest tip other than src
        return t, dd[t], prev

    t1, _, _ = far(0)
    t2, total, prev = far(t1)
    # walk from t2 back to t1 to find the midpoint edge
    path = [t2]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    path.reverse()  # t1 ... t2
    half, acc, root_edge = total / 2.0, 0.0, None
    for x, y in zip(path, path[1:]):
        w = adj[x][y]
        if acc + w >= half:
            root_edge = (x, y)
            break
        acc += w
    if root_edge is None:
        root_edge = (path[0], path[1])
    rx, ry = root_edge
    root = nxt + 1
    wxy = adj[rx].pop(ry)
    adj[ry].pop(rx)
    adj[root] = {rx: half - acc, ry: wxy - (half - acc)}
    adj[rx][root] = half - acc
    adj[ry][root] = wxy - (half - acc)
    # rooted tree: collect child tip sets per internal node (iterative post-order)
    out, tips = {}, {}
    order, st, parent = [], [root], {root: None}
    while st:
        x = st.pop()
        order.append(x)
        for y in adj[x]:
            if y != parent[x]:
                parent[y] = x
                st.append(y)
    for x in reversed(order):
        if x < n:
            tips[x] = [x]
            continue
        kids = [y for y in adj[x] if y != parent[x]]
        sets = [tips[k] for k in kids]
        label = f"n{x - n}"
        for p in range(len(sets)):
            for q in range(p + 1, len(sets)):
                for i in sets[p]:
                    for j in sets[q]:
                        out[(min(i, j), max(i, j))] = label
        tips[x] = [i for sset in sets for i in sset]
    return out


def assign_nodes(fam_rows):
    """fam_rows: list of row dicts of one family. Sets row['node'] (NJ on identity distance)."""
    import numpy as np
    for r in fam_rows:
        if r["status"] == "identical":
            r["node"] = "i_" + r["g1"]  # identical copies form their own node (dS = 0)
    prs = [r for r in fam_rows if r["status"] != "identical"]
    if not prs:
        return
    fam_rows = prs
    genes = sorted({r["g1"] for r in fam_rows} | {r["g2"] for r in fam_rows})
    ix = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    d = np.full((n, n), np.nan)
    for r in fam_rows:
        idn = r["alignmentidentity"]
        if idn == idn:
            d[ix[r["g1"]], ix[r["g2"]]] = d[ix[r["g2"]], ix[r["g1"]]] = 1.0 - idn
    fill = np.nanmax(d) if np.isfinite(d).any() else 1.0
    d[np.isnan(d)] = fill
    np.fill_diagonal(d, 0.0)
    lab = nj_midpoint_nodes(genes, d)
    for r in fam_rows:
        i, j = ix[r["g1"]], ix[r["g2"]]
        r["node"] = lab.get((min(i, j), max(i, j)), "")


GROUPS = []  # (family, member genes) of collapsed groups, filled by jobs()


def jobs():
    n = 0
    with open(A.families) as fh:
        hdr = next(fh, None)  # header: <TAB>genome file name
        if hdr is None:
            print("ERROR: families file is empty: " + A.families, file=sys.stderr)
            sys.exit(2)
        for line in fh:
            line = line.rstrip("\n")
            if "\t" not in line:
                continue
            fam, genes = line.split("\t", 1)
            genes = [g.strip() for g in genes.split(",") if g.strip()]
            if A.no_collapse:
                groups = [[g] for g in genes]
            else:
                by, groups = {}, []
                for g in genes:
                    key = cds.get(g)
                    if key is None:
                        groups.append([g])  # missing gene: its own group, reported as "missing"
                    elif key in by:
                        by[key].append(g)
                    else:
                        by[key] = [g]
                        groups.append(by[key])
            for a, b in itertools.combinations(groups, 2):
                yield (fam, a[0], b[0], len(a), len(b), "pair")
                n += 1
                if A.max_pairs and n >= A.max_pairs:
                    return
            for grp in groups:
                if len(grp) > 1:
                    GROUPS.append((fam, grp))
                    yield (fam, grp[0], grp[1], len(grp), len(grp), "ident")
                    n += 1


def flush(rows, w):
    if not rows:
        return
    assign_nodes(rows)
    for r in rows:
        w.write("\t".join(fmt(r[k]) for k in HDR) + "\n")



def phase1(pair):
    """Cheap pass for node sampling: protein identity and the number of usable codon columns."""
    g1, g2 = pair
    s1, s2 = cds[g1], cds[g2]
    p1 = str(Seq(s1).translate(table=A.table))
    p2 = str(Seq(s2).translate(table=A.table))
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    p1 = p1.replace("*", "X")
    p2 = p2.replace("*", "X")
    if not p1 or not p2:
        return float("nan"), 0
    aln = aligner.align(p1, p2)[0]
    ncols = nident = nc = 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0):
            i, j = a0 + k, b0 + k
            ncols += 1
            if p1[i] == p2[j]:
                nident += 1
            x = s1[3 * i: 3 * i + 3]
            y = s2[3 * j: 3 * j + 3]
            if len(x) == 3 == len(y) and set(x) <= set(NT) and set(y) <= set(NT) and x not in STOPS and y not in STOPS:
                nc += 1
    return (nident / ncols if ncols else float("nan")), nc


def nj_label_matrix(D):
    """Neighbour-joining tree from the distance matrix D (numpy, n x n), midpoint-rooted.
    Returns an n x n int32 matrix: entry (i, j) is the id (0..n-2) of the last common ancestor of tips
    i and j (diagonal -1). Vectorised with numpy, so a family of 3,000 genes takes seconds."""
    import numpy as np
    n = D.shape[0]
    lab = np.full((n, n), -1, dtype=np.int32)
    if n < 2:
        return lab
    if n == 2:
        lab[0, 1] = lab[1, 0] = 0
        return lab
    M = np.array(D, dtype=np.float64)
    ids = list(range(n))
    adj = {i: {} for i in range(n)}
    nxt = n
    while len(ids) > 3:
        m = len(ids)
        r = M.sum(axis=1)
        Q = (m - 2) * M - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = divmod(int(np.argmin(Q)), m)
        dij = M[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2.0 * (m - 2))
        lj = dij - li
        u = nxt
        nxt += 1
        adj[u] = {ids[i]: li, ids[j]: lj}
        adj[ids[i]][u] = li
        adj[ids[j]][u] = lj
        newrow = 0.5 * (M[i] + M[j] - dij)
        keep = np.ones(m, dtype=bool)
        keep[[i, j]] = False
        M2 = np.zeros((m - 1, m - 1))
        M2[:-1, :-1] = M[np.ix_(keep, keep)]
        M2[-1, :-1] = newrow[keep]
        M2[:-1, -1] = newrow[keep]
        M = M2
        ids = [ids[x] for x in range(m) if keep[x]] + [u]
    a, b, c = ids
    la = 0.5 * (M[0, 1] + M[0, 2] - M[1, 2])
    lb = 0.5 * (M[0, 1] + M[1, 2] - M[0, 2])
    lc = 0.5 * (M[0, 2] + M[1, 2] - M[0, 1])
    center = nxt
    adj[center] = {a: la, b: lb, c: lc}
    adj[a][center] = la
    adj[b][center] = lb
    adj[c][center] = lc

    def far(src):
        dd, prev, st = {src: 0.0}, {src: None}, [src]
        while st:
            x = st.pop()
            for y, w in adj[x].items():
                if y not in dd:
                    dd[y] = dd[x] + w
                    prev[y] = x
                    st.append(y)
        t = max((i for i in range(n) if i != src), key=lambda i: dd[i])
        return t, dd[t], prev

    t1, _, _ = far(0)
    t2, total, prev = far(t1)
    path = [t2]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    path.reverse()
    half, acc, root_edge = total / 2.0, 0.0, None
    for x, y in zip(path, path[1:]):
        w = adj[x][y]
        if acc + w >= half:
            root_edge = (x, y)
            break
        acc += w
    if root_edge is None:
        root_edge = (path[0], path[1])
    rx, ry = root_edge
    root = center + 1
    wxy = adj[rx].pop(ry)
    adj[ry].pop(rx)
    adj[root] = {rx: half - acc, ry: wxy - (half - acc)}
    adj[rx][root] = half - acc
    adj[ry][root] = wxy - (half - acc)
    order, st, parent = [], [root], {root: None}
    while st:
        x = st.pop()
        order.append(x)
        for y in adj[x]:
            if y != parent[x]:
                parent[y] = x
                st.append(y)
    tips, nid = {}, 0
    for x in reversed(order):
        if x < n:
            tips[x] = np.array([x], dtype=np.int64)
            continue
        kids = [y for y in adj[x] if y != parent[x]]
        sets = [tips.pop(k) for k in kids]
        for p in range(len(sets)):
            for q in range(p + 1, len(sets)):
                lab[np.ix_(sets[p], sets[q])] = nid
                lab[np.ix_(sets[q], sets[p])] = nid
        tips[x] = np.concatenate(sets)
        nid += 1
    return lab


def group_genes(genes):
    """Collapse genes with an identical CDS. Returns a list of member lists (first = representative)."""
    if A.no_collapse:
        return [[g] for g in genes]
    by, groups = {}, []
    for g in genes:
        key = cds.get(g)
        if key is None:
            groups.append([g])  # missing gene: its own group
        elif key in by:
            by[key].append(g)
        else:
            by[key] = [g]
            groups.append(by[key])
    return groups


BLOCK_PAIRS = 200000  # families are processed in blocks of about this many pairs so the pool stays busy


def run_sampled():
    import numpy as np, zlib
    K = A.node_sample
    if A.max_ngds > 0:
        print("ERROR: --node-sample does not support --max-ngds", file=sys.stderr)
        sys.exit(2)
    oopen = gzip.open if A.out.endswith(".gz") else open
    counts, gp = {}, {}
    nodes_hdr = ["family", "node", "kind", "node_pairs", "node_rep_pairs", "node_pass", "node_sampled", "node_gated"]

    def bump(status, rows, pairs):
        counts[status] = counts.get(status, 0) + rows
        gp[status] = gp.get(status, 0) + pairs

    def families(fh):
        for line in fh:
            line = line.rstrip("\n")
            if "\t" not in line:
                continue
            fam, glist = line.split("\t", 1)
            genes = [g.strip() for g in glist.split(",") if g.strip()]
            groups = group_genes(genes)
            for grp in groups:
                if len(grp) > 1:
                    GROUPS.append((fam, grp))
                    r = work((fam, grp[0], grp[1], len(grp), len(grp), "ident"))
                    r["node"] = "i_" + grp[0]
                    r["node_pairs"] = r["npairs"]
                    r["node_sampled"] = 0.0
                    bump("identical", 1, int(r["npairs"]))
                    out.write("\t".join(fmt(r[k]) for k in HDR) + "\n")
                    nout.write(f"{fam}\ti_{grp[0]}\tidentical\t{int(r['npairs'])}\t0\t0\t0\t0\n")
            reps = [g[0] for g in groups]
            copies = np.array([len(g) for g in groups], dtype=np.int64)
            have = np.array([x in cds for x in reps])
            if (~have).any():
                mp = int(sum(copies[a_] * copies[b_] for a_ in range(len(reps)) for b_ in range(a_ + 1, len(reps))
                             if not (have[a_] and have[b_])))
                bump("missing", mp, mp)
            sel = np.nonzero(have)[0]
            if len(sel) >= 2:
                yield fam, [reps[i] for i in sel], copies[sel]

    def process_block(block):
        # 1. cheap pass over every pair of the block (identity, usable codons)
        def gen():
            for fam, R, C in block:
                n = len(R)
                for a_ in range(n):
                    for b_ in range(a_ + 1, n):
                        yield R[a_], R[b_]
        res = pool.imap(phase1, gen(), chunksize=500)
        states = []
        for fam, R, C in block:
            n = len(R)
            ia, ib = np.triu_indices(n, 1)
            P = len(ia)
            ident = np.empty(P)
            ncod = np.empty(P, dtype=np.int32)
            for t in range(P):
                ident[t], ncod[t] = next(res)
            passed = (ident >= A.min_identity) & (ncod >= A.min_ng_codons)
            # 2. identity-distance NJ tree -> node of every pair
            D = np.zeros((n, n))
            dist = 1.0 - ident
            fill = np.nanmax(dist) if np.isfinite(dist).any() else 1.0
            dist = np.where(np.isnan(dist), fill, dist)
            D[ia, ib] = dist
            D[ib, ia] = dist
            lab = nj_label_matrix(D)
            nid = lab[ia, ib].astype(np.int64)
            nnode = int(nid.max()) + 1
            w = C[ia] * C[ib]
            # 3. pass-pairs of every node in a seeded random order; rank = position within the node
            pidx = np.nonzero(passed)[0]
            rng = np.random.default_rng((A.seed * 1000003 + zlib.crc32(fam.encode())) & 0xFFFFFFFF)
            order = np.lexsort((rng.random(len(pidx)), nid[pidx]))
            cand = pidx[order]
            cn = nid[cand]
            start = np.r_[True, cn[1:] != cn[:-1]] if len(cn) else np.zeros(0, dtype=bool)
            first = np.maximum.accumulate(np.where(start, np.arange(len(cn)), 0)) if len(cn) else np.zeros(0, dtype=np.int64)
            rank = np.arange(len(cn)) - first
            states.append({"fam": fam, "R": R, "C": C, "ia": ia, "ib": ib, "nid": nid, "nnode": nnode, "w": w,
                           "passed": passed, "cand": cand, "cn": cn, "rank": rank,
                           "npass": np.bincount(nid[passed], minlength=nnode),
                           "nrep": np.bincount(nid, minlength=nnode),
                           "npairs": np.bincount(nid, weights=w, minlength=nnode).astype(np.int64),
                           "drawn": np.zeros(nnode, dtype=np.int64), "gated": np.zeros(nnode, dtype=np.int64),
                           "rows": []})
        # 4. codeml in rounds: every node draws pairs until it has K gated pairs or MAXD draws
        MAXD = max(A.node_max_draws, K)
        while True:
            jobs_r, where = [], []
            for st in states:
                need = np.clip(K - st["gated"], 0, None)
                take = np.minimum(np.minimum(need, MAXD - st["drawn"]), st["npass"] - st["drawn"])
                take = np.clip(take, 0, None)
                if len(st["cand"]) == 0 or take.max() == 0:
                    continue
                sel = (st["rank"] >= st["drawn"][st["cn"]]) & (st["rank"] < (st["drawn"] + take)[st["cn"]])
                for t in st["cand"][sel].tolist():
                    jobs_r.append((st["fam"], st["R"][st["ia"][t]], st["R"][st["ib"][t]],
                                   int(st["C"][st["ia"][t]]), int(st["C"][st["ib"][t]]), "pair"))
                    where.append((st, t))
                st["drawn"] = st["drawn"] + take
            if not jobs_r:
                break
            for r, (st, t) in zip(pool.imap(work, jobs_r, chunksize=2), where):
                k_ = int(st["nid"][t])
                r["node"] = f"n{k_}"
                st["rows"].append((r, k_))
                ok_gate = (r["status"] == "ok" and r["dS"] > 0 and (A.gate_max_ds <= 0 or r["dS"] <= A.gate_max_ds)
                           and (r["N"] + r["S"]) / 3.0 >= A.gate_min_codons)
                if ok_gate:
                    st["gated"][k_] += 1
        # 5. write rows and the node table
        for st in states:
            fam = st["fam"]
            for r, k_ in st["rows"]:
                r["node_pairs"] = float(st["npairs"][k_])
                r["node_sampled"] = float(st["drawn"][k_])
                bump(r["status"], 1, int(r["npairs"]))
                out.write("\t".join(fmt(r[k]) for k in HDR) + "\n")
            sw = int(sum(r["npairs"] for r, _ in st["rows"]))
            bump("prefiltered", int((~st["passed"]).sum()), int(st["w"][~st["passed"]].sum()))
            bump("not_sampled", int(st["passed"].sum()) - len(st["rows"]), int(st["w"][st["passed"]].sum()) - sw)
            for k_ in range(st["nnode"]):
                nout.write(f"{fam}\tn{k_}\ttree\t{int(st['npairs'][k_])}\t{int(st['nrep'][k_])}\t{int(st['npass'][k_])}\t{int(st['drawn'][k_])}\t{int(st['gated'][k_])}\n")

    with open(A.families) as fh, Pool(A.threads) as pool, oopen(A.out, "wt") as out, \
            gzip.open(A.out + ".nodes.tsv.gz", "wt") as nout:
        hdr = next(fh, None)
        if hdr is None:
            print("ERROR: families file is empty: " + A.families, file=sys.stderr)
            sys.exit(2)
        out.write("\t".join(HDR) + "\n")
        nout.write("\t".join(nodes_hdr) + "\n")
        block, bp = [], 0
        for item in families(fh):
            n = len(item[1])
            block.append(item)
            bp += n * (n - 1) // 2
            if bp >= BLOCK_PAIRS:
                process_block(block)
                block, bp = [], 0
        if block:
            process_block(block)
    return counts, gp


def run_all_pairs():
    oopen = gzip.open if A.out.endswith(".gz") else open
    counts, gp = {}, {}
    cur, buf = None, []
    with Pool(A.threads) as pool, oopen(A.out, "wt") as out:
        out.write("\t".join(HDR) + "\n")
        for r in pool.imap(work, jobs(), chunksize=50):
            counts[r["status"]] = counts.get(r["status"], 0) + 1
            gp[r["status"]] = gp.get(r["status"], 0) + int(r["npairs"])
            if r["family"] != cur:
                flush(buf, out)
                cur, buf = r["family"], []
            buf.append(r)
        flush(buf, out)
    return counts, gp


counts, gp = run_sampled() if A.node_sample > 0 else run_all_pairs()
total = sum(counts.values())
print("rows by status:", counts, file=sys.stderr)
print("gene pairs represented, by status:", gp, "total", sum(gp.values()), file=sys.stderr)
with open(A.out + ".status.txt", "w") as fh:
    for k in sorted(counts):
        fh.write(f"{k}\t{counts[k]}\t{gp.get(k, 0)}\n")
if GROUPS:
    with gzip.open(A.out + ".groups.tsv.gz", "wt") as fh:
        fh.write("family\trepresentative\tn_copies\tmembers\n")
        for fam, grp in GROUPS:
            fh.write(f"{fam}\t{grp[0]}\t{len(grp)}\t{','.join(grp)}\n")
if total == 0:
    print("WARNING: no pairs (no family with 2 or more genes)", file=sys.stderr)
miss = counts.get("missing", 0)
if total and miss / total > 0.01:
    print(f"ERROR: {miss} of {total} pairs have genes missing from the CDS file", file=sys.stderr)
    sys.exit(3)
if miss:
    print(f"WARNING: {miss} pairs have genes missing from the CDS file", file=sys.stderr)
att = sum(counts.get(k, 0) for k in ("ok", "timeout", "no_out", "no_estimate"))
bad = sum(counts.get(k, 0) for k in ("timeout", "no_out", "no_estimate"))
if att and bad / att > 0.05:
    print(f"ERROR: {bad} of {att} codeml runs failed", file=sys.stderr)
    sys.exit(4)
