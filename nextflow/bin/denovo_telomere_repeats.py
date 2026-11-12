#!/usr/bin/env python3
"""De novo terminal tandem-repeat finder (no motif list).

For each scaffold end, read the terminal window inward from the edge (the 3'
end is reverse-complemented, so every end is read from the edge inward). For
each period 4..40 bp, find the longest run that starts within --tolerance bp of
the edge and keeps at least --purity match to its own sequence shifted by the
period. Keep the period with the longest run (smallest period on near-ties),
infer the repeat unit by majority vote per phase, and report the unit in a
canonical form (G-rich strand, minimum rotation). No motif is assumed.

A background rate is estimated by running the same detector on random internal
windows. A unit counts as a telomere-repeat candidate when it is found at
several scaffold ends and is rare at internal positions.

Outputs (the first is gzip-compressed when the name ends in .gz):
  --ends-out     one row per scaffold end examined that has a tandem tract
  --summary-out  one row per candidate unit with end counts, enrichment and
                 polarity (share of ends in one orientation), plus
                 a per-genome status row (unit = '.')

Irregular repeats (for example Saccharomyces TG1-3) have no fixed period. They
are flagged by strand composition (``ca_frac_100``) and reported in the
summary row as ``strand_biased_ends``.

Input may be plain or gzip FASTA.
"""

import argparse
import csv
import gzip
import random
import re
import sys
from collections import Counter, defaultdict

import numpy as np

COMP = str.maketrans("ACGTN", "TGCAN")
CODE = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate("ACGT"):
    CODE[ord(c)] = i
    CODE[ord(c.lower())] = i


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def primitive(u: str) -> str:
    for k in range(1, len(u)):
        if len(u) % k == 0 and u[:k] * (len(u) // k) == u:
            return u[:k]
    return u


def canonical(u: str) -> str:
    """G-rich strand, minimum rotation. TTAGGG and CCCTAA give the same label."""
    def gt_minus_ca(x):
        return sum(c in "GT" for c in x) - sum(c in "CA" for c in x)
    s = rc(u) if gt_minus_ca(rc(u)) > gt_minus_ca(u) else u
    return min(s[i:] + s[:i] for i in range(len(s)))


def _rotations(u: str) -> set:
    return {u[i:] + u[:i] for i in range(len(u))}


def orientation(observed: str, canon: str) -> str:
    """'fwd' if the observed unit is a rotation of canon, 'rev' if of its
    reverse complement, 'amb' if both, 'none' if neither."""
    f = canon in _rotations(observed)
    r = canon in _rotations(rc(observed))
    return "amb" if (f and r) else "fwd" if f else "rev" if r else "none"


def polarity(observed_units, canon: str):
    """Share of ends in the majority orientation, or None if no end has a
    defined orientation.

    Every end is read from the edge inward (the 3' end is reverse-complemented),
    so a telomere shows one orientation at all ends: ~1.0. Sequence that is not
    strand-polar (for example AT-rich low-complexity repeats at contig ends)
    gives ~0.5. With few ends the value is weakly informative.
    """
    c = Counter(orientation(o, canon) for o in observed_units)
    n = c["fwd"] + c["rev"]
    return max(c["fwd"], c["rev"]) / n if n else None


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    name, buf = None, []
    with open_text(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def detect(seq: str, pmin: int, pmax: int, purity: float, tol: int,
           min_len: int, kmer_fallback: bool = True):
    """Return (period, start, end, unit, purity, method) or None.

    ``seq`` is read from the edge inward, upper case. ``method`` is
    ``periodic`` (fixed-period match) or ``kmer_runs`` (fallback for repeats
    broken by insertions or deletions; ``purity`` is then the covered fraction
    of the first 300 bp).
    """
    r = _detect_periodic(seq, pmin, pmax, purity, tol, min_len)
    if r:
        return r + ("periodic",)
    if kmer_fallback:
        r = detect_kmer_runs(seq, pmin, min(pmax, 12), tol)
        if r:
            return r + ("kmer_runs",)
    return None


def _detect_periodic(seq: str, pmin: int, pmax: int, purity: float, tol: int,
                     min_len: int):
    a = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a)
    best = None
    for p in range(pmin, pmax + 1):
        if n < min_len + p:
            break
        eq = ((a[:-p] == a[p:]) & (a[:-p] < 4)).astype(np.int32)
        w = max(30, 2 * p)
        if len(eq) < w:
            continue
        cs = np.concatenate(([0], np.cumsum(eq)))
        frac = (cs[w:] - cs[:-w]) / w
        ok = frac >= purity
        if not ok[: tol + 1].any():
            continue
        i0 = int(np.argmax(ok[: tol + 1]))
        # Extend while windows stay pure; tolerate no gaps.
        bad = np.flatnonzero(~ok[i0:])
        run_end_win = i0 + (int(bad[0]) if len(bad) else len(ok) - i0)
        end = run_end_win + w + p - 1
        end = min(end, n)
        length = end - i0
        if length < min_len:
            continue
        # Prefer the longest run; on near-ties keep the smaller period.
        if best is None or length > best[0] * 1.05:
            best = (length, p, i0, end)
    if best is None:
        return None
    length, p, i0, end = best
    tract = seq[i0:end]
    unit = "".join(
        Counter(tract[k::p]).most_common(1)[0][0] for k in range(p)
    )
    unit = primitive(unit)
    if len(set(unit)) < 2:
        return None
    # Purity of the final tract against the inferred unit.
    best_match = 0.0
    for rot in range(len(unit)):
        r = (unit[rot:] + unit[:rot]) * (len(tract) // len(unit) + 2)
        m = sum(x == y for x, y in zip(tract, r)) / len(tract)
        best_match = max(best_match, m)
    return len(unit), i0, end, unit, round(best_match, 3)


def detect_kmer_runs(seq: str, kmin: int, kmax: int, tol: int,
                     span: int = 300, min_cover: float = 0.3,
                     min_copies: int = 3, min_unit: int = 4):
    """Fallback for repeats broken by insertions or deletions.

    Count the bases of the first ``span`` bp covered by exact tandem runs
    (>= ``min_copies`` copies) of one unit and its rotations/strand, and return
    the unit with the highest coverage when it is >= ``min_cover`` and the
    first run starts within ``tol`` bp of the edge.
    Returns (period, start, end, unit, coverage) or None.
    """
    head = seq[:span + tol]
    best = None
    for k in range(kmin, kmax + 1):
        pat = re.compile(r"([ACGT]{%d})\1{%d,}" % (k, min_copies - 1))
        cover = defaultdict(int)
        first = {}
        last = {}
        for m in pat.finditer(head):
            unit = primitive(m.group(1))
            if len(unit) < min_unit or len(set(unit)) < 2:
                continue
            cu = canonical(unit)
            lo, hi = m.start(), min(m.end(), span)
            if hi > lo:
                cover[cu] += hi - lo
            first.setdefault(cu, m.start())
            last[cu] = m.end()
        for cu, c in cover.items():
            if first[cu] <= tol and (best is None or c > best[0]):
                best = (c, cu, first[cu], last[cu])
    if best is None or best[0] / span < min_cover:
        return None
    c, cu, start, end = best
    return len(cu), start, end, cu, round(c / span, 3)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("--ends-out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--asmid", default=None)
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--tolerance", type=int, default=50,
                    help="max distance of tract start from the scaffold edge")
    ap.add_argument("--min-period", type=int, default=4)
    ap.add_argument("--max-period", type=int, default=40)
    ap.add_argument("--purity", type=float, default=0.85)
    ap.add_argument("--min-length", type=int, default=60)
    ap.add_argument("--min-scaffold", type=int, default=5000)
    ap.add_argument("--background", type=int, default=300,
                    help="number of random internal windows for the background")
    ap.add_argument("--min-ends", type=int, default=2,
                    help="ends needed before a unit is called a candidate")
    ap.add_argument("--min-enrichment", type=float, default=5.0)
    ap.add_argument("--no-kmer-fallback", action="store_true",
                    help="disable the indel-tolerant k-mer run fallback")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)

    asmid = args.asmid or args.genome.split("/")[-1].split(".")[0]
    W = args.window
    rng = random.Random(args.seed)

    scaffolds = []  # (name, seq)
    for name, seq in read_fasta(args.genome):
        if len(seq) >= args.min_scaffold:
            scaffolds.append((name, seq.upper()))

    det = dict(pmin=args.min_period, pmax=args.max_period,
               purity=args.purity, tol=args.tolerance, min_len=args.min_length,
               kmer_fallback=not args.no_kmer_fallback)

    ends_rows = []
    unit_ends = defaultdict(list)
    unit_obs = defaultdict(list)
    n_ends = 0
    biased = 0
    for name, seq in scaffolds:
        for end_name, win in (("5prime", seq[:W]), ("3prime", rc(seq[-W:]))):
            n_ends += 1
            head = win[:100]
            ca = (head.count("C") + head.count("A")) / max(1, len(head))
            gt = (head.count("G") + head.count("T")) / max(1, len(head))
            if len(head) == 100 and max(ca, gt) >= 0.9 and head.count("N") == 0:
                biased += 1
            r = detect(win, **det)
            if r:
                period, i0, end, unit, pur, method = r
                cu = canonical(unit)
                unit_ends[cu].append(end - i0)
                unit_obs[cu].append(unit)
                ends_rows.append([asmid, name, len(seq), end_name, cu, unit,
                                  period, i0, end, end - i0, pur,
                                  round(ca, 3), method])

    # Background: random internal windows read in a random orientation.
    long_scaf = [(n, s) for n, s in scaffolds if len(s) >= 4 * W]
    bg_hits = Counter()
    bg_n = 0
    if long_scaf:
        for _ in range(args.background):
            name, seq = rng.choice(long_scaf)
            pos = rng.randint(W, len(seq) - 2 * W)
            win = seq[pos:pos + W]
            if rng.random() < 0.5:
                win = rc(win)
            bg_n += 1
            r = detect(win, **det)
            if r:
                bg_hits[canonical(r[3])] += 1

    open_out = gzip.open if args.ends_out.endswith(".gz") else open
    with open_out(args.ends_out, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["ASMID", "scaffold", "scaffold_len", "end", "unit_canonical",
                    "unit_observed", "period", "start", "end_coord",
                    "tract_length", "purity", "ca_frac_100", "method"])
        w.writerows(ends_rows)

    cols = ["ASMID", "unit", "period", "n_ends", "n_ends_examined", "frac_ends",
            "median_tract_len", "bg_windows", "bg_hits", "enrichment",
            "polarity", "strand_biased_ends", "status"]
    rows = []
    for cu, lens in sorted(unit_ends.items(), key=lambda kv: -len(kv[1])):
        k = len(lens)
        bg = bg_hits.get(cu, 0)
        # Ratio of the end-window rate to the internal-window rate. With no
        # internal hit the ratio is undefined; report "bg0" and call on end count.
        end_rate = k / max(1, n_ends)
        if bg_n == 0 or bg == 0:
            enrich = None
            called = k >= args.min_ends
        else:
            enrich = end_rate / (bg / bg_n)
            called = k >= args.min_ends and enrich >= args.min_enrichment
        med = sorted(lens)[k // 2]
        pol = polarity(unit_obs[cu], cu)
        rows.append([asmid, cu, len(cu), k, n_ends, round(end_rate, 3), med,
                     bg_n, bg, "bg0" if enrich is None else round(enrich, 1),
                     "NA" if pol is None else round(pol, 3), biased, "candidate" if called else "weak"])
    called_any = any(r[-1] == "candidate" for r in rows)
    if called_any:
        status = "repeat_found"
    elif biased >= args.min_ends:
        status = "strand_biased_only"
    elif n_ends == 0:
        status = "no_scaffolds"
    else:
        status = "no_terminal_repeat"
    rows.append([asmid, ".", 0, 0, n_ends, 0, 0, bg_n, 0, "NA", "NA", biased, status])
    with open(args.summary_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)
    print(f"# denovo_telomere_repeats {asmid}: {n_ends} ends, "
          f"{len(ends_rows)} with tandem tract, status={status}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
