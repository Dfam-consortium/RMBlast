#!/usr/bin/env python3
"""Compare a faithful rmblastn tabular output against a --prelim-cull run.

Usage: cull_compare.py faithful.tsv culled.tsv

Both files must use an outfmt whose first four fields are
"score qseqid qstart qend" (further fields are included in exact-hit
comparison).  Reports exact-tuple hit differences and the base-pair
delta of the union query coverage (the masking-sensitivity metric).
"""
import sys


def load(path):
    hits = set()
    spans = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        hits.add(tuple(f))
        spans.append((int(f[2]) - 1, int(f[3])))
    return hits, spans


def union_bp(spans):
    total, end = 0, -1
    for s, e in sorted(spans):
        if s > end:
            total += e - s
            end = e
        elif e > end:
            total += e - end
            end = e
    return total


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    base_hits, base_spans = load(sys.argv[1])
    hits, spans = load(sys.argv[2])
    base_bp = union_bp(base_spans)
    bp = union_bp(spans)
    lost = base_hits - hits
    gained = hits - base_hits
    lost_scores = sorted(int(t[0]) for t in lost)
    med = lost_scores[len(lost_scores) // 2] if lost_scores else 0
    mx = lost_scores[-1] if lost_scores else 0
    print(f"faithful: {len(base_hits)} hits, {base_bp} bp query coverage")
    print(f"culled:   {len(hits)} hits, {bp} bp query coverage")
    print(f"lost={len(lost)} ({100*len(lost)/max(len(base_hits),1):.2f}%)  "
          f"gained={len(gained)}  "
          f"coverage_delta={base_bp-bp} bp ({100*(base_bp-bp)/max(base_bp,1):.3f}%)  "
          f"lost_score med/max={med}/{mx}")


if __name__ == "__main__":
    main()
