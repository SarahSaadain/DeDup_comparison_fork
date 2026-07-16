#!/usr/bin/env python3
"""Generate a synthetic single-locus "big pileup" BAM for benchmarking DeDup's
merged-mode fast path: many reads sharing one start position but with enough
distinct lengths to form many distinct (start,end) duplicate groups, plus
several duplicate reads per group.

This is the shape that made the old (pre-bucket-map) merged-mode duplicate
resolution go O(depth^2): the resolve trigger never fires while same-start
reads stream in (a later read's start never exceeds an earlier same-start
read's end), so the whole pileup piles into RMDupper's sliding-window buffer
at once and is only drained group-by-group at end of stream — the old code
rescanned the full buffer on every group resolution instead of doing a
bucket-map lookup.

Mirrors DeDup_fork's
MergedHighDepthLocusTest#mergedFastPath_bigLocalPileup_correctnessAndTiming,
scaled up further here for a real wall-clock comparison between the original
jar and the fork.

Usage:
    python3 generate_big_pileup_bam.py -o merged_big_pileup.bam \
        --groups 1000 --reads-per-group 40 --seed 42
"""
import argparse
import array
import random

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="output BAM path")
    ap.add_argument("--start", type=int, default=10_001,
                     help="shared 1-based start position for every read")
    ap.add_argument("--groups", type=int, default=1000,
                     help="number of distinct (start,end) groups (distinct lengths)")
    ap.add_argument("--reads-per-group", type=int, default=40,
                     help="duplicate reads per group")
    ap.add_argument("--min-len", type=int, default=30,
                     help="length of group 0; group i has length min_len + i")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.reads_per_group > 93:
        raise SystemExit("--reads-per-group must be <= 93 (distinct Phred qualities per group)")

    rng = random.Random(args.seed)
    max_len = args.min_len + args.groups - 1
    contig_len = args.start + max_len + 1000

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": contig_len}],
    }

    with pysam.AlignmentFile(args.out, "wb", header=header) as af:
        for g in range(args.groups):
            length = args.min_len + g
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            # Distinct qualities per group so every group has one unmistakable
            # highest-quality winner (matches the JUnit fixture this mirrors).
            quals = rng.sample(range(1, 94), args.reads_per_group)
            for r in range(args.reads_per_group):
                rec = pysam.AlignedSegment()
                rec.query_name = f"M_g{g}_r{r}"
                rec.query_sequence = seq
                rec.flag = 0
                rec.reference_id = 0
                rec.reference_start = args.start - 1  # pysam is 0-based
                rec.mapping_quality = 37
                rec.cigartuples = [(0, length)]
                rec.query_qualities = array.array("B", [quals[r]] * length)
                af.write(rec)

    total = args.groups * args.reads_per_group
    print(f"Done: {total:,} reads written to {args.out} "
          f"({args.groups:,} groups x {args.reads_per_group}/group, all at chr1:{args.start})")


if __name__ == "__main__":
    main()
