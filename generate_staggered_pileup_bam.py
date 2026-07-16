#!/usr/bin/env python3
"""Generate a synthetic staggered, non-merged (mixed M_/F_/R_) pileup BAM for benchmarking
DeDup's non-merged-mode candidate-narrowing fix: many distinct start positions across a
narrow band, each with several reads piled on (a few distinct lengths per position), so the
whole band stays buffered in RMDupper's sliding window at once (the first read's end exceeds
the last read's start throughout) -- the non-merged analogue of generate_big_pileup_bam.py's
single-locus shape, but with real position staggering, which is what the old full-buffer
scan in checkForDuplication() (before RecordBufferHeap-based candidate narrowing) paid
O(depth) for on every anchor resolution.

Read length is fixed and decoupled from --reads (depth is scaled via reads-per-position with
--positions/--length held constant) -- this matters for benchmarking: an earlier attempt
scaled read length together with read count to keep the window open, which made total BAM
I/O volume scale as O(N^2) and completely masked the dedup algorithm's own complexity. Fixing
read length and scaling depth via reads-per-position keeps I/O volume O(N), so wall-clock
time actually isolates the algorithm being benchmarked.

Usage:
    python3 generate_staggered_pileup_bam.py -o staggered_pileup.bam \
        --reads 160000 --positions 300 --length 400 --seed 7
"""
import argparse
import array
import random

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="output BAM path")
    ap.add_argument("--reads", type=int, default=20000,
                     help="approximate total reads (rounded down to a multiple of --positions)")
    ap.add_argument("--positions", type=int, default=300,
                     help="number of distinct start positions spanning the band")
    ap.add_argument("--length", type=int, default=400,
                     help="base read length; must exceed --positions to keep the whole band buffered")
    ap.add_argument("--start", type=int, default=10_001,
                     help="1-based start of the first position")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    if args.length <= args.positions:
        raise SystemExit("--length must exceed --positions to keep the whole band buffered")

    rng = random.Random(args.seed)
    contig_len = args.start + args.positions + args.length + 1000
    prefixes = ["M_", "F_", "R_"]

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": contig_len}],
    }

    reads_per_pos = max(1, args.reads // args.positions)
    total = reads_per_pos * args.positions

    with pysam.AlignmentFile(args.out, "wb", header=header) as af:
        idx = 0
        for p in range(args.positions):
            start0 = args.start - 1 + p  # pysam is 0-based
            for r in range(reads_per_pos):
                length = args.length - (r % 5)  # a few distinct lengths per position
                seq = "".join(rng.choices("ACGT", k=length))
                rec = pysam.AlignedSegment()
                rec.query_name = f"{prefixes[idx % 3]}p{p}_r{r}"
                rec.query_sequence = seq
                rec.flag = 16 if (idx % 2) else 0
                rec.reference_id = 0
                rec.reference_start = start0
                rec.mapping_quality = 37
                rec.cigartuples = [(0, length)]
                rec.query_qualities = array.array("B", [rng.randint(2, 40) for _ in range(length)])
                af.write(rec)
                idx += 1

    print(f"Done: {total:,} reads written to {args.out} "
          f"({args.positions:,} positions x ~{reads_per_pos}/position, fixed length~{args.length}, "
          f"chr1:{args.start}-{args.start + args.positions})")


if __name__ == "__main__":
    main()
