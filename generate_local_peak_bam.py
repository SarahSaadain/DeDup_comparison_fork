#!/usr/bin/env python3
"""Generate a synthetic single-locus "local peak" BAM for benchmarking DeDup's
*non-merged* (default, prefix-based) duplicate resolution: many distinct,
1bp-apart start positions inside one narrow window, each with several
duplicate F_ (forward) reads.

This is the non-merged-mode analogue of generate_big_pileup_bam.py's
merged-mode shape. In default mode, an F_-F_ duplicate match requires only
equal alignStart (see DuplicateLogic's buffer_read_one/maybed_read_one +
equal_alignment_start entries) — read length is irrelevant to F_-F_ matching,
so unlike merged mode (which keys off exact (start,end) and therefore needs
same-start-different-length reads to form distinct groups), non-merged
forward duplicate groups are already distinct per start position on their
own. That also means length is the only lever left to control *when* the
window resolves: if reads were short relative to the position spacing, each
position's group would resolve (and drain from recordBuffer) as soon as the
next position streamed in, and old and new code would both be trivially
O(total reads).

To reproduce the true O(depth^2) shape, --read-len defaults to
`positions + 100` — long enough that the very first read (at the first
position) still overlaps the very last position in the peak. Nothing can
resolve while positions keep streaming in, so the whole peak (all
positions * reads_per_position reads) piles into recordBuffer at once, same
as generate_big_pileup_bam.py's single shared start, and only drains
group-by-group at end of stream/chromosome. The old code rescanned the
*entire remaining buffer* on every one of those `positions` group
resolutions (O(depth) each, depth ~ positions * reads_per_position) instead
of narrowing to the ~reads_per_position reads that actually share the
resolving group's start — true O(positions^2 * reads_per_position). Pass an
explicit smaller --read-len to instead model a "reads resolve as they go"
shallow peak, where both old and new code are already linear.

Mirrors DeDup_fork's RecordBufferHeapDifferentialTest / the fork's
startIndex/endIndex narrowing added alongside RecordBufferHeap.

Usage:
    python3 generate_local_peak_bam.py -o local_peak.bam \
        --positions 5000 --reads-per-position 20 --read-len 100 --seed 42
"""
import argparse
import array
import random

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="output BAM path")
    ap.add_argument("--start", type=int, default=10_001,
                     help="1-based start position of the first read")
    ap.add_argument("--positions", type=int, default=5000,
                     help="number of distinct 1bp-apart start positions (= duplicate groups)")
    ap.add_argument("--reads-per-position", type=int, default=20,
                     help="duplicate F_ reads sharing each start position")
    ap.add_argument("--read-len", type=int, default=None,
                     help="length of every read; default positions+100 so the first "
                          "read still overlaps the last position (nothing resolves "
                          "until end of stream — the O(depth^2) shape). Pass a smaller "
                          "value to model a shallow peak that resolves as it streams.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.reads_per_position > 93:
        raise SystemExit("--reads-per-position must be <= 93 (distinct Phred qualities per group)")

    read_len = args.read_len if args.read_len is not None else args.positions + 100
    rng = random.Random(args.seed)
    contig_len = args.start + args.positions + read_len + 1000

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": contig_len}],
    }

    with pysam.AlignmentFile(args.out, "wb", header=header) as af:
        for p in range(args.positions):
            pos = args.start - 1 + p  # pysam is 0-based
            seq = "".join(rng.choice("ACGT") for _ in range(read_len))
            # Distinct qualities per group so every group has one unmistakable
            # highest-quality winner (matches generate_big_pileup_bam.py's convention).
            quals = rng.sample(range(1, 94), args.reads_per_position)
            for r in range(args.reads_per_position):
                rec = pysam.AlignedSegment()
                rec.query_name = f"F_p{p}_r{r}"
                rec.query_sequence = seq
                rec.flag = 0  # forward strand
                rec.reference_id = 0
                rec.reference_start = pos
                rec.mapping_quality = 37
                rec.cigartuples = [(0, read_len)]
                rec.query_qualities = array.array("B", [quals[r]] * read_len)
                af.write(rec)

    total = args.positions * args.reads_per_position
    print(f"Done: {total:,} reads written to {args.out} "
          f"({args.positions:,} positions x {args.reads_per_position}/position, "
          f"read_len={read_len}, chr1:{args.start}-{args.start + args.positions})")


if __name__ == "__main__":
    main()
