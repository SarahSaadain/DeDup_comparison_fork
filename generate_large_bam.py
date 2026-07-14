#!/usr/bin/env python3
"""Generate a synthetic coordinate-sorted BAM for large-scale DeDup benchmarking.

DeDup only cares about read-name prefix (M_/F_/R_), reference position,
strand, CIGAR-derived length, and mapping quality — it never inspects the
sequence or mate fields. So instead of simulating a real reference and
aligner, this script fabricates a header with N contigs and writes reads
directly, calibrated to hit a target read count and PCR-duplicate rate.

Reads are generated as "clusters": one genomic locus shared by 1+ reads
(a duplicate stack), sorted by ascending position within each contig and by
contig order in the header, so the output is already coordinate-sorted
without a separate samtools sort pass.

Usage:
    python3 generate_large_bam.py -o large_10M_15contigs.bam \
        --reads 10000000 --contigs 15 --seed 42
"""
import argparse
import array
import random
import sys

import pysam

# Real human chr1-15 lengths (GRCh38) — used only as relative proportions so
# contigs vary in size the way real chromosomes do; the genome is otherwise
# synthetic (no real sequence).
CHR_RATIOS = [
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
]

BASES = "ACGT"


def contig_lengths(n_contigs, genome_length):
    ratios = (CHR_RATIOS * ((n_contigs // len(CHR_RATIOS)) + 1))[:n_contigs]
    ratio_sum = sum(ratios)
    return [max(1_000_000, round(genome_length * r / ratio_sum)) for r in ratios]


def build_header(n_contigs, genome_length):
    lengths = contig_lengths(n_contigs, genome_length)
    sq = [{"SN": f"chr{i + 1}", "LN": ln} for i, ln in enumerate(lengths)]
    return {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}, lengths


def _geometric(rng, p):
    # random.random() < p support starting at 1 (mean = 1/p)
    k = 1
    while rng.random() >= p:
        k += 1
    return k


def make_read(name, tid, pos, length, reverse, mapq, seq_pool, qual_pool, seq_off):
    r = pysam.AlignedSegment()
    r.query_name = name
    r.query_sequence = seq_pool[seq_off:seq_off + length]
    r.flag = 16 if reverse else 0
    r.reference_id = tid
    r.reference_start = pos
    r.mapping_quality = mapq
    r.cigartuples = [(0, length)]
    r.query_qualities = qual_pool[:length]
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", required=True, help="output BAM path")
    ap.add_argument("--reads", type=int, default=10_000_000, help="approximate total read count")
    ap.add_argument("--contigs", type=int, default=15, help="number of contigs")
    ap.add_argument("--genome-length", type=int, default=300_000_000,
                     help="total synthetic genome length in bp, split across contigs")
    ap.add_argument("--avg-advance", type=float, default=150.0,
                     help="target average bp advance per duplicate cluster; "
                          "lower values push the duplication rate up to fit "
                          "the read count into the genome")
    ap.add_argument("--p-merged", type=float, default=0.5,
                     help="fraction of clusters that are M_ (merged) vs F_/R_ paired")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    header, lengths = build_header(args.contigs, args.genome_length)
    total_length = sum(lengths)

    seq_pool = "".join(rng.choice(BASES) for _ in range(2000))
    qual_pool = array.array("B", (rng.randint(20, 40) for _ in range(200)))

    stats = {"M": 0, "F": 0, "R": 0, "clusters": 0}
    read_id = 0

    with pysam.AlignmentFile(args.out, "wb", header=header) as af:
        for tid, length in enumerate(lengths):
            target_reads_for_contig = round(args.reads * length / total_length)
            # mean_reads_per_step = p_merged * E[cluster] + (1-p_merged) * 2 * E[cluster]
            factor = args.p_merged + (1 - args.p_merged) * 2

            cursor = 0
            reads_in_contig = 0
            while cursor < length - 300 and reads_in_contig < target_reads_for_contig:
                # Recompute the required duplication density from what's actually left
                # to place, so the calibration self-corrects instead of drifting off
                # the target read count as read_len/gap variance accumulates.
                remaining_reads = target_reads_for_contig - reads_in_contig
                remaining_length = max(1, length - cursor)
                mean_reads_per_step = (remaining_reads / remaining_length) * args.avg_advance
                mean_cluster = max(1.01, mean_reads_per_step / factor)
                cluster_size = min(2000, _geometric(rng, min(0.98, max(0.02, 1.0 / mean_cluster))))
                stats["clusters"] += 1

                if rng.random() < args.p_merged:
                    read_len = rng.randint(30, 150)
                    read_len = min(read_len, length - cursor - 1)
                    reverse_cluster = rng.random() < 0.5
                    for _ in range(cluster_size):
                        read_id += 1
                        name = f"M_r{read_id}_c{stats['clusters']}"
                        mapq = rng.randint(1, 60)
                        seq_off = rng.randint(0, len(seq_pool) - read_len)
                        r = make_read(name, tid, cursor, read_len,
                                      reverse_cluster, mapq, seq_pool, qual_pool, seq_off)
                        af.write(r)
                        stats["M"] += 1
                        reads_in_contig += 1
                    advance = read_len
                else:
                    read_len_f = rng.randint(30, 150)
                    insert = read_len_f + rng.randint(20, 150)
                    insert = min(insert, length - cursor - 1)
                    read_len_r = min(rng.randint(30, 150), max(1, insert - 1))
                    r_pos = cursor + insert - read_len_r
                    # Emit all F_ reads (pos=cursor) before any R_ reads (pos=r_pos >=
                    # cursor) so the file stays coordinate-sorted within the cluster.
                    for _ in range(cluster_size):
                        read_id += 1
                        seq_off_f = rng.randint(0, len(seq_pool) - read_len_f)
                        rf = make_read(f"F_r{read_id}_c{stats['clusters']}", tid,
                                       cursor, read_len_f, False, rng.randint(1, 60),
                                       seq_pool, qual_pool, seq_off_f)
                        af.write(rf)
                        stats["F"] += 1
                        reads_in_contig += 1
                    for _ in range(cluster_size):
                        read_id += 1
                        seq_off_r = rng.randint(0, len(seq_pool) - read_len_r)
                        rr = make_read(f"R_r{read_id}_c{stats['clusters']}", tid,
                                       r_pos, read_len_r, True, rng.randint(1, 60),
                                       seq_pool, qual_pool, seq_off_r)
                        af.write(rr)
                        stats["R"] += 1
                        reads_in_contig += 1
                    advance = insert

                gap = rng.randint(0, max(1, int(args.avg_advance)))
                cursor += advance + gap

            print(f"  chr{tid + 1}: length={length:,} target={target_reads_for_contig:,} "
                  f"written={reads_in_contig:,}", file=sys.stderr)

    total = stats["M"] + stats["F"] + stats["R"]
    print(f"\nDone: {total:,} reads written to {args.out}", file=sys.stderr)
    print(f"  M_={stats['M']:,} F_={stats['F']:,} R_={stats['R']:,} "
          f"clusters={stats['clusters']:,} "
          f"mean_cluster_size={total / stats['clusters']:.2f}", file=sys.stderr)


if __name__ == "__main__":
    main()
