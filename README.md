# DeDup fork vs. original — comparison & benchmarks

Infrastructure for verifying that the fork ([`DeDup_fork/`](https://github.com/SarahSaadain/DeDup)) produces identical output to the
original ([`DeDup_original/`](https://github.com/apeltzer/DeDup)) and for measuring how much faster it is.

## Scripts

- `run_comparison.sh` — runs both jars against the fixture BAMs plus a generated large BAM and
  a generated big local pileup, diffs kept-read-name sets, and writes a report to
  `results/report_<timestamp>.md`. Pass `--no-large` to skip the large-BAM section, `--no-pileup`
  to skip the big-pileup section.
- `generate_large_bam.py` — fabricates a coordinate-sorted synthetic BAM directly via `pysam`
  (no real reference; DeDup only reads position/strand/length/MAPQ/name-prefix).
  `--reads` / `--contigs` / `--genome-length` control scale; duplicate-cluster density is
  self-calibrated per contig to hit the target read count (see caveat below).
- `generate_big_pileup_bam.py` — fabricates a single-locus "big pileup" BAM: many reads sharing
  one start position with enough distinct lengths to form many distinct (start,end) duplicate
  groups (`--groups`, default 1000), each with several duplicate reads (`--reads-per-group`,
  default 40). This is the O(depth²) shape merged mode's bucket-map fast path targets — the
  resolve trigger never fires while same-start reads stream in, so the whole pileup buffers at
  once and only drains group-by-group at end of stream, where the old code rescanned the full
  buffer per group. Mirrors `DeDup_fork`'s `MergedHighDepthLocusTest` big-pileup unit test, scaled
  up for a real wall-clock comparison. `run_comparison.sh` always runs this scenario with `-m`.
- `scale_benchmark.sh [scale ...]` — runs orig / fork-seq / fork-par once each across a list of
  read-count scales, sharing one cached genome per scale. Default scales:
  `100000 1000000 5000000 10000000 25000000 50000000 100000000`. Each scale is benchmarked in
  **both dedup modes** — `default` (prefix-based M_/F_/R_) and `merged` (`-m`, both ends for
  every read) — since merged mode is the path the fork optimizes (O(depth²) → O(group)). The
  TSV carries a `mode` column; override the set of modes with `DEDUP_MODES` (e.g.
  `DEDUP_MODES="merged" ./scale_benchmark.sh 5000000`). Results go to
  `results/scale_benchmark_<timestamp>.tsv` (+ `.log`).

Generated BAMs live in `large_bam/` (gitignored, multi-GB) — delete freely and regenerate on
demand; generation is deterministic via `--seed`.

## Correctness

Fork output is **byte-identical** to the original (kept reads, `.log` stats, `.hist`) across:
- all 16 fixture-BAM test scenarios (sequential + parallel fork modes)
- every scale in the benchmark below, 100K → ~90M reads
- the big local pileup scenario (below), merged mode

## Scaling benchmark (2026-07-12, default mode)

Same 15-contig / 300 Mbp synthetic genome at every scale, so duplicate density grows naturally
with depth (mirroring deeper sequencing of one library). 8 logical CPUs, single run per config.
These numbers are **default (prefix-based) mode**; the benchmark now also runs **merged mode**
(`-m`) at every scale — the path the fork actually optimizes — recorded under the TSV's `mode`
column.

| reads | orig | fork -t1 | fork -t8 | seq × | par × | match |
|---|---|---|---|---|---|---|
| 100K | 1.5s | 1.5s | 1.7s | 1.01× | 0.90× | MATCH |
| 1M | 6.4s | 5.4s | 7.8s | 1.18× | 0.81× | MATCH |
| 5M | 19.6s | 13.2s | 13.7s | 1.48× | 1.43× | MATCH |
| 10M | 31.3s | 18.3s | 15.6s | 1.71× | 2.01× | MATCH |
| 25M | 65.7s | 30.5s | 19.0s | 2.16× | 3.45× | MATCH |
| 50M | 142.6s | 52.8s | 27.2s | 2.70× | 5.25× | MATCH |
| ~90M (100M target) | 273.4s | 89.4s | 35.4s | 3.06× | **7.72×** | MATCH |

Raw data: `results/scale_benchmark_20260712_140426.tsv` (100K/1M/10M/100M) and
`results/scale_benchmark_20260712_143458.tsv` (5M/25M/50M).

**Crossover:** `-t8` parallel mode trails `-t1` sequential below ~5M reads — thread-pool
startup overhead dominates a short run. It overtakes between 5M and 10M reads, and the gap
widens steadily through 90M, reaching 7.7× over the original at the top end.

**Caveat:** the generator's duplicate-density calibration undershoots at very high density — a
100,000,000-read request produced 90,395,866 reads (contigs ran out of `avg_advance` budget
before hitting target). Fine for benchmarking, not exact for anything needing a precise read
count. Shown as "~90M" above; timings are for the file actually run.

## Big local pileup benchmark (2026-07-15, merged mode)

One locus, 1000 distinct duplicate groups sharing a start position (40 duplicate reads each,
40,000 reads total) — the O(depth²) shape merged mode's bucket-map fast path targets. Both jars
run with `-m`; single run per config.

| reads | orig | fork -t1 | fork -t4 | seq × | par × | match |
|---|---|---|---|---|---|---|
| 40,000 | 3.6s | 0.6s | 0.8s | **6.11×** | 4.66× | MATCH |

Raw data: `results/report_20260715_233339.md` (section 7).

**Parallel trails sequential here:** the fork's parallelism splits work by chromosome, and this
scenario is a single contig, so `-t4` can't split the work — it only adds thread-pool overhead
on top of the same single-threaded bucket-map resolution. This is expected and mirrors the
small-file crossover seen in the scaling benchmark above.

## Reproducing

```bash
cd DeDup_comparison

# Scaling benchmark — both modes (default + merged) at every scale:
./scale_benchmark.sh 100000 1000000 5000000 10000000 25000000 50000000 100000000
# Just one mode:
DEDUP_MODES="merged" ./scale_benchmark.sh 100000 1000000 5000000 10000000 25000000 50000000 100000000

# Full comparison report, including the big local pileup section:
./run_comparison.sh
# Skip the 10M large-scale section and/or the big-pileup section:
./run_comparison.sh --no-large --no-pileup
```

Requires `DeDup-0.12.9.jar` (original) and `DeDup-0.13.0.jar` (fork) built
