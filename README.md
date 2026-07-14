# DeDup fork vs. original — comparison & benchmarks

Infrastructure for verifying that the fork (`DeDup_fork/`) produces identical output to the
original (`DeDup_original/`) and for measuring how much faster it is.

## Scripts

- `run_comparison.sh` — runs both jars against the fixture BAMs plus a generated large BAM,
  diffs kept-read-name sets, and writes a report to `results/report_<timestamp>.md`.
  Pass `--no-large` to skip the large-BAM section.
- `generate_large_bam.py` — fabricates a coordinate-sorted synthetic BAM directly via `pysam`
  (no real reference; DeDup only reads position/strand/length/MAPQ/name-prefix).
  `--reads` / `--contigs` / `--genome-length` control scale; duplicate-cluster density is
  self-calibrated per contig to hit the target read count (see caveat below).
- `scale_benchmark.sh [scale ...]` — runs orig / fork-seq / fork-par once each across a list of
  read-count scales, sharing one cached genome per scale. Default scales:
  `100000 1000000 5000000 10000000 25000000 50000000 100000000`. Results go to
  `results/scale_benchmark_<timestamp>.tsv` (+ `.log`).

Generated BAMs live in `large_bam/` (gitignored, multi-GB) — delete freely and regenerate on
demand; generation is deterministic via `--seed`.

## Correctness

Fork output is **byte-identical** to the original (kept reads, `.log` stats, `.hist`) across:
- all 16 fixture-BAM test scenarios (sequential + parallel fork modes)
- every scale in the benchmark below, 100K → ~90M reads

## Scaling benchmark (2026-07-12)

Same 15-contig / 300 Mbp synthetic genome at every scale, so duplicate density grows naturally
with depth (mirroring deeper sequencing of one library). 8 logical CPUs, single run per config.

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

## Reproducing

```bash
cd DeDup_comparison
./scale_benchmark.sh 100000 1000000 5000000 10000000 25000000 50000000 100000000
```

Requires `DeDup-0.12.9.jar` (original) and `DeDup-0.13.0.jar` (fork) built — see the top-level
`CLAUDE.md` for build commands.
