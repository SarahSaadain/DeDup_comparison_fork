#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# DeDup local-peak scaling benchmark — original vs fork across duplicate-group
# counts at a single high-depth locus (the O(depth^2) shape, not genome-wide
# throughput — see scale_benchmark.sh for that).
#
# Runs TWO scenarios at each requested group count:
#   merged  — generate_big_pileup_bam.py, one locus, all reads share the exact
#             same start, distinct lengths give `groups` distinct (start,end)
#             duplicate groups of `--reads-per-group` reads each. Run with -m.
#             This is the fork's merged-mode bucket-map fast path (see
#             dedup-merged-oN2-fix memory).
#   default — generate_local_peak_bam.py, one locus, `groups` distinct 1bp-apart
#             F_ start positions with `--reads-per-position` duplicate reads
#             each, read-len sized so the first read still overlaps the last
#             position (nothing resolves until end of stream). Run without -m.
#             This is the fork's non-merged startIndex/endIndex narrowing,
#             added alongside RecordBufferHeap.
#
# Both scenarios keep reads-per-group/-position FIXED and scale the number of
# distinct groups, since that's the dimension the old O(depth^2) buffer
# rescan is quadratic in (see each generator's docstring).
#
# Usage:
#   ./local_peak_benchmark.sh [group-count ...]
#   # default group counts: 500 1000 2000 4000 8000
#   # scenarios default to "merged default"; override via env, e.g.
#   PEAK_SCENARIOS="default" ./local_peak_benchmark.sh 2000 4000
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
POC_DIR="$(dirname "$SCRIPT_DIR")"
ORIG_DIR="$POC_DIR/DeDup_original"
FORK_DIR="$POC_DIR/DeDup_fork"
LARGE_BAM_DIR="$SCRIPT_DIR/large_bam"
RESULTS_DIR="$SCRIPT_DIR/results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_TSV="$RESULTS_DIR/local_peak_benchmark_$TIMESTAMP.tsv"
OUT_LOG="$RESULTS_DIR/local_peak_benchmark_$TIMESTAMP.log"

READS_PER_GROUP=40     # merged scenario: fixed duplicate reads per (start,end) group
READS_PER_POSITION=10  # default scenario: fixed duplicate reads per start position

mkdir -p "$LARGE_BAM_DIR" "$RESULTS_DIR"

GROUP_COUNTS=("$@")
[[ ${#GROUP_COUNTS[@]} -eq 0 ]] && GROUP_COUNTS=(500 1000 2000 4000 8000)

read -r -a SCENARIOS <<< "${PEAK_SCENARIOS:-merged default}"

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$OUT_LOG" >&2; }

JAVA=""
for candidate in /usr/local/opt/java/bin/java /opt/homebrew/opt/java/bin/java java; do
  if command -v "$candidate" >/dev/null 2>&1; then
    ver=$("$candidate" -version 2>&1 | awk -F'"' '/version/{print $2}' | cut -d. -f1)
    if [[ "$ver" -ge 11 ]] 2>/dev/null; then JAVA="$candidate"; break; fi
  fi
done
[[ -z "$JAVA" ]] && { log "ERROR: No Java 11+ found"; exit 1; }

ORIG_JAR=$(ls "$ORIG_DIR"/build/libs/DeDup-*.jar 2>/dev/null | grep -v sources | head -1 || true)
FORK_JAR=$(ls "$FORK_DIR"/build/libs/DeDup-*.jar 2>/dev/null | grep -v sources | head -1 || true)
[[ -z "$ORIG_JAR" || -z "$FORK_JAR" ]] && { log "ERROR: build JARs first (gradle jar in both dirs)"; exit 1; }
log "Original JAR: $(basename "$ORIG_JAR")"
log "Fork JAR:     $(basename "$FORK_JAR")"

N_CPUS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)

time_ms() {
  local start end
  start=$(python3 -c "import time; print(int(time.time()*1000))")
  "$@" >/dev/null 2>&1
  end=$(python3 -c "import time; print(int(time.time()*1000))")
  echo $((end - start))
}

bam_read_count() { samtools view -c -F 4 "$1"; }
bam_read_names() { samtools view "$1" | cut -f1 | sort; }

printf 'groups\tscenario\tbam\treads\torig_ms\tfork_seq_ms\tfork_par_ms\tseq_speedup\tpar_speedup\torig_kept\tfork_kept\tmatch\n' > "$OUT_TSV"

for groups in "${GROUP_COUNTS[@]}"; do
  for scenario in "${SCENARIOS[@]}"; do
    log "── groups=$groups scenario=$scenario ──"

    if [[ "$scenario" == "merged" ]]; then
      bam="$LARGE_BAM_DIR/peak_merged_${groups}g${READS_PER_GROUP}.bam"
      if [[ ! -f "$bam" ]]; then
        log "  generating $bam ..."
        python3 "$SCRIPT_DIR/generate_big_pileup_bam.py" -o "$bam" \
          --groups "$groups" --reads-per-group "$READS_PER_GROUP" --seed 42 2>>"$OUT_LOG"
      else
        log "  reusing cached $bam"
      fi
      mflag=(-m)
    else
      bam="$LARGE_BAM_DIR/peak_local_${groups}p${READS_PER_POSITION}.bam"
      if [[ ! -f "$bam" ]]; then
        log "  generating $bam ..."
        python3 "$SCRIPT_DIR/generate_local_peak_bam.py" -o "$bam" \
          --positions "$groups" --reads-per-position "$READS_PER_POSITION" --seed 42 2>>"$OUT_LOG"
      else
        log "  reusing cached $bam"
      fi
      mflag=()
    fi
    [[ -f "${bam}.bai" ]] || { log "  indexing..."; samtools index "$bam"; }

    reads=$(bam_read_count "$bam")

    od=$(mktemp -d); fsd=$(mktemp -d); fpd=$(mktemp -d)

    log "  running original ($scenario)..."
    orig_ms=$(time_ms "$JAVA" -jar "$ORIG_JAR" -i "$bam" -o "$od" "${mflag[@]+"${mflag[@]}"}")
    log "  running fork -t1 ($scenario)..."
    fseq_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$bam" -o "$fsd" -t 1 "${mflag[@]+"${mflag[@]}"}")
    log "  running fork -t$N_CPUS ($scenario)..."
    fpar_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$bam" -o "$fpd" -t "$N_CPUS" "${mflag[@]+"${mflag[@]}"}")

    orig_bam=$(ls "$od"/*_rmdup.bam 2>/dev/null | head -1 || true)
    fork_bam=$(ls "$fsd"/*_rmdup.bam 2>/dev/null | head -1 || true)
    orig_kept="?"; fork_kept="?"; match="?"
    if [[ -f "$orig_bam" && -f "$fork_bam" ]]; then
      orig_kept=$(bam_read_count "$orig_bam")
      fork_kept=$(bam_read_count "$fork_bam")
      if [[ "$(bam_read_names "$orig_bam")" == "$(bam_read_names "$fork_bam")" ]]; then
        match="MATCH"
      else
        match="DIFF"
      fi
    fi

    seq_speedup=$(python3 -c "print(f'{$orig_ms/$fseq_ms:.2f}x')" 2>/dev/null || echo "—")
    par_speedup=$(python3 -c "print(f'{$orig_ms/$fpar_ms:.2f}x')" 2>/dev/null || echo "—")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$groups" "$scenario" "$(basename "$bam")" "$reads" "$orig_ms" "$fseq_ms" "$fpar_ms" \
      "$seq_speedup" "$par_speedup" "$orig_kept" "$fork_kept" "$match" >> "$OUT_TSV"

    log "  [$scenario] groups=$groups reads=$reads orig=${orig_ms}ms fork-seq=${fseq_ms}ms(${seq_speedup}) fork-par=${fpar_ms}ms(${par_speedup}) match=$match"

    rm -rf "$od" "$fsd" "$fpd"
  done
done

log "─────────────────────────────────────────────────────"
log "Results: $OUT_TSV"
column -t -s $'\t' "$OUT_TSV" | tee -a "$OUT_LOG"
