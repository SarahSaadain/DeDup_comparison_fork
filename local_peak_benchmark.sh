#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# DeDup local-peak scaling benchmark — original vs fork at a single high-depth
# locus (the O(depth^2)-shaped bug both fixes target — not genome-wide
# throughput, see scale_benchmark.sh for that).
#
# Runs TWO scenarios, each swept across a --scale list of total read counts:
#   merged  — generate_big_pileup_bam.py: one locus, all reads share the exact
#             same start, distinct lengths give `groups` distinct (start,end)
#             duplicate groups of READS_PER_GROUP reads each. Run with -m.
#             Targets the merged-mode bucket-map fast path (see
#             dedup-merged-oN2-fix memory / checkForDuplicationMerged).
#   default — generate_staggered_pileup_bam.py: POSITIONS distinct start
#             positions and LENGTH held fixed (LENGTH > POSITIONS so the whole
#             band stays buffered at once); each scale's --reads controls
#             reads-per-position, i.e. depth. Read length is deliberately
#             decoupled from --reads here (see that generator's docstring) —
#             scaling length together with read count would make total BAM
#             I/O volume O(N^2) and mask the dedup algorithm's own cost. Run
#             without -m. Targets checkForDuplication()'s non-merged path
#             (RecordBufferHeap + startIndex/endIndex candidate narrowing).
#
# Usage:
#   ./local_peak_benchmark.sh [scale ...]
#   # default scales: 20000 50000 100000 200000 400000 800000
#   # scenarios default to "merged default"; override via env, e.g.
#   PEAK_SCENARIOS="default" ./local_peak_benchmark.sh 100000 400000
#   # jar paths default to sibling checkouts; override for portability:
#   ./local_peak_benchmark.sh --orig-jar PATH --fork-jar PATH 100000
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

READS_PER_GROUP=40   # merged scenario: fixed duplicate reads per (start,end) group; scale via --groups
STAGGER_POSITIONS=300 # default scenario: fixed distinct start positions (generate_staggered_pileup_bam.py)
STAGGER_LENGTH=400     # default scenario: fixed read length (must exceed STAGGER_POSITIONS)

mkdir -p "$LARGE_BAM_DIR" "$RESULTS_DIR"

# Pull --orig-jar/--fork-jar out of the argument list before what's left is
# treated as the scale list, so both styles can be combined freely.
ORIG_JAR_ARG=""; FORK_JAR_ARG=""; POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --orig-jar) ORIG_JAR_ARG="$2"; shift 2 ;;
    --orig-jar=*) ORIG_JAR_ARG="${1#*=}"; shift ;;
    --fork-jar) FORK_JAR_ARG="$2"; shift 2 ;;
    --fork-jar=*) FORK_JAR_ARG="${1#*=}"; shift ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

SCALES=("$@")
[[ ${#SCALES[@]} -eq 0 ]] && SCALES=(20000 50000 100000 200000 400000 800000)

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

if [[ -n "$ORIG_JAR_ARG" ]]; then
  ORIG_JAR="$ORIG_JAR_ARG"
  [[ -f "$ORIG_JAR" ]] || { log "ERROR: --orig-jar not found: $ORIG_JAR"; exit 1; }
else
  ORIG_JAR=$(ls "$ORIG_DIR"/build/libs/DeDup-*.jar 2>/dev/null | grep -v sources | head -1 || true)
fi
if [[ -n "$FORK_JAR_ARG" ]]; then
  FORK_JAR="$FORK_JAR_ARG"
  [[ -f "$FORK_JAR" ]] || { log "ERROR: --fork-jar not found: $FORK_JAR"; exit 1; }
else
  FORK_JAR=$(ls "$FORK_DIR"/build/libs/DeDup-*.jar 2>/dev/null | grep -v sources | head -1 || true)
fi
[[ -z "$ORIG_JAR" || -z "$FORK_JAR" ]] && { log "ERROR: JARs not found — build them first (gradle jar in both dirs) or pass --orig-jar/--fork-jar"; exit 1; }
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

fmt_name() {
  local scale="$1"
  if   (( scale >= 1000000 )); then printf '%sM' "$((scale / 1000000))"
  elif (( scale >= 1000 ));    then printf '%sK' "$((scale / 1000))"
  else printf '%s' "$scale"
  fi
}

printf 'scale\tscenario\tbam\treads\torig_ms\tfork_seq_ms\tfork_par_ms\tseq_speedup\tpar_speedup\torig_kept\tfork_kept\tmatch\n' > "$OUT_TSV"

for scale in "${SCALES[@]}"; do
  name=$(fmt_name "$scale")
  for scenario in "${SCENARIOS[@]}"; do
    log "── scale=$scale ($name) scenario=$scenario ──"

    if [[ "$scenario" == "merged" ]]; then
      groups=$(( scale / READS_PER_GROUP ))
      bam="$LARGE_BAM_DIR/peak_merged_${name}.bam"
      if [[ ! -f "$bam" ]]; then
        log "  generating $bam (groups=$groups, reads-per-group=$READS_PER_GROUP) ..."
        python3 "$SCRIPT_DIR/generate_big_pileup_bam.py" -o "$bam" \
          --groups "$groups" --reads-per-group "$READS_PER_GROUP" --seed 42 2>>"$OUT_LOG"
      else
        log "  reusing cached $bam"
      fi
      mflag=(-m)
    else
      bam="$LARGE_BAM_DIR/peak_staggered_${name}.bam"
      if [[ ! -f "$bam" ]]; then
        log "  generating $bam (positions=$STAGGER_POSITIONS, length=$STAGGER_LENGTH, reads=$scale) ..."
        python3 "$SCRIPT_DIR/generate_staggered_pileup_bam.py" -o "$bam" \
          --reads "$scale" --positions "$STAGGER_POSITIONS" --length "$STAGGER_LENGTH" --seed 7 2>>"$OUT_LOG"
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
      "$name" "$scenario" "$(basename "$bam")" "$reads" "$orig_ms" "$fseq_ms" "$fpar_ms" \
      "$seq_speedup" "$par_speedup" "$orig_kept" "$fork_kept" "$match" >> "$OUT_TSV"

    log "  [$scenario] scale=$name reads=$reads orig=${orig_ms}ms fork-seq=${fseq_ms}ms(${seq_speedup}) fork-par=${fpar_ms}ms(${par_speedup}) match=$match"

    rm -rf "$od" "$fsd" "$fpd"
  done
done

log "─────────────────────────────────────────────────────"
log "Results: $OUT_TSV"
column -t -s $'\t' "$OUT_TSV" | tee -a "$OUT_LOG"
