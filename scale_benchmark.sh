#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# DeDup scaling benchmark — original vs fork across read-count scales.
#
# Generates (or reuses cached) synthetic BAMs at each requested scale, all
# sharing the same 15-contig / 300Mbp synthetic genome, so duplication rate
# grows naturally with depth (mirrors sequencing a library deeper). For each
# scale AND each dedup mode (default prefix-based / merged -m), runs original,
# fork -t1, fork -t<ncpu> once each and records wall time + correctness
# (kept-read-name-set match).
#
# Merged mode (-m) is the path the fork actually optimizes (O(depth^2) →
# O(group) bucket resolution), so its speedups are the interesting ones;
# default mode is kept for coverage / regression.
#
# Usage:
#   ./scale_benchmark.sh [--orig-jar PATH] [--fork-jar PATH] [scale ...]
#   # default scales: 100000 1000000 5000000 10000000 25000000 50000000 100000000
#   # modes default to "default merged"; override via env, e.g.
#   DEDUP_MODES="merged" ./scale_benchmark.sh 1000000 10000000
#
#   By default the jars are auto-discovered as sibling checkouts
#   (../DeDup_original, ../DeDup_fork relative to this script) built via
#   `gradle jar`. Pass --orig-jar/--fork-jar to point at jars anywhere else —
#   e.g. on a machine that only has the built jars, not the source checkouts:
#   ./scale_benchmark.sh --orig-jar ~/DeDup-0.12.9.jar --fork-jar ~/DeDup-0.13.0.jar 1000000
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
POC_DIR="$(dirname "$SCRIPT_DIR")"
ORIG_DIR="$POC_DIR/DeDup_original"
FORK_DIR="$POC_DIR/DeDup_fork"
LARGE_BAM_DIR="$SCRIPT_DIR/large_bam"
RESULTS_DIR="$SCRIPT_DIR/results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUT_TSV="$RESULTS_DIR/scale_benchmark_$TIMESTAMP.tsv"
OUT_LOG="$RESULTS_DIR/scale_benchmark_$TIMESTAMP.log"
GENOME_LENGTH=300000000

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
[[ ${#SCALES[@]} -eq 0 ]] && SCALES=(100000 1000000 5000000 10000000 25000000 50000000 100000000)

# Dedup modes to run per scale. "default" = prefix-based (M_/F_/R_); "merged" = -m
# (both ends for every read, the fork's optimized path). Override via DEDUP_MODES.
read -r -a MODES <<< "${DEDUP_MODES:-default merged}"

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
  if   (( scale >= 1000000000 )); then printf '%sB' "$((scale / 1000000000))"
  elif (( scale >= 1000000 ));    then printf '%sM' "$((scale / 1000000))"
  elif (( scale >= 1000 ));       then printf '%sK' "$((scale / 1000))"
  else printf '%s' "$scale"
  fi
}

printf 'scale\tmode\tbam\treads\torig_ms\tfork_seq_ms\tfork_par_ms\tseq_speedup\tpar_speedup\torig_kept\tfork_kept\tmatch\n' > "$OUT_TSV"

for scale in "${SCALES[@]}"; do
  name=$(fmt_name "$scale")
  bam="$LARGE_BAM_DIR/scale_${name}_15contigs.bam"
  log "── scale=$scale ($name reads) ──"

  if [[ ! -f "$bam" ]]; then
    log "  generating $bam ..."
    python3 "$SCRIPT_DIR/generate_large_bam.py" -o "$bam" --reads "$scale" \
      --contigs 15 --genome-length "$GENOME_LENGTH" --seed 42 2>>"$OUT_LOG"
  else
    log "  reusing cached $bam"
  fi
  [[ -f "${bam}.bai" ]] || { log "  indexing..."; samtools index "$bam"; }

  reads=$(bam_read_count "$bam")

  for mode in "${MODES[@]}"; do
    # Merged mode adds -m to every invocation; default mode adds nothing. The
    # ${mflag[@]+...} guard expands a possibly-empty array safely under set -u
    # (portable back to bash 3.2 on macOS).
    mflag=()
    [[ "$mode" == "merged" ]] && mflag=(-m)
    log "  ── mode=$mode ──"

    od=$(mktemp -d); fsd=$(mktemp -d); fpd=$(mktemp -d)

    log "  running original ($mode)..."
    orig_ms=$(time_ms "$JAVA" -jar "$ORIG_JAR" -i "$bam" -o "$od" "${mflag[@]+"${mflag[@]}"}")
    log "  running fork -t1 ($mode)..."
    fseq_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$bam" -o "$fsd" -t 1 "${mflag[@]+"${mflag[@]}"}")
    log "  running fork -t$N_CPUS ($mode)..."
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
      "$name" "$mode" "$(basename "$bam")" "$reads" "$orig_ms" "$fseq_ms" "$fpar_ms" \
      "$seq_speedup" "$par_speedup" "$orig_kept" "$fork_kept" "$match" >> "$OUT_TSV"

    log "  [$mode] orig=${orig_ms}ms fork-seq=${fseq_ms}ms(${seq_speedup}) fork-par=${fpar_ms}ms(${par_speedup}) match=$match"

    rm -rf "$od" "$fsd" "$fpd"
  done
done

log "─────────────────────────────────────────────────────"
log "Results: $OUT_TSV"
column -t -s $'\t' "$OUT_TSV" | tee -a "$OUT_LOG"
