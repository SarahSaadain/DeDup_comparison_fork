#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# DeDup scaling benchmark — original vs fork across read-count scales.
#
# Generates (or reuses cached) synthetic BAMs at each requested scale, all
# sharing the same 15-contig / 300Mbp synthetic genome, so duplication rate
# grows naturally with depth (mirrors sequencing a library deeper). Runs
# original, fork -t1, fork -t<ncpu> once each per scale and records wall time
# + correctness (kept-read-name-set match).
#
# Usage:
#   ./scale_benchmark.sh [scale ...]
#   # default: 100000 1000000 5000000 10000000 25000000 50000000 100000000
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

SCALES=("$@")
[[ ${#SCALES[@]} -eq 0 ]] && SCALES=(100000 1000000 5000000 10000000 25000000 50000000 100000000)

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

N_CPUS=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

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

printf 'scale\tbam\treads\torig_ms\tfork_seq_ms\tfork_par_ms\tseq_speedup\tpar_speedup\torig_kept\tfork_kept\tmatch\n' > "$OUT_TSV"

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

  od=$(mktemp -d); fsd=$(mktemp -d); fpd=$(mktemp -d)

  log "  running original..."
  orig_ms=$(time_ms "$JAVA" -jar "$ORIG_JAR" -i "$bam" -o "$od")
  log "  running fork -t1..."
  fseq_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$bam" -o "$fsd" -t 1)
  log "  running fork -t$N_CPUS..."
  fpar_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$bam" -o "$fpd" -t "$N_CPUS")

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

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$name" "$(basename "$bam")" "$reads" "$orig_ms" "$fseq_ms" "$fpar_ms" \
    "$seq_speedup" "$par_speedup" "$orig_kept" "$fork_kept" "$match" >> "$OUT_TSV"

  log "  orig=${orig_ms}ms fork-seq=${fseq_ms}ms(${seq_speedup}) fork-par=${fpar_ms}ms(${par_speedup}) match=$match"

  rm -rf "$od" "$fsd" "$fpd"
done

log "─────────────────────────────────────────────────────"
log "Results: $OUT_TSV"
column -t -s $'\t' "$OUT_TSV" | tee -a "$OUT_LOG"
