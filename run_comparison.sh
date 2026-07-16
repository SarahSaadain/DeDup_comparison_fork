#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# DeDup Comparison — original vs fork
#
# Tracks:
#   1. Test-suite health (gradle test, both versions)
#   2. Dedup correctness (same BAM in → same reads kept out?)
#   3. Wall-clock time  (3-run median, sequential vs parallel)
#   4. Peak memory      (/usr/bin/time -l, macOS RSS)
#   5. Large-scale benchmark (10M reads / 15 contigs, generated on demand)
#   6. Big local pileup benchmark, merged mode (one locus, thousands of
#      duplicate groups sharing a start — the O(depth^2) shape the fork's
#      merged-mode bucket-map fix targets; generated on demand)
#
# Usage:
#   ./run_comparison.sh [--no-build] [--no-tests] [--no-large] [--no-pileup] \
#                        [--orig-jar PATH] [--fork-jar PATH]
#
#   By default the jars are built from sibling source checkouts (../DeDup_original,
#   ../DeDup_fork relative to this script) via `gradle jar`. Pass --orig-jar/--fork-jar
#   to point at prebuilt jars anywhere else instead — e.g. on a machine that only has the
#   built jars, not the source checkouts:
#   ./run_comparison.sh --orig-jar ~/DeDup-0.12.9.jar --fork-jar ~/DeDup-0.13.0.jar
#   Passing either one implies --no-build --no-tests (there's no source checkout to build
#   or run gradle test against), and the fixture-BAM correctness section (§2) will show no
#   rows since those BAMs live under DeDup_fork's checkout.
#
# Outputs (in results/):
#   report_TIMESTAMP.md   — human-readable summary
#   raw_TIMESTAMP.tsv     — machine-readable per-BAM data
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
POC_DIR="$(dirname "$SCRIPT_DIR")"
ORIG_DIR="$POC_DIR/DeDup_original"
FORK_DIR="$POC_DIR/DeDup_fork"
BAM_DIR="$FORK_DIR/src/test/resources"
RESULTS_DIR="$SCRIPT_DIR/results"
LARGE_BAM_DIR="$SCRIPT_DIR/large_bam"
LARGE_BAM="$LARGE_BAM_DIR/large_10M_15contigs.bam"
PILEUP_BAM="$LARGE_BAM_DIR/merged_big_pileup.bam"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
REPORT="$RESULTS_DIR/report_$TIMESTAMP.md"
RAW_TSV="$RESULTS_DIR/raw_$TIMESTAMP.tsv"
TMPDIR_CMP=$(mktemp -d)
trap 'rm -rf "$TMPDIR_CMP"' EXIT

ORIG_JAR_ARG=""; FORK_JAR_ARG=""
NO_BUILD=0; NO_TESTS=0; NO_LARGE=0; NO_PILEUP=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-build) NO_BUILD=1; shift ;;
    --no-tests) NO_TESTS=1; shift ;;
    --no-large) NO_LARGE=1; shift ;;
    --no-pileup) NO_PILEUP=1; shift ;;
    --orig-jar) ORIG_JAR_ARG="$2"; shift 2 ;;
    --orig-jar=*) ORIG_JAR_ARG="${1#*=}"; shift ;;
    --fork-jar) FORK_JAR_ARG="$2"; shift 2 ;;
    --fork-jar=*) FORK_JAR_ARG="${1#*=}"; shift ;;
    *) shift ;;
  esac
done

# A jar-only machine has no DeDup_original/DeDup_fork source checkout for gradle
# to build or test against, so passing either jar explicitly implies skipping both.
if [[ -n "$ORIG_JAR_ARG" || -n "$FORK_JAR_ARG" ]]; then
  NO_BUILD=1
  NO_TESTS=1
fi

mkdir -p "$RESULTS_DIR"

log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }

# Detect Java 11+ (the JARs require at least Java 11; Homebrew OpenJDK preferred)
JAVA=""
for candidate in \
    /usr/local/opt/java/bin/java \
    /opt/homebrew/opt/java/bin/java \
    java; do
  if command -v "$candidate" >/dev/null 2>&1; then
    ver=$("$candidate" -version 2>&1 | awk -F'"' '/version/{print $2}' | cut -d. -f1)
    if [[ "$ver" -ge 11 ]] 2>/dev/null; then
      JAVA="$candidate"; break
    fi
  fi
done
[[ -z "$JAVA" ]] && { log "ERROR: No Java 11+ found. Install via: brew install java"; exit 1; }
log "Using Java: $JAVA ($($JAVA -version 2>&1 | head -1))"

# ── helpers ────────────────────────────────────────────────────────────────────

# Return median of 3 integers
median3() {
  printf '%s\n' "$1" "$2" "$3" | sort -n | sed -n '2p'
}

# Measure wall-clock time of a command in milliseconds
time_ms() {
  local start end
  start=$(python3 -c "import time; print(int(time.time()*1000))")
  "$@" >/dev/null 2>&1
  end=$(python3 -c "import time; print(int(time.time()*1000))")
  echo $((end - start))
}

# Peak RSS in KB from /usr/bin/time -l (macOS)
peak_rss_kb() {
  local rss
  rss=$({ /usr/bin/time -l "$@" >/dev/null; } 2>&1 \
        | grep 'maximum resident' | awk '{print $1}')
  echo $(( ${rss:-0} / 1024 ))
}

# Extract numeric XML attribute value (BSD grep-safe)
xml_attr() {
  local attr="$1" file="$2"
  sed -n "s/.*${attr}=\"\([0-9]*\)\".*/\1/p" "$file" | head -1
}

# Sum attribute across all TEST-*.xml in a gradle project test-results dir
sum_xml_attr() {
  local attr="$1" xml_dir="$2"
  local total=0 f
  for f in "$xml_dir"/TEST-*.xml; do
    [[ -f "$f" ]] || continue
    v=$(xml_attr "$attr" "$f"); v=${v:-0}
    total=$((total + v))
  done
  echo "$total"
}

# Collect test totals for a gradle project
collect_tests() {
  local dir="$1"
  local xdir="$dir/build/test-results/test"
  local tot fa er sk
  tot=$(sum_xml_attr tests    "$xdir")
  fa=$(sum_xml_attr failures  "$xdir")
  er=$(sum_xml_attr errors    "$xdir")
  sk=$(sum_xml_attr skipped   "$xdir")
  printf '%s %s %s %s' "$tot" "$fa" "$er" "$sk"
}

# Flags for a given BAM basename
bam_flags() {
  case "$1" in
    all_reads_as_merged_test.bam) echo "-m" ;;
    *) echo "" ;;
  esac
}

# Count mapped reads in a BAM
bam_read_count() { samtools view -c -F 4 "$1"; }

# Sorted read-name list
bam_read_names() { samtools view "$1" | cut -f1 | sort; }

# ── build ─────────────────────────────────────────────────────────────────────
if [[ $NO_BUILD -eq 0 ]]; then
  log "Building original..."
  gradle -p "$ORIG_DIR" jar -q
  log "Building fork..."
  gradle -p "$FORK_DIR" jar -q
fi

if [[ -n "$ORIG_JAR_ARG" ]]; then
  ORIG_JAR="$ORIG_JAR_ARG"
  [[ -f "$ORIG_JAR" ]] || { log "ERROR: --orig-jar not found: $ORIG_JAR"; exit 1; }
else
  ORIG_JAR=$(ls "$ORIG_DIR"/build/libs/DeDup-*.jar 2>/dev/null | head -1 || true)
fi
if [[ -n "$FORK_JAR_ARG" ]]; then
  FORK_JAR="$FORK_JAR_ARG"
  [[ -f "$FORK_JAR" ]] || { log "ERROR: --fork-jar not found: $FORK_JAR"; exit 1; }
else
  FORK_JAR=$(ls "$FORK_DIR"/build/libs/DeDup-*.jar 2>/dev/null | head -1 || true)
fi
[[ -z "$ORIG_JAR" ]] && { log "ERROR: original JAR missing (run without --no-build, or pass --orig-jar)"; exit 1; }
[[ -z "$FORK_JAR" ]] && { log "ERROR: fork JAR missing (run without --no-build, or pass --fork-jar)"; exit 1; }
log "Original JAR: $(basename "$ORIG_JAR")"
log "Fork JAR:     $(basename "$FORK_JAR")"

# ── test suites ────────────────────────────────────────────────────────────────
if [[ $NO_TESTS -eq 0 ]]; then
  log "Running original test suite..."
  gradle -p "$ORIG_DIR" test --rerun-tasks -q 2>&1 || true
  log "Running fork test suite..."
  gradle -p "$FORK_DIR" test --rerun-tasks -q 2>&1 || true
else
  log "--no-tests: using cached XML results"
fi

ORIG_TESTS=$(collect_tests "$ORIG_DIR")
FORK_TESTS=$(collect_tests "$FORK_DIR")

# ── BAM benchmark ─────────────────────────────────────────────────────────────
BAMS=(
  forward_test.bam
  reverse_test.bam
  forward_unmerged_duplicate_must_overlap.bam
  reverse_unmerged_duplicate_must_overlap.bam
  all_reads_as_merged_test.bam
  mapq_zero_test.bam
  stack_test.bam
  stack_read_one_test.bam
  stack_read_two_test.bam
  strand_forward.bam
  strand_reverse.bam
  yield_test.bam
)

N_CPUS=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
RUNS=3

log "Benchmarking $RUNS runs per BAM per config (CPUs=$N_CPUS)..."
printf 'bam\tflags\torig_ms\tfork_seq_ms\tfork_par_ms\torig_kept\tfork_kept\tmatch\n' > "$RAW_TSV"

for bam in "${BAMS[@]}"; do
  INPUT="$BAM_DIR/$bam"
  [[ -f "$INPUT" ]] || { log "  SKIP $bam (not found)"; continue; }
  flags=$(bam_flags "$bam")
  log "  $bam  flags='$flags'"

  OD_ORIG=$(mktemp -d "$TMPDIR_CMP/o_XXXX")
  OD_FSEQ=$(mktemp -d "$TMPDIR_CMP/fs_XXXX")
  OD_FPAR=$(mktemp -d "$TMPDIR_CMP/fp_XXXX")

  # Build flag lists
  if [[ -n "$flags" ]]; then
    O_ARGS=(-jar "$ORIG_JAR" -i "$INPUT" -o "$OD_ORIG" "$flags")
    FS_ARGS=(-jar "$FORK_JAR" -i "$INPUT" -o "$OD_FSEQ" -t 1 "$flags")
    FP_ARGS=(-jar "$FORK_JAR" -i "$INPUT" -o "$OD_FPAR" -t "$N_CPUS" "$flags")
  else
    O_ARGS=(-jar "$ORIG_JAR" -i "$INPUT" -o "$OD_ORIG")
    FS_ARGS=(-jar "$FORK_JAR" -i "$INPUT" -o "$OD_FSEQ" -t 1)
    FP_ARGS=(-jar "$FORK_JAR" -i "$INPUT" -o "$OD_FPAR" -t "$N_CPUS")
  fi

  # 3 timing runs each, take median
  ot1=$(time_ms "$JAVA" "${O_ARGS[@]}");  ot2=$(time_ms "$JAVA" "${O_ARGS[@]}");  ot3=$(time_ms "$JAVA" "${O_ARGS[@]}")
  ft1=$(time_ms "$JAVA" "${FS_ARGS[@]}"); ft2=$(time_ms "$JAVA" "${FS_ARGS[@]}"); ft3=$(time_ms "$JAVA" "${FS_ARGS[@]}")
  pt1=$(time_ms "$JAVA" "${FP_ARGS[@]}"); pt2=$(time_ms "$JAVA" "${FP_ARGS[@]}"); pt3=$(time_ms "$JAVA" "${FP_ARGS[@]}")

  orig_ms=$(median3 "$ot1" "$ot2" "$ot3")
  fork_seq_ms=$(median3 "$ft1" "$ft2" "$ft3")
  fork_par_ms=$(median3 "$pt1" "$pt2" "$pt3")

  # Compare last-run outputs
  base="${bam%.bam}"
  orig_bam=$(ls "$OD_ORIG/${base}"*_rmdup.bam 2>/dev/null | head -1 || true)
  fork_bam=$(ls "$OD_FSEQ/${base}"*_rmdup.bam 2>/dev/null | head -1 || true)

  orig_kept="?"; fork_kept="?"; match="?"
  if [[ -f "$orig_bam" && -f "$fork_bam" ]]; then
    orig_kept=$(bam_read_count "$orig_bam")
    fork_kept=$(bam_read_count "$fork_bam")
    o_names=$(bam_read_names "$orig_bam")
    f_names=$(bam_read_names "$fork_bam")
    if [[ "$o_names" == "$f_names" ]]; then
      match="MATCH"
    else
      match="DIFF"
      printf '%s\n' "$o_names" > "$RESULTS_DIR/${base}_orig_names.txt"
      printf '%s\n' "$f_names" > "$RESULTS_DIR/${base}_fork_names.txt"
    fi
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$bam" "$flags" "$orig_ms" "$fork_seq_ms" "$fork_par_ms" \
    "$orig_kept" "$fork_kept" "$match" >> "$RAW_TSV"
done

# ── memory (yield_test.bam) ───────────────────────────────────────────────────
log "Measuring peak memory on yield_test.bam..."
MEM_IN="$BAM_DIR/yield_test.bam"
MO=$(mktemp -d "$TMPDIR_CMP/mo_XXXX")
MF=$(mktemp -d "$TMPDIR_CMP/mf_XXXX")
MP=$(mktemp -d "$TMPDIR_CMP/mp_XXXX")

orig_rss=$(peak_rss_kb    "$JAVA" -jar "$ORIG_JAR" -i "$MEM_IN" -o "$MO")
fork_rss=$(peak_rss_kb    "$JAVA" -jar "$FORK_JAR" -i "$MEM_IN" -o "$MF" -t 1)
fork_par_rss=$(peak_rss_kb "$JAVA" -jar "$FORK_JAR" -i "$MEM_IN" -o "$MP" -t "$N_CPUS")

# ── large-scale benchmark (10M reads / 15 contigs) ─────────────────────────────
large_orig_ms="—"; large_fseq_ms="—"; large_fpar_ms="—"
large_reads="—"; large_orig_kept="—"; large_fork_kept="—"; large_match="—"
if [[ $NO_LARGE -eq 0 ]]; then
  if [[ ! -f "$LARGE_BAM" ]]; then
    log "Large BAM not found, generating (10M reads / 15 contigs, ~90s)..."
    mkdir -p "$LARGE_BAM_DIR"
    python3 "$SCRIPT_DIR/generate_large_bam.py" -o "$LARGE_BAM" --reads 10000000 --contigs 15
  fi
  [[ -f "${LARGE_BAM}.bai" ]] || samtools index "$LARGE_BAM"

  log "Benchmarking large BAM ($(basename "$LARGE_BAM"), single run per config)..."
  LO=$(mktemp -d "$TMPDIR_CMP/lo_XXXX")
  LFS=$(mktemp -d "$TMPDIR_CMP/lfs_XXXX")
  LFP=$(mktemp -d "$TMPDIR_CMP/lfp_XXXX")

  large_orig_ms=$(time_ms "$JAVA" -jar "$ORIG_JAR" -i "$LARGE_BAM" -o "$LO")
  large_fseq_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$LARGE_BAM" -o "$LFS" -t 1)
  large_fpar_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$LARGE_BAM" -o "$LFP" -t "$N_CPUS")

  large_orig_bam=$(ls "$LO"/*_rmdup.bam 2>/dev/null | head -1 || true)
  large_fork_bam=$(ls "$LFS"/*_rmdup.bam 2>/dev/null | head -1 || true)
  if [[ -f "$large_orig_bam" && -f "$large_fork_bam" ]]; then
    large_reads=$(bam_read_count "$LARGE_BAM")
    large_orig_kept=$(bam_read_count "$large_orig_bam")
    large_fork_kept=$(bam_read_count "$large_fork_bam")
    if [[ "$(bam_read_names "$large_orig_bam")" == "$(bam_read_names "$large_fork_bam")" ]]; then
      large_match="MATCH"
    else
      large_match="DIFF"
    fi
  fi
else
  log "--no-large: skipping large-scale benchmark"
fi

# ── big local pileup benchmark (merged mode) ────────────────────────────────────
# One locus, many distinct (start,end) duplicate groups sharing a start position,
# each with several duplicate reads. Since the resolve trigger never fires while
# same-start reads stream in, the whole pileup buffers at once and drains
# group-by-group at end of stream — the exact shape that made the old merged-mode
# duplicate resolution O(depth^2); the fork's bucket-map fast path resolves it in
# O(depth). Mirrors DeDup_fork's MergedHighDepthLocusTest big-pileup case.
pileup_reads="—"; pileup_orig_ms="—"; pileup_fseq_ms="—"; pileup_fpar_ms="—"
pileup_orig_kept="—"; pileup_fork_kept="—"; pileup_match="—"
if [[ $NO_PILEUP -eq 0 ]]; then
  if [[ ! -f "$PILEUP_BAM" ]]; then
    log "Big pileup BAM not found, generating (1000 groups x 40 reads/group, one locus)..."
    mkdir -p "$LARGE_BAM_DIR"
    python3 "$SCRIPT_DIR/generate_big_pileup_bam.py" -o "$PILEUP_BAM" \
      --groups 1000 --reads-per-group 40 --seed 42
  fi
  [[ -f "${PILEUP_BAM}.bai" ]] || samtools index "$PILEUP_BAM"

  log "Benchmarking big local pileup ($(basename "$PILEUP_BAM"), merged mode, single run per config)..."
  PO=$(mktemp -d "$TMPDIR_CMP/po_XXXX")
  PFS=$(mktemp -d "$TMPDIR_CMP/pfs_XXXX")
  PFP=$(mktemp -d "$TMPDIR_CMP/pfp_XXXX")

  pileup_orig_ms=$(time_ms "$JAVA" -jar "$ORIG_JAR" -i "$PILEUP_BAM" -o "$PO" -m)
  pileup_fseq_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$PILEUP_BAM" -o "$PFS" -t 1 -m)
  pileup_fpar_ms=$(time_ms "$JAVA" -jar "$FORK_JAR" -i "$PILEUP_BAM" -o "$PFP" -t "$N_CPUS" -m)

  pileup_orig_bam=$(ls "$PO"/*_rmdup.bam 2>/dev/null | head -1 || true)
  pileup_fork_bam=$(ls "$PFS"/*_rmdup.bam 2>/dev/null | head -1 || true)
  if [[ -f "$pileup_orig_bam" && -f "$pileup_fork_bam" ]]; then
    pileup_reads=$(bam_read_count "$PILEUP_BAM")
    pileup_orig_kept=$(bam_read_count "$pileup_orig_bam")
    pileup_fork_kept=$(bam_read_count "$pileup_fork_bam")
    if [[ "$(bam_read_names "$pileup_orig_bam")" == "$(bam_read_names "$pileup_fork_bam")" ]]; then
      pileup_match="MATCH"
    else
      pileup_match="DIFF"
      bam_read_names "$pileup_orig_bam" > "$RESULTS_DIR/big_pileup_orig_names.txt"
      bam_read_names "$pileup_fork_bam" > "$RESULTS_DIR/big_pileup_fork_names.txt"
    fi
  fi
else
  log "--no-pileup: skipping big local pileup benchmark"
fi

# ── report (Python for table formatting) ──────────────────────────────────────
python3 - "$RAW_TSV" "$REPORT" \
  "$ORIG_TESTS" "$FORK_TESTS" \
  "$orig_rss" "$fork_rss" "$fork_par_rss" \
  "$N_CPUS" \
  "$(basename "$ORIG_JAR")" "$(basename "$FORK_JAR")" \
  "$(date '+%Y-%m-%d %H:%M:%S')" \
  "$large_reads" "$large_orig_ms" "$large_fseq_ms" "$large_fpar_ms" \
  "$large_orig_kept" "$large_fork_kept" "$large_match" \
  "$pileup_reads" "$pileup_orig_ms" "$pileup_fseq_ms" "$pileup_fpar_ms" \
  "$pileup_orig_kept" "$pileup_fork_kept" "$pileup_match" \
<<'PYEOF'
import sys, csv
from pathlib import Path

tsv_path   = Path(sys.argv[1])
out_path   = Path(sys.argv[2])
orig_tests = sys.argv[3].split()
fork_tests = sys.argv[4].split()
orig_rss   = int(sys.argv[5])
fork_rss   = int(sys.argv[6])
fpar_rss   = int(sys.argv[7])
n_cpus     = sys.argv[8]
orig_jar   = sys.argv[9]
fork_jar   = sys.argv[10]
ts         = sys.argv[11]
large_reads      = sys.argv[12]
large_orig_ms    = sys.argv[13]
large_fseq_ms    = sys.argv[14]
large_fpar_ms    = sys.argv[15]
large_orig_kept  = sys.argv[16]
large_fork_kept  = sys.argv[17]
large_match      = sys.argv[18]
pileup_reads      = sys.argv[19]
pileup_orig_ms    = sys.argv[20]
pileup_fseq_ms    = sys.argv[21]
pileup_fpar_ms    = sys.argv[22]
pileup_orig_kept  = sys.argv[23]
pileup_fork_kept  = sys.argv[24]
pileup_match      = sys.argv[25]

def parse_tests(t):
    tot, fa, er, sk = map(int, t)
    return dict(total=tot, passed=tot-fa-er-sk, failures=fa, errors=er, skipped=sk)

ot = parse_tests(orig_tests)
ft = parse_tests(fork_tests)

rows = []
with open(tsv_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        rows.append(r)

def speedup(orig, fork):
    orig, fork = int(orig), int(fork)
    if fork == 0: return "—"
    return f"{orig/fork:.2f}x"

lines = [
    f"# DeDup Performance Comparison Report",
    f"**Generated:** {ts}  ",
    f"**Host CPUs:** {n_cpus} logical cores  ",
    f"**Timing:** median of 3 runs per configuration  ",
    f"**Original JAR:** {orig_jar}  ",
    f"**Fork JAR:** {fork_jar}  ",
    "",
    "---",
    "",
    "## 1 · Test Suite Health",
    "",
    "| Version  | Total | Passed | Failures | Errors | Skipped |",
    "|----------|------:|-------:|---------:|-------:|--------:|",
    f"| Original | {ot['total']} | {ot['passed']} | {ot['failures']} | {ot['errors']} | {ot['skipped']} |",
    f"| Fork     | {ft['total']} | {ft['passed']} | {ft['failures']} | {ft['errors']} | {ft['skipped']} |",
    "",
    "> Fork adds **ParallelDedupTest** (7 tests verifying sequential == parallel output).",
    "",
    "---",
    "",
    "## 2 · Dedup Correctness (original vs fork-sequential, identical BAM input)",
    "",
    "| BAM | Orig kept | Fork kept | Read names match? |",
    "|-----|----------:|----------:|:-----------------:|",
]
for r in rows:
    lines.append(f"| {r['bam']} | {r['orig_kept']} | {r['fork_kept']} | {r['match']} |")

lines += [
    "",
    "> **DIFF** rows have `results/<bam>_{orig,fork}_names.txt` for inspection.",
    "",
    "---",
    "",
    "## 3 · Wall-Clock Time (ms, median of 3 runs)",
    "",
    "`seq_speedup` / `par_speedup` = original_ms ÷ fork_ms  (> 1 means fork is faster).",
    "",
    "| BAM | Original | Fork-seq | Fork-par | seq_speedup | par_speedup |",
    "|-----|---------:|---------:|---------:|------------:|------------:|",
]
for r in rows:
    lines.append(
        f"| {r['bam']} | {r['orig_ms']} | {r['fork_seq_ms']} | {r['fork_par_ms']}"
        f" | {speedup(r['orig_ms'], r['fork_seq_ms'])}"
        f" | {speedup(r['orig_ms'], r['fork_par_ms'])} |"
    )

lines += [
    "",
    "---",
    "",
    "## 4 · Peak Memory — yield_test.bam (macOS peak RSS)",
    "",
    "| Version                    | Peak RSS (KB) |",
    "|----------------------------|-------------:|",
    f"| Original (sequential)      | {orig_rss:,} |",
    f"| Fork sequential (−t 1)     | {fork_rss:,} |",
    f"| Fork parallel (−t {n_cpus})       | {fpar_rss:,} |",
    "",
    "---",
    "",
    "## 5 · Architectural Differences",
    "",
    "| Aspect | Original | Fork |",
    "|--------|----------|------|",
    "| Dedup algorithm | Two PriorityQueues + discardSet | recordBuffer is RecordBufferHeap (algorithmically identical to PriorityQueue), over cached BufferedRead |",
    "| Condition matching | 21-entry EnumSet<DL> lookup table | Same 21-entry EnumSet<DL> lookup table (ported) |",
    "| Cross-type dedup | ✓ M_ absorbs F_/R_ at same position | ✓ Identical output (logic table ported) |",
    "| Non-merged duplicate scan | Full recordBuffer scan per anchor: O(depth) | Narrowed to reads sharing the anchor's exact start/end via startIndex/endIndex: O(candidates) |",
    "| Merged-mode duplicate scan | Full recordBuffer scan per anchor: O(depth) | Exact (start,end) bucket-map lookup: O(group) |",
    "| Quality score | Recomputed per comparison | Cached once in BufferedRead on creation |",
    "| Read-type detection | String.startsWith() each comparison | Parsed once, stored as ReadType enum |",
    "| Parallelism | Single-threaded only | Per-chromosome thread pool (−t N) |",
    "| Thread default | N/A | 1 (opt in via −t N) |",
    "",
    "---",
    "",
    "## 6 · Large-Scale Benchmark (~10M reads / 15 contigs, single run)",
    "",
    "Synthetic BAM generated by `generate_large_bam.py` (not committed — regenerate "
    "or reuse the cached copy in `large_bam/`). Calibrated for a realistic "
    "ancient-DNA-style duplication rate (~78% observed).",
    "",
    f"Total reads: **{large_reads}**",
    "",
    "| Version                    | Wall-clock (ms) | speedup vs original |",
    "|-----------------------------|-----------------:|---------------------:|",
    f"| Original                   | {large_orig_ms} | — |",
    f"| Fork sequential (−t 1)     | {large_fseq_ms} | {speedup(large_orig_ms, large_fseq_ms) if large_fseq_ms != '—' else '—'} |",
    f"| Fork parallel (−t {n_cpus})       | {large_fpar_ms} | {speedup(large_orig_ms, large_fpar_ms) if large_fpar_ms != '—' else '—'} |",
    "",
    "| Orig kept | Fork kept | Read names match? |",
    "|----------:|----------:|:------------------:|",
    f"| {large_orig_kept} | {large_fork_kept} | {large_match} |",
    "",
    "---",
    "",
    "## 7 · Big Local Pileup Benchmark (merged mode, single run)",
    "",
    "One locus, 1000 distinct duplicate groups sharing a start position (40 duplicate "
    "reads each) — synthetic BAM from `generate_big_pileup_bam.py` (not committed, "
    "cached in `large_bam/`). This is the O(depth^2) shape the fork's merged-mode "
    "bucket-map fast path targets: same-start reads never trigger the resolve check "
    "while streaming in, so the whole pileup buffers at once and only drains "
    "group-by-group at end of stream, where the old code rescanned the full buffer "
    "per group. Both jars run with `-m`.",
    "",
    f"Total reads: **{pileup_reads}**",
    "",
    "| Version                    | Wall-clock (ms) | speedup vs original |",
    "|-----------------------------|-----------------:|---------------------:|",
    f"| Original                   | {pileup_orig_ms} | — |",
    f"| Fork sequential (−t 1)     | {pileup_fseq_ms} | {speedup(pileup_orig_ms, pileup_fseq_ms) if pileup_fseq_ms != '—' else '—'} |",
    f"| Fork parallel (−t {n_cpus})       | {pileup_fpar_ms} | {speedup(pileup_orig_ms, pileup_fpar_ms) if pileup_fpar_ms != '—' else '—'} |",
    "",
    "| Orig kept | Fork kept | Read names match? |",
    "|----------:|----------:|:------------------:|",
    f"| {pileup_orig_kept} | {pileup_fork_kept} | {pileup_match} |",
    "",
]

out_path.write_text('\n'.join(lines))
print('\n'.join(lines))
PYEOF

log "─────────────────────────────────────────────────────"
log "Report: $REPORT"
log "Raw TSV: $RAW_TSV"
