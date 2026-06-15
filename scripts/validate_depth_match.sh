#!/usr/bin/env bash
#
# validate_depth_match.sh — Spot-check `impg depth` output rows against the
# `impg query -x -m 0` gold standard.
#
# Plan reference: PLAN_depth_100pct.md §2.2 — large-scale acceptance harness.
#
# What it does:
#   1. Streams a depth.tsv file (4 cols: id, length, depth, positions).
#   2. Keeps rows with length > MIN_LEN (default 200), reservoir-samples N of them.
#   3. For each sampled row, takes positions[0] as the anchor and runs
#        impg query --alignment-list ALIST -t THREADS -r ANCHOR -x -m 0
#      then compares the unique sample prefixes to the row's reported set.
#   4. Emits a TSV with per-row diagnostics and a one-line summary.
#
# Per-row TSV columns (written to OUT_TSV):
#   id len reported_depth reported_samples anchor q_count q_samples \
#   subset_ok over_extra under_missing direction match
#
#   subset_ok       Y if reported ⊆ query (depth claim is consistent with truth)
#   direction       match | over | under | both | query_failed
#   match           Y iff sample sets are equal AND depth == |sample set|
#
# Usage:
#   validate_depth_match.sh <depth.tsv> <alist> <bin> <N> \
#       [--datadir DIR] [--min-len LEN] [--threads T] \
#       [--index INDEX_FILE] [--out OUT_TSV] [--seed SEED]
#
# Required positional args:
#   depth.tsv   header `#id\tlength\tdepth\tpositions` then data rows
#   alist       passed verbatim to impg --alignment-list
#   bin         path to the impg binary
#   N           number of rows to sample (must be > 0)
#
# Common flags:
#   --datadir DIR     cwd to invoke `bin` in (so relative paths in `alist` resolve)
#   --min-len LEN     length threshold for sampling (default: 200)
#   --threads T       --threads passed to impg query (default: 4)
#   --index FILE      passed verbatim to `impg query -i` (cached IMPG index, optional)
#   --out FILE        path to write per-row TSV (default: <depth.tsv>.match.tsv)
#   --seed SEED       seed for awk reservoir sampling (default: epoch + PID)
#
# Exit codes:
#   0  match_pct == 100%
#   1  match_pct  <  threshold (default 99% or via --threshold)
#   2  CLI / IO error
#
# Acceptance gates per plan:
#   ≥ 95%   fix PR allowed to merge
#   ≥ 99%   issue can be closed
#   ≥ 100%  release-tag eligibility
#

set -euo pipefail

usage() {
    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
}

# ---- Defaults ----
MIN_LEN=200
THREADS=4
INDEX=""
OUT_TSV=""
SEED="$(($(date +%s) ^ $$))"
THRESHOLD_PCT=99
DATADIR="."

# ---- Parse positional + named args ----
POS=()
while (($# > 0)); do
    case "$1" in
        --datadir) DATADIR="$2"; shift 2 ;;
        --min-len) MIN_LEN="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --index)   INDEX="$2"; shift 2 ;;
        --out)     OUT_TSV="$2"; shift 2 ;;
        --seed)    SEED="$2"; shift 2 ;;
        --threshold) THRESHOLD_PCT="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; POS+=("$@"); break ;;
        -*) echo "unknown flag: $1" >&2; usage ;;
        *)  POS+=("$1"); shift ;;
    esac
done
((${#POS[@]} == 4)) || { echo "expected 4 positional args, got ${#POS[@]}" >&2; usage; }

DEPTH_TSV="${POS[0]}"
ALIST="${POS[1]}"
BIN="${POS[2]}"
N="${POS[3]}"

[[ -r "$DEPTH_TSV" ]] || { echo "depth.tsv unreadable: $DEPTH_TSV" >&2; exit 2; }
[[ -x "$BIN" ]]       || { echo "impg binary not executable: $BIN" >&2; exit 2; }
[[ "$N" =~ ^[0-9]+$ ]] && (( N > 0 )) || { echo "N must be a positive integer, got: $N" >&2; exit 2; }

if [[ -z "$OUT_TSV" ]]; then
    OUT_TSV="${DEPTH_TSV}.match.tsv"
fi

# ---- Reservoir-sample N rows with length > MIN_LEN ----
# Streams depth.tsv (skips header), keeps length-filtered rows,
# Vitter's algorithm R via awk for unbiased reservoir sampling.
SAMPLE_TSV=$(mktemp)
trap 'rm -f "$SAMPLE_TSV"' EXIT

awk -F'\t' -v n="$N" -v min_len="$MIN_LEN" -v seed="$SEED" '
    BEGIN { srand(seed); count = 0 }
    /^#/ { next }
    NF < 4 { next }
    ($2 + 0) <= min_len { next }
    {
        count++
        if (count <= n) {
            buf[count] = $0
        } else {
            r = int(rand() * count) + 1
            if (r <= n) buf[r] = $0
        }
    }
    END {
        kept = (count < n) ? count : n
        for (i = 1; i <= kept; i++) print buf[i]
    }
' "$DEPTH_TSV" > "$SAMPLE_TSV"

SAMPLED=$(wc -l < "$SAMPLE_TSV")
if (( SAMPLED == 0 )); then
    echo "no rows with length > $MIN_LEN found in $DEPTH_TSV" >&2
    exit 2
fi

# ---- Header for the per-row diagnostics TSV ----
{
    printf 'id\tlen\treported_depth\treported_samples\tanchor\tq_count\tq_samples\tsubset_ok\tover_extra\tunder_missing\tdirection\tmatch\n'
} > "$OUT_TSV"

# ---- Per-row probe ----
ok_match=0
ok_subset=0
over=0
under=0
both=0
query_failed=0

while IFS=$'\t' read -r id length depth positions; do
    [[ -z "$id" ]] && continue
    anchor="${positions%%;*}"

    reported_samples=$(echo "$positions" | tr ';' '\n' | awk -F'#' 'NF{print $1}' | sort -u | tr '\n' ',' | sed 's/,$//')

    declare -a query_args=(
        query
        --alignment-list "$ALIST"
        -t "$THREADS"
        -r "$anchor"
        -x -m 0
        -o bed
    )
    if [[ -n "$INDEX" ]]; then
        query_args+=( -i "$INDEX" )
    fi

    if q_output=$(cd "$DATADIR" && "$BIN" "${query_args[@]}" 2>/dev/null); then
        :
    else
        # Capture the failure but keep going
        q_output=""
        query_failed=$((query_failed + 1))
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$id" "$length" "$depth" "$reported_samples" "$anchor" \
            "" "" "N" "" "" "query_failed" "N" >> "$OUT_TSV"
        continue
    fi

    q_samples_list=$(printf '%s\n' "$q_output" | awk -F'\t' 'NF>=3{split($1,a,"#"); print a[1]}' | sort -u | tr '\n' ',' | sed 's/,$//')
    q_count=$(printf '%s\n' "$q_output" | awk -F'\t' 'NF>=3{split($1,a,"#"); print a[1]}' | sort -u | grep -c .)

    # Set diff: reported \ query  (over-reports)
    over_extra=$(awk -v rep="$reported_samples" -v qry=",${q_samples_list}," '
        BEGIN {
            n = split(rep, a, ",")
            out = ""; sep = ""
            for (i = 1; i <= n; i++) {
                if (a[i] == "") continue
                if (index(qry, "," a[i] ",") == 0) { out = out sep a[i]; sep = "," }
            }
            print out
        }
    ')

    # Set diff: query \ reported (under-reports)
    under_missing=$(awk -v qry="$q_samples_list" -v rep=",${reported_samples}," '
        BEGIN {
            n = split(qry, a, ",")
            out = ""; sep = ""
            for (i = 1; i <= n; i++) {
                if (a[i] == "") continue
                if (index(rep, "," a[i] ",") == 0) { out = out sep a[i]; sep = "," }
            }
            print out
        }
    ')

    if [[ -z "$over_extra" ]]; then subset_ok="Y"; else subset_ok="N"; fi

    if [[ -z "$over_extra" && -z "$under_missing" && "$q_count" == "$depth" ]]; then
        match="Y"; direction="match"
        ok_match=$((ok_match + 1))
    else
        match="N"
        if [[ -n "$over_extra" && -n "$under_missing" ]]; then
            direction="both"; both=$((both + 1))
        elif [[ -n "$over_extra" ]]; then
            direction="over"; over=$((over + 1))
        elif [[ -n "$under_missing" ]]; then
            direction="under"; under=$((under + 1))
        else
            # Sample sets agree but depth count differs.
            direction="depth_count_mismatch"
        fi
    fi

    [[ "$subset_ok" == "Y" ]] && ok_subset=$((ok_subset + 1))

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$id" "$length" "$depth" "$reported_samples" "$anchor" \
        "$q_count" "$q_samples_list" "$subset_ok" "$over_extra" \
        "$under_missing" "$direction" "$match" >> "$OUT_TSV"
done < "$SAMPLE_TSV"

# ---- Summary ----
total=$(($(wc -l < "$OUT_TSV") - 1))
match_pct=$(awk -v m="$ok_match" -v t="$total" 'BEGIN{ if (t == 0) print "0.0"; else printf "%.2f", 100.0 * m / t }')
subset_pct=$(awk -v m="$ok_subset" -v t="$total" 'BEGIN{ if (t == 0) print "0.0"; else printf "%.2f", 100.0 * m / t }')

cat <<EOF >&2

==================== validate_depth_match.sh summary ====================
depth.tsv       : $DEPTH_TSV
alist           : $ALIST
binary          : $BIN
sampled rows    : $total  (min length > $MIN_LEN, requested N=$N, seed=$SEED)
out             : $OUT_TSV

match           : $ok_match  ($match_pct%)
subset_ok       : $ok_subset  ($subset_pct%)
over-report     : $over     (reported ⊋ query — depth fabricated samples)
under-report    : $under    (reported ⊊ query — depth missed samples)
both directions : $both
query_failed    : $query_failed

acceptance gate : ≥ ${THRESHOLD_PCT}% match required to exit 0
=========================================================================
EOF

awk -v p="$match_pct" -v g="$THRESHOLD_PCT" 'BEGIN{ exit (p+0 < g+0) ? 1 : 0 }'
