#!/bin/sh
#SBATCH -A MST109178
#SBATCH -J verify549
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o verify549_out.log
#SBATCH -e verify549_err.log
#SBATCH --mail-user=yhlin.md05@nycu.edu.tw
#SBATCH --mail-type=BEGIN,END,FAIL
#
# verify_fix549.sh
# ================
# Confirms that fixing star_caller.py:549 (get_snp_position argument order)
# does not change any existing genotype call.
#
# WHY ONLY 10 SAMPLES
#   The fix makes call_var42127803hap functional. Its only consumer is the
#   *119/*2 vs *1/*41 tie-break in match_star_allele. The risk is that a sample
#   currently called *1/*41 flips to *119/*2 (wrong: no *119 exists in either
#   cohort). Every sample whose gold genotype contains *41, *119, *32 or *27 is
#   therefore tested. Nothing else can reach that branch.
#
# WHAT IT DOES
#   Runs the 10 samples twice - once with the current (fixed) star_caller.py and
#   once with the pre-fix backup - into a scratch directory, then diffs the TSVs.
#   Nothing is written into the benchmark results tree.
#
# USAGE
#   sbatch verify_fix549.sh
#   (or run interactively: sh verify_fix549.sh)

set -eu

SRC=/work/u7715055/staging/biology/u7715055/CyriPanel
BACKUP=$SRC/star_caller.py.bak_fix549_20260902
SCRATCH=/work/u7715055/staging/biology/u7715055/verify_fix549
BWA=/work/u7715055/staging/biology/u7715055/bwa_results
VENV=/home/u7715055/cyrius_env/bin/activate

C1_BAM=$BWA/HPRC_PGpanel_altaware
C2_BAM=$BWA/PG_20250711_altaware

# Samples whose gold genotype can reach the affected tie-breaks
C1_SAMPLES="HG00733 HG01071 HG01109 HG01175 HG03492"
C2_SAMPLES="HG00280 HG01167 HG01192 HG02559 HG04184"

echo "=== verify_fix549: start $(date) ==="

if [ ! -f "$BACKUP" ]; then
    echo "FATAL: pre-fix backup not found at $BACKUP" >&2
    exit 1
fi

# Sanity: the two versions must actually differ, and only at line 549
echo "--- difference between fixed and pre-fix star_caller.py ---"
diff "$BACKUP" "$SRC/star_caller.py" || true
echo "-----------------------------------------------------------"

rm -rf "$SCRATCH"
mkdir -p "$SCRATCH"/work "$SCRATCH"/out

# A working copy so the live tree is never modified
cp -r "$SRC"/caller "$SRC"/depth_calling "$SRC"/data "$SCRATCH"/work/
cp "$SRC"/star_caller.py "$SCRATCH"/work/star_caller_fixed.py
cp "$BACKUP"            "$SCRATCH"/work/star_caller_prefix.py
mkdir -p "$SCRATCH"/work/ref_dir
rm -rf "$SCRATCH"/work/depth_calling/cnv_panelizer_results
mkdir -p "$SCRATCH"/work/depth_calling/cnv_panelizer_results

# shellcheck disable=SC1090
. "$VENV"

run_cohort () {
    cohort=$1; bamdir=$2; samples=$3; dipcsv=$4

    # Reference panel: first 5 known diploids of this cohort.
    # The choice is arbitrary - it only has to be IDENTICAL across the two runs.
    rm -f "$SCRATCH"/work/ref_dir/*.bam "$SCRATCH"/work/ref_dir/*.bai
    refs=$(tail -n +2 "$dipcsv" | head -5 | tr -d '\r')
    echo "  reference panel ($cohort): $(echo "$refs" | tr '\n' ' ')"
    for r in $refs; do
        ln -sf "$bamdir/$r/$r.markdup.bam"     "$SCRATCH/work/ref_dir/$r.bam"
        ln -sf "$bamdir/$r/$r.markdup.bam.bai" "$SCRATCH/work/ref_dir/$r.bam.bai"
    done

    man=$SCRATCH/${cohort}_manifest.txt
    : > "$man"
    for s in $samples; do
        echo "$bamdir/$s/$s.markdup.bam" >> "$man"
    done

    for variant in fixed prefix; do
        echo "--- $cohort / $variant ---"
        cp "$SCRATCH/work/star_caller_${variant}.py" "$SCRATCH/work/star_caller.py"
        rm -f "$SCRATCH"/work/depth_calling/cnv_panelizer_results/*.csv
        find "$SCRATCH/work" -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        ( cd "$SCRATCH/work" && python3 star_caller.py \
              -m "$man" --genome 38 \
              -o "$SCRATCH/out" -p "${cohort}_${variant}" -t 14 )
    done
}

run_cohort cohort1 "$C1_BAM" "$C1_SAMPLES" "$SRC/data/CYP2D6_cohort1_33diploids.csv"
run_cohort cohort2 "$C2_BAM" "$C2_SAMPLES" "$SRC/data/CYP2D6_cohort2_28diploids.csv"

echo
echo "==================== RESULT ===================="
status=0
for cohort in cohort1 cohort2; do
    a=$SCRATCH/out/${cohort}_prefix.tsv
    b=$SCRATCH/out/${cohort}_fixed.tsv
    if [ ! -f "$a" ] || [ ! -f "$b" ]; then
        echo "$cohort: MISSING OUTPUT - check the logs"; status=1; continue
    fi
    if diff -q "$a" "$b" >/dev/null; then
        echo "$cohort: IDENTICAL - the fix changed no call"
    else
        echo "$cohort: *** CALLS CHANGED *** - inspect before trusting the results"
        diff "$a" "$b" || true
        status=1
    fi
done
echo "==============================================="
echo "Outputs: $SCRATCH/out"
echo "=== verify_fix549: end $(date) ==="
exit $status
