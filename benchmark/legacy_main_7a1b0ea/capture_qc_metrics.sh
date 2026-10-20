#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J cyri_qc
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o qc_out.log
#SBATCH -e qc_err.log
#SBATCH --mail-user=yhlin.md05@nycu.edu.tw
#SBATCH --mail-type=BEGIN,END,FAIL
#
# capture_qc_metrics.sh
# ---------------------
# Cohort-wide sequencing QC for Table 1: on-target fraction, duplicate rate,
# insert size and CYP2D-locus depth, for every BAM in all four benchmark arms.
#
# All counts are index-assisted and region-restricted, so this is fast: it never
# streams a whole BAM. Output: one TSV, one row per (arm, sample).
#
# Usage:  sbatch capture_qc_metrics.sh
#     or: bash capture_qc_metrics.sh          (login node; ~20-40 min)

set -u

ROOT=/work/u7715055/staging/biology/u7715055

# --- locate samtools -------------------------------------------------------
# Errors are shown, not swallowed: the first run of this script failed silently
# because a module name was wrong and the failure went to /dev/null.
if ! command -v module >/dev/null 2>&1; then
  for init in /usr/share/Modules/init/bash /etc/profile.d/modules.sh \
              /opt/ohpc/admin/lmod/lmod/init/bash /usr/share/lmod/lmod/init/bash; do
    [ -f "$init" ] && . "$init" && break
  done
fi
if command -v module >/dev/null 2>&1; then
  module load biology         || echo "note: 'module load biology' failed"
  module load Samtools/1.15.1 \
    || module load SAMtools/1.15.1 \
    || module load samtools    \
    || echo "note: no samtools module loaded"
else
  echo "note: no module command available; relying on PATH"
fi
# Fall back to the project virtualenvs, which is where the other scripts get it.
if ! command -v samtools >/dev/null 2>&1; then
  for env in "$ROOT/analysis_env" /home/u7715055/cyrius_env /home/u7715055/aldy_new_env; do
    [ -x "$env/bin/samtools" ] && PATH="$env/bin:$PATH" && export PATH && break
  done
fi
SAMTOOLS=${SAMTOOLS_BIN:-$(command -v samtools 2>/dev/null || true)}
if [ -z "$SAMTOOLS" ]; then
  {
    echo "ERROR: samtools not found - nothing was computed."
    echo "  PATH = $PATH"
    echo
    echo "  Find it with:"
    echo "     module load biology && module avail 2>&1 | grep -i samtool"
    echo "     ls /home/u7715055/*/bin/samtools $ROOT/*/bin/samtools 2>/dev/null"
    echo "  then re-run with the right module name, or pass the path directly:"
    echo "     SAMTOOLS_BIN=/full/path/to/samtools sbatch capture_qc_metrics.sh"
  } >&2
  exit 1
fi
echo "using samtools: $SAMTOOLS"
"$SAMTOOLS" --version | head -1
echo
BED=$ROOT/CyriPanel/data/PGxProbe_region_hg38.bed
LOCUS=chr22:42121996-42155994          # whole CYP2D locus
D6=chr22:42123192-42132032             # CYP2D6plusREP6_hapcn2
OUT=${QC_OUT:-$ROOT/CyriPanel_seq_qc.tsv}
echo "writing $OUT"

printf 'arm\tsample\ttotal_mapped\ton_target\tpct_on_target\tdup_on_target\tpct_dup_on_target\tinsert_median\tinsert_mean\tinsert_sd\tdepth_locus\tdepth_CYP2D6\tcov_pct_CYP2D6\n' > "$OUT"

run_arm () {
  arm=$1; dir=$2; suffix=$3; n=0
  for d in "$dir"/*/; do
    s=$(basename "$d")
    bam="$d$s$suffix"
    [ -f "$bam" ] || { bam=$(ls "$d"*.bam 2>/dev/null | head -1); }
    [ -f "$bam" ] || { echo "  skip $arm/$s (no bam)" >&2; continue; }
    [ -f "$bam.bai" ] || [ -f "${bam%.bam}.bai" ] || { echo "  skip $arm/$s (no index)" >&2; continue; }

    # total primary mapped reads, from the index (instant)
    tot=$("$SAMTOOLS" idxstats "$bam" | awk '{s+=$3} END{print s+0}')
    # on-target primary reads
    ont=$("$SAMTOOLS" view -c -F 0x904 -L "$BED" "$bam")
    # duplicates among on-target reads
    dup=$("$SAMTOOLS" view -c -f 0x400 -F 0x904 -L "$BED" "$bam")
    # insert size over the CYP2D locus, from samtools stats SN lines
    read -r imed imean isd <<EOF
$("$SAMTOOLS" stats -@ 1 "$bam" "$LOCUS" 2>/dev/null | awk -F'\t' '
  /^SN\tinsert size average:/  {mean=$3}
  /^SN\tinsert size standard deviation:/ {sd=$3}
  /^IS\t/ {c[$2]=$3; n+=$3}
  END{ t=0; for(i=0;i<=2000;i++){ if(i in c){ t+=c[i]; if(t>=n/2 && med==""){med=i} } }
       printf "%s %s %s", (med==""?"NA":med), (mean==""?"NA":mean), (sd==""?"NA":sd) }')
EOF
    # mean depth over the locus and over the CYP2D6 interval
    read -r dl <<EOF
$("$SAMTOOLS" coverage -r "$LOCUS" "$bam" 2>/dev/null | awk 'NR==2{print $7}')
EOF
    read -r d6d d6c <<EOF
$("$SAMTOOLS" coverage -r "$D6" "$bam" 2>/dev/null | awk 'NR==2{print $7, $6}')
EOF
    pct=$(awk -v a="$ont" -v b="$tot" 'BEGIN{printf "%.2f", (b>0? 100*a/b : 0)}')
    pdup=$(awk -v a="$dup" -v b="$ont" 'BEGIN{printf "%.2f", (b>0? 100*a/b : 0)}')
    n=$((n+1)); [ $((n % 10)) -eq 0 ] && echo "    ... $arm $n samples" >&2
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$arm" "$s" "$tot" "$ont" "$pct" "$dup" "$pdup" \
      "${imed:-NA}" "${imean:-NA}" "${isd:-NA}" "${dl:-NA}" "${d6d:-NA}" "${d6c:-NA}" >> "$OUT"
  done
}

echo "cohort1  BWA-MEM"; run_arm cohort1        "$ROOT/bwa_results/HPRC_PGpanel_altaware" ".markdup.bam"
echo "cohort2  BWA-MEM"; run_arm cohort2        "$ROOT/bwa_results/PG_20250711_altaware"  ".markdup.bam"
echo "cohort1  DRAGEN";  run_arm cohort1_dragen "$ROOT/dragen_markdup/cohort1"            ".bam"
echo "cohort2  DRAGEN";  run_arm cohort2_dragen "$ROOT/dragen_markdup/cohort2"            ".bam"

echo
echo "Wrote $OUT"
echo "Now run:  python3 inspect_seq_qc.py"
