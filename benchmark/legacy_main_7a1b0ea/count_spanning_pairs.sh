#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J span_pairs
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -t 01:00:00
#SBATCH -o span_pairs_out.log
#SBATCH -e span_pairs_err.log
#SBATCH --mail-user=yhlin.md05@nycu.edu.tw
#SBATCH --mail-type=BEGIN,END,FAIL
set -u
module load biology >/dev/null 2>&1; module load Samtools/1.15.1 >/dev/null 2>&1
command -v samtools >/dev/null || { echo "samtools not found: run  module load biology; module load Samtools/1.15.1"; exit 1; }
B=${SPAN_BASE:-/work/u7715055/staging/biology/u7715055}
ARMS="BWA-MEM cohort1 $B/bwa_results/HPRC_PGpanel_altaware .markdup.bam
BWA-MEM cohort2 $B/bwa_results/PG_20250711_altaware .markdup.bam
DRAGEN cohort1 $B/dragen_markdup/cohort1 .bam
DRAGEN cohort2 $B/dragen_markdup/cohort2 .bam"
span () {
  comm -12 <(samtools view -F 0xF04 "$1" "chr22:$2-$2" | cut -f1 | sort -u) \
           <(samtools view -F 0xF04 "$1" "chr22:$3-$3" | cut -f1 | sort -u) | wc -l | tr -d ' '
}
OUT=spanning_pairs_327bp.tsv
printf 'aligner\tcohort\tsample\tpairs\n' > "$OUT"
echo "$ARMS" | while read -r aln coh dir suf; do
  for d in "$dir"/*/; do
    s=$(basename "$d"); bam="$d$s$suf"; [ -f "$bam" ] || continue
    printf '%s\t%s\t%s\t%s\n' "$aln" "$coh" "$s" "$(span "$bam" 42126611 42126938)" >> "$OUT"
  done
done
echo "== 327-bp site pair in every sample: median [IQR], range"
python3 - "$OUT" <<'PY'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
def q(v, p):
    v = sorted(v); h = (len(v) - 1) * p; lo = int(h); hi = min(lo + 1, len(v) - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])
for aln in ('BWA-MEM', 'DRAGEN'):
    for coh in ('cohort1', 'cohort2'):
        v = [int(r['pairs']) for r in rows if r['aligner'] == aln and r['cohort'] == coh]
        if v:
            print('  %-8s %s  n=%-3d median %g [IQR %g-%g]  range %d-%d' % (aln, coh, len(v), st.median(v), q(v, .25), q(v, .75), min(v), max(v)))
    print('  %-8s HG02559 = %s' % (aln, ', '.join(r['pairs'] for r in rows if r['aligner'] == aln and r['sample'] == 'HG02559')))
PY
echo "== Site pairs of the multiple-candidate samples (Table S5a: 4, 10, 0, 0, 0)"
echo "$ARMS" | while read -r aln coh dir suf; do
  for spec in "HG01361 42128945 42129130 185" "HG01952 42128945 42129130 185" \
              "HG01361 42128945 42130692 1747" "HG01952 42128945 42130692 1747" \
              "HG02630 42129075 42130715 1640"; do
    set -- $spec; bam="$dir/$1/$1$suf"; [ -f "$bam" ] || continue
    printf '  %-8s %-8s %5s bp  %s\n' "$aln" "$1" "$4" "$(span "$bam" "$2" "$3")"
  done
done
echo "Per-sample counts: $OUT"
