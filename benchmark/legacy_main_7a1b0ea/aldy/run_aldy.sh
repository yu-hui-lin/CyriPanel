#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J aldy_bench
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --time=24:00:00
#SBATCH --array=1-122%60
#SBATCH -o /work/u7715055/staging/biology/u7715055/aldy_bench/logs/task%a_%A.out
#SBATCH -e /work/u7715055/staging/biology/u7715055/aldy_bench/logs/task%a_%A.err
set -uo pipefail

module purge
module load biology
module load Samtools/1.15.1
module load python/3.12.2
source /home/u7715055/aldy_new_env/bin/activate

NEW=/work/u7715055/staging/biology/u7715055
WD=$NEW/aldy_bench
CN="chr22:42151472-42152258"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $WD/manifest_refs.txt)
ARM=$(echo "$LINE"  | cut -f1)
REF=$(echo "$LINE"  | cut -f2)
REFBAM=$(echo "$LINE" | cut -f3)

OUTDIR=$WD/$ARM/ref_$REF
mkdir -p "$OUTDIR"

# 該 arm 所有 test BAM（排除 reference 自己）
BAMDIR=$(dirname "$(dirname "$REFBAM")")
SUF=$(basename "$REFBAM" | sed "s/^${REF}//")

echo "=== task $SLURM_ARRAY_TASK_ID | $ARM | ref=$REF | $(date) ==="
echo "    bamdir=$BAMDIR  suffix=$SUF"

run_one () {
  local BAM="$1" S
  S=$(basename "$BAM" | sed "s/${SUF}\$//")
  [ "$S" = "$REF" ] && return 0
  [ -s "$OUTDIR/$S.cyp2d6.aldy" ] && return 0
  aldy genotype -p "$REFBAM" -g cyp2d6 -n "$CN" \
      -o "$OUTDIR/$S.cyp2d6.aldy" -l "$OUTDIR/$S.cyp2d6.log" \
      "$BAM" > /dev/null 2>&1 \
    || echo "    FAILED: $S"
}
export -f run_one
export REF REFBAM OUTDIR CN SUF

# Pass 1: 並行 3（峰值可達 ~26 GB/程序）
ls "$BAMDIR"/*/*"$SUF" | xargs -P 3 -I{} bash -c 'run_one "$@"' _ {}

# 清掉 OOM 產生的 0-byte 檔
find "$OUTDIR" -name "*.cyp2d6.aldy" -size 0 -delete

# Pass 2: 序列重跑失敗者（獨佔記憶體）
echo "--- retry pass ---"
for BAM in "$BAMDIR"/*/*"$SUF"; do
  S=$(basename "$BAM" | sed "s/${SUF}\$//")
  [ "$S" = "$REF" ] && continue
  [ -s "$OUTDIR/$S.cyp2d6.aldy" ] && continue
  echo "    retrying $S"
  aldy genotype -p "$REFBAM" -g cyp2d6 -n "$CN" \
      -o "$OUTDIR/$S.cyp2d6.aldy" -l "$OUTDIR/$S.cyp2d6.log" \
      "$BAM" > /dev/null 2>&1 || echo "    STILL FAILED: $S"
done
find "$OUTDIR" -name "*.cyp2d6.aldy" -size 0 -delete

N=$(ls "$OUTDIR"/*.cyp2d6.aldy 2>/dev/null | wc -l)
echo "=== done | $N results | $(date) ==="
