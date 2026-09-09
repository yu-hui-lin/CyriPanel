#!/bin/bash
#SBATCH --job-name=cyri_cohort2_dragen_n20
#SBATCH --partition=ngs92G
#SBATCH --account=MST109178
#SBATCH --cpus-per-task=14
#SBATCH --mem=92g
#SBATCH --time=24:13:00
#SBATCH --array=0-29%10
#SBATCH --output=/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_logs/cohort2_dragen_n20_iter%a_%A.out
#SBATCH --error=/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_logs/cohort2_dragen_n20_iter%a_%A.err

set -euo pipefail
source /home/u7715055/cyrius_env/bin/activate

ITER=$SLURM_ARRAY_TASK_ID
SEED=$(( 42 + ITER ))
ITER_PAD=$(printf '%03d' "$ITER")

OUT_DIR="/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/cohort2_dragen/panel20/iter_${ITER_PAD}"
WORK_DIR="$OUT_DIR/work"
mkdir -p "$OUT_DIR"

if [ -f "$OUT_DIR/iteration_results.json" ]; then
    echo "[SKIP] iter $ITER already complete"; exit 0
fi

echo "=== task $ITER | cohort2_dragen | panel=20 | $(date) ==="
python /work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/run_single_iteration.py \
    --cohort cohort2_dragen \
    --panel-size 20 \
    --iteration $ITER \
    --seed $SEED \
    --work-dir "$WORK_DIR" \
    --out-dir  "$OUT_DIR"
echo "=== finished $ITER | $(date) ==="
