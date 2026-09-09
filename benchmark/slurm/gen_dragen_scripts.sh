NEW=/work/u7715055/staging/biology/u7715055
SD=$NEW/cyripanel_benchmark_results/slurm_scripts
mkdir -p $NEW/cyripanel_benchmark_results/slurm_logs

for COHORT in cohort1_dragen cohort2_dragen; do
  for P in 5 10 20; do
    TAG="${COHORT}_n${P}"
    cat > $SD/cyri_${TAG}.sh << INNER
#!/bin/bash
#SBATCH --job-name=cyri_${TAG}
#SBATCH --partition=ngs92G
#SBATCH --account=MST109178
#SBATCH --cpus-per-task=14
#SBATCH --mem=92g
#SBATCH --time=24:13:00
#SBATCH --array=0-29%10
#SBATCH --output=${NEW}/cyripanel_benchmark_results/slurm_logs/${TAG}_iter%a_%A.out
#SBATCH --error=${NEW}/cyripanel_benchmark_results/slurm_logs/${TAG}_iter%a_%A.err

set -euo pipefail
source /home/u7715055/cyrius_env/bin/activate

ITER=\$SLURM_ARRAY_TASK_ID
SEED=\$(( 42 + ITER ))
ITER_PAD=\$(printf '%03d' "\$ITER")

OUT_DIR="${NEW}/cyripanel_benchmark_results/${COHORT}/panel${P}/iter_\${ITER_PAD}"
WORK_DIR="\$OUT_DIR/work"
mkdir -p "\$OUT_DIR"

if [ -f "\$OUT_DIR/iteration_results.json" ]; then
    echo "[SKIP] iter \$ITER already complete"; exit 0
fi

echo "=== task \$ITER | ${COHORT} | panel=${P} | \$(date) ==="
python ${NEW}/CyriPanel/benchmark/run_single_iteration.py \\
    --cohort ${COHORT} \\
    --panel-size ${P} \\
    --iteration \$ITER \\
    --seed \$SEED \\
    --work-dir "\$WORK_DIR" \\
    --out-dir  "\$OUT_DIR"
echo "=== finished \$ITER | \$(date) ==="
INNER
    echo "generated: cyri_${TAG}.sh"
  done
done
