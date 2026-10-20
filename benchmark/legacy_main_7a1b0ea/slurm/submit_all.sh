#!/bin/bash
# Submit all CyriPanel benchmark array jobs.
# Run from any directory.

sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort1_n5.sh
sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort1_n10.sh
sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort1_n20.sh
sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort2_n5.sh
sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort2_n10.sh
sbatch /work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results/slurm_scripts/cyri_cohort2_n20.sh
