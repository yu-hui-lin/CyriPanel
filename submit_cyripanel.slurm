#!/bin/bash
#SBATCH -J CyriPanel
#SBATCH -A MST109178
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --output=/staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel/logs/cyripanel_%j.out
#SBATCH --error=/staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel/logs/cyripanel_%j.err

# ============================================
# CyriPanel SLURM Submission Script
# ============================================

# Load conda
source /home/u7715055/.conda/envs/cyripanel/bin/activate
conda activate cyripanel

# Get arguments
MANIFEST=$1
PREFIX=$2

# Check arguments
if [ -z "$MANIFEST" ] || [ -z "$PREFIX" ]; then
    echo "Usage: sbatch submit_cyripanel.slurm <manifest_file> <output_prefix>"
    exit 1
fi

# Change to CyriPanel directory
cd /staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel

# Run CyriPanel
bash run_cyripanel.sh --manifest "$MANIFEST" "$PREFIX"

echo "CyriPanel job completed at $(date)"
