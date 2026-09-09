#!/usr/bin/env python3
"""
CyriPanel Benchmark Configuration
===================================
Central configuration for all benchmark scripts.
Edit paths and parameters here before running.
"""

# ──────────────────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────────────────

# Root of the CyriPanel SOURCE tree (where star_caller.py lives).
# This is NEVER modified by the benchmark — scripts only read from here.
CYRIPANEL_SRC_DIR = "/work/u7715055/staging/biology/u7715055/CyriPanel"

# Top-level directory for all benchmark outputs (results + per-iteration work dirs).
# Kept OUTSIDE CyriPanel to avoid polluting the git repo and to prevent
# concurrent array jobs from interfering with other users of CyriPanel.
BENCHMARK_ROOT = "/work/u7715055/staging/biology/u7715055/cyripanel_benchmark_results"

# ──────────────────────────────────────────────────────────────────────
# COHORT DEFINITIONS
# ──────────────────────────────────────────────────────────────────────
# Each cohort dict must contain:
#   bam_dir        – directory holding all BAM files for the cohort
#                    (subdirectory layout: bam_dir/{sample_id}/{sample_id}{bam_suffix})
#   diploid_csv    – CSV listing known-diploid sample IDs (col: Sample_ID)
#   gold_csv       – CSV listing gold-standard genotypes for ALL samples
#                    (cols: Sample_ID, Genotype)
#   total_samples  – expected number of samples (for sanity check)
#   bam_suffix     – file extension to match (.markdup.bam, .sorted.bam, etc.)

COHORTS = {
    "cohort1": {
        "bam_dir": "/work/u7715055/staging/biology/u7715055/bwa_results/HPRC_PGpanel_altaware",
        "diploid_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_33diploids.csv",
        "gold_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_gold_standard_v2.csv",
        "total_samples": 44,
        "bam_suffix": ".markdup.bam",
    },
    "cohort2": {
        "bam_dir": "/work/u7715055/staging/biology/u7715055/bwa_results/PG_20250711_altaware",
        "diploid_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_28diploids.csv",
        "gold_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_gold_standard_v2.csv",
        "total_samples": 72,
        "bam_suffix": ".markdup.bam",
    },
    "cohort1_dragen": {
        "bam_dir": "/work/u7715055/staging/biology/u7715055/dragen_markdup/cohort1",
        "diploid_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_33diploids.csv",
        "gold_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort1_gold_standard_v2.csv",
        "total_samples": 44,
        "bam_suffix": ".bam",
    },
    "cohort2_dragen": {
        "bam_dir": "/work/u7715055/staging/biology/u7715055/dragen_markdup/cohort2",
        "diploid_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_28diploids.csv",
        "gold_csv": "/work/u7715055/staging/biology/u7715055/CyriPanel/benchmark/data/CYP2D6_cohort2_gold_standard_v2.csv",
        "total_samples": 72,
        "bam_suffix": ".bam",
    },
}

# ──────────────────────────────────────────────────────────────────────
# BENCHMARK PARAMETERS
# ──────────────────────────────────────────────────────────────────────

# Reference panel sizes to test. Start with just 5 — extend later.
REFERENCE_PANEL_SIZES = [5, 10, 20]

# Number of random iterations per (cohort, panel_size) combination.
N_ITERATIONS = 30

# Random seed base. Iteration i uses seed = SEED_BASE + i.
SEED_BASE = 42

# ──────────────────────────────────────────────────────────────────────
# CSV COLUMN NAMES
# ──────────────────────────────────────────────────────────────────────
# diploid_csv: single-column  "Sample_ID"
# gold_csv:    two-column      "Sample_ID,Genotype"
COL_SAMPLE_ID = "Sample_ID"
COL_GENOTYPE = "Genotype"

# ──────────────────────────────────────────────────────────────────────
# SLURM SETTINGS
# ──────────────────────────────────────────────────────────────────────
SLURM_PARTITION = "ngs92G"
SLURM_ACCOUNT = "MST109178"
SLURM_CPUS_PER_TASK = 14
SLURM_MEM = "92g"

# Estimated wall-time per iteration (minutes).
# Measured: ~300 min for cohort1 panel20 with 24 test samples.
# panel_size=5 has more test samples (39 for cohort1, 67 for cohort2),
# so we allocate more time. 360 min = 6h gives ~10% headroom.
SLURM_MINUTES_PER_ITER = 1320

# Maximum concurrent tasks in the SLURM array (--array=0-29%N)
SLURM_ARRAY_MAX_CONCURRENT = 10

# Virtual environment to activate in SLURM jobs (for scipy, pysam, etc.)
VENV_ACTIVATE = "/home/u7715055/cyrius_env/bin/activate"
