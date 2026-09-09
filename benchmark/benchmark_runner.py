#!/usr/bin/env python3
"""
benchmark_runner.py
====================
Orchestrates benchmark execution across cohorts and reference-panel sizes.

Three modes:
  dry-run  – print the plan, do nothing
  local    – run iterations sequentially on the current node (for smoke tests)
  slurm    – generate SLURM ARRAY job scripts; submit separately via sbatch

Architecture:
  For SLURM mode, one array job is generated per (cohort, panel_size)
  combination. The array task ID == iteration number; seed = SEED_BASE + iter.
  Each task runs in its own isolated work_dir so array tasks can execute
  concurrently on different nodes without clobbering each other.
"""

import os
import sys
import argparse
import subprocess
import logging

# Resolve benchmark_config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import benchmark_config as cfg


def parse_args():
    p = argparse.ArgumentParser(
        description="Run CyriPanel benchmark iterations (local or SLURM array)"
    )
    p.add_argument("--mode", required=True, choices=["dry-run", "local", "slurm"],
                   help="dry-run: print plan; local: sequential; slurm: generate array scripts")
    p.add_argument("--cohort", default=None,
                   help="Run only this cohort (default: all cohorts in config)")
    p.add_argument("--panel-size", type=int, default=None,
                   help="Run only this panel size (default: all from config)")
    p.add_argument("--n-iter", type=int, default=None,
                   help="Override number of iterations (default: N_ITERATIONS from config)")
    p.add_argument("--start-iter", type=int, default=0,
                   help="Iteration to start from (default 0) — useful for resuming")
    p.add_argument("--skip-done", action="store_true", default=True,
                   help="Skip iterations where iteration_results.json already exists")
    return p.parse_args()


def iter_out_dir(cohort, panel_size, iteration):
    """Path: {BENCHMARK_ROOT}/{cohort}/panel{N}/iter_{III}/"""
    return os.path.join(
        cfg.BENCHMARK_ROOT,
        cohort,
        "panel%d" % panel_size,
        "iter_%03d" % iteration,
    )


def iter_work_dir(cohort, panel_size, iteration):
    """Path: {iter_out_dir}/work/"""
    return os.path.join(iter_out_dir(cohort, panel_size, iteration), "work")


def is_iteration_done(cohort, panel_size, iteration):
    """True if iteration_results.json already exists."""
    return os.path.exists(
        os.path.join(iter_out_dir(cohort, panel_size, iteration),
                     "iteration_results.json")
    )


def compute_wall_time(n_iter):
    """
    For SLURM array mode, each task runs ONE iteration, so the per-task
    wall-time is SLURM_MINUTES_PER_ITER with a 10% buffer.
    Returns 'HH:MM:SS' string.
    """
    minutes = int(cfg.SLURM_MINUTES_PER_ITER * 1.1) + 1
    h = minutes // 60
    m = minutes % 60
    return "%02d:%02d:00" % (h, m)


# ──────────────────────────────────────────────────────────────────────
# LOCAL mode (sequential, for smoke tests)
# ──────────────────────────────────────────────────────────────────────

def run_local(cohort, panel_size, n_iter, start_iter, skip_done):
    logging.info("=" * 60)
    logging.info("  %s | panel_size=%d | iterations %d..%d (LOCAL)",
                 cohort, panel_size, start_iter, start_iter + n_iter - 1)
    logging.info("=" * 60)

    n_ok = n_skip = n_fail = 0

    for i in range(n_iter):
        it = start_iter + i
        seed = cfg.SEED_BASE + it

        if skip_done and is_iteration_done(cohort, panel_size, it):
            logging.info("[SKIP] iter %d already complete", it)
            n_skip += 1
            continue

        out_dir = iter_out_dir(cohort, panel_size, it)
        work_dir = iter_work_dir(cohort, panel_size, it)
        os.makedirs(out_dir, exist_ok=True)

        logging.info("[RUN] iter %d | %s | n=%d", it, cohort, panel_size)
        cmd = [
            sys.executable,
            os.path.join(SCRIPT_DIR, "run_single_iteration.py"),
            "--cohort", cohort,
            "--panel-size", str(panel_size),
            "--iteration", str(it),
            "--seed", str(seed),
            "--work-dir", work_dir,
            "--out-dir", out_dir,
        ]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            n_ok += 1
        else:
            n_fail += 1
            logging.error("[FAIL] iter %d exited with code %d",
                          it, result.returncode)

    logging.info("-- Done: %d OK, %d skipped, %d failed --",
                 n_ok, n_skip, n_fail)


# ──────────────────────────────────────────────────────────────────────
# SLURM ARRAY mode
# ──────────────────────────────────────────────────────────────────────

_ARRAY_SCRIPT_TEMPLATE = '''#!/bin/bash
#SBATCH --job-name=cyri_{cohort}_n{panel_size}
#SBATCH --partition={partition}
#SBATCH --account={account}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={wall_time}
#SBATCH --array={array_range}%{max_concurrent}
#SBATCH --output={log_dir}/{cohort}_n{panel_size}_iter%a_%A.out
#SBATCH --error={log_dir}/{cohort}_n{panel_size}_iter%a_%A.err

set -euo pipefail

# Activate virtual environment
source {venv_activate}

ITER=$SLURM_ARRAY_TASK_ID
SEED=$(( {seed_base} + ITER ))
ITER_PAD=$(printf '%03d' "$ITER")

OUT_DIR="{benchmark_root}/{cohort}/panel{panel_size}/iter_${{ITER_PAD}}"
WORK_DIR="$OUT_DIR/work"

mkdir -p "$OUT_DIR"

# Skip if already complete (resumable arrays)
if [ -f "$OUT_DIR/iteration_results.json" ]; then
    echo "[SKIP] iter $ITER already complete"
    exit 0
fi

echo "=== Array task $ITER | {cohort} | panel_size={panel_size} | $(date) ==="

python {script_dir}/run_single_iteration.py \\
    --cohort {cohort} \\
    --panel-size {panel_size} \\
    --iteration $ITER \\
    --seed $SEED \\
    --work-dir "$WORK_DIR" \\
    --out-dir  "$OUT_DIR"

echo "=== Array task $ITER finished | $(date) ==="
'''


def generate_slurm_scripts(cohorts_to_run, panel_sizes_to_run, n_iter, start_iter):
    """
    Generate one array job script per (cohort, panel_size) combination.
    Writes a submit_all.sh convenience wrapper.
    """
    scripts_dir = os.path.join(cfg.BENCHMARK_ROOT, "slurm_scripts")
    log_dir = os.path.join(cfg.BENCHMARK_ROOT, "slurm_logs")
    os.makedirs(scripts_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    wall_time = compute_wall_time(n_iter)
    array_range = "%d-%d" % (start_iter, start_iter + n_iter - 1)

    generated = []
    for cohort in cohorts_to_run:
        for panel_size in panel_sizes_to_run:
            script_name = "cyri_%s_n%d.sh" % (cohort, panel_size)
            script_path = os.path.join(scripts_dir, script_name)

            content = _ARRAY_SCRIPT_TEMPLATE.format(
                cohort=cohort,
                panel_size=panel_size,
                partition=cfg.SLURM_PARTITION,
                account=cfg.SLURM_ACCOUNT,
                cpus=cfg.SLURM_CPUS_PER_TASK,
                mem=cfg.SLURM_MEM,
                wall_time=wall_time,
                array_range=array_range,
                max_concurrent=cfg.SLURM_ARRAY_MAX_CONCURRENT,
                log_dir=log_dir,
                venv_activate=cfg.VENV_ACTIVATE,
                seed_base=cfg.SEED_BASE,
                benchmark_root=cfg.BENCHMARK_ROOT,
                script_dir=SCRIPT_DIR,
            )
            with open(script_path, "w") as f:
                f.write(content)
            os.chmod(script_path, 0o755)
            generated.append(script_path)
            logging.info("  Generated: %s", script_path)
            logging.info("    array=%s  max_concurrent=%d  wall/task=%s",
                         array_range, cfg.SLURM_ARRAY_MAX_CONCURRENT, wall_time)

    # submit_all.sh
    submit_all = os.path.join(scripts_dir, "submit_all.sh")
    with open(submit_all, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all CyriPanel benchmark array jobs.\n")
        f.write("# Run from any directory.\n\n")
        for script in generated:
            f.write("sbatch %s\n" % script)
    os.chmod(submit_all, 0o755)
    logging.info("  Submit wrapper: %s", submit_all)

    return generated


# ──────────────────────────────────────────────────────────────────────
# DRY-RUN mode
# ──────────────────────────────────────────────────────────────────────

def dry_run(cohorts_to_run, panel_sizes_to_run, n_iter, start_iter):
    logging.info("=" * 70)
    logging.info("  CyriPanel Benchmark -- DRY RUN")
    logging.info("=" * 70)
    logging.info("  CyriPanel src:   %s", cfg.CYRIPANEL_SRC_DIR)
    logging.info("  Benchmark root:  %s", cfg.BENCHMARK_ROOT)
    logging.info("  Iterations:      %d (starting from %d)", n_iter, start_iter)
    logging.info("  Seed base:       %d", cfg.SEED_BASE)
    logging.info("  Wall time/task:  %s  (%d min/iter x 1.1 buffer)",
                 compute_wall_time(n_iter), cfg.SLURM_MINUTES_PER_ITER)
    logging.info("  Array concurr.:  %d", cfg.SLURM_ARRAY_MAX_CONCURRENT)

    total_runs = 0
    total_jobs = 0
    for cohort in cohorts_to_run:
        c = cfg.COHORTS[cohort]
        logging.info("  -- %s --", cohort)
        logging.info("     BAM dir:       %s", c["bam_dir"])
        logging.info("     Total samples: %d", c["total_samples"])
        logging.info("     Diploid CSV:   %s", c["diploid_csv"])
        logging.info("     Gold CSV:      %s", c["gold_csv"])
        for ps in panel_sizes_to_run:
            n_test = c["total_samples"] - ps
            logging.info("     panel_size=%-2d:  %d iters x %d test samples  "
                         "=> 1 array job w/ %d tasks",
                         ps, n_iter, n_test, n_iter)
            total_runs += n_iter * n_test
            total_jobs += 1

    logging.info("  Total CyriPanel runs:   %d", total_runs)
    logging.info("  Total SLURM array jobs: %d", total_jobs)
    logging.info("=" * 70)


# ──────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
    )

    cohorts_to_run = [args.cohort] if args.cohort else list(cfg.COHORTS.keys())
    panel_sizes_to_run = [args.panel_size] if args.panel_size \
        else cfg.REFERENCE_PANEL_SIZES
    n_iter = args.n_iter if args.n_iter is not None else cfg.N_ITERATIONS

    # Sanity checks
    for c in cohorts_to_run:
        if c not in cfg.COHORTS:
            logging.error("Unknown cohort: %s", c)
            sys.exit(1)

    if args.mode == "dry-run":
        dry_run(cohorts_to_run, panel_sizes_to_run, n_iter, args.start_iter)

    elif args.mode == "local":
        os.makedirs(cfg.BENCHMARK_ROOT, exist_ok=True)
        for cohort in cohorts_to_run:
            for panel_size in panel_sizes_to_run:
                run_local(cohort, panel_size, n_iter,
                          args.start_iter, args.skip_done)

    elif args.mode == "slurm":
        os.makedirs(cfg.BENCHMARK_ROOT, exist_ok=True)
        logging.info("Generating SLURM array scripts...")
        scripts = generate_slurm_scripts(
            cohorts_to_run, panel_sizes_to_run, n_iter, args.start_iter
        )
        logging.info("")
        logging.info("Generated %d array job script(s).", len(scripts))
        logging.info("To submit:")
        logging.info("  bash %s/slurm_scripts/submit_all.sh",
                     cfg.BENCHMARK_ROOT)

    return 0


if __name__ == "__main__":
    sys.exit(main())
