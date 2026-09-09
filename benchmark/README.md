# Benchmark

Everything needed to reproduce the evaluation in the CyriPanel manuscript: CyriPanel
against Aldy 4, on two cohorts, under two aligners, with reference panels of 5, 10 and 20
samples drawn at random over 30 iterations each. Each script carries a docstring
describing what it does; this file says how they fit together.

## Which script produced which table

| Manuscript table | Script | Writes |
|---|---|---|
| Table 1, Suppl. S5 | `capture_qc_metrics.sh`, then `inspect_seq_qc.py` | `published_results/CyriPanel_seq_qc.tsv` |
| Table 2 | `recompute_concordance.py` | `published_results/iteration_summary_v2.csv` and `per_sample_v2.csv` |
| Table 3a-c, Suppl. S1 and S3 | `regenerate_table3.py` | CSVs plus a markdown summary |
| Panel-size comparison | `common_sample_analysis.py` | per-arm and per-sample summaries |

`run_single_iteration.py` defines `normalize_genotype()` and `compare_genotype()`, which
score **both** tools, so neither is advantaged by the notation handling.

## Regenerating the tables without re-running anything

`published_results/` and `aldy/aldy_results_long.csv` hold the per-evaluation records:

    python3 regenerate_table3.py --selftest
    python3 regenerate_table3.py

The self-test checks the Newcombe method-10 interval against a worked example before any
data is read, and prints SELFTEST PASS.

## Re-running from BAMs

`benchmark_runner.py` generates one SLURM array job per (cohort, panel size); `slurm/` holds
the scripts as submitted. Seeding is deterministic - iteration i uses seed 42 + i - so any
draw can be regenerated exactly.

`aldy/` holds the Aldy 4 side. Aldy was run as an exhaustive single-reference benchmark:
every known-diploid sample serves once as the sole reference profile, genotyping all others
(6,814 runs, none failed). `aldy/METHODS_provenance.txt` has the command line and versions.

`verify_fix549.sh` and `inspect_fix549.py` are the evidence that fixing star_caller.py line
549 changed no existing call: 10 samples run under both versions and diffed.

## Paths

Every script hard-codes the paths of the environment the published results were produced in
(Taiwania 3, NCHC, Taiwan). To run elsewhere, edit the `bam_dir`, `gold_csv` and
`diploid_csv` entries in `benchmark_config.py` and the `ROOT` variable at the top of each
analysis script. They are left as they were rather than parameterised, so that what is
published is what was run.

## Not included

- The 180 per-iteration output directories and Aldy's 6,814 raw outputs: large, and
  regenerable from the scripts here.
- BAM files; see the manuscript's data-availability statement.
- The v1 gold standard: superseded, and the released version of deCYPher is v2.
- Development and smoke-test scripts, and the supplementary-workbook generators, which
  predate the corrected tables and would produce numbers inconsistent with them.
