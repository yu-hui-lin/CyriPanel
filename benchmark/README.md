# Benchmark

Scripts and records for the evaluation in the CyriPanel preprint: CyriPanel against
Aldy 4, on two cohorts, under two aligners, with reference panels of 5, 10 and 20 samples
drawn at random over 30 iterations each. Each script carries a docstring describing what
it does; this file says how they fit together.

## Which script produced which table

Table numbers follow the preprint.

| Preprint table or analysis | Script | Reads |
|---|---|---|
| Table 2, Supplemental Table S2 | `capture_qc_metrics.sh`, then `inspect_seq_qc.py` | BAM files |
| Tables 3 and 4, single-candidate sensitivity analysis | `caller_comparison_n5.py` | `published_results/`, `aldy/`, `data/` |
| Table 5, concordance | `recalculate.py` | `published_results/per_sample_v2.csv` |
| Table 5, all columns | `table5.py` | `published_results/per_sample_v2.csv` |
| Reference-panel size, per-sample means and paired differences | `recalculate.py` | `published_results/per_sample_v2.csv` |
| Reference-panel size, Friedman and Wilcoxon tests | `common_sample_analysis.py` | per-iteration outputs |
| Discordant samples and sequencing quality | `inspect_seq_qc.py` | `published_results/CyriPanel_seq_qc.tsv`, sample calls from `regenerate_table3.py` |
| Supplemental Table S1 | none | `data/` |
| Supplemental Tables S3 and S4 | `supplemental_s3_s4.py` | `published_results/`, `aldy/`, `data/` |
| Supplemental Table S5 | `count_spanning_pairs.sh` | BAM files |
| Supplemental Table S6 | none | per-iteration outputs |
| Supplemental Table S7 | `mem_audit.py` for the memory columns | SLURM accounting records |

`recompute_concordance.py` read the per-iteration outputs and wrote `per_sample_v2.csv`
and `iteration_summary_v2.csv` in `published_results/`. It also reported the PASS, no-call
and ambiguous-call rates in Table 5, which `table5.py` derives again from those records.

`run_single_iteration.py` defines `normalize_genotype()` and `compare_genotype()`, which
score **both** tools, so neither is advantaged by the notation handling.

`regenerate_table3.py` belongs to an earlier draft. Its caller comparison used one modal
call per sample, McNemar tests and Newcombe intervals, and the preprint does not report
it. `inspect_seq_qc.py` still reads the per-sample calls it writes. On these data, its
modal-call rule flags the same discordant samples as the rule in Supplemental Methods S5.

`recalculate.py` is the script used for the preprint, adapted to this folder as its
docstring describes. Its outputs go to `recalculate_output/`. `table5.py` and
`supplemental_s3_s4.py` were written for release v1.0.2. They reproduce every value in
Table 5 and in Supplemental Tables S3 and S4 from the records in this folder.
Supplemental Table S6 and the runtime columns of Supplemental Table S7 have no script in
this release. They were tabulated on the cluster from the per-iteration outputs and the
SLURM accounting records.

`cache_test.sh` is the control job in Supplemental Methods S4. It shows that the memory
accounted on this cluster includes file cache.

## Regenerating the tables without re-running anything

`published_results/`, `aldy/aldy_results_long.csv` and `data/` hold the per-evaluation
records. From this directory:

    python3 caller_comparison_n5.py
    python3 recalculate.py
    python3 table5.py
    python3 supplemental_s3_s4.py

All four need only the Python standard library. `caller_comparison_n5.py` and `table5.py`
print their input checks first and end with the list of failed checks, which should read
`none`. `recalculate.py` stops with an error if `per_sample_v2.csv` differs from the
archived file.

## Columns in `published_results/`

`per_sample_v2.csv` has one row per CyriPanel evaluation. `gold` is the reference-standard
diplotype from `benchmark/data/`. `comparison_single` is `concordant` only when the call is a
single diplotype equal to `gold`, so a call that lists two candidates counts as discordant.
`iteration_summary_v2.csv` sums these per iteration (`n_concordant_single`,
`concordance_single`). The candidate-inclusive concordance reported in the preprint, which
accepts a call when any listed candidate equals `gold`, is computed from the same rows by
`recalculate.py` and `caller_comparison_n5.py`.

## Re-running from BAMs

`benchmark_runner.py` generates one SLURM array job per (cohort, panel size); `slurm/` holds
the scripts as submitted. Seeding is deterministic - iteration i uses seed 42 + i - so any
draw can be regenerated exactly.

`aldy/` holds the Aldy 4 side. Aldy was run as an exhaustive single-reference benchmark:
every known-diploid sample serves once as the sole reference profile, genotyping all others
(6,814 runs, none failed). `aldy/METHODS_provenance.txt` has the command line and versions.

`verify_fix549.sh` and `inspect_fix549.py` are the evidence that fixing star_caller.py line
549 changed no existing call: 10 samples run under both versions and diffed.

## Software versions

CyriPanel ran in the environment recorded in `requirements-validated.txt` at the top of the
repository. The analysis scripts behind the preprint ran under Python 3.12.2 with
numpy 2.5.2, scipy 1.18.1 and pandas 3.0.5.

## Paths

Scripts that read the BAM files, the per-iteration outputs or the SLURM accounting records
hard-code the paths of the environment the results were produced in (Taiwania 3, NCHC,
Taiwan). To run elsewhere, edit the `bam_dir`, `gold_csv` and `diploid_csv` entries in
`benchmark_config.py` and the path variables at the top of each script. They are left as
they were rather than parameterised, so that what is archived is what was run. The one
exception is release v1.0.2, which removed the superseded reference standard from these
scripts and renamed the columns that held it. The scoring is unchanged.

## Not included

- The 360 per-iteration output directories and Aldy's 6,814 raw outputs: large, and
  regenerable from the scripts here.
- BAM files; see the preprint's data-availability statement.
- SLURM accounting records.
- Development and smoke-test scripts, and the supplementary-workbook generators, which
  predate the corrected tables and would produce numbers inconsistent with them.
