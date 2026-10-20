# Reproduce the benchmark analyses

This directory contains archived evaluation records and three entry points for
reconstructing manuscript Tables 3-5, exploratory panel-size tests, and sensitivity
analyses. These entry points use repository-relative inputs and require no BAMs,
Slurm account, or access to the original cluster.

The six input CSVs are pinned to release `v1.0.1`, commit
`5d33b820dd2ba8b874f1bbea51f411b5f081856f`. The three reconstruction scripts and
`reproduced_*` outputs were added after that release. They document new
reconstructions from those archived inputs, not a rerun of the original full
benchmark. Original genotyping and the current example validation have separate
provenance.

## Inputs and entry points

| Input, relative to this directory | Records |
| --- | ---: |
| `published_results/per_sample_v2.csv` | 16,680 CyriPanel evaluations, 360 iterations |
| `aldy/aldy_results_long.csv` | 6,814 retained Aldy major-call records |
| `data/CYP2D6_cohort1_gold_standard_v2.csv` | 44 samples |
| `data/CYP2D6_cohort2_gold_standard_v2.csv` | 72 samples |
| `data/CYP2D6_cohort1_33diploids.csv` | 33 eligible references |
| `data/CYP2D6_cohort2_28diploids.csv` | 28 eligible references |

| Entry point | Analysis | Recorded outputs |
| --- | --- | --- |
| [reproduce_tables.py](reproduce_tables.py) | Tables 3, 4 and 5: 8, 28 and 6 display rows | [Table results](reproduced_tables/tables.md) |
| [panel_size_tests.py](panel_size_tests.py) | Friedman and Wilcoxon tests; minimum 1 and 5 evaluations | [Panel-size results](reproduced_panel_tests/summary.md) |
| [sensitivity_analyses.py](sensitivity_analyses.py) | Single-candidate scoring, hypothetical *49 rescoring, discordance classification | [Sensitivity methods](reproduced_sensitivity/summary.md), [eight result rows](reproduced_sensitivity/sensitivity.csv) |

All three scripts check the input files they use against recorded SHA256 hashes.
The latter two load a hash-pinned copy of `reproduce_tables.py` for common scoring
functions. Keep these scripts together. Expected table display-row hashes and
analysis checks were independently verified during preparation; the scripts do
not read the manuscript on each run.

## Run from the repository root

First create and activate the environment following the [installation instructions](../README.md).
The reconstruction was tested with Python 3.10.22. Table reproduction and
sensitivity analysis use the Python standard library; panel-size tests require
NumPy 1.26.4 and SciPy 1.15.2 and check those versions explicitly.

The following commands create a fresh temporary result directory. Each script
refuses an existing output directory, so recorded evidence is not overwritten.

```bash
(
set -eu
repro_dir=$(mktemp -d "${TMPDIR:-/tmp}/cyripanel-benchmark.XXXXXX")
python -I benchmark/reproduce_tables.py --outdir "$repro_dir/tables"
python -I benchmark/panel_size_tests.py --outdir "$repro_dir/panel-tests"
python -I benchmark/sensitivity_analyses.py --outdir "$repro_dir/sensitivity"
for part in tables panel-tests sensitivity; do
    (cd "$repro_dir/$part" && sha256sum -c SHA256SUMS)
done
printf 'Results: %s\n' "$repro_dir"
)
```

Success markers are `TABLES_3_4_5_REPRODUCED`, `PANEL_SIZE_TESTS_REPRODUCED`, and
`SENSITIVITY_ANALYSES_REPRODUCED`. Each output directory includes `validation.txt`,
`INPUT_SHA256SUMS`, and `SHA256SUMS`, together with its analysis results. The last
manifest validates files within that output directory; input manifest paths are
relative to `benchmark/` in the layout above. Validation metadata may differ with
the Python version even when result tables agree.

The recorded outputs in this checkout can be checked without rerunning analyses:

```bash
(
set -eu
for folder in reproduced_tables reproduced_panel_tests reproduced_sensitivity; do
    (cd "benchmark/$folder" && sha256sum -c SHA256SUMS)
done
)
```

## Scoring and interpretation

CyriPanel scoring accepts a matching candidate among semicolon-separated calls.
No-calls remain in evaluation denominators. Aldy scoring re-parses the retained
`raw_major` field into core star alleles. Its archived `n_solutions=1` values do
not establish that original Aldy solution files contained no other candidates.

Tables 3 and 4 use CyriPanel panels of five references. Table 3 averages the
per-sample concordance proportions with equal sample weights. Table 4 reports
correct/evaluated counts for each hybrid sample. Table 5 uses per-iteration
concordance over 30 iterations, summarized as mean and sample SD (denominator
29); failures include no-calls. These aggregation levels are different.

Panel-size tests compare sample-level rates across 5, 10 and 20 references within
each cohort and aligner. Each included sample must have at least 1 evaluation at
every panel size, or at least 5 in the minimum-evaluation sensitivity analysis.
The latter excludes HG01530 in both
cohort-2 arms because its evaluation counts are 24, 17 and 4. Wilcoxon differences
are calculated as exact fractions before conversion to floating point to preserve
true ties. The p-values are asymptotic, two-sided, with zero differences removed
and no continuity correction; Holm adjustment covers the three panel pairs within
each arm and minimum-evaluation setting. All-zero differences are assigned p=1
with an undefined statistic; all-tied Friedman results are reported as NA.

Correct tie handling changes the cohort-2 BWA-MEM 5-versus-10 raw p-value from
0.345231 to 0.278517 and its Holm p-value from 1.000000 to 0.835551 under both
minimum-evaluation settings. The rounded
minimum p-value summaries in the manuscript remain consistent. The tests are
exploratory: they do not account for shared-reference dependence or all uncertainty
from reference selection. Nonsignificance does not establish equivalence.

The *49 sensitivity scenario makes all evaluations of HG00597, HG00609 and NA24631
incorrect, holding denominators fixed. It is hypothetical rescoring, not an
experiment removing the g.42129180A>T exemption or a proven bound on its effects.
No *49 call was observed in 15,960 archived evaluations of the 112 samples whose
reference-standard diplotypes lacked *49, across all panel sizes and both aligners;
no-calls are included in this evaluation count.

Across all 16,680 CyriPanel evaluations, 986 are discordant under candidate-inclusive
scoring: 106 no-calls, 430 phased copy-distribution differences, 172 unphased
total-copy differences, 270 phased normal/hybrid-count differences, and 8 unresolved
unphased evaluations of HG00597. The classifier treats *5 as zero copies and collapses *13, *36,
*68 and *83 into one hybrid class. Ordinary star identities are ignored. These
coarse categories do not establish biological equivalence or the cause of an error;
the eight unresolved calls prevent attributing all discordance to structure.

## Relationship to BAM processing and historical scripts

The [executable example](../example/README.md) separately validates five reference
BAMs and two test BAMs, including archive restoration and environment rebuilding.
Its slices preserve panel-wide counts; they are not CYP2D-locus-only slices.
That two-sample example does not establish full-cohort equivalence after code changes.

The CSV reconstructions do not recover all CNVPanelizer matrices, site-level counts,
or original CyriPanel/Aldy output directories. The complete original BAM benchmark
and Aldy candidate-selection procedure are not rerun by the commands above.
Additional diagnostic claims, such as read-pair spanning and per-iteration depth
mechanisms, require their own supporting inputs.

Older benchmark runners, configuration, analysis scripts and Slurm files remain
as historical material. Some contain original site paths and resource settings;
they are not prerequisites for these three reconstruction entry points. Adapting
them for a new full-BAM benchmark requires separate configuration and validation.
Use the table mapping above for the current manuscript, rather than older script
names or table numbering. `v1.0.1` remains the source of the pinned input CSVs;
new software/data release identifiers must be recorded when those releases are published.
