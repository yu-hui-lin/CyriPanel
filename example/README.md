# CyriPanel reproducibility example

This example runs CyriPanel with five assay-matched diploid references and two test samples, then compares the results with frozen observed outputs.

**Release status:** The data archive is available at https://doi.org/10.6084/m9.figshare.34724976.v1. Use a checkout containing this example/ directory and the matching archive. The original v1.0.1 tag predates this example.

## Data source

The seven example BAMs were prepared from targeted-capture sequencing data generated for the CyriPanel study using genomic DNA purchased from the [Coriell Institute for Medical Research](https://www.coriell.org/). All seven samples belong to study cohort 2 and are represented in [HPRC Release 2](https://humanpangenome.org/hprc-data-release-2/). Sample identifiers and their reference or test roles are listed in [manifest.tsv](manifest.tsv). The raw reads for these samples are in the European Nucleotide Archive under study PRJEB128569, as runs ERR18060456 (NA18522), ERR18060366 (HG00280), ERR18060363 (HG00133), ERR18060464 (NA19043), ERR18060396 (HG01530), ERR18060377 (HG00639) and ERR18060416 (HG02129).

Reads were aligned to GRCh38 using ALT-aware BWA-MEM, and duplicates were marked with Picard. The example BAMs retain capture-panel targets and their mates to preserve the inputs used for panel-wide normalization. HPRC haplotype-resolved assemblies were used to derive the study's reference-standard diplotypes.

## Data license

The study-generated BAM files and accompanying indexes, cleaned headers, and metadata in `cyripanel-example-data.tar.gz` are licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

These data may be reused, modified, and redistributed for any purpose, including commercial use. Users must credit the dataset creators, include a link to the license, and identify any changes they make.

The CyriPanel software is separately licensed under the GNU General Public License, version 3 or later, as described in [LICENSE](../LICENSE).

## Data and expected results

The reference samples are **NA18522, HG00280, HG00133, NA19043, and HG01530**. Sample roles and paths are recorded in [manifest.tsv](manifest.tsv); its paths are relative to this directory.

| Test sample | Expected observed genotype | Filter |
| --- | --- | --- |
| HG00639 | `*1/*1` | PASS |
| HG02129 | `*1/*36+*10` | PASS |

| Test sample | CYP2D6plusREP6_hapcn2 MeanRatio | CYP2D7 MeanRatio | D7spacer_hapcn1 MeanRatio |
| --- | --- | --- | --- |
| HG00639 | 1.02 | 1.05 | 1.03 |
| HG02129 | 1.34 | 1.01 | 1.43 |

These are observed software outputs, with genotype strings preserved exactly as emitted. `PASS` is the caller's filter label. The example checks computational reproducibility for these two samples; independent biological accuracy and the full-cohort benchmark require separate evaluation.

## Environment

Run the commands below from the repository root on **Linux x86_64**. Conda and Bash are required.

Create a new environment using the tested YAML installation workflow:

```bash
bash setup_environment.sh cyripanel-example
conda activate cyripanel-example
```

If the validated `cyripanel-example` environment already exists, activate it instead of creating it again. To choose a new name, use that name in both commands.

The reconstruction test used the name `cyripanel-rebuild`. Its key versions were Python 3.10.22, R 4.2.3, CNVPanelizer 1.30.0, pysam 0.23.0, and samtools 1.24. [environment.yml](../environment.yml) contains the resolved package versions used to create it; [environment.example.yml](../environment.example.yml) contains the direct dependencies used for the initial setup.

The earlier archive-installation test used the original environment. The rebuilt environment then reproduced all 10 primary comparisons and both complete parsed JSON records using the previously restored code and data, with the existing package cache available. See the [environment rebuild evidence](validation/environment-rebuild/evidence.json).

The [original explicit package list](../reproducibility/conda-linux-64.explicit.txt) and [rebuilt environment package list](validation/environment-rebuild/conda-linux-64.explicit.txt) record installed package builds and URLs. Direct recreation from these explicit lists has not been separately tested; optional instructions are in the [root README](../README.md#installation).

## Install the data archive

Public download: https://doi.org/10.6084/m9.figshare.34724976.v1. Obtain the archive matching these metadata:

- Filename: `cyripanel-example-data.tar.gz`
- Size: **2,691,360,951 bytes** (approximately 2.507 GiB)
- SHA256: `194110431cc3244b533d151563ba1e62b2da71435d3bfc0a4a3d34b7581cacf6`

Keep the downloaded archive outside the repository. Replace the path below with its location:

```bash
archive="/path/to/cyripanel-example-data.tar.gz"
python -I example/install_data.py --archive "$archive" --check
```

This lightweight check verifies the small metadata files and archive size. It does not hash or extract the archive, or check whether the installation destination is available.

Install the data:

```bash
python -I example/install_data.py --archive "$archive"
```

The installer verifies the complete archive SHA256, all 24 member checksums, and gzip integrity. It installs seven BAMs, seven BAI indexes, seven cleaned SAM headers, and supporting metadata beneath `example/data/`. It preserves the existing matching metadata files.

Installation requires `example/data/bams/`, `example/data/headers/`, and `example/data/header-preparation.json` to be absent. Existing files, directories, or symbolic links are refused. Use a fresh code copy for another installation test.

On an HPC system, perform full installation and analysis in an allocated compute job. Activate the environment and set the archive path inside that job. The lightweight `--check` commands do not read entire BAM files.

## Run and compare

After installation, check the example inputs and frozen outputs:

```bash
bash example/run_example.sh --check
```

Then run the example, choosing an output directory that does not already exist:

```bash
bash example/run_example.sh --out example/runs/run-01
```

The runner verifies all 14 BAM/index checksums, runs `samtools quickcheck` on the BAMs, creates an isolated working copy, runs CNVPanelizer and CyriPanel, and compares the resulting outputs. Input BAMs and frozen expected files are retained unchanged.

Successful analysis prints:

```text
HG00639 *1/*1 PASS
HG02129 *1/*36+*10 PASS
EXAMPLE_OUTPUTS_MATCH: 10 metrics and 2 complete parsed JSON records
```

The 10 comparisons are genotype, filter, and three stored MeanRatio values for each sample. Complete parsed JSON records are also compared. MeanRatio equality is checked at the precision stored in the reports; the check does not compare every column in the CNV CSV files.

The chosen output directory contains:

- `results/cyripanel.tsv` and `results/cyripanel.json`: caller outputs.
- `results/cnv_reports/`: the two CNVPanelizer reports.
- `comparison.tsv` and `json-comparison.json`: comparison results.
- `input-SHA256SUMS` and `source-sha256.json`: input and runtime provenance.
- `caller.stdout.log`, `caller.stderr.log`, and environment logs: execution details.
- `EXAMPLE_COMPLETE`: written only after all comparisons pass.

If execution fails, inspect the error message and output logs. Use a new output directory for a subsequent run.

## Data preparation and validation

The seven GRCh38 BAM slices retain capture-panel targets and fetched mates so CNVPanelizer can perform panel-wide normalization. The [panel BED](../data/PGxProbe_region_hg38.bed) contains 784 data rows without a header; the tested R wrapper reads it with `skip=0`.

All **5,488 region counts** (seven samples multiplied by 784 regions) matched between the full BAMs and panel slices under the tested counting settings. Both test samples also matched the full-BAM run for the 10 primary comparisons and complete parsed JSON records. These findings apply to these samples and settings.

After slicing, command-line fields (`@PG.CL`) were removed from BAM headers. Alignment-record bytes and the binary reference dictionary were verified unchanged by that header edit, and new BAI indexes were generated. Details are recorded in [data/bam-validation.json](data/bam-validation.json); the archive also contains cleaned headers and `data/header-preparation.json`.

Validation records:

- [expected/](expected/): frozen observed outputs and full-versus-sliced comparison records.
- [validation/evidence.json](validation/evidence.json): the completed public-runner test.
- [validation/archive-installation/evidence.json](validation/archive-installation/evidence.json): archive installation into a clean code copy followed by a successful two-sample run.
- [validation/environment-rebuild/evidence.json](validation/environment-rebuild/evidence.json): a new environment created from the exported YAML followed by a successful run using the previously restored code and data.
- [data/archive-report.json](data/archive-report.json): archive identity, size, and membership summary.

Checksums are provided in `expected/SHA256SUMS`, `validation/SHA256SUMS`, `validation/archive-installation/SHA256SUMS`, and `validation/environment-rebuild/SHA256SUMS`. The archive-installation snapshot manifest describes the historical tested code copy; later documentation and ignore-rule edits may differ from that snapshot.

Keep the BAM archive outside Git. Code, expected results, and small verification records belong in the repository; `example/data/bams/` remains ignored by Git.
