# CyriPanel

CYP2D6 star-allele genotyping from laboratory-developed targeted sequencing panels.

CyriPanel is a computational research tool developed for methodological evaluation of targeted capture sequencing data. It is **not** intended for clinical diagnostics, patient management, or prescribing decisions.

CyriPanel is derived from [Cyrius](https://github.com/Illumina/Cyrius) (Illumina) and
incorporates the updated star-allele definitions of
[BCyrius](https://gitlab.com/andreassh/bcyrius). It keeps Cyrius's paralogue-aware
variant calling and star-allele matching, but replaces the copy-number step: instead of a
Gaussian mixture model fitted to whole-genome coverage, copy number is estimated with
[CNVPanelizer](https://bioconductor.org/packages/CNVPanelizer/) against a small cohort of
known-diploid samples sequenced on the *same* panel. That makes CYP2D6 genotyping possible
on capture data, where WGS-derived coverage expectations do not transfer.

CyriPanel resolves single-nucleotide and indel variants, copy-number changes, and
CYP2D6–CYP2D7 hybrid alleles, and reports a diplotype with a filter status.

**GRCh38 only.**

## Reproducibility example

Start with the [worked example](example/README.md): five assay-matched diploid
references and two test samples, with frozen observed outputs and verification
records. Archive installation into a clean code copy was verified with the
original Conda environment. A subsequent environment rebuilt from
`environment.yml` reproduced both sample outputs using that restored copy.

These checks assess whether the observed software outputs can be reproduced for
two samples. They do not establish independent biological accuracy or reproduce
the full benchmark.
The example data sources and CC BY 4.0 license are documented in the [example README](example/README.md#data-source). The public download link will be added after the dataset is published on Zenodo. Use a checkout containing `example/`. The original `v1.0.1` tag predates this example.

---

## Requirements

### Tested example environment

| Component | Version or requirement |
|---|---|
| Platform | Linux x86_64 with Conda and Bash |
| Python | 3.10.22 |
| R / CNVPanelizer | R 4.2.3 / CNVPanelizer 1.30.0 |
| BAM tools | pysam 0.23.0, samtools 1.24, bedtools 2.31.1 |
| Example input | GRCh38 panel-sliced BAMs with matching BAI indexes |

These versions describe the new reproducibility example. The earlier Python
dependency records in `requirements*.txt` are retained separately. They do not
describe this example's environment. Full dependencies are recorded in the
Conda files linked below.

For other datasets, use coordinate-sorted, duplicate-marked BAM files aligned to
GRCh38. Every test and reference BAM needs a matching index. The validated
workflow uses BAM input. See the CRAM guidance in [Usage](#usage).

### Alignment requirements

CYP2D6 and CYP2D7 share extensive sequence similarity. Reads from these regions
can receive low mapping quality, which affects CNVPanelizer counts because the
primary workflow requires MAPQ ≥20. Use the same reference assembly and alignment
workflow for test samples and reference samples.

The benchmark evaluated ALT-aware BWA-MEM and Illumina DRAGEN. Alignment and
duplicate-marking tools are separate prerequisites for generating input BAMs.
The BWA-MEM settings were:

```bash
# Replace the FASTQ paths with your paired-end reads.
bwa mem -Y -K 100000000 GRCh38_full_analysis_set_plus_decoy_hla.fa \
    /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz > sample.sam
```

This command assumes that the BWA index and accompanying `.alt` file are
available beside the reference. Sort the resulting alignments, mark duplicates,
and create a BAM index before genotyping.

The DRAGEN benchmark used version 4.2.4 with the `hg38_graph_v4.2` reference and
ALT-aware mapping enabled. Duplicates were marked with Picard after alignment in
both benchmark workflows.

Retain duplicate-marked reads in the BAMs. The primary CNVPanelizer workflow
includes them in depth counts. The read-counting fallback and ordinary variant
counting apply different duplicate filters, as described below.

---

## Installation

From the root of a checkout containing `example/`, create a new environment on
Linux x86_64 using the tested YAML installation workflow:

```bash
bash setup_environment.sh cyripanel-example
conda activate cyripanel-example
```

The setup script reads [environment.yml](environment.yml). If this environment
already exists from the tested setup, activate it instead of creating it again.
To choose a new name, use that name in both commands.

The reconstruction test used the name `cyripanel-rebuild`. Its two-sample run
matched all 10 primary comparisons and both complete parsed JSON records.
See the [environment rebuild evidence](example/validation/environment-rebuild/evidence.json).
It reused the previously restored code and data, with the existing package
cache available.

Dependency records:

- [environment.example.yml](environment.example.yml): direct dependencies used
  for the initial environment creation.
- [Original explicit package list](reproducibility/conda-linux-64.explicit.txt):
  package builds and URLs from the original environment.
- [Rebuilt environment package list](example/validation/environment-rebuild/conda-linux-64.explicit.txt):
  package builds and URLs recorded in the successful rebuild run.

An optional alternative is to recreate the original recorded package builds
directly from their explicit list. This alternative has not been separately
tested. Use a new environment name:

```bash
conda create --name cyripanel-explicit \
  --file reproducibility/conda-linux-64.explicit.txt
conda activate cyripanel-explicit
```

Continue with the [data installation and example run](example/README.md).
The verified analysis entry point is `example/run_example.sh`.
The general `run_cyripanel.sh` wrapper and site deployment scripts remain
outside this example's verification.

---

## Two things you must provide

### 1. A panel BED file in `data/`

Keep exactly one headerless, four-column BED file (`chrom  start  end  name`) in
`data/`. Chromosome names must match those in the BAMs. The included
`PGxProbe_region_hg38.bed` is the capture panel used in the study and the verified
example.

When using another panel, work in a separate checkout and move the bundled BED
outside `data/` before adding its replacement. If multiple BED files are present,
the program warns and selects the first file it finds. The R wrapper uses
`skip=0`, so the first line must contain a target interval rather than column names.

**Your BED must contain three intervals with these exact names**, because the
copy-number step looks them up by name. The start and end values below reproduce
the corresponding fields in the bundled BED:

| Name | Chromosome | BED start | BED end | Purpose |
|---|---|---|---|---|
| `CYP2D6plusREP6_hapcn2` | chr22 | 42123192 | 42132032 | CYP2D6 plus REP6 |
| `CYP2D7` | chr22 | 42139676 | 42145745 | CYP2D7 |
| `D7spacer_hapcn1` | chr22 | 42138124 | 42139676 | CYP2D7 downstream spacer, used to resolve hybrids |

Include the capture-panel targets needed for normalization, rather than reducing
the BED to these three intervals. The assay also needs coverage of the
allele-defining and paralog-discriminating sites used by the caller. Matching
region names alone does not establish that another capture design will work.

### 2. Reference BAMs in `ref_dir/`

Between 5 and 20 samples that are **CYP2D6 diploid with no structural variant**, sequenced
with the *same* wet-lab protocol and processed with the *same* alignment pipeline as the
samples you want to genotype. Each needs its matching `.bai` index. Do not include
a test sample in its own reference panel.

In the benchmark, concordance did not differ detectably between panels of 5, 10 and 20
references. Larger panels took longer to run.

```
ref_dir/
├── SAMPLE1.bam
├── SAMPLE1.bam.bai
└── ...
```

---

## Usage

Run these commands from the repository root with the Conda environment active.
Prepare the BED and reference BAMs described above, then create the working
directories and a manifest containing one absolute test-BAM path per line:

```bash
mkdir -p manifests results
ls /absolute/path/to/samples/*.bam > manifests/my_batch.manifest
```

Replace `/absolute/path/to/samples/` with your input directory and check that the
manifest lists the intended test samples. Each BAM must have a matching index.
The sample identifier is the BAM filename without `.bam`, so filenames must be
unique within a run even when the BAMs are stored in different directories.

Call the tool directly:

```bash
python3 star_caller.py \
    --manifest manifests/my_batch.manifest \
    --genome 38 \
    --outDir results \
    --prefix my_batch \
    --threads 14
```

Alternatively, use the wrapper, which creates the output directory and saves a log:

```bash
bash run_cyripanel.sh --manifest manifests/my_batch.manifest my_batch
```

For CRAM inputs, first convert them to BAM using the matching reference FASTA,
ensure coordinate sorting, and index the resulting BAMs. End-to-end CRAM
processing has not been validated. Although `star_caller.py` accepts `-r`, this
argument is not passed to the CNVPanelizer or read-counting fallback stages.

### Output

`results/my_batch.tsv` contains one row per processed sample. Check that all
intended sample identifiers are present. Missing input paths are skipped with a
warning. For example, the two tested example samples produced:

```text
Sample      Genotype        Filter
HG00639     *1/*1           PASS
HG02129     *1/*36+*10      PASS
```

`results/my_batch.json` contains the calls and available supporting information,
including total and spacer copy number, the CNV/hybrid configuration, per-site
copy-number calls, called variants, and raw allele counts. The available fields
depend on how far processing progressed for that sample.

### Filter values

The main workflow can produce these final filter values:

| Filter | Meaning |
|---|---|
| `PASS` | A single diplotype was assigned under the implemented acceptance rules |
| `More_than_one_possible_genotype` | More than one candidate diplotype remains after the matching rules. Candidates are separated by semicolons |
| `Not_assigned_to_haplotypes` | Alleles were identified but could not be split into two haplotypes |
| `CN_Determination_Failed` | The copy-number step did not produce a usable estimate |
| `None` in TSV / `null` in JSON | No usable genotype was found after the fallback attempts |

`Invalid_Total_CN` and `CNV_Group_Not_Accepted` are statuses returned by individual
calling attempts. The main workflow can replace them during fallback processing.

Treat anything other than `PASS` as requiring review or orthogonal testing. A
`PASS` status means only that the call met the implemented acceptance rules. In
the benchmark, one sample received `PASS` in every evaluation although most of
its calls missed one `*36` copy.

---

## How it works

1. **Copy number.** CNVPanelizer applies TMM normalization and runs 10,000 bootstrap
   replicates. CyriPanel multiplies each relevant MeanRatio by two and rounds it
   to an integer. The total CYP2D6 plus CYP2D7 copy number is the sum of the two
   rounded gene estimates. The spacer copy number is obtained from the spacer
   ratio. Bootstrap confidence bounds are not used in this conversion. If
   CNVPanelizer returns an execution error or both gene ratios are below 0.3, a
   fallback estimator counts non-duplicate primary alignments at MAPQ ≥0 and
   divides sample counts by the corresponding reference means without TMM
   normalization.
2. **Structural configuration.** Patterns across 117 paralog-discriminating sites
   are interpreted as deletions, duplications, or hybrid alleles. The caller
   checks each inferred configuration for consistency with the regional
   copy-number estimates.
3. **Sequence variants.** The caller evaluates 159 variant targets in the bundled
   tables against the copy number expected under the accepted structural
   configuration. Some positions have multiple alternate variants. Ordinary
   variant counting excludes duplicate-marked reads. The workflow includes
   adjustments for high-depth noise and CYP2D7 reads that cross-map into CYP2D6.
4. **Star alleles.** The called variant set is matched against the BCyrius
   definition table. If the initial configuration fails the acceptance rules,
   alternative copy-number configurations are evaluated in a deterministic order.
   The first result with `PASS` or `More_than_one_possible_genotype` is accepted.

---

## Benchmark

CyriPanel was evaluated on 116 Human Pangenome Reference Consortium cell lines sequenced
on a laboratory-developed pharmacogenomic panel, using assembly-derived diplotypes as the
reference standard, with reference panels of 5, 10 and 20 samples drawn at random over 30
iterations each. Details, including the comparison with Aldy 4, are in the accompanying
manuscript.

The reference-standard diplotypes, the per-evaluation results of both callers and the
analysis scripts are in [`benchmark/`](benchmark/README.md), which lists the script
behind each table of the manuscript.

### Known limitations

- **`*1/*32` versus `*27/*41`** cannot be distinguished by the current caller.
  Separating them requires phasing g.42126938C>T with g.42126611C>G, which are
  327 bp apart. The study identified read pairs spanning these positions, but
  positional coverage alone does not establish reliable allele phasing. Read-pair
  phasing has not been implemented or validated, and no `*32`-carrying sample was
  available for validation.
- **`*15.003`** is not resolved as a sub-allele. `*15` itself is called through the
  standard variant path.
- Other capture designs require adequate coverage of the three copy-number
  intervals and the allele-defining and paralog-discriminating sites, together
  with assay-matched diploid reference samples. Performance on another capture
  design requires separate evaluation.

---

## Running on an HPC cluster

Run full data installation, BAM checksum verification, and analysis in an
allocated compute job. Activate the Conda environment inside the job and follow
the [example instructions](example/README.md) from the repository root.
Use the account, partition, and resource allocation required by your cluster.

The general `deployment/submit_cyripanel.slurm` template inherits the Conda
environment activated before submission. Supply your cluster's resource options
to `sbatch`, followed by the script, an absolute repository path, and the
arguments for `run_cyripanel.sh`. Its default scheduler logs are written in the
submission directory. This template is separate from the verified example.

Scripts under `benchmark/slurm/` retain site-specific configuration and need
separate updates before external use.

---

## Citation

Manuscript and software-release citation details will be added when available.

CyriPanel builds on:

- Chen X, Shen F, Gonzaludo N, et al. Cyrius: accurate CYP2D6 genotyping using whole-genome
  sequencing data. *Pharmacogenomics J.* 2021, 21:251–261.
- Halman A, Conyers R. BCyrius: an upgraded version of Cyrius for accurate CYP2D6 genotyping
  from short-read sequencing data. *Pharmacol Res Perspect.* 2025, 13:e70065.
- Oliveira C, Wolf T. CNVPanelizer: reliable CNV detection in targeted sequencing
  applications. Bioconductor.

## License

Copyright © Yu-Hui Lin (林育慧).

CyriPanel is distributed under the GNU General Public License, version 3 or later. See
[`LICENSE`](LICENSE). Source files derived from Cyrius and BCyrius retain their original
copyright notices. BCyrius is distributed under the Apache License 2.0, and a copy of that
license is in [`LICENSE-APACHE.txt`](LICENSE-APACHE.txt).

The example BAM data archive is licensed separately under CC BY 4.0. See the
[data license](example/README.md#data-license) for its scope and attribution terms.

## Contact

Yu-Hui Lin: yhlin.md05@nycu.edu.tw
