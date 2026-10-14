# CyriPanel

CYP2D6 star-allele genotyping from laboratory-developed targeted sequencing panels.

⚠️ **Disclaimer: Research Use Only (RUO)**
CyriPanel is a computational research tool developed for methodological evaluation of targeted capture sequencing data. It is **not** intended for clinical diagnostics, patient management, or prescribing decisions.

CyriPanel is derived from [Cyrius](https://github.com/Illumina/Cyrius) (Illumina) and
incorporates the updated star-allele definitions of
[BCyrius](https://github.com/Halman-Lab/BCyrius). It keeps Cyrius's paralogue-aware
variant calling and star-allele matching, but replaces the copy-number step: instead of a
Gaussian mixture model fitted to whole-genome coverage, copy number is estimated with
[CNVPanelizer](https://bioconductor.org/packages/CNVPanelizer/) against a small cohort of
known-diploid samples sequenced on the *same* panel. That makes CYP2D6 genotyping possible
on capture data, where WGS-derived coverage expectations do not transfer.

CyriPanel resolves single-nucleotide and indel variants, copy-number changes, and
CYP2D6–CYP2D7 hybrid alleles, and reports a diplotype with a filter status.

**GRCh38 only.**

---

## Requirements

|        |                                                              |
| ------ | ------------------------------------------------------------ |
| Python | 3.6.8 (the validated interpreter). Install with `pip install -r requirements.txt`; the exact versions used for the benchmark in the accompanying preprint are in `requirements-validated.txt` |
| R      | ≥ 4.2 with [CNVPanelizer](https://bioconductor.org/packages/CNVPanelizer/) (validated with v1.30.0, Bioconductor 3.16) |
| Other  | `samtools`; an aligner (see below)                           |
| Input  | Coordinate-sorted, duplicate-marked BAM/CRAM aligned to GRCh38 |

### Alignment matters more than usual

CYP2D6 and CYP2D7 are ~94 % identical, so a read that maps equally well to both is given
mapping quality 0. With a plain BWA-MEM alignment, **80–90 % of reads across the CYP2D
locus carry MAPQ 0**, and because CNVPanelizer counts only reads at MAPQ ≥ 20 the
copy-number step is then starved of data at exactly the locus of interest.

Two alignment routes were evaluated in the benchmark:

```bash
# ALT-aware BWA-MEM
bwa mem -Y -K 100000000 GRCh38_full_analysis_set_plus_decoy_hla.fa ...
# the accompanying .alt file must be present beside the reference
```

or **Illumina DRAGEN**, which retains > 85 % of reads at MAPQ ≥ 20 natively and needs no
ALT-aware workaround.

Duplicates should be marked but **not removed** — capture duplicates carry depth signal.

---

## Installation

```bash
git clone https://github.com/yu-hui-lin/CyriPanel.git
cd CyriPanel

# Python dependencies
pip install -r requirements.txt

# R dependency (pinned to the validated release)
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("CNVPanelizer", version="3.16", update=FALSE)'
```

`setup_environment.sh` automates the same steps into a conda environment and is convenient
on an HPC where `module load miniconda3` is available.

---

## Two things you must provide

### 1. A panel BED file in `data/`

A standard 4-column BED (`chrom  start  end  name`) describing your capture targets.
CyriPanel picks up the single `.bed` file it finds in `data/`; the panel used in the
benchmark, `PGxProbe_region_hg38.bed`, is included as a reference.

**Your BED must contain three intervals with these exact names**, because the copy-number
step looks them up by name:

| Name                    | GRCh38 interval             | Purpose                                           |
| ----------------------- | --------------------------- | ------------------------------------------------- |
| `CYP2D6plusREP6_hapcn2` | chr22:42,123,192-42,132,032 | CYP2D6 plus REP6                                  |
| `CYP2D7`                | chr22:42,139,676-42,145,745 | CYP2D7                                            |
| `D7spacer_hapcn1`       | chr22:42,138,124-42,139,676 | CYP2D7 downstream spacer, used to resolve hybrids |

A panel that does not capture these regions cannot be genotyped by CyriPanel. If yours
covers them under different names, rename them in your BED.

### 2. Reference BAMs in `ref_dir/`

Between 5 and 20 samples that are **CYP2D6 diploid with no structural variant**, sequenced
with the *same* wet-lab protocol and processed with the *same* alignment pipeline as the
samples you want to genotype. Each needs its `.bai` index.

In the benchmark, concordance did not differ detectably between panels of 5, 10 and 20
references. Larger panels only lengthened runtime.

```
ref_dir/
├── SAMPLE1.bam
├── SAMPLE1.bam.bai
└── ...
```

---

## Usage

Create a manifest listing one BAM path per line:

```bash
ls /path/to/samples/*.bam > manifests/my_batch.manifest
```

Then either call the tool directly:

```bash
python3 star_caller.py \
    --manifest manifests/my_batch.manifest \
    --genome 38 \
    --outDir results \
    --prefix my_batch \
    --threads 14
```

or use the wrapper, which creates `results/` and handles logging:

```bash
bash run_cyripanel.sh --manifest manifests/my_batch.manifest my_batch
```

For CRAM input, add `-r /path/to/GRCh38.fa`.

### Output

`results/my_batch.tsv` — one row per sample:

```
Sample      Genotype        Filter
HG02129     *1/*10+*36      PASS
HG00436     *2/*41          PASS
```

`results/my_batch.json` — the same calls with the evidence behind them: total and spacer
copy number, the CNV/hybrid configuration, per-site copy-number calls, called variants and
raw allele counts.

### Filter values

| Filter                                         | Meaning                                                      |
| ---------------------------------------------- | ------------------------------------------------------------ |
| `PASS`                                         | A single diplotype was assigned                              |
| `More_than_one_possible_genotype`              | Several diplotypes explain the data equally well; all are listed, separated by `;` |
| `Not_assigned_to_haplotypes`                   | Alleles were identified but could not be split into two haplotypes |
| `CNV_Group_Not_Accepted`                       | No consistent structural configuration could be inferred     |
| `Invalid_Total_CN` / `CN_Determination_Failed` | The copy-number step did not produce a usable estimate       |

Treat anything other than `PASS` as requiring review or orthogonal testing. A `PASS`
status means only that the call met the implemented acceptance rules. In the benchmark,
one sample received `PASS` in every evaluation although most of its calls missed one
`*36` copy.

---

## How it works

1. **Copy number.** CNVPanelizer normalises read counts (TMM) against the reference panel
   and bootstraps 10,000 replicates to obtain the total CYP2D6 + CYP2D7 copy number and the
   CYP2D7 spacer copy number. If normalisation fails — both ratios implausibly low, which
   happens when MAPQ filtering has starved it — a fallback estimator counts reads directly
   at MAPQ ≥ 0.
2. **Structural configuration.** Copy number is called at 117 paralogous CYP2D6/CYP2D7
   sites; transitions between segments are translated into deletions, duplications and
   hybrid alleles, and accepted only if they reconstruct the observed total copy number.
3. **Sequence variants.** 159 target sites are called against the copy number expected at
   each site under that configuration, with panel-specific adjustments for high-depth noise
   and for CYP2D7 reads that cross-map into CYP2D6.
4. **Star alleles.** The called variant set is matched against the BCyrius definition table.
   If no configuration yields a consistent diplotype, alternative copy-number
   configurations are tried in a fixed, deterministic order.

---

## Benchmark

CyriPanel was evaluated on 116 Human Pangenome Reference Consortium cell lines sequenced
on a laboratory-developed pharmacogenomic panel, using assembly-derived diplotypes as the
reference standard, with reference panels of 5, 10 and 20 samples drawn at random over 30
iterations each. Details, including the comparison with Aldy 4, are in the accompanying
preprint.

The reference-standard diplotypes, the per-evaluation results of both callers and the
analysis scripts are in [`benchmark/`](benchmark/README.md), which lists the script
behind each table of the preprint.

### Known limitations

- **`*1/*32` versus `*27/*41`** cannot be distinguished. Separating them requires phasing
  g.42126938C>T against g.42126611C>G, 327 bp apart — beyond read length. Paired-end
  fragments covering both positions are plentiful (median 120 per sample, range 6-434 across the 116 evaluated samples; none had fewer than six), so read-pair
  phasing is feasible, but it is not implemented and no `*32`-carrying sample was available
  to validate it against.
- **`*15.003`** is not resolved as a sub-allele. `*15` itself is called through the standard
  variant path.
- Requires a panel that captures the three intervals listed above.

---

## Citation

If you use CyriPanel, please cite the accompanying preprint:

- Lin Y-H, Chang T-Y, Huang K-T, et al. CyriPanel for *CYP2D6* genotyping from targeted
  capture panels using assay-matched diploid references. *bioRxiv* 2026.
  doi: (added when the preprint is posted)

CyriPanel builds on:

- Chen X, Shen F, Gonzaludo N, et al. Cyrius: accurate CYP2D6 genotyping using whole-genome
  sequencing data. *Pharmacogenomics J.* 2021;21:251–261.
- Halman A, Conyers R. BCyrius: an upgraded version of Cyrius for accurate CYP2D6 genotyping
  from short-read sequencing data. *Pharmacol Res Perspect.* 2025;13:e70065.
- Oliveira C, Wolf T. CNVPanelizer: reliable CNV detection in targeted sequencing
  applications. Bioconductor.

## License

See [`LICENSE`](LICENSE). Source files derived from Cyrius and BCyrius retain their
original copyright notices.

## Contact

Yu-Hui Lin — yhlin.md05@nycu.edu.tw
