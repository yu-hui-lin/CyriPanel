# CyriPanel
An accurate CYP2D6 genotyper for laboratory-developed targeted panels.

## Introduction
Cyrius, an open-source tool by Illumina, excels at CYP2D6 genotyping from whole-genome sequencing (WGS) data. **However, its reliance on a Gaussian Mixture Model (GMM) trained on WGS coverage patterns from 2,504 samples in the 1000 Genomes Project creates significant limitations for laboratory-developed targeted sequencing panels, whose coverage patterns differ substantially from WGS data.**

To overcome this limitation, we developed CyriPanel, an integrated pipeline that combines CNVPanelizer with a modified Cyrius framework. **CyriPanel replaces the WGS-dependent GMM with CNVPanelizer's bootstrap-based approach, generating synthetic reference sets directly from panel data. We further optimized variant calling algorithms for high-depth sequencing noise and implemented fallback mechanisms to resolve ambiguous CNV or SV calls.** CyriPanel enables accurate CYP2D6 genotyping from targeted sequencing without requiring extensive reference cohorts, making comprehensive pharmacogenomic testing more accessible for clinical laboratories using custom panels.

---

## Key Modifications
### 1.  Copy Number Determination
CyriPanel replaces WGS-dependent GMM with CNVPanelizer, which uses non-parametric bootstrap subsampling to calculate normalized read count ratios without requiring extensive reference cohorts, providing accurate CYP2D6 copy number estimates for targeted sequencing panels.

### 2.  Algorithm Optimization
- To accommodate targeted sequencing variability, the strand-bias filters are relaxed for known noisy variants, and position-specific read thresholds are implemented.
-  A dual-validation system compares SNP-based and spacer-based copy number estimates to ensure high-confidence calls of structural variants.
-  For ambiguous scenarios, the algorithm preferentially assigns CN=2 over CN=3 when confidence is below the threshold.
-  To prevent false positives, CN=1 calling is enhanced with allele fraction validation and increased error parameters.

### 3.  Fallback Mechanisms
When initial estimations fail to produce high-confidence genotypes, the system automatically tests other copy number configurations, prioritized by proximity to initial estimates, ensuring robust genotyping in challenging cases.

---

## Installation Guide
CyriPanel requires both Python (≥3.6) and R (≥4.0) environments with specific packages.

### Step 1: Clone Repository
```
git clone https://github.com/yourusername/CyriPanel.git
cd CyriPanel
```

### Step 2: Install Python Dependencies
```
# Using pip
pip install -r requirements.txt

# Or using conda
conda env create -f environment.yml
conda activate cyripanel
```

### Step 3: Install R and CNVPanelizer
Prerequisites:
- R version ≥4.0.0 (check with `R --version`)
- Bioconductor package manager

```
# Install R (if not already installed)
# On Ubuntu/Debian:
sudo apt-get install r-base r-base-dev

# On macOS:
brew install r

# Install Bioconductor package manager
Rscript -e "install.packages('BiocManager')"

# Install CNVPanelizer
Rscript -e "BiocManager::install('CNVPanelizer')"
```

### Step 4: Verify Installation
```
# Test Python dependencies
python3 -c "import numpy, scipy, pysam, statsmodels, pandas; print('✓ Python packages OK')"

# Test R dependencies
Rscript -e "library(CNVPanelizer); cat('✓ R packages OK\n')"

# Run test suite (optional)
pytest tests/ -v
```

### Step 5: Prepare BED file defining targeted panel regions 
CyriPanel requires a BED file defining the genomic regions captured by your panel:
1. Create a BED file in 0-based (standard BED format), **GRCh38 coordinates with 4-column format**.
   Columns: chromosome, start, end, region name
   Example BED file format (PGxProbe_region_hg38.bed):
   ```
   chr22	42121996	42155994	CYP2D6/CYP2D7/CYP2D8atReverseStrand
   chr22	42123192	42132032	CYP2D6plusREP6_hapcn2
   chr22	42135344	42138124	REP7
   chr22	42138124	42139676	D7spacer_hapcn1
   chr22	42139676	42145745	CYP2D7
   ```
   
2. The BED file **must contain other genomic regions beyond CYP2D6 locus that are captured by your panel**. This is essential for CNVPanelizer normalization across samples.
3. **Place the BED file in the /data directory**:
   ```
   CyriPanel/
   ├── star_caller.py
   ├── depth_calling/
   ├── caller/
   ├── ref_dir
   └── data/
       ├── PGxProbe_region_hg38.bed
       ├── star_table.txt
       ├── CYP2D6_SNP_38.txt
       ├── CYP2D6_haplotype_38.txt
       ├── CYP2D6_target_variant_38.txt
       └── CYP2D6_target_variant_homology_region_38.txt
   ```

### Step 6: Prepare Reference Panel
CyriPanel requires a reference panel of diploid samples (CN=2 at CYP2D6):
1. **Collect reference BAM files from samples with known diploid CYP2D6.**
   **- 20+ samples recommended**
   **- Same sequencing platform and capture kit as test samples**
3. **Ensure all BAM files have index files(.bai).** 
4. **Place these BAM and BAI files in the /ref_dir directory:**
   ```
   CyriPanel/
   ├── data/
   ├── star_caller.py
   ├── depth_calling/
   ├── caller/
   └── ref_dir/
       ├── sample01.bam
       ├── sample01.bam.bai
       ├── sample02.bam
       ├── sample02.bam.bai
       └── ...
   ```
   ```
   # Organize reference BAMs (20+ samples with known CN=2 at CYP2D6)
   mkdir reference_panel
   # Copy your reference BAM files here
   cp /path/to/reference*.bam reference_panel/
   cp /path/to/reference*.bam.bai reference_panel/
   ```


---

## Running the program
This Python3 program can be run as follows:
```
python3 star_caller.py \
    --manifest sample01.manifest \        # Path to manifest file
    --genome 38 \                         # Reference genome (only GRCh38 supported)
    --prefix sample01 \                   # Output file prefix
    --outDir results/ \                   # Output directory
    --threads 8                           # Number of threads
```

1. CyriPanel only supports **GRCh38**. If --genome input is not 38, CyriPanel will terminate with a warning.
2. The manifest is a text file in which each line should list the absolute path to an input BAM file.
```
for bam in /BAM_dir/*.bam; do 
    echo $bam > /manifest_dir/$(basename $bam .bam).manifest
done
```

---

## Output Format
### TSV File

File: `{prefix}.tsv`

Columns:
- `Sample` - Sample identifier
- `Genotype` - Final star allele call (e.g., *1/*2, *2x2/*4)
- `Filter` - Filters on the genotype call
  | Filter Value | Description |
  |-------------|-------------|
  | `PASS` | A passing, confident call |
  | `More_than_one_possible_genotype` | In rare cases, Cyrius reports two possible genotypes that cannot be distinguished from each other. These are different sets of star alleles that result in the same set of variants that cannot be phased with short reads (e.g., *1/*46 and *43/*45). The two possible genotypes are reported together, separated by a semicolon |
  | `Not_assigned_to_haplotypes` | In a very small portion of samples with more than two copies of CYP2D6, Cyrius calls a set of star alleles but they can be assigned to haplotypes in more than one way. Cyrius reports the star alleles joined by underscores (e.g., *1_*2_*68 could be *1+*68/*2, *2+*68/*1, or *1+*2/*68) |
  | `LowQ_high_CN` | In rare cases at high copy number (≥6 copies of CYP2D6), Cyrius uses less strict approximation in calling copy numbers to account for higher noise in depth. The genotype call could be lower confidence than usual |
  **Note:** A genotype of "None" indicates a no-call.

Example:
```tsv
Sample	   Genotype      Filter
NA18747	   *1/*36+*10	   PASS
```

### JSON File
File: `{prefix}.json`

Columns:
- `Sample` - Sample identifier
- `Genotype` - Final star allele call (e.g., *1/*2, *2x2/*4)
- `Filter` - PASS or reason for failure

Example:
```tsv
Sample	   Genotype      Filter
HG00544	   *1/*36+*10	   PASS
```







---

## Troubleshooting

### Issue 1: R Package Installation Fails
BiocManager::install() fails with compilation errors

Solutions:
```
# On Ubuntu/Debian
sudo apt-get install build-essential r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev

# On macOS
xcode-select --install
brew install libxml2
```

### Issue 2: pysam Installation Fails
Error compiling pysam C extensions

**Solutions:**
```
# On Ubuntu/Debian
sudo apt-get install libbz2-dev liblzma-dev

# On macOS
brew install xz bzip2

# Then retry pip install
pip install --upgrade pysam
```

### Issue 3: CNVPanelizer Not Found in Subprocess
**Symptoms:** Python script can't find Rscript or CNVPanelizer

**Solution:**
```
# Verify R is in PATH
which Rscript
R --version

# Set R library path explicitly (add to ~/.bashrc or ~/.zshrc)
export R_LIBS_USER="$HOME/R/library"

# Or specify full Rscript path in panel_cn.py
```

### Issue 4: Personal Library Path Issues
**Symptoms:** Permission denied when installing R packages

**Solution:** The run_CNVPanelizer.R script automatically handles this by setting:
```
personal_lib_path <- file.path(Sys.getenv("HOME"), "R", "library")
.libPaths(c(personal_lib_path, .libPaths()))
```

---

## Quick Start Verification Script

Save as `test_installation.sh`:
```
#!/bin/bash
echo "Testing CyriPanel Installation..."
echo ""

echo "1. Testing Python environment..."
python3 -c "import numpy, scipy, pysam, statsmodels, pandas; print('✓ Python packages OK')"

echo "2. Testing R environment..."
Rscript -e 'library(CNVPanelizer); cat("✓ R packages OK\n")'

echo "3. Testing Rscript in PATH..."
which Rscript && echo "✓ Rscript found in PATH"

echo ""
echo "Installation verification complete!"
```

Run with: `bash test_installation.sh`

---

## Support and Contact
Yu-Hui Lin, M.D.
yhlin.md05@nycu.edu.tw
