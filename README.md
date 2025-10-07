# CyriPanel
CYP2D6 genotyper for targeted panels with modified Cyrius and CNVPanelizer

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
1. BED file in GRCh38 and 4-column format
```
chr22	42121996	42155994	CYP2D6/CYP2D7/CYP2D8atReverseStrand
      chr22	42123192	42132032	CYP2D6plusREP6_hapcn2
      chr22	42135344	42138124	REP7
      chr22	42138124	42139676	D7spacer_hapcn1
      chr22	42139676	42145745	CYP2D7
```
3. Place the 

### Step 6: Prepare Reference Panel
CyriPanel requires a reference panel of diploid samples (CN=2 at CYP2D6):
1. **Collect reference BAM files from samples with known diploid CYP2D6**
   **- 20+ samples recommended**
   **- Same sequencing platform and capture kit as test samples**
3. **Ensure all BAM files have index files** (.bai)
4. **Place these BAM and BAI files in the ref_dir:**
  ```
   CyriPanel/
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

## Usage

```
python3 star_caller.py \
    --manifest sample01.manifest \
    --genome 38 \
    --prefix output_prefix \
    --outDir results/ \
    --threads 4
```

---

## Input Files


---

## Output Format









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
