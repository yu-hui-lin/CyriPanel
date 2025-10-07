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
  | PASS | A passing, confident call |
  | More_than_one_possible_genotype | In rare cases, Cyrius reports two possible genotypes that cannot be distinguished from each other. These are different sets of star alleles that result in the same set of variants that cannot be phased with short reads (e.g., *1/*46 and *43/*45). The two possible genotypes are reported together, separated by a semicolon |
  | Not_assigned_to_haplotypes | In a very small portion of samples with more than two copies of CYP2D6, Cyrius calls a set of star alleles but they can be assigned to haplotypes in more than one way. Cyrius reports the star alleles joined by underscores (e.g., *1_*2_*68 could be *1+*68/*2, *2+*68/*1, or *1+*2/*68) |
  | LowQ_high_CN | In rare cases at high copy number (≥6 copies of CYP2D6), Cyrius uses less strict approximation in calling copy numbers to account for higher noise in depth. The genotype call could be lower confidence than usual |
  
  **Note:** A genotype of "None" indicates a no-call.

Example:
```tsv
Sample	   Genotype      Filter
NA18747	   *1/*36+*10	   PASS
```

### JSON File

File: `{prefix}.json`

Structure: Nested JSON object containing detailed calling information.

Fields for each sample:
- `Coverage_MAD` - Median absolute deviation of depth **(0.0 as CNVPanelizer skips GMM-based calculation)**
- `Median_depth` - Median sequencing depth across the region **(0.0 as CNVPanelizer skips GMM-based calculation)**
- `Total_CN` - Total copy number call of CYP2D6+CYP2D7
- `Spacer_CN` - Copy number of CYP2D7 spacer region
- `Total_CN_raw` - Raw normalized depth of CYP2D6+CYP2D7
- `Spacer_CN_raw` - Raw normalized depth of CYP2D7 spacer region
- `Variants_called` - Targeted variants called in CYP2D6
- `CNV_group` - An identifier for the sample's CNV/fusion status
- `Genotype` - Final star allele call
- `Filter` - Quality filter on the genotype call (see filter values table above)
- `Raw_star_allele` - Raw star allele call
- `Call_info` - Information about the matching confidence (e.g., unique_match)
- `Exon9_CN` - Copy number of exon 9
- `CNV_consensus` - Consensus copy number calls across regions
- `d67_snp_call` - CYP2D6 copy number call at CYP2D6/7 differentiating sites
- `d67_snp_raw` - Raw CYP2D6 copy number at CYP2D6/7 differentiating sites
- `Variant_raw_count` - Supporting reads for each variant

Example:

```json
<details>
<summary><b>Click to expand full field descriptions</b></summary>
{
    "NA18747": {
        "Coverage_MAD": 0.0,
        "Median_depth": 0.0,
        "Total_CN": 5,
        "Spacer_CN": 3,
        "Total_CN_raw": 5.0,
        "Spacer_CN_raw": 0.0,
        "Variants_called": [
            "g.42126611C>G",
            "g.42129130C>G",
            "g.42129130C>G",
            "g.42129754G>A",
            "g.42129754G>A",
            "g.42130692G>A",
            "g.42130692G>A",
            "exon9gc",
            "g.42126611C>G"
        ],
        "CNV_group": "exon9hyb",
        "Genotype": "*1/*36+*10",
        "Filter": "PASS",
        "Raw_star_allele": [
            "*1_*10.002_*36",
            "*1.011_*10_*36"
        ],
        "Call_info": "unique_star",
        "Exon9_CN": 2,
        "CNV_consensus": "2,2,3,3,3",
        "d67_snp_call": "2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3",
        "d67_snp_raw": "2.272,2.161,2.255,2.218,2.503,2.517,2.3,1.947,1.891,3.093,2.983,3.252,3.195,3.111,3.061,3.016,3.052,3.171,3.071,2.955,3.018,2.979,3.143,3.1,3.187,3.125,2.976,3.099,2.955,2.898,2.977,2.92,3.028,3.147,3.058,2.976,3.143,3.218,2.874,2.958,2.952,3.131,2.928,3.175,3.025,4.057,2.995,2.957,2.868,3.046,3.13,3.166,3.105,3.166,3.277,3.148,3.131,3.305,3.241,3.218,3.368,3.558,3.554,3.029,3.449,3.128,3.315,3.141,3.239,3.113,3.083,3.309,3.156,3.071,3.076,3.023,3.101,3.074,2.955,3.05,2.887,2.957,2.741,3.043,3.054,3.002,3.03,3.084,3.119,3.092,3.006,3.05,3.015,3.109,2.965,3.053,2.931,2.835,2.855,2.863,2.922,2.928,3.073,2.927,3.173,3.003,2.977,2.927,2.945,2.987,3.001,3.057,2.804,2.652,2.669,2.827,2.559",
        "Variant_raw_count": {
            "g.42126611C>G": "1793,1850",
            "g.42126656-42126657insAGTGGGCAC": "0,3806",
            "g.42126715C>T": "0,3681",
            "g.42126758G>A": "0,3031",
            "g.42126877G>A": "0(0:0),2337",
            "g.42126914C>G": "0,2465",
            "g.42126914C>T": "0,2465",
            "g.42126926G>A": "0,2585",
            "g.42127457C>T": "1,4480",
            "g.42127473C>T": "0(0:0),4373",
            "g.42127512C>T": "0,4251",
            "g.42127514A>G": "0,4591",
            "g.42127523A>G": "1,4420",
            "g.42127533-42127534insAC": "0,4579",
            "g.42127565T>C": "0,3667",
            "g.42127589C>T": "0,6356",
            "g.42127590G>A": "3,6428",
            "g.42127593G>C": "0,6487",
            "g.42127593G>A": "2,6485",
            "g.42127602C>T": "0,6672",
            "g.42127605T>G": "9,6407",
            "g.42127608C>T": "1,6813",
            "g.42127610T>C": "0,6415",
            "g.42127611C>T": "2,6797",
            "g.42127619T>G": "11,6984",
            "g.42127631C>G": "1,7410",
            "g.42127634C>A": "0,7408",
            "g.42127803C>T": "0,7703",
            "g.42127841C>G": "1,7181",
            "g.42127842G>A": "5,7129",
            "g.42127846delCACATCCGGATGTAGGATC": "2,7000",
            "g.42127852C>T": "0,7207",
            "g.42127853G>A": "0,7206",
            "g.42127856T>G": "30,6686",
            "g.42127899T>C": "1,7070",
            "g.42127922A>G": "0,6656",
            "g.42127938T>G": "17,6244",
            "g.42127941G>A": "7,6695",
            "g.42127963delG": "3,6486",
            "g.42127973T>C": "2(2:0),4763",
            "g.42128071C>G": "0,6356",
            "g.42128160G>T": "8,6203",
            "g.42128174delCTT": "1,5328",
            "g.42128181A>T": "3,5273",
            "g.42128185C>T": "0,6275",
            "g.42128201-42128202insA": "0,5346",
            "g.42128203delTCAG": "0,5249",
            "g.42128218-42128219insG": "57,5703",
            "g.42128212G>A": "0,6465",
            "g.42128215G>A": "0,6438",
            "g.42128217G>T": "9,6454",
            "g.42128235G>A": "0,6509",
            "g.42128242delT": "0,6052",
            "g.42128251delAGTT": "1,5922",
            "g.42128272T>G": "2,6474",
            "g.42128308C>A": "5,6461",
            "g.42128325A>G": "0(0:0),5950",
            "g.42128326delGGA": "5,5719",
            "g.42128329G>A": "0,6336",
            "g.42128332C>T": "5,6327",
            "g.42128340G>A": "9,6264",
            "g.42128351C>T": "0,6209",
            "g.42128445delCTGCC": "0,0",
            "g.42128796delG": "5,5995",
            "g.42128807C>T": "0,6424",
            "g.42128817star20complex": "0,5723",
            "g.42128812A>G": "0,5961",
            "g.42128813G>A": "1,6232",
            "g.42128817-42128818insC": "0,5783",
            "g.42128815C>T": "0,6329",
            "g.42128820-42128821insC": "0,5768",
            "g.42128848C>T": "5,6020",
            "g.42128864C>T": "0,6103",
            "g.42128876C>T": "2,6145",
            "g.42128878C>A": "3,6105",
            "g.42128879A>G": "1,5683",
            "g.42128903-42128904insTA": "4,5848",
            "g.42128926T>C": "0,5834",
            "g.42128936-42128937insGGGGCGAAAGGGGCGAAA": "0,5333",
            "g.42128936-42128937insGGGGCGAAA": "0,5333",
            "g.42128932C>T": "7,6386",
            "g.42128933G>A": "0,6440",
            "g.42128945C>T": "0,6461",
            "g.42129033C>A": "3,6335",
            "g.42129033C>T": "0,6338",
            "g.42129036A>C": "0,6018",
            "g.42129037G>T": "1,6281",
            "g.42129042T>C": "0(0:0),6040",
            "g.42129044G>A": "0,6336",
            "g.42129056C>G": "0,6567",
            "g.42129071T>A": "2,6235",
            "g.42129071T>G": "13,6224",
            "g.42129075C>T": "0,6789",
            "g.42129084delA": "0,6397",
            "g.42129087G>C": "0,6654",
            "g.42129098T>C": "3,6020",
            "g.42129113A>G": "0,6147",
            "g.42129130C>G": "4448,2339",
            "g.42129132C>T": "7,6792",
            "g.42129134G>A": "1,6739",
            "g.42129155C>T": "0,6773",
            "g.42129156A>G": "0,6388",
            "g.42129166A>C": "29,6191",
            "g.42129174C>A": "0(0:0),6718",
            "g.42129180A>T": "5(1:4),5927",
            "g.42129183C>A": "0,6648",
            "g.42129183C>T": "0,6648",
            "g.42129740T>C": "0,5150",
            "g.42129753C>G": "0,5527",
            "g.42129754G>A": "3509,1748",
            "g.42129759C>T": "0,5367",
            "g.42129765T>C": "0,4609",
            "g.42129770G>A": "0,5460",
            "g.42129770G>T": "6,5454",
            "g.42129771T>A": "5,4954",
            "g.42129779A>G": "1,5039",
            "g.42129780C>T": "0,5380",
            "g.42129799G>C": "1,5331",
            "g.42129809T>C": "0,5115",
            "g.42129819G>T": "0,5516",
            "g.42129821G>A": "0,5511",
            "g.42129827C>G": "1,5540",
            "g.42129836G>A": "2,5434",
            "g.42129880G>T": "2,4588",
            "g.42129899C>T": "3,4478",
            "g.42129906G>A": "6,4375",
            "g.42130638G>A": "0,7082",
            "g.42130655-42130656insA": "0,6739",
            "g.42130667C>T": "0,7305",
            "g.42130668C>T": "0,7351",
            "g.42130670G>A": "4,7169",
            "g.42130692G>A": "4656,2435",
            "g.42130710G>A": "6,6877",
            "g.42130715C>T": "11,6829",
            "g.42130719G>A": "3,6813",
            "g.42130728delG": "6,6314",
            "g.42130761C>T": "12,6900",
            "g.42130773C>T": "0,6852",
            "g.42130778G>A": "3,6877",
            "g.42130783-42130784insT": "0,6324",
            "g.42126578C>T": "15,10825",
            "g.42126605G>A": "0,10613",
            "g.42126634A>C": "22,10720",
            "g.42126647C>T": "1,11458",
            "g.42126681G>C": "1,11222",
            "g.42126697G>T": "14,11336",
            "g.42126719C>T": "12,11365",
            "g.42126735C>T": "0,11862",
            "g.42126746C>T": "3,11636",
            "g.42126747G>A": "5,11516",
            "g.42126749C>T": "12,11706",
            "g.42126752C>T": "3,11487",
            "g.42126896G>A": "0,11818",
            "g.42126904A>G": "2,11457",
            "g.42126956T>G": "16,12589",
            "g.42126962C>G": "1,12842",
            "g.42126982-42126983delGA": "62,11708",
            "g.42129887A>C": "19,7844",
            "g.42129910C>G": "2,8613",
            "g.42126938C>T": "0,0",
            "g.42127526C>T": "0,0",
            "g.42127556T>C": "0,0"
        }
    }
}
</details>
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
