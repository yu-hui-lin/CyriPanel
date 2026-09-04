# CyriPanel - CYP2D6 Genotyper for Targeted Sequencing Panels

##  國網位址 (Location)

```bash
/staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/CyriPanel
```

## 快速開始 (Quick Start)

### 0. 放置 5-20 個已知為 diploid、無 CNV/SV 的檢體至
```bash
/staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel/ref_dir
```

-  **！！！ BAM 檔需由 Illumina Dragen 作為 aligner 生成 ！！！**
- 須為與待檢檢體相同 wet lab 流程產出
- 每個 BAM 檔需有對應的 .bai 索引檔

預設為 NA21102, HG01940, HG00639, HG00350, HG00232 五個檢體

### 1. 啟動環境

```bash
# 載入miniconda
module load biology
module load miniconda3

# 載入並啟動 CyriPanel 環境
cd /staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel
source activate_cyripanel.sh
```

### 2. 驗證安裝 (Validation Test)

執行測試樣本 HG02129 以確認安裝正確：

```bash
cd /staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel
bash run_test_HG02129.sh
```

預期結果：

```bash
[INFO] Called Genotype: *1/*36+*10
[INFO] Filter Status: PASS
[INFO] Expected Genotype: *1/*36+*10
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   ✓ TEST PASSED - EXACT MATCH                               ║
║                                                              ║
║   CyriPanel is correctly installed and configured!          ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
```

### 3. 建立 .manifest 

假設要分析 /sample_directory 資料夾內所有 .bam

```bash
# 自動列出所有 BAM 檔案
cd /staging/reserve/paylong_ntu/Clinical/__Clinical_pipeline/PG_pipeline/software/CyriPanel
ls /sample_directory/*.bam > manifests/my_samples.manifest

# 檢查內容
cat manifests/my_samples.manifest
```

### 4. 執行批次分析 

##### 方式 A：直接在登入節點執行（小量批次樣本）

```bash
bash run_cyripanel.sh --manifest manifests/my_samples.manifest batch_name
```

batch_name: 就是一個輸出檔案的前綴名稱，你可以自己取名，例如：PG_batch20260129
輸出檔案會變成：results/PG_batch20260129_cyp2d6.tsv

##### 方式 B：提交 SLURM 工作 (把工作提交到計算節點，不會佔用登入節點資源)

```bash
sbatch submit_cyripanel.slurm manifests/my_samples.manifest batch_name
```

batch_name: 就是一個輸出檔案的前綴名稱，你可以自己取名，例如：PG_batch20260129
輸出檔案會變成：results/PG_batch20260129_cyp2d6.tsv

## 輸出檔案 (Output Files)

| 檔案              | 說明                             |
| --------------- | ------------------------------ |
| `{prefix}.tsv`  | 簡化結果（Sample, Genotype, Filter） |
| `{prefix}.json` | 完整結果（包含 CNV 資訊、變異位點等）          |

### TSV 輸出範例

```
Sample      Genotype        Filter
HG02129     *1/*10+*36      PASS
HG00436     *2/*41          PASS
```

## 目錄結構 (Directory Structure)

```
CyriPanel/
├── star_caller.py           # 主程式
├── run_cyripanel.sh         # 執行腳本
├── run_test_HG02129.sh      # 測試腳本
├── activate_cyripanel.sh    # 環境啟動腳本
├── submit_cyripanel.slurm   # SLURM 工作腳本
├── caller/                  # Star allele 分析模組
├── depth_calling/           # CNV 分析模組
├── data/                    # 參考資料檔
├── ref_dir/                 # CNVPanelizer 參考樣本 (20 個 diploid BAMs)
├── manifests/               # .manifest files
├── test_data/               # 測試樣本
├── results/                 # 輸出結果
└── logs/                    # 日誌檔案
```

---
Copyright © Yu-Hui Lin (林育慧) yhlin.md05@nycu.edu.tw
Github: https://github.com/yu-hui-lin/CyriPanel/tree/main
*CyriPanel is licensed under GPL-3.0. See LICENSE file for details.*