# 01b_QC_DenseMatrix.R — Dense matrix 文字檔的 QC

`01_QC_PerSample.R` 的姊妹版。處理 GEO 上那些**不是 10x 三件套、而是 dense matrix 純文字檔**
(gene × cell 大表)的資料集,做和 `01` **完全相同**的 QC。

## 為什麼需要它

`01` 只認 10x 的 `barcodes` / `features` / `matrix`。但很多 GEO 資料集給的是這種:

```
            AAACCTGAGC   AAACGGGAGT   ...
hg19_SAMD11      0            5
hg19_NOC2L       1            0
...
```

格式還每家不一樣(tab/逗號、基因有沒有前綴、一檔一樣本還是全部合併、混進 RPKM)。
這支腳本把這些雜亂格式**自動讀進來、統一處理**,然後套用 `01` 的 QC 邏輯。

## 設計原則:換讀取層,QC 不變

```
各種 dense matrix (txt/tsv/csv)
        │
        ▼  自動偵測(分隔符 / 方向 / 前綴 / RPKM / 合併拆樣本)
   統一成 gene × cell sparse matrix
        │
        ▼  ← 與 01 完全相同的 QC
   門檻 200~5000、MT<30、物種偵測、pre/post RDS、stats、plots、per-sample row
```

輸出格式和 `01` 一模一樣,所以 `02_Aggregate_Summary.R` 可以同時吃 `01` 和 `01b` 的結果。

## 自動偵測能力(全部經測試驗證)

| 偵測項 | 說明 | 例子 |
|---|---|---|
| **分隔符** | Tab 還是逗號 | `.tsv`/`.txt` → tab；`.csv` → 逗號 |
| **矩陣方向** | gene×cell 或 cell×gene(後者自動轉置) | 用「哪一維像 barcode/gene」判斷 |
| **基因名前綴** | 偵測並剝除 | `hg19_SAMD11` → `SAMD11` |
| **非 raw counts** | 檔名含這些字就跳過 | `rpkm`/`fpkm`/`tpm`/`cpm`/`normalized` |
| **合併大檔** | barcode 帶樣本前綴 → 依前綴拆成多樣本各自 QC | `KL14_ACGT...` / `NS5_ACGT...` → 拆成 2 樣本 |

> 偵測邏輯在交付前已用 7 種真實格式案例(含 GSE129611/137897/155816/165816/191067 + 轉置 + RPKM)
> 驗證全部正確,特別確認「純 barcode 單樣本檔不會被誤拆」這個關鍵邊界。

## 對應你目前手上的 5 個 GSE

| GSE | 這支腳本怎麼處理 |
|---|---|
| GSE129611 | tab、gene×cell、自動去 `hg19_` 前綴 |
| GSE137897 | 只讀 `_read.txt`(counts),**自動跳過 `_rpkm.txt`** |
| GSE155816 | tab、多個 GSM 各自一檔,逐一處理 |
| GSE165816 | 逗號 csv、55 個 GSM 逐一處理 |
| GSE191067 | 單一合併大檔 → **依 barcode 前綴拆樣本**後各自 QC |

## 用法

```r
# 1. 改路徑(腳本最上方)
SRC_ROOT <- "Z:/Dataset_Online/scRNA-seq/##_Keloid"   # 放 dense matrix 的根目錄
OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"  # 與 01 共用輸出

# 2. 直接 source 執行
source("01b_QC_DenseMatrix.R")
```

資料夾結構(和 `01` 一樣是 `<GSE>/` 之下放檔):

```
##_Keloid/
├── #GSE129611/
│   ├── GSM..._10X1data.tsv
│   └── GSM..._comb_clean.dge.txt
├── #GSE137897/
│   ├── GSE137897_read.txt      ← 處理這個
│   └── GSE137897_rpkm.txt      ← 自動跳過
└── #GSE191067/
    └── GSE191067_all.UMI.matrix.csv   ← 自動拆樣本
```

## 建議先裝 data.table

大檔(像 GSE191067 那個 412 MB)用 base R 的 `read.table` 會很慢。裝了 `data.table`
會自動改用 `fread`,快很多:

```r
install.packages("data.table")
```

沒裝也能跑,只是慢。

## 輸出

和 `01` 完全相同:

```
#Keloid_QC/
├── <GSE>/<sample>/
│   ├── <sample>_pre_QC.rds
│   ├── <sample>_post_QC.rds
│   ├── <sample>_QC_stats.csv
│   └── plots/
└── _run_log/
    ├── per_sample/<GSE>__<sample>.csv      ← 02 腳本彙整用(欄位與 01 一致)
    ├── structure_scan_densematrix.csv      ← 哪些檔讀了/跳過/失敗
    └── run_log_densematrix.txt
```

## 可調參數(腳本上方)

和 `01` 對齊的:`MIN_FEATURE` / `MAX_FEATURE` / `MAX_MT_PCT` / `RUN_DOUBLETFINDER` / `IGNORE_FOLDERS`。

dense matrix 專屬的:
- `SKIP_PATTERNS` — 視為正規化檔要跳過的關鍵字
- `GENE_PREFIXES` — 要偵測剝除的基因名前綴
- `MIN_CELLS_PER_SPLIT` — 合併檔中,一個前綴至少幾個 barcode 才算一個樣本(預設 20,濾雜訊)

## 注意事項

- **方向判斷靠啟發式**:絕大多數情況正確,但若遇到極不典型的矩陣(基因名和 barcode 都很怪),
  可能誤判。腳本的 `structure_scan_densematrix.csv` 會記錄每個檔判到的方向(`orient=...`),
  跑完掃一眼確認即可。
- **condition 自動推測**:從樣本名/前綴的 `KL`/`Keloid`/`NS`/`Normal` 等關鍵字猜,猜不到記 NA。
  日後可在 02 階段用對照表補正。
- 這支處理「**已經是 count 表**」的資料;若 GSE 只提供 raw FASTQ/SRA,那要走另一條
  `sra_to_cellranger.sh` → Cell Ranger 的路。
