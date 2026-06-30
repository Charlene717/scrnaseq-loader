# 01d_QC_H5.R — 10x .h5 (HDF5) 的 per-sample QC

`01` 系列第四支。處理 **10x Cell Ranger 打包的 `.h5` 檔**
(`filtered_feature_bc_matrix.h5`),每個 `.h5` 各自跑與 `01` 完全相同的 QC。

格式家族到此補齊:

| 腳本 | 吃什麼格式 |
|---|---|
| `01_QC_PerSample.R` | 10x 三件套 (barcodes/features/matrix) |
| `01b_QC_DenseMatrix.R` | dense matrix 文字檔 (txt/tsv/csv) |
| `01c_QC_FromSeuratRDS.R` | 已整合的 Seurat RDS (按 sample 拆) |
| **`01d_QC_H5.R`** | **10x .h5 (HDF5)** |

四支輸出格式一致,`02_Aggregate_Summary.R` 可一起彙整。

## 為什麼 .h5 最單純

`.h5` 是把 barcodes/features/matrix 三件套打包成單一 HDF5 檔。Seurat 內建
`Read10X_h5()` 直接讀,**不必像 `01b` 那樣猜分隔符/方向/前綴**。所以這支最乾淨。

唯一要自動判斷的:`Read10X_h5()` 可能回傳兩種東西

- **單一矩陣** → 純單樣本(最常見,你截圖的 `filtered_feature_bc_matrix.h5` 就是)
- **一個 list** → 檔案含多個 feature type(如 `Gene Expression` + `Antibody Capture`/CITE-seq,
  或多 genome)→ 自動取 `Gene Expression` 那塊做 scRNA QC

腳本自動偵測回傳型別並處理。**這兩種情況我都用模擬的 10x v3 .h5 驗證過**:單樣本檔正確讀成
單一矩陣;多 feature type 檔正確抽出 `Gene Expression`(含 `MT-` 基因供物種偵測),矩陣值
能正確還原成 gene × cell。

## 「一個 .h5 = 一個樣本」嗎?自動判斷

你不確定,所以腳本這樣處理:**預設一個 `.h5` 當一個樣本**(用檔名/GSM 當樣本名,
如截圖的 `RUN_1` / `RUN_2` 各一檔 = 兩個樣本)。但如果某個 `.h5` 內部其實含多 feature
type,腳本會偵測到並只取 `Gene Expression`,不會把抗體那類混進 scRNA QC。

## 前置需求:hdf5r 套件

讀 `.h5` 需要 `hdf5r`,沒裝會直接報錯提示:

```r
install.packages("hdf5r")
```

## 用法

```r
# 1. 改路徑(腳本上方)
SRC_ROOT <- "Z:/Dataset_Online/scRNA-seq/h5"    # 放 .h5 的根目錄
OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"

# 2. source 執行
source("01d_QC_H5.R")
```

資料夾結構(和 `01` 一樣 `<GSE>/` 之下放檔):

```
h5/
└── ##GSE166950/
    ├── GSM5089673_RUN_1_filtered_feature_bc_matrix.h5   ← 樣本 1
    └── GSM5089674_RUN_2_filtered_feature_bc_matrix.h5   ← 樣本 2
```

## 輸出(與 01/01b/01c 完全相同)

```
#Keloid_QC/
├── <GSE>/<sample>/
│   ├── <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv
│   └── plots/  (violin pre/post, scatter pre, before/after)
└── _run_log/
    ├── per_sample/<GSE>__<sample>.csv      ← 02 彙整用,欄位與 01 一致
    ├── structure_scan_h5.csv               ← 每個 .h5 讀取狀態 + 偵測到的型別
    └── run_log_h5.txt
```

樣本名會自動把 `_filtered_feature_bc_matrix` 這類後綴去掉
(`GSM5089673_RUN_1_filtered_feature_bc_matrix` → `GSM5089673_RUN_1`)。

## 可調參數(腳本上方)

和 `01` 對齊:`MIN_FEATURE` / `MAX_FEATURE` / `MAX_MT_PCT` / `RUN_DOUBLETFINDER` / `IGNORE_FOLDERS`。

`.h5` 專屬:`GEX_FEATURE_TYPE`(多 feature type 時取哪個,預設 `"Gene Expression"`)。

## 小提醒

- `structure_scan_h5.csv` 會記錄每個 `.h5` 讀到的型別(`single matrix` 或
  `list[...] -> 取 'Gene Expression'`),跑完掃一眼就知道有沒有 CITE-seq 之類混進來。
- 這支處理「已經是 count 的 `.h5`」(filtered 或 raw matrix 都可)。若是 `.h5ad`
  (AnnData / scanpy 的格式,不是 10x HDF5),那是**不同格式**,要另外用 `SeuratDisk` 或
  `anndata` 套件轉,不能用 `Read10X_h5()`。你倉庫規劃裡有列 h5ad,若之後遇到再做一支 `01e`。
