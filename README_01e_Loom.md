# 01e_QC_Loom.R — .loom (HDF5-based) 的 per-sample QC

`01` 系列第五支。處理 **`.loom` 檔**(含 gzip 壓縮的 `.loom.gz`),每個 `.loom` 各自跑
與 `01` 完全相同的 QC。

格式家族現在覆蓋:

| 腳本 | 格式 |
|---|---|
| `01_QC_PerSample.R` | 10x 三件套 |
| `01b_QC_DenseMatrix.R` | dense matrix 文字檔 |
| `01c_QC_FromSeuratRDS.R` | 已整合 Seurat RDS |
| `01d_QC_H5.R` | 10x .h5 (HDF5) |
| **`01e_QC_Loom.R`** | **.loom / .loom.gz** |

五支輸出格式一致,`02_Aggregate_Summary.R` 全部吃得下。

## .loom 的兩個坑(都自動處理)

### 1. 常以 gzip 壓縮(.loom.gz)
你截圖的 `GSM4647785_P2_1.loom.gz` 就是——WinRAR 顯示裡面才是 `.loom`。腳本:**有 `.gz`
就先解壓到暫存再讀,讀完刪暫存;沒有 `.gz` 就直讀**。用 R 內建 `gzfile` 串流解壓,不需
額外裝工具。

### 2. 內部欄位命名 / 矩陣方向不統一 ← 這是 .loom 真正的坑
`.loom` 是 HDF5,結構是 `/matrix` + `/row_attrs` + `/col_attrs`,但**不同工具產生的 loom,
gene/cell 名的 attribute 叫法都不同**:

| 工具 | gene 名欄位 | cell 名欄位 |
|---|---|---|
| loompy 標準 | `Gene` | `CellID` |
| velocyto | `Gene`(+ spliced/unspliced layers) | `CellID` |
| scanpy 匯出 | `var_names` | `obs_names` |
| 其他 | `gene_names` | `barcode` |

而且矩陣**理論上** gene×cell,但偶有工具存成 cell×gene(要轉置)。

腳本的做法:**以「內容」判斷方向**——不管 attribute 叫什麼,先看哪一維的值像 gene symbol、
哪一維像 barcode,用內容決定 gene/cell 維度與是否轉置;attribute 名稱只當內容不明確時的
後備。**這套邏輯用 5 種工具風格的 loom 驗證過**(loompy / velocyto / scanpy / 別名命名 /
轉置),全部正確,特別是轉置的 loom 也能正確判出來。

## 用法

```r
# 1. 改路徑
SRC_ROOT <- "Z:/Dataset_Online/scRNA-seq/loom"
OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"

# 2. source 執行
source("01e_QC_Loom.R")
```

資料夾結構(`<GSE>/` 之下放檔,每檔一樣本):

```
loom/
└── #GSE153596/
    ├── GSM4647785_P2_1.loom.gz      ← 自動解壓後讀
    ├── GSM4647791_P21Wound.loom.gz
    └── ...
```

## 前置需求

- **hdf5r**(讀 loom):`install.packages("hdf5r")`,沒裝會直接報錯提示。
- 不需要 `SeuratDisk`。刻意不用 `SeuratDisk::LoadLoom`,因為它對「非標準命名」的 loom
  常直接失敗;這裡用 `hdf5r` 直接讀 + 自己偵測,較穩健。

## 輸出(與 01~01d 完全相同)

```
#Keloid_QC/
├── <GSE>/<sample>/
│   ├── <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv
│   └── plots/  (violin pre/post, scatter pre, before/after)
└── _run_log/
    ├── per_sample/<GSE>__<sample>.csv      ← 02 彙整用,欄位與 01 一致
    ├── structure_scan_loom.csv             ← 每個 loom 讀到的方向/欄位(load_note)
    └── run_log_loom.txt
```

樣本名自動去 `.loom` / `.loom.gz`(`GSM4647785_P2_1.loom.gz` → `GSM4647785_P2_1`)。

## 可調參數(腳本上方)

和 `01` 對齊:`MIN_FEATURE` / `MAX_FEATURE` / `MAX_MT_PCT` / `RUN_DOUBLETFINDER` / `IGNORE_FOLDERS`。

`.loom` 專屬:
- `GZ_TMP_DIR` — 解壓 `.loom.gz` 的暫存目錄(預設系統暫存;**大檔請確認空間足夠**,
  你截圖有些 loom 解壓後上百 MB)。
- `LOOM_MATRIX_LAYER` — 主矩陣用哪個 layer(預設 `"matrix"`;velocyto loom 若想只用
  剪接後的可改 `"spliced"`)。

## 小提醒(重要)

- **spliced/unspliced**:velocyto 產生的 loom 帶 `spliced`/`unspliced` layers,`/matrix`
  通常是總表現量。做一般 QC 用 `/matrix` 即可(預設)。若你要做 RNA velocity 分析那是另一回事。
- **這批可能是 spatial**:你截圖分頁標題有 "spatial and single-cell"、圖裡有組織影像。
  若這個 GSE 其實是**空間轉錄體**(Visium 等),`.loom` 裡的「cell」可能是 spot,QC 門檻
  (200~5000 / mt<30)不一定適用 spatial —— spatial 的 per-spot QC 慣例不同。先跑沒問題,
  但**留意結果**,必要時 spatial 資料另訂門檻。
- 跑完看 `structure_scan_loom.csv` 的 `load_note` 欄,確認每個 loom 判到的方向與 gene/cell
  欄位是否合理。
