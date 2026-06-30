# 01c_QC_FromSeuratRDS.R — 從「已整合 Seurat RDS」按 sample 做 QC

`01` / `01b` 的第三個姊妹版。處理「**已經整合 + 註解好的單一大 Seurat 物件 (.rds)**」,
把它依 `meta.data` 的某個欄位(預設 `sample`)拆開,每個樣本各自跑與 `01` 相同的 QC。

適用你手上的 GSE307504:13.3 GB、58,064 cells、11 個 sample
(Skin01–05 / Ke01–02 / Hyper01–02 / NScar01–02)。

## 三個這個物件特有、會影響正確性的處理(務必了解)

### 1. percent.mito 單位問題 → 一律重算 percent.mt
你的物件 `meta.data` 已有 `percent.mito`,但值是 `0.027`、`0.010` 這種範圍 —— 這是
**比例 (0~1)**,不是百分比 (0~100)。若直接拿它套 `01` 的 `MAX_MT_PCT = 30` 門檻,
**永遠不會濾掉任何細胞**(因為最大值才 1)。

所以本腳本**不沿用舊欄位**,一律用 `PercentageFeatureSet` 重新算 `percent.mt`(0~100),
和 `01` 完全一致。這是刻意的選擇,避免單位不一致導致過濾失效。

### 2. 這是「已 QC + 已整合」資料 → 過濾用開關控制
RNA assay 的 counts 是原始處理時**就已過濾過**的 counts。再套 200~5000 / mt<30 等於
**第二次過濾**。所以過濾做成開關:

```r
APPLY_FILTER <- TRUE   # 與 01 同理,會依門檻砍細胞(這是第二次過濾)
APPLY_FILTER <- FALSE  # 只算 QC 統計 + 畫圖,不砍任何細胞(post = pre)
```

如果你的目的只是「**看每個 sample 的品質分佈 / 出 QC 報告**」,而不想動已經整合好的細胞組成,
建議設 `FALSE`。如果你要的是「**比照 01 重跑一次過濾**」,維持 `TRUE`。

### 3. 記憶體 + Seurat 版本相容
- 13.3 GB 物件 + 你環境已用 ~12.85 GB。本腳本用「**切 counts 矩陣**」而非 `SplitObject`
  (後者會複製整個物件結構),逐 sample 處理並 `gc()`,盡量省記憶體。
- 即便如此,擷取 RNA counts 仍會多佔幾 GB。**若記憶體吃緊**,跑前可手動釋放整合 assay:
  ```r
  DefaultAssay(SeuratObject) <- "RNA"
  SeuratObject[["integrated"]] <- NULL
  gc()
  ```
- 相容 **Seurat v4 與 v5**(counts 取法不同;v5 若 layer 被拆分會先 `JoinLayers`)。

## 用法

```r
# 1. 你已經在 RStudio 把物件讀進來了(變數名 SeuratObject)
#    -> 本腳本偵測到環境已有 SeuratObject 就直接用,不會再讀一次 13GB

# 2. 確認/修改腳本上方設定
SAMPLE_COL <- "sample"        # 用哪個欄位分樣本(你的是 sample)
GSE_LABEL  <- "GSE307504"     # 輸出資料夾名
APPLY_FILTER <- TRUE          # 或 FALSE(只出統計+圖)
OUT_ROOT   <- "Z:/.../#Keloid_QC"

# 3. source 執行
source("01c_QC_FromSeuratRDS.R")
```

若物件還沒載入,把上方 `RDS_PATH` 設成你的 `.rds` 路徑,腳本會自己讀。

## 條件 (condition) 自動標註

依 sample 名前綴推測(此資料集):

| sample 前綴 | condition |
|---|---|
| `Ke...` | Keloid |
| `Hyper...` | HypertrophicScar |
| `NScar...` | NormalScar |
| `Skin...` | NormalSkin |

推測規則在腳本的 `guess_condition()`,可自行調整。

## 輸出(與 01 / 01b 完全相同)

```
#Keloid_QC/
├── GSE307504/<sample>/
│   ├── <sample>_pre_QC.rds      ← 該樣本乾淨的 per-sample 物件(過濾前)
│   ├── <sample>_post_QC.rds     ← 過濾後(APPLY_FILTER=FALSE 時 = pre)
│   ├── <sample>_QC_stats.csv
│   └── plots/  (violin pre/post, scatter pre, before/after)
└── _run_log/
    ├── per_sample/GSE307504__<sample>.csv   ← 02 彙整用,欄位與 01 一致
    └── run_log_fromrds.txt
```

每個 per-sample 物件是**從 RNA counts 重建的乾淨物件**(不帶整合結構,輕量),
但會把既有的 `celltype` 等註解帶進 `meta.data` 作參照(可在腳本 `CARRY_META` 調整)。

## 重要說明(請讀)

- **我無法在這裡對你的真實 13GB 物件實際執行**(這個環境沒有 R,也沒有你的物件)。
  腳本的結構、括號配對、函式定義順序、Seurat API 用法我都靜態檢查過,但首次執行時
  請先確認:① `Assays(SeuratObject)` 裡有 `RNA`;② `SeuratObject$sample` 是你要的分組。
- 建議**先用一個樣本試跑**:把 `samples` 那行暫時改成 `samples <- samples[1]`,確認流程
  與輸出無誤,再跑全部 11 個。
- 若 `RNA` assay 的 counts 其實不是原始 counts(有些物件只留 integrated),`percent.mt`
  和過濾會失真 —— 腳本會在 log 警告,留意。
