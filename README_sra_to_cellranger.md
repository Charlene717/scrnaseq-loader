# SRA → Cell Ranger 流程腳本

把從 SRA 下載的 10x scRNA-seq 原始資料(`.sra`),轉成可直接接 `01_QC_PerSample.R` 的
count matrix(`barcodes` / `features` / `matrix`)。

兩支腳本:

| 腳本 | 做什麼 | 何時用 |
|---|---|---|
| `inspect_sra.sh` | 唯讀偵測:是不是 10x、化學版本(v2/v3)、物種、read 配置 | **開跑前先確認**,只轉前 N 條 reads,不動原始檔 |
| `sra_to_cellranger.sh` | 一條龍:`.sra` → fastq → 改名 → `cellranger count` → 整理 matrix | 確認後正式跑,支援批次與續跑 |

---

## 為什麼需要這兩步

`.sra` 不是分析能直接吃的格式。10x 資料的完整鏈是:

```
.sra  →(fasterq-dump)→  10x FASTQ  →(Cell Ranger)→  barcodes+features+matrix  →  01_QC_PerSample.R
```

兩個容易踩的坑,腳本都自動處理了:

1. **fastq 命名**：Cell Ranger 要求 `<Sample>_S1_L001_R1_001.fastq.gz` 這種格式,但 SRA 轉出來是 `SRR..._1/_2/_3`,必須改名。
2. **R1/R2/I1 對應**：SRA 轉出的 `_1/_2/_3` 不保證對應 barcode/cDNA/index。腳本用 **read 長度** 自動判斷：
   - cDNA(R2)= 最長(90~150bp)
   - barcode+UMI(R1)= 26bp(v2)或 28bp(v3)
   - sample index(I1)= 8bp
   - ⚠️ 注意:index(8bp)比 barcode(28bp)**還短**,所以不能用「最短=barcode」,腳本已正確處理。

---

## 依賴安裝

```bash
# 必要
conda install -c bioconda sra-tools          # fasterq-dump
# Cell Ranger 需從 10x 官網下載後加入 PATH(需註冊):
#   https://www.10xgenomics.com/support/software/cell-ranger/downloads

# 選用(加速 / 自動物種偵測)
conda install -c bioconda pigz minimap2 seqkit
```

參考基因組(Cell Ranger 官方預建):

```bash
# human (GRCh38)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
# mouse 同理下載 refdata-gex-GRCm39-2024-A
```

---

## 用法

### 第一步:先偵測(建議)

```bash
./inspect_sra.sh /data/PRJNA772373/SRR16475068
```

會印出是不是 10x、化學版本、R1/R2 對應。想順便自動判物種,給它 transcriptome:

```bash
HUMAN_REF=/ref/human_cdna.fa.gz MOUSE_REF=/ref/mouse_cdna.fa.gz \
  ./inspect_sra.sh SRR16475068
```

### 第二步:正式跑

**單筆,手動指定物種:**

```bash
./sra_to_cellranger.sh \
    --species human \
    --ref /ref/refdata-gex-GRCh38-2020-A \
    /data/PRJNA772373/SRR16475068
```

**整個 PRJNA 資料夾(每個 SRR 一個樣本),自動偵測物種:**

```bash
./sra_to_cellranger.sh \
    --human-ref /ref/refdata-gex-GRCh38-2020-A \
    --mouse-ref /ref/refdata-gex-GRCm39-2024-A \
    --auto-species \
    --threads 16 --localmem 64 \
    /data/PRJNA772373/
```

### 主要參數

| 參數 | 說明 |
|---|---|
| `--species human\|mouse` | 手動指定物種 |
| `--ref PATH` | 直接指定 Cell Ranger 參考(優先於 species 對應) |
| `--human-ref` / `--mouse-ref` | 各物種的參考路徑(搭配 `--auto-species`) |
| `--auto-species` | 自動偵測物種(需 `inspect_sra.sh` + minimap2 + transcriptome) |
| `--out DIR` | 輸出根目錄(預設 `./cellranger_out`) |
| `--threads N` | 執行緒(預設 8) |
| `--localmem GB` | Cell Ranger 記憶體上限(預設 32) |
| `--chemistry auto` | 化學版本(預設 auto,讓 Cell Ranger 自己判斷) |
| `--expect-cells N` | 預期細胞數(選用) |
| `--keep-fastq` | 保留中間 fastq(預設刪除省空間) |

---

## 輸出

```
cellranger_out/
├── SRR16475068/
│   ├── SRR16475068_cellranger/outs/   # Cell Ranger 完整輸出
│   └── filtered_feature_bc_matrix/    # ← 接 01_QC 的三個檔在這
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── _logs/
    ├── run_summary.csv                # 每筆狀態彙整(OK/FAILED/SKIPPED + 細胞數)
    ├── SRR16475068_fasterqdump.log
    └── SRR16475068_cellranger.log
```

把 `filtered_feature_bc_matrix/` 放進你的 `#Keloid/<GSE>/<GSM>/` 結構,就能接 `01_QC_PerSample.R`。

---

## 續跑

兩支腳本都可中斷後重跑:已經有 Cell Ranger `outs/` 的樣本會自動跳過,
失敗的樣本記錄在 `run_summary.csv` 後跳過、不中斷整批。

---

## 關於那兩個下載檔(`SRR...` vs `SRR....lite.1`)

`prefetch` 有時會同時抓「標準版」和「lite 版」。兩者序列相同,差別只在品質分數的精細度,
對 Cell Ranger 流程**結果幾乎無差異**。**留一個就好**(建議留標準版),刪掉另一個省空間。
本腳本對兩種都能處理。
