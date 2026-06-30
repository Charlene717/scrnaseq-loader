###############################################################################
## 01_QC_PerSample.R
## ---------------------------------------------------------------------------
## scRNA-seq (10X Genomics) 廣泛品質評估 + 初步 QC  (per-sample)
##
## 資料結構 (來源):
##   Z:\Dataset_Online\scRNA-seq\#Keloid\<GSE...>\<GSM..._XX>\{barcodes,features,matrix}
##   - GSE 資料夾可能有 "#" 前綴 (視為「不符規定 / 未完成」-> 記錄但仍嘗試)
##   - 也可能有非資料夾項目 (如 GSE293834_Ori.txt) -> 跳過並記錄
##   - GSM 內 10X 三檔可能為 .gz 或未壓縮; matrix 可能是 .mtx / .mtx.gz
##
## 輸出 (目的地):
##   Z:\Dataset_Online\scRNA-seq\#Keloid_QC\<GSE...>\<GSM...>\
##     - <GSM>_pre_QC.rds          (QC 前 Seurat 物件)
##     - <GSM>_post_QC.rds         (QC 後 Seurat 物件)
##     - <GSM>_QC_stats.csv        (該樣本前後狀態詳盡紀錄)
##     - plots/  (violin pre/post, scatter pre/post, before-after, doublet)
##   ...以及 _run_log/ 內每個樣本的 per-sample metrics row (供 02 腳本彙整)
##
## 特性:
##   * 可續跑: 若該樣本所有輸出 (2 RDS + stats + 預期圖) 都齊全則跳過
##   * 物種自動偵測 (^MT- 人類 / ^mt- 小鼠)
##   * 過濾門檻採實驗室慣例 (固定): nFeature 200~5000, percent.mt < 30
##     (同時計算 MAD adaptive 門檻僅作參考, 不套用)
##   * DoubletFinder 偵測 (若安裝)
##   * 跑不出東西的樣本 -> 記錄錯誤, 跳過, 繼續下一個
###############################################################################

## ----------------------------- 0. 路徑設定 -------------------------------- ##
## 在 RStudio (Windows 端) 直接跑時，用 Windows 路徑；若在 Linux 端跑請改成掛載路徑。
SRC_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid"      # 來源根目錄
OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"   # QC 輸出根目錄

## ------ 過濾門檻 (實驗室慣例，固定值) ------ ##
MIN_FEATURE <- 200      # nFeature_RNA 下限
MAX_FEATURE <- 5000     # nFeature_RNA 上限
MAX_MT_PCT  <- 30       # percent.mt 上限 (%)
MIN_CELLS_GENE <- 3     # CreateSeuratObject min.cells (基因至少在幾顆細胞表現)
MIN_FEATURES_CELL <- 0  # CreateSeuratObject min.features (建物件時先不濾，QC 步驟再濾)

## ====== 使用者可調參數 ====== ##
## 是否執行 DoubletFinder 偵測 (TRUE/FALSE)。
##   TRUE  -> 每樣本跑 doublet 偵測 (需標準化/PCA/UMAP，大樣本較耗時)；沒裝套件會自動略過並記錄
##   FALSE -> 完全跳過 doublet 偵測，大表中 n_doublets 記為 NA、doublet_method 記為 "skipped(user)"
RUN_DOUBLETFINDER <- FALSE

## 要忽略 (不處理) 的頂層資料夾名稱清單；會記錄為 IGNORED_BY_USER 但不做 QC。
IGNORE_FOLDERS <- c("GEO_10X_auto")
## =============================== ##

DEFAULT_DOUBLET_RATE <- 0.075  # 預設 doublet 率 (~7.5%, 10X 約每 1000 cells 0.8%); 下面會依細胞數估算
SEED <- 42

## ----------------------------- 1. 套件 ------------------------------------ ##
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("需要 Seurat 套件，請先 install.packages('Seurat')")
  library(Seurat)
  library(Matrix)
})
## 選用套件 (沒裝不致命)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAVE_GGPLOT  <- has_pkg("ggplot2");  if (HAVE_GGPLOT) suppressPackageStartupMessages(library(ggplot2))
HAVE_PATCH   <- has_pkg("patchwork");if (HAVE_PATCH)  suppressPackageStartupMessages(library(patchwork))
HAVE_DF      <- has_pkg("DoubletFinder")
HAVE_SCDBL   <- has_pkg("scDblFinder")  # 備援 doublet 方法
set.seed(SEED)

## ----------------------------- 2. 工具函式 -------------------------------- ##

## 安全建立資料夾
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

## 時間戳記字串 (寫 log 用)
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

## 跨平台穩健的 UTF-8 (帶 BOM) CSV 寫出，避免 Windows/Linux locale 轉換破壞中文，
## 且 Windows Excel 直接開啟即正確顯示中文。
write_csv_utf8 <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_excel_csv(df, path)            # 自帶 UTF-8 BOM
    return(invisible(path))
  }
  lines <- utils::capture.output(utils::write.csv(df, file = "", row.names = FALSE))
  con <- file(path, open = "wb")                 # binary，不做 locale 轉換
  writeBin(charToRaw("\xEF\xBB\xBF"), con)        # UTF-8 BOM
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
  close(con)
  invisible(path)
}

## 寫一行到全域 run log
LOG_DIR <- file.path(OUT_ROOT, "_run_log")
log_line <- function(msg) {
  ensure_dir(LOG_DIR)
  cat(sprintf("[%s] %s\n", ts(), msg))
  cat(sprintf("[%s] %s\n", ts(), msg),
      file = file.path(LOG_DIR, "run_log.txt"), append = TRUE)
}

## 在 GSM 資料夾中尋找 10X 三個檔案 (容忍 .gz / 命名差異)，回傳 list 或 NULL
## 回傳: list(dir=, barcodes=, features=, matrix=, gzipped=, note=)
locate_10x <- function(gsm_dir) {
  files <- list.files(gsm_dir, full.names = FALSE)
  lf <- tolower(files)

  find_one <- function(patterns) {
    for (pat in patterns) {
      hit <- files[grepl(pat, lf)]
      if (length(hit) >= 1) return(hit[1])
    }
    return(NA_character_)
  }

  bc  <- find_one(c("^barcodes\\.tsv\\.gz$", "^barcodes\\.tsv$",
                    "barcodes.*\\.tsv\\.gz$", "barcodes.*\\.tsv$"))
  ft  <- find_one(c("^features\\.tsv\\.gz$", "^features\\.tsv$",
                    "^genes\\.tsv\\.gz$",    "^genes\\.tsv$",
                    "features.*\\.tsv\\.gz$", "features.*\\.tsv$",
                    "genes.*\\.tsv\\.gz$",   "genes.*\\.tsv$"))
  mx  <- find_one(c("^matrix\\.mtx\\.gz$", "^matrix\\.mtx$",
                    "matrix.*\\.mtx\\.gz$", "matrix.*\\.mtx$"))

  note <- character(0)
  if (is.na(bc)) note <- c(note, "缺 barcodes")
  if (is.na(ft)) note <- c(note, "缺 features/genes")
  if (is.na(mx)) note <- c(note, "缺 matrix.mtx")

  ok <- !any(is.na(c(bc, ft, mx)))
  list(
    ok       = ok,
    dir      = gsm_dir,
    barcodes = if (!is.na(bc)) file.path(gsm_dir, bc) else NA,
    features = if (!is.na(ft)) file.path(gsm_dir, ft) else NA,
    matrix   = if (!is.na(mx)) file.path(gsm_dir, mx) else NA,
    note     = if (length(note)) paste(note, collapse = "; ") else "OK"
  )
}

## 讀取 10X 三檔成 sparse matrix (處理壓縮與否、features 欄位數)
read_10x_manual <- function(loc) {
  read_lines_any <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
    on.exit(close(con))
    readLines(con)
  }
  barcodes <- read_lines_any(loc$barcodes)
  feat_raw <- read_lines_any(loc$features)
  feat_tab <- do.call(rbind, strsplit(feat_raw, "\t", fixed = TRUE))
  ## features.tsv 第2欄通常是 gene symbol；只有1欄就用第1欄
  gene_names <- if (ncol(feat_tab) >= 2) feat_tab[, 2] else feat_tab[, 1]
  gene_names <- make.unique(as.character(gene_names))

  mat <- Matrix::readMM(if (grepl("\\.gz$", loc$matrix)) gzfile(loc$matrix) else loc$matrix)
  mat <- as(mat, "CsparseMatrix")
  rownames(mat) <- gene_names
  colnames(mat) <- make.unique(as.character(barcodes))
  mat
}

## 偵測物種 / 粒線體 pattern；回傳 list(pattern, species, n_mt, ribo_pattern)
detect_mt_pattern <- function(genes) {
  n_human <- sum(grepl("^MT-",  genes))
  n_mouse <- sum(grepl("^mt-",  genes))
  if (n_human >= n_mouse && n_human > 0) {
    list(pattern = "^MT-", species = "human", n_mt = n_human, ribo = "^RP[SL]")
  } else if (n_mouse > 0) {
    list(pattern = "^mt-", species = "mouse", n_mt = n_mouse, ribo = "^Rp[sl]")
  } else {
    ## 找不到任何 MT 基因 -> 仍回傳人類 pattern，但標記 0
    list(pattern = "^MT-", species = "unknown(no MT genes)", n_mt = 0, ribo = "^RP[SL]")
  }
}

## 由一個 Seurat 物件擷取摘要統計 (用於前後比較與大表)
summarise_obj <- function(obj, label) {
  md <- obj@meta.data
  n_cells <- ncol(obj)
  n_genes_detected <- sum(Matrix::rowSums(GetAssayData(obj, slot = "counts") > 0) > 0)
  safe_mean   <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
  safe_median <- function(x) if (length(x)) median(x, na.rm = TRUE) else NA_real_
  data.frame(
    state                 = label,
    n_cells               = n_cells,
    n_genes_detected      = n_genes_detected,
    mean_nCount           = safe_mean(md$nCount_RNA),
    median_nCount         = safe_median(md$nCount_RNA),
    mean_nFeature         = safe_mean(md$nFeature_RNA),
    median_nFeature       = safe_median(md$nFeature_RNA),
    mean_percent_mt       = safe_mean(md$percent.mt),
    median_percent_mt     = safe_median(md$percent.mt),
    mean_percent_ribo     = if (!is.null(md$percent.ribo)) safe_mean(md$percent.ribo) else NA_real_,
    stringsAsFactors = FALSE
  )
}

## MAD adaptive 門檻 (僅報告，不套用)
mad_bounds <- function(x, nmads = 3, log = FALSE) {
  v <- if (log) log1p(x) else x
  med <- median(v, na.rm = TRUE); m <- mad(v, na.rm = TRUE)
  lo <- med - nmads * m; hi <- med + nmads * m
  if (log) { lo <- expm1(lo); hi <- expm1(hi) }
  c(lower = lo, upper = hi)
}

## ----------------------------- 3. 掃描資料夾結構 -------------------------- ##
log_line("==== QC run 開始 ====")
log_line(sprintf("來源: %s", SRC_ROOT))
log_line(sprintf("輸出: %s", OUT_ROOT))
ensure_dir(OUT_ROOT); ensure_dir(LOG_DIR)

if (!dir.exists(SRC_ROOT)) stop(sprintf("找不到來源資料夾: %s (請確認 Z: 已掛載)", SRC_ROOT))

## 用於記錄「結構掃描」結果 (含不符規定者) -> 給 02 腳本與 manifest
structure_rows <- list()
add_struct <- function(gse_raw, gse_clean, gsm, status, note,
                       has_hash = FALSE, src_path = "") {
  structure_rows[[length(structure_rows) + 1]] <<- data.frame(
    GSE_folder_raw = gse_raw, GSE = gse_clean, GSM = gsm,
    structure_status = status, structure_note = note,
    gse_has_hash_prefix = has_hash, source_path = src_path,
    stringsAsFactors = FALSE
  )
}

## 列出 SRC_ROOT 下所有項目
top_items <- list.files(SRC_ROOT, full.names = TRUE)
gse_dirs  <- character(0)

for (it in top_items) {
  base <- basename(it)
  ## 使用者指定忽略的資料夾 -> 記錄但不處理
  if (base %in% IGNORE_FOLDERS) {
    add_struct(base, sub("^#", "", base), NA, "IGNORED_BY_USER",
               sprintf("使用者指定忽略 (%s)，未做 QC", base),
               has_hash = grepl("^#", base), src_path = it)
    log_line(sprintf("[ignore] 依使用者設定跳過資料夾: %s", base))
    next
  }
  if (!dir.exists(it)) {
    ## 非資料夾 (如 .txt) -> 記錄並跳過
    add_struct(base, sub("^#", "", tools::file_path_sans_ext(base)),
               NA, "NON_FOLDER_SKIPPED",
               sprintf("頂層非資料夾項目 (%s), 已跳過", base),
               has_hash = grepl("^#", base), src_path = it)
    log_line(sprintf("跳過非資料夾頂層項目: %s", base))
    next
  }
  gse_dirs <- c(gse_dirs, it)
}

log_line(sprintf("偵測到 %d 個 GSE 資料夾", length(gse_dirs)))

## ----------------------------- 4. 逐樣本主迴圈 ---------------------------- ##
## 每個樣本完成後寫一個 per-sample metrics row 到 _run_log/per_sample/<GSE>__<GSM>.csv
## (02 腳本會把這些 row 合併成大表，避免中途中斷遺失進度)
PERSAMPLE_DIR <- file.path(LOG_DIR, "per_sample")
ensure_dir(PERSAMPLE_DIR)

## 寫單一樣本 metrics row 的函式 (欄位即大表欄位)
write_metrics_row <- function(row) {
  fn <- file.path(PERSAMPLE_DIR,
                  sprintf("%s__%s.csv", row$GSE, ifelse(is.na(row$GSM), "NA", row$GSM)))
  write_csv_utf8(row, fn)
}

for (gse_dir in gse_dirs) {
  gse_raw   <- basename(gse_dir)
  gse_hash  <- grepl("^#", gse_raw)
  gse_clean <- sub("^#", "", gse_raw)        # 去掉 # 前綴後的 GSE 編號
  if (gse_hash)
    log_line(sprintf("[!] GSE '%s' 有 # 前綴 (視為未完成/不符規定)，仍嘗試處理", gse_raw))

  gsm_items <- list.files(gse_dir, full.names = TRUE)
  gsm_dirs  <- gsm_items[dir.exists(gsm_items)]

  if (length(gsm_dirs) == 0) {
    add_struct(gse_raw, gse_clean, NA, "EMPTY_GSE",
               "GSE 內無 GSM 子資料夾", has_hash = gse_hash, src_path = gse_dir)
    log_line(sprintf("[!] GSE '%s' 內無 GSM 子資料夾，跳過", gse_raw))
    write_metrics_row(data.frame(
      GSE = gse_clean, GSE_folder_raw = gse_raw, GSM = NA, condition = NA,
      status = "NO_GSM_SUBFOLDER", error_note = "GSE 內無 GSM 子資料夾",
      gse_has_hash_prefix = gse_hash,
      species = NA, mt_pattern = NA,
      pre_n_cells = NA, pre_n_genes = NA, pre_mean_nCount = NA, pre_mean_nFeature = NA,
      pre_mean_pct_mt = NA, pre_median_pct_mt = NA,
      post_n_cells = NA, post_n_genes = NA, post_mean_nCount = NA, post_mean_nFeature = NA,
      post_mean_pct_mt = NA, post_median_pct_mt = NA,
      cells_removed = NA, pct_cells_removed = NA,
      n_doublets = NA, doublet_method = NA,
      filter_used = NA, source_path = gse_dir, output_path = NA,
      processed_time = ts(), stringsAsFactors = FALSE))
    next
  }

  for (gsm_dir in gsm_dirs) {
    gsm <- basename(gsm_dir)
    ## 由資料夾名推測 condition (KL=Keloid / NS=Normal Skin 等後綴)
    cond <- sub("^GSM[0-9]+_?", "", gsm); if (cond == gsm || cond == "") cond <- NA

    out_dir   <- file.path(OUT_ROOT, gse_clean, gsm)
    plots_dir <- file.path(out_dir, "plots")
    pre_rds   <- file.path(out_dir, paste0(gsm, "_pre_QC.rds"))
    post_rds  <- file.path(out_dir, paste0(gsm, "_post_QC.rds"))
    stats_csv <- file.path(out_dir, paste0(gsm, "_QC_stats.csv"))

    ## ---- 續跑檢查: 全部齊全 -> 跳過 ----
    expected_plots <- c(
      file.path(plots_dir, paste0(gsm, "_violin_preQC.pdf")),
      file.path(plots_dir, paste0(gsm, "_violin_postQC.pdf")),
      file.path(plots_dir, paste0(gsm, "_scatter_preQC.pdf")),
      file.path(plots_dir, paste0(gsm, "_beforeafter.pdf"))
    )
    all_exist <- file.exists(pre_rds) && file.exists(post_rds) &&
      file.exists(stats_csv) && all(file.exists(expected_plots))
    if (all_exist) {
      log_line(sprintf("[skip] %s/%s 已有完整輸出，跳過", gse_clean, gsm))
      add_struct(gse_raw, gse_clean, gsm, "DONE_SKIPPED",
                 "已有完整輸出，續跑跳過", has_hash = gse_hash, src_path = gsm_dir)
      next
    }

    log_line(sprintf("==> 處理 %s / %s", gse_clean, gsm))

    ## ---- 用 tryCatch 包住整個樣本，失敗就記錄並繼續 ----
    res <- tryCatch({
      ## (a) 定位 10X 檔案
      loc <- locate_10x(gsm_dir)
      if (!loc$ok) {
        add_struct(gse_raw, gse_clean, gsm, "MISSING_10X_FILES", loc$note,
                   has_hash = gse_hash, src_path = gsm_dir)
        stop(sprintf("10X 檔案不齊全: %s", loc$note))
      }
      add_struct(gse_raw, gse_clean, gsm, "OK", loc$note,
                 has_hash = gse_hash, src_path = gsm_dir)

      ensure_dir(out_dir); ensure_dir(plots_dir)

      ## (b) 讀取 -> 先試 Seurat::Read10X (整資料夾)，失敗改手動
      mat <- tryCatch({
        ## Read10X 需要標準命名；若檔名非標準會失敗，故包 try
        Read10X(data.dir = gsm_dir)
      }, error = function(e) {
        log_line(sprintf("    Read10X 失敗(%s)，改用手動讀取", conditionMessage(e)))
        read_10x_manual(loc)
      })
      ## 若回傳的是 list (多 assay)，取 Gene Expression
      if (is.list(mat) && !inherits(mat, "dgCMatrix")) {
        if ("Gene Expression" %in% names(mat)) mat <- mat[["Gene Expression"]]
        else mat <- mat[[1]]
      }

      ## (c) 建 Seurat 物件
      obj <- CreateSeuratObject(counts = mat, project = gsm,
                                min.cells = MIN_CELLS_GENE,
                                min.features = MIN_FEATURES_CELL)
      ## 依實驗室慣例加 metadata
      obj$orig.ident1 <- gse_clean
      obj$orig.ident2 <- gsm
      if (!is.na(cond)) obj$condition <- cond

      ## (d) 物種 / MT pattern 偵測
      det <- detect_mt_pattern(rownames(obj))
      obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = det$pattern)
      obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = det$ribo)
      log_line(sprintf("    物種偵測: %s | MT pattern=%s | MT基因數=%d",
                       det$species, det$pattern, det$n_mt))

      ## (e) QC 前摘要 + MAD adaptive 門檻 (僅報告)
      pre_sum <- summarise_obj(obj, "pre_QC")
      mad_feat  <- mad_bounds(obj$nFeature_RNA, 3, log = TRUE)
      mad_count <- mad_bounds(obj$nCount_RNA,   3, log = TRUE)
      mad_mt    <- mad_bounds(obj$percent.mt,   3, log = FALSE)

      ## (f) 存 QC 前 RDS
      saveRDS(obj, pre_rds)

      ## (g) 繪圖: 小提琴 (前)
      if (HAVE_GGPLOT) {
        feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
        feats <- feats[feats %in% colnames(obj@meta.data)]
        p_vln_pre <- VlnPlot(obj, features = feats, ncol = length(feats), pt.size = 0.1)
        pdf(expected_plots[1], width = 4 * length(feats), height = 5)
        print(p_vln_pre); dev.off()

        ## 散點圖 (前): count vs feature, count vs mt
        s1 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + ggtitle("nCount vs nFeature")
        s2 <- FeatureScatter(obj, "nCount_RNA", "percent.mt")   + ggtitle("nCount vs percent.mt")
        pdf(expected_plots[3], width = 12, height = 5)
        if (HAVE_PATCH) print(s1 + s2) else { print(s1); print(s2) }
        dev.off()
      }

      ## (h) 套用過濾 (固定門檻, 實驗室慣例)
      filter_expr <- sprintf("nFeature_RNA > %d & nFeature_RNA < %d & percent.mt < %d",
                             MIN_FEATURE, MAX_FEATURE, MAX_MT_PCT)
      n_before <- ncol(obj)
      obj_post <- subset(obj, subset = nFeature_RNA > MIN_FEATURE &
                                       nFeature_RNA < MAX_FEATURE &
                                       percent.mt   < MAX_MT_PCT)
      n_after <- ncol(obj_post)

      ## (i) DoubletFinder (在過濾後、標準化後跑)；失敗或未裝就記錄
      n_doublets <- NA_integer_
      dbl_method <- if (!RUN_DOUBLETFINDER) "skipped(user)" else "none"
      if (RUN_DOUBLETFINDER && n_after > 50) {
        dbl <- tryCatch({
          run_doublet(obj_post, det)
        }, error = function(e) {
          log_line(sprintf("    Doublet 偵測失敗: %s", conditionMessage(e)))
          NULL
        })
        if (!is.null(dbl)) {
          obj_post   <- dbl$obj
          n_doublets <- dbl$n_doublets
          dbl_method <- dbl$method
          ## doublet 視覺化
          if (HAVE_GGPLOT && !is.null(dbl$plot)) {
            pdf(file.path(plots_dir, paste0(gsm, "_doublet.pdf")), width = 7, height = 6)
            print(dbl$plot); dev.off()
          }
        }
      }

      ## (j) QC 後摘要
      post_sum <- summarise_obj(obj_post, "post_QC")

      ## (k) 小提琴 (後) + before/after 對比
      if (HAVE_GGPLOT) {
        feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
        feats <- feats[feats %in% colnames(obj_post@meta.data)]
        p_vln_post <- VlnPlot(obj_post, features = feats, ncol = length(feats), pt.size = 0.1)
        pdf(expected_plots[2], width = 4 * length(feats), height = 5)
        print(p_vln_post); dev.off()

        ## before/after 對比 (同一張圖比 nFeature/nCount/percent.mt 分佈)
        ba <- rbind(
          data.frame(stage = "pre",  nFeature = obj$nFeature_RNA,
                     nCount = obj$nCount_RNA, percent.mt = obj$percent.mt),
          data.frame(stage = "post", nFeature = obj_post$nFeature_RNA,
                     nCount = obj_post$nCount_RNA, percent.mt = obj_post$percent.mt)
        )
        g1 <- ggplot(ba, aes(stage, nFeature, fill = stage)) + geom_violin() +
              geom_boxplot(width = .1, outlier.size = .3) + theme_bw() + ggtitle("nFeature_RNA")
        g2 <- ggplot(ba, aes(stage, nCount, fill = stage)) + geom_violin() +
              geom_boxplot(width = .1, outlier.size = .3) + theme_bw() +
              scale_y_log10() + ggtitle("nCount_RNA (log10)")
        g3 <- ggplot(ba, aes(stage, percent.mt, fill = stage)) + geom_violin() +
              geom_boxplot(width = .1, outlier.size = .3) + theme_bw() + ggtitle("percent.mt")
        pdf(expected_plots[4], width = 12, height = 5)
        if (HAVE_PATCH) print(g1 + g2 + g3) else { print(g1); print(g2); print(g3) }
        dev.off()
      }

      ## (l) 存 QC 後 RDS
      saveRDS(obj_post, post_rds)

      ## (m) 寫該樣本詳盡 stats CSV (前後 + MAD 門檻 + 過濾資訊)
      stats_df <- rbind(pre_sum, post_sum)
      stats_df$GSE <- gse_clean; stats_df$GSM <- gsm; stats_df$condition <- cond
      stats_df$species <- det$species; stats_df$mt_pattern <- det$pattern
      stats_df$n_mt_genes <- det$n_mt
      stats_df$filter_used <- filter_expr
      stats_df$MAD_nFeature_lower <- round(mad_feat["lower"], 1)
      stats_df$MAD_nFeature_upper <- round(mad_feat["upper"], 1)
      stats_df$MAD_nCount_lower   <- round(mad_count["lower"], 1)
      stats_df$MAD_nCount_upper   <- round(mad_count["upper"], 1)
      stats_df$MAD_percent_mt_upper <- round(mad_mt["upper"], 2)
      stats_df$cells_removed <- n_before - n_after
      stats_df$pct_cells_removed <- round(100 * (n_before - n_after) / n_before, 2)
      stats_df$n_doublets <- n_doublets
      stats_df$doublet_method <- dbl_method
      stats_df$processed_time <- ts()
      write_csv_utf8(stats_df, stats_csv)

      ## (n) 寫 per-sample metrics row (大表用)
      write_metrics_row(data.frame(
        GSE = gse_clean, GSE_folder_raw = gse_raw, GSM = gsm, condition = cond,
        status = "OK", error_note = "",
        gse_has_hash_prefix = gse_hash,
        species = det$species, mt_pattern = det$pattern,
        pre_n_cells = pre_sum$n_cells, pre_n_genes = pre_sum$n_genes_detected,
        pre_mean_nCount = round(pre_sum$mean_nCount, 1),
        pre_mean_nFeature = round(pre_sum$mean_nFeature, 1),
        pre_mean_pct_mt = round(pre_sum$mean_percent_mt, 3),
        pre_median_pct_mt = round(pre_sum$median_percent_mt, 3),
        post_n_cells = post_sum$n_cells, post_n_genes = post_sum$n_genes_detected,
        post_mean_nCount = round(post_sum$mean_nCount, 1),
        post_mean_nFeature = round(post_sum$mean_nFeature, 1),
        post_mean_pct_mt = round(post_sum$mean_percent_mt, 3),
        post_median_pct_mt = round(post_sum$median_percent_mt, 3),
        cells_removed = n_before - n_after,
        pct_cells_removed = round(100 * (n_before - n_after) / n_before, 2),
        n_doublets = n_doublets, doublet_method = dbl_method,
        filter_used = filter_expr,
        source_path = gsm_dir, output_path = out_dir,
        processed_time = ts(), stringsAsFactors = FALSE))

      log_line(sprintf("    完成 %s: 細胞 %d -> %d (移除 %.1f%%), doublets=%s",
                       gsm, n_before, n_after,
                       100 * (n_before - n_after) / n_before,
                       ifelse(is.na(n_doublets), "NA", n_doublets)))
      "OK"
    },
    error = function(e) {
      ## 失敗: 記錄並寫一個 FAILED row，然後繼續下一個樣本
      msg <- conditionMessage(e)
      log_line(sprintf("    [FAILED] %s/%s: %s", gse_clean, gsm, msg))
      write_metrics_row(data.frame(
        GSE = gse_clean, GSE_folder_raw = gse_raw, GSM = gsm, condition = cond,
        status = "FAILED", error_note = msg,
        gse_has_hash_prefix = gse_hash,
        species = NA, mt_pattern = NA,
        pre_n_cells = NA, pre_n_genes = NA, pre_mean_nCount = NA, pre_mean_nFeature = NA,
        pre_mean_pct_mt = NA, pre_median_pct_mt = NA,
        post_n_cells = NA, post_n_genes = NA, post_mean_nCount = NA, post_mean_nFeature = NA,
        post_mean_pct_mt = NA, post_median_pct_mt = NA,
        cells_removed = NA, pct_cells_removed = NA,
        n_doublets = NA, doublet_method = NA,
        filter_used = NA, source_path = gsm_dir, output_path = out_dir,
        processed_time = ts(), stringsAsFactors = FALSE))
      "FAILED"
    })
    ## 釋放記憶體
    suppressWarnings(rm(list = intersect(c("obj","obj_post","mat"), ls())))
    gc(verbose = FALSE)
  } # end GSM loop
} # end GSE loop

## ----------------------------- 5. 寫結構掃描表 ---------------------------- ##
if (length(structure_rows)) {
  struct_df <- do.call(rbind, structure_rows)
  write_csv_utf8(struct_df, file.path(LOG_DIR, "structure_scan.csv"))
  log_line(sprintf("結構掃描表已寫出 (%d 列): %s",
                   nrow(struct_df), file.path(LOG_DIR, "structure_scan.csv")))
}

log_line("==== QC run 全部完成。請接著執行 02_Aggregate_Summary.R 產生大表 ====")


###############################################################################
## DoubletFinder / scDblFinder 包裝函式
## 放最後，因為主迴圈以 source 後函式已定義的順序執行沒問題 (R 會先 parse 全檔)
###############################################################################
run_doublet <- function(obj_post, det) {
  ## 標準前處理 (DoubletFinder 需要)
  obj <- NormalizeData(obj_post, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)

  n_cells <- ncol(obj)
  ## 依細胞數估 doublet 率 (10X: 約 0.8% / 1000 cells)
  exp_rate <- min(0.008 * (n_cells / 1000), 0.30)
  if (!is.finite(exp_rate) || exp_rate <= 0) exp_rate <- DEFAULT_DOUBLET_RATE

  if (HAVE_DF) {
    library(DoubletFinder)
    ## pK 掃描
    sweep <- paramSweep(obj, PCs = 1:20, sct = FALSE)
    sweep_stats <- summarizeSweep(sweep, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    pk_opt <- as.numeric(as.character(
      bcmvn$pK[which.max(bcmvn$BCmetric)]))
    ## 估 homotypic 並調整
    homotypic <- modelHomotypic(obj$seurat_clusters)
    nExp <- round(exp_rate * n_cells)
    nExp_adj <- round(nExp * (1 - homotypic))
    obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = pk_opt,
                         nExp = nExp_adj, sct = FALSE)
    df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
    obj$doublet_class <- obj@meta.data[[df_col]]
    n_db <- sum(obj$doublet_class == "Doublet", na.rm = TRUE)
    plt <- if (HAVE_GGPLOT)
      DimPlot(obj, group.by = "doublet_class") + ggtitle("DoubletFinder") else NULL
    return(list(obj = obj, n_doublets = n_db, method = "DoubletFinder", plot = plt))

  } else if (HAVE_SCDBL) {
    ## 備援: scDblFinder
    suppressPackageStartupMessages(library(scDblFinder))
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, dbr = exp_rate)
    obj$doublet_class <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
    n_db <- sum(obj$doublet_class == "Doublet", na.rm = TRUE)
    plt <- if (HAVE_GGPLOT)
      DimPlot(obj, group.by = "doublet_class") + ggtitle("scDblFinder") else NULL
    return(list(obj = obj, n_doublets = n_db, method = "scDblFinder", plot = plt))

  } else {
    stop("未安裝 DoubletFinder 或 scDblFinder，略過 doublet 偵測")
  }
}
