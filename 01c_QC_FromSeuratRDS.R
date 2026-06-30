###############################################################################
## 01c_QC_FromSeuratRDS.R
## ---------------------------------------------------------------------------
## scRNA-seq 廣泛品質評估 + 初步 QC —— 「已整合 Seurat RDS」版本 (per-sample)
##
## 為什麼有這支:
##   01_QC_PerSample.R 吃 10X 三件套;01b 吃 dense matrix 文字檔。
##   但有些資料集你拿到的是「已經整合 + 註解好的單一大 Seurat 物件 (.rds)」,
##   裡面用一個 meta.data 欄位 (預設 'sample') 區分多個樣本。
##   這支把該物件「依 sample 欄拆開」,每個 sample 各自跑與 01 相同的 QC,
##   輸出格式與 01/01b 一致 -> 02_Aggregate_Summary.R 可一起彙整。
##
## ★ 這個物件的三個特殊處理 (重要):
##   (1) percent.mito 可能已存在但為「比例(0~1)」而非百分比 -> 本腳本一律
##       用 PercentageFeatureSet 重算 percent.mt (0~100),與 01 一致,不沿用舊欄。
##   (2) 這是「已 QC + 已整合」的資料,RNA counts 是先前處理過的 counts。
##       再套門檻 = 第二次過濾。故過濾以開關 APPLY_FILTER 控制 (預設 TRUE)。
##       若你只想要 per-sample 的 QC 統計 + 圖、不再砍細胞 -> 設 FALSE。
##   (3) 記憶體:用「切 counts 矩陣」而非 SplitObject (省記憶體),逐 sample + gc。
##       並相容 Seurat v4/v5 (counts 取法不同)。
##
## 輸出 (與 01 完全相同):
##   <OUT_ROOT>/<GSE>/<sample>/
##     - <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv / plots/
##   _run_log/per_sample/<GSE>__<sample>.csv  (欄位與 01 一致)
###############################################################################

## ----------------------------- 0. 路徑與輸入 ------------------------------ ##
## 來源 RDS (已整合的大物件)。若該物件已經在環境中 (變數名 SeuratObject),
## 本腳本不會重複讀取,直接用既有的,省下再讀 13GB 的時間。
RDS_PATH   <- "X:/Dataset_Online/##_Keloid/scRNA-seq/rds/GSE307504/Integration_0427_seurat_annotated.rds"
OUT_ROOT   <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"   # 與 01 共用輸出
GSE_LABEL  <- "GSE307504"        # 輸出資料夾用的 GSE 名 (會去掉開頭非英數字元)
SAMPLE_COL <- "sample"           # 用哪個 meta.data 欄位區分樣本
RNA_ASSAY  <- "RNA"              # QC 要用的 assay (原始 counts 所在)

## ------ 過濾門檻 (與 01 完全相同) ------ ##
MIN_FEATURE <- 200
MAX_FEATURE <- 5000
MAX_MT_PCT  <- 30
MIN_CELLS_GENE    <- 3
MIN_FEATURES_CELL <- 0

## ====== 使用者可調參數 ====== ##
## 是否實際套用過濾 (砍細胞)。
##   TRUE  -> 與 01 同理,產生 pre/post 並依門檻過濾 (注意這是對已處理資料的「第二次過濾」)
##   FALSE -> 只算 QC 統計 + 畫圖,不砍任何細胞 (post = pre);適合只想看 per-sample 品質分佈
APPLY_FILTER <- TRUE

## 是否把既有註解 (如 celltype) 帶進每個 per-sample 物件 (方便日後參照;不影響 QC)
CARRY_META <- c("celltype", "condition", "subj", "age")

## 條件 (condition) 由 sample 名推測的規則 (此資料集: Skin/Ke/Hyper/NScar)
guess_condition <- function(s) {
  if (grepl("^Ke",    s, ignore.case = TRUE)) return("Keloid")
  if (grepl("^Hyper", s, ignore.case = TRUE)) return("HypertrophicScar")
  if (grepl("^NScar", s, ignore.case = TRUE)) return("NormalScar")
  if (grepl("^Skin",  s, ignore.case = TRUE)) return("NormalSkin")
  NA_character_
}
## ============================ ##

SEED <- 42; set.seed(SEED)

## ----------------------------- 1. 套件 ------------------------------------ ##
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("需要 Seurat 套件")
  library(Seurat); library(Matrix)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAVE_GGPLOT <- has_pkg("ggplot2"); if (HAVE_GGPLOT) suppressPackageStartupMessages(library(ggplot2))
HAVE_PATCH  <- has_pkg("patchwork"); if (HAVE_PATCH) suppressPackageStartupMessages(library(patchwork))

## ----------------------------- 2. 工具函式 (與 01/01b 一致) --------------- ##
ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
write_csv_utf8 <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) { readr::write_excel_csv(df, path); return(invisible(path)) }
  lines <- utils::capture.output(utils::write.csv(df, file = "", row.names = FALSE))
  con <- file(path, open = "wb"); writeBin(charToRaw("\xEF\xBB\xBF"), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con); close(con); invisible(path)
}
LOG_DIR <- file.path(OUT_ROOT, "_run_log")
log_line <- function(msg) {
  ensure_dir(LOG_DIR); cat(sprintf("[%s] %s\n", ts(), msg))
  cat(sprintf("[%s] %s\n", ts(), msg), file = file.path(LOG_DIR, "run_log_fromrds.txt"), append = TRUE)
}
PERSAMPLE_DIR <- file.path(LOG_DIR, "per_sample"); ensure_dir(PERSAMPLE_DIR)
write_metrics_row <- function(row) {
  fn <- file.path(PERSAMPLE_DIR, sprintf("%s__%s.csv", row$GSE,
                  gsub("[^A-Za-z0-9_.-]", "_", as.character(row$GSM))))
  write_csv_utf8(row, fn)
}

detect_mt_pattern <- function(genes) {
  n_h <- sum(grepl("^MT-", genes)); n_m <- sum(grepl("^mt-", genes))
  if (n_h >= n_m && n_h > 0) list(pattern="^MT-", species="human", n_mt=n_h, ribo="^RP[SL]")
  else if (n_m > 0)          list(pattern="^mt-", species="mouse", n_mt=n_m, ribo="^Rp[sl]")
  else                       list(pattern="^MT-", species="unknown(no MT genes)", n_mt=0, ribo="^RP[SL]")
}
summarise_obj <- function(obj, label) {
  md <- obj@meta.data
  sm <- function(x) if (length(x)) mean(x, na.rm=TRUE) else NA_real_
  sd_ <- function(x) if (length(x)) median(x, na.rm=TRUE) else NA_real_
  data.frame(state=label, n_cells=ncol(obj),
    n_genes_detected=sum(Matrix::rowSums(get_counts_obj(obj) > 0) > 0),
    mean_nCount=sm(md$nCount_RNA), median_nCount=sd_(md$nCount_RNA),
    mean_nFeature=sm(md$nFeature_RNA), median_nFeature=sd_(md$nFeature_RNA),
    mean_percent_mt=sm(md$percent.mt), median_percent_mt=sd_(md$percent.mt),
    mean_percent_ribo=if(!is.null(md$percent.ribo)) sm(md$percent.ribo) else NA_real_,
    stringsAsFactors=FALSE)
}
mad_bounds <- function(x, nmads=3, log=FALSE) {
  v <- if (log) log1p(x) else x
  med <- median(v, na.rm=TRUE); m <- mad(v, na.rm=TRUE)
  lo <- med-nmads*m; hi <- med+nmads*m
  if (log) { lo <- expm1(lo); hi <- expm1(hi) }; c(lower=lo, upper=hi)
}

## --- Seurat v4/v5 相容: 取某物件的 counts (sparse) ---
get_counts_obj <- function(obj, assay = NULL) {
  a <- if (is.null(assay)) DefaultAssay(obj) else assay
  out <- tryCatch(SeuratObject::GetAssayData(obj, assay = a, layer = "counts"),
                  error = function(e) NULL)                 # v5
  if (is.null(out))
    out <- tryCatch(SeuratObject::GetAssayData(obj, assay = a, slot = "counts"),
                    error = function(e) NULL)               # v4
  if (is.null(out)) stop("無法取得 counts (assay=", a, ")")
  out
}

## ----------------------------- 3. 載入大物件 ------------------------------ ##
log_line("==== From-RDS per-sample QC 開始 ====")
ensure_dir(OUT_ROOT); ensure_dir(LOG_DIR)

if (!exists("SeuratObject")) {
  log_line(sprintf("讀取 RDS: %s", RDS_PATH))
  if (!file.exists(RDS_PATH)) stop(sprintf("找不到 RDS: %s", RDS_PATH))
  SeuratObject <- readRDS(RDS_PATH)
} else {
  log_line("偵測到環境中已有 SeuratObject,直接使用 (不重讀)")
}
obj_all <- SeuratObject

## 檢查 assay 與 sample 欄
if (!(RNA_ASSAY %in% Assays(obj_all))) {
  log_line(sprintf("[!] 找不到 assay '%s',改用 DefaultAssay='%s' (注意:若非 RNA 原始 counts,QC 結果可能失真)",
                   RNA_ASSAY, DefaultAssay(obj_all)))
  RNA_ASSAY <- DefaultAssay(obj_all)
}
if (!(SAMPLE_COL %in% colnames(obj_all@meta.data)))
  stop(sprintf("meta.data 沒有欄位 '%s';現有欄位: %s",
               SAMPLE_COL, paste(colnames(obj_all@meta.data), collapse=", ")))

## v5 若 RNA assay 的 layer 被拆分 (counts.1/counts.2...),先 JoinLayers 合併
obj_all <- tryCatch({
  if (inherits(obj_all[[RNA_ASSAY]], "Assay5")) JoinLayers(obj_all, assay = RNA_ASSAY) else obj_all
}, error = function(e) { log_line(sprintf("JoinLayers 略過: %s", conditionMessage(e))); obj_all })

## --- 一次取出 RNA counts (sparse) + sample 向量 + 要攜帶的 meta ---
log_line(sprintf("擷取 assay '%s' 的 counts ...", RNA_ASSAY))
counts_all <- get_counts_obj(obj_all, assay = RNA_ASSAY)
sample_vec <- as.character(obj_all@meta.data[[SAMPLE_COL]])
names(sample_vec) <- rownames(obj_all@meta.data)

carry_df <- NULL
carry_cols <- intersect(CARRY_META, colnames(obj_all@meta.data))
if (length(carry_cols)) carry_df <- obj_all@meta.data[, carry_cols, drop = FALSE]

GSE_CLEAN <- sub("^[^A-Za-z0-9]+", "", GSE_LABEL)
samples <- sort(unique(sample_vec[!is.na(sample_vec) & sample_vec != ""]))
log_line(sprintf("偵測到 %d 個樣本: %s", length(samples), paste(samples, collapse=", ")))
log_line(sprintf("APPLY_FILTER = %s (門檻 nFeature %d~%d, percent.mt < %d)",
                 APPLY_FILTER, MIN_FEATURE, MAX_FEATURE, MAX_MT_PCT))

## 記憶體提示:大物件仍在環境中。若記憶體吃緊,可在跑前手動釋放 integrated assay:
##   DefaultAssay(SeuratObject) <- "RNA"; SeuratObject[["integrated"]] <- NULL; gc()

## ----------------------------- 4. 逐樣本 QC ------------------------------- ##
structure_rows <- list()

for (sname in samples) {
  out_dir   <- file.path(OUT_ROOT, GSE_CLEAN, sname)
  plots_dir <- file.path(out_dir, "plots")
  pre_rds   <- file.path(out_dir, paste0(sname, "_pre_QC.rds"))
  post_rds  <- file.path(out_dir, paste0(sname, "_post_QC.rds"))
  stats_csv <- file.path(out_dir, paste0(sname, "_QC_stats.csv"))
  expected_plots <- c(
    file.path(plots_dir, paste0(sname, "_violin_preQC.pdf")),
    file.path(plots_dir, paste0(sname, "_violin_postQC.pdf")),
    file.path(plots_dir, paste0(sname, "_scatter_preQC.pdf")),
    file.path(plots_dir, paste0(sname, "_beforeafter.pdf")))

  ## 續跑檢查 (與 01 相同)
  if (file.exists(pre_rds) && file.exists(post_rds) && file.exists(stats_csv) &&
      all(file.exists(expected_plots))) {
    log_line(sprintf("[skip] %s/%s 已完整,跳過", GSE_CLEAN, sname)); next
  }
  ensure_dir(out_dir); ensure_dir(plots_dir)
  cond <- guess_condition(sname)
  log_line(sprintf("==> 處理 %s / %s (condition=%s)", GSE_CLEAN, sname, ifelse(is.na(cond),"NA",cond)))

  res <- tryCatch({
    ## (a) 切出此樣本的 counts 欄
    cells <- names(sample_vec)[which(sample_vec == sname)]
    submat <- counts_all[, cells, drop = FALSE]
    ## 丟掉全 0 的基因 (此樣本沒表現的),避免 PercentageFeatureSet 受影響
    keep_g <- Matrix::rowSums(submat) > 0
    submat <- submat[keep_g, , drop = FALSE]

    ## (b) 建立乾淨的 per-sample 物件 (不帶整合結構,輕量)
    obj <- CreateSeuratObject(counts = submat, project = sname,
                              min.cells = MIN_CELLS_GENE, min.features = MIN_FEATURES_CELL)
    obj$orig.ident1 <- GSE_CLEAN; obj$orig.ident2 <- sname
    if (!is.na(cond)) obj$condition <- cond
    ## 攜帶既有 meta (celltype 等),僅作參照
    if (!is.null(carry_df)) {
      cc <- intersect(colnames(obj), rownames(carry_df))
      for (col in colnames(carry_df)) {
        if (col == "condition" && !is.na(cond)) next   # 已用推測值
        obj@meta.data[cc, col] <- carry_df[cc, col]
      }
    }

    ## (c) 物種 / MT pattern 偵測 + 重算 percent.mt (0~100),與 01 一致
    det <- detect_mt_pattern(rownames(obj))
    obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = det$pattern)
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = det$ribo)
    log_line(sprintf("    物種=%s MT=%s(%d) | 細胞=%d 基因=%d | mean percent.mt=%.2f%%",
                     det$species, det$pattern, det$n_mt, ncol(obj), nrow(obj),
                     mean(obj$percent.mt, na.rm=TRUE)))

    ## (d) QC 前摘要 + MAD (僅報告)
    pre_sum  <- summarise_obj(obj, "pre_QC")
    mad_feat <- mad_bounds(obj$nFeature_RNA, 3, TRUE)
    mad_cnt  <- mad_bounds(obj$nCount_RNA,   3, TRUE)
    mad_mt   <- mad_bounds(obj$percent.mt,   3, FALSE)
    saveRDS(obj, pre_rds)

    ## (e) 繪圖 (前)
    if (HAVE_GGPLOT) {
      feats <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), colnames(obj@meta.data))
      pdf(expected_plots[1], width=4*length(feats), height=5); print(VlnPlot(obj, features=feats, ncol=length(feats), pt.size=0.1)); dev.off()
      s1 <- FeatureScatter(obj,"nCount_RNA","nFeature_RNA")+ggtitle("nCount vs nFeature")
      s2 <- FeatureScatter(obj,"nCount_RNA","percent.mt")+ggtitle("nCount vs percent.mt")
      pdf(expected_plots[3], width=12, height=5); if (HAVE_PATCH) print(s1+s2) else { print(s1); print(s2) }; dev.off()
    }

    ## (f) 過濾 (依 APPLY_FILTER)
    n_before <- ncol(obj)
    if (APPLY_FILTER) {
      filter_expr <- sprintf("nFeature_RNA > %d & nFeature_RNA < %d & percent.mt < %d",
                             MIN_FEATURE, MAX_FEATURE, MAX_MT_PCT)
      obj_post <- subset(obj, subset = nFeature_RNA > MIN_FEATURE &
                                       nFeature_RNA < MAX_FEATURE & percent.mt < MAX_MT_PCT)
    } else {
      filter_expr <- "none (stats only, APPLY_FILTER=FALSE)"
      obj_post <- obj
    }
    n_after <- ncol(obj_post)
    post_sum <- summarise_obj(obj_post, "post_QC")

    ## (g) 繪圖 (後 + before/after)
    if (HAVE_GGPLOT) {
      feats <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), colnames(obj_post@meta.data))
      pdf(expected_plots[2], width=4*length(feats), height=5); print(VlnPlot(obj_post, features=feats, ncol=length(feats), pt.size=0.1)); dev.off()
      ba <- rbind(
        data.frame(stage="pre",  nFeature=obj$nFeature_RNA, nCount=obj$nCount_RNA, percent.mt=obj$percent.mt),
        data.frame(stage="post", nFeature=obj_post$nFeature_RNA, nCount=obj_post$nCount_RNA, percent.mt=obj_post$percent.mt))
      g1 <- ggplot(ba,aes(stage,nFeature,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+ggtitle("nFeature_RNA")
      g2 <- ggplot(ba,aes(stage,nCount,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+scale_y_log10()+ggtitle("nCount_RNA (log10)")
      g3 <- ggplot(ba,aes(stage,percent.mt,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+ggtitle("percent.mt")
      pdf(expected_plots[4], width=12, height=5); if (HAVE_PATCH) print(g1+g2+g3) else { print(g1);print(g2);print(g3) }; dev.off()
    }
    saveRDS(obj_post, post_rds)

    ## (h) stats CSV (與 01 相同欄位 + load_note)
    stats_df <- rbind(pre_sum, post_sum)
    stats_df$GSE <- GSE_CLEAN; stats_df$GSM <- sname; stats_df$condition <- cond
    stats_df$species <- det$species; stats_df$mt_pattern <- det$pattern; stats_df$n_mt_genes <- det$n_mt
    stats_df$filter_used <- filter_expr
    stats_df$MAD_nFeature_lower <- round(mad_feat["lower"],1); stats_df$MAD_nFeature_upper <- round(mad_feat["upper"],1)
    stats_df$MAD_nCount_lower <- round(mad_cnt["lower"],1);    stats_df$MAD_nCount_upper <- round(mad_cnt["upper"],1)
    stats_df$MAD_percent_mt_upper <- round(mad_mt["upper"],2)
    stats_df$cells_removed <- n_before-n_after
    stats_df$pct_cells_removed <- if (n_before>0) round(100*(n_before-n_after)/n_before,2) else NA
    stats_df$n_doublets <- NA_integer_; stats_df$doublet_method <- "skipped(from-rds)"
    stats_df$source_assay <- RNA_ASSAY; stats_df$processed_time <- ts()
    write_csv_utf8(stats_df, stats_csv)

    ## (i) per-sample metrics row (大表用,欄位與 01 一致)
    write_metrics_row(data.frame(
      GSE=GSE_CLEAN, GSE_folder_raw=GSE_LABEL, GSM=sname, condition=cond,
      status="OK", error_note="", gse_has_hash_prefix=grepl("^#", GSE_LABEL),
      species=det$species, mt_pattern=det$pattern,
      pre_n_cells=pre_sum$n_cells, pre_n_genes=pre_sum$n_genes_detected,
      pre_mean_nCount=round(pre_sum$mean_nCount,1), pre_mean_nFeature=round(pre_sum$mean_nFeature,1),
      pre_mean_pct_mt=round(pre_sum$mean_percent_mt,3), pre_median_pct_mt=round(pre_sum$median_percent_mt,3),
      post_n_cells=post_sum$n_cells, post_n_genes=post_sum$n_genes_detected,
      post_mean_nCount=round(post_sum$mean_nCount,1), post_mean_nFeature=round(post_sum$mean_nFeature,1),
      post_mean_pct_mt=round(post_sum$mean_percent_mt,3), post_median_pct_mt=round(post_sum$median_percent_mt,3),
      cells_removed=n_before-n_after,
      pct_cells_removed=if (n_before>0) round(100*(n_before-n_after)/n_before,2) else NA,
      n_doublets=NA, doublet_method="skipped(from-rds)", filter_used=filter_expr,
      source_path=RDS_PATH, output_path=out_dir, processed_time=ts(), stringsAsFactors=FALSE))

    log_line(sprintf("    完成 %s: 細胞 %d -> %d (移除 %.1f%%)",
                     sname, n_before, n_after,
                     if (n_before>0) 100*(n_before-n_after)/n_before else 0))
    "OK"
  }, error = function(e) {
    msg <- conditionMessage(e)
    log_line(sprintf("    [FAILED] %s/%s: %s", GSE_CLEAN, sname, msg))
    write_metrics_row(data.frame(
      GSE=GSE_CLEAN, GSE_folder_raw=GSE_LABEL, GSM=sname, condition=guess_condition(sname),
      status="FAILED", error_note=msg, gse_has_hash_prefix=grepl("^#", GSE_LABEL),
      species=NA, mt_pattern=NA,
      pre_n_cells=NA,pre_n_genes=NA,pre_mean_nCount=NA,pre_mean_nFeature=NA,pre_mean_pct_mt=NA,pre_median_pct_mt=NA,
      post_n_cells=NA,post_n_genes=NA,post_mean_nCount=NA,post_mean_nFeature=NA,post_mean_pct_mt=NA,post_median_pct_mt=NA,
      cells_removed=NA,pct_cells_removed=NA,n_doublets=NA,doublet_method=NA,filter_used=NA,
      source_path=RDS_PATH, output_path=out_dir, processed_time=ts(), stringsAsFactors=FALSE))
    "FAILED"
  })

  ## 釋放此樣本暫存
  suppressWarnings(rm(list = intersect(c("obj","obj_post","submat"), ls())))
  gc(verbose = FALSE)
}

log_line("==== From-RDS per-sample QC 全部完成。可接 02_Aggregate_Summary.R ====")
