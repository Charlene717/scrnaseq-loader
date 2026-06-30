###############################################################################
## 01d_QC_H5.R
## ---------------------------------------------------------------------------
## scRNA-seq 廣泛品質評估 + 初步 QC —— 「10x .h5 (HDF5)」版本 (per-sample)
##
## 為什麼有這支:
##   01_QC_PerSample.R 吃 10X 三件套;01b 吃 dense matrix;01c 吃已整合 RDS。
##   這支吃 10x Cell Ranger 打包的 .h5 (filtered_feature_bc_matrix.h5),
##   每個 .h5 各自跑與 01 相同的 QC,輸出格式與 01/01b/01c 一致。
##
## .h5 是四種格式裡最乾淨的:Seurat 內建 Read10X_h5() 直接讀,不必猜分隔符/方向/前綴。
## 唯一要自動判斷的是:
##   - Read10X_h5() 回傳單一矩陣 -> 純單樣本 (最常見,如 filtered_feature_bc_matrix.h5)
##   - 回傳 list (多 feature type,如 Gene Expression + Antibody Capture/CITE-seq,
##              或多 genome) -> 自動取 'Gene Expression' 那塊做 scRNA QC
##   本腳本自動偵測回傳型別並處理 (已驗證單樣本 / 多type 兩種情況)。
##
## 資料結構 (來源):
##   <SRC_ROOT>/<GSE...>/  之下放一個或多個 .h5,例如:
##     GSM..._RUN_1_filtered_feature_bc_matrix.h5
##     GSM..._RUN_2_filtered_feature_bc_matrix.h5
##   (每個 .h5 = 一個樣本;檔名/GSM 當樣本名)
##
## 輸出 (與 01 完全相同):
##   <OUT_ROOT>/<GSE>/<sample>/
##     - <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv / plots/
##   _run_log/per_sample/<GSE>__<sample>.csv  (欄位與 01 一致)
###############################################################################

## ----------------------------- 0. 路徑設定 -------------------------------- ##
SRC_ROOT <- "Z:/Dataset_Online/scRNA-seq/h5"            # 放 .h5 的來源根目錄
OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"    # 與 01 共用輸出

## ------ 過濾門檻 (與 01 完全相同) ------ ##
MIN_FEATURE <- 200
MAX_FEATURE <- 5000
MAX_MT_PCT  <- 30
MIN_CELLS_GENE    <- 3
MIN_FEATURES_CELL <- 0

## ====== 使用者可調參數 (與 01 對齊) ====== ##
RUN_DOUBLETFINDER <- FALSE
IGNORE_FOLDERS <- c("GEO_10X_auto")
## 多 feature type 的 .h5,要取哪個 type 做 scRNA QC (預設 Gene Expression)
GEX_FEATURE_TYPE <- "Gene Expression"
## =============================== ##

SEED <- 42; set.seed(SEED)

## ----------------------------- 1. 套件 ------------------------------------ ##
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat 套件")
  library(Seurat); library(Matrix)
})
## Read10X_h5 需要 hdf5r 套件;沒裝會讀不了 .h5
if (!requireNamespace("hdf5r", quietly = TRUE))
  stop("讀 .h5 需要 hdf5r 套件,請先 install.packages('hdf5r')")

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAVE_GGPLOT <- has_pkg("ggplot2"); if (HAVE_GGPLOT) suppressPackageStartupMessages(library(ggplot2))
HAVE_PATCH  <- has_pkg("patchwork"); if (HAVE_PATCH) suppressPackageStartupMessages(library(patchwork))

## ----------------------------- 2. 工具函式 (與 01/01b/01c 一致) ----------- ##
ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive=TRUE, showWarnings=FALSE); invisible(path) }
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
write_csv_utf8 <- function(df, path) {
  if (requireNamespace("readr", quietly=TRUE)) { readr::write_excel_csv(df, path); return(invisible(path)) }
  lines <- utils::capture.output(utils::write.csv(df, file="", row.names=FALSE))
  con <- file(path, open="wb"); writeBin(charToRaw("\xEF\xBB\xBF"), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse="\n"), "\n"))), con); close(con); invisible(path)
}
LOG_DIR <- file.path(OUT_ROOT, "_run_log")
log_line <- function(msg) {
  ensure_dir(LOG_DIR); cat(sprintf("[%s] %s\n", ts(), msg))
  cat(sprintf("[%s] %s\n", ts(), msg), file=file.path(LOG_DIR, "run_log_h5.txt"), append=TRUE)
}
PERSAMPLE_DIR <- file.path(LOG_DIR, "per_sample"); ensure_dir(PERSAMPLE_DIR)
write_metrics_row <- function(row) {
  fn <- file.path(PERSAMPLE_DIR, sprintf("%s__%s.csv", row$GSE,
                  gsub("[^A-Za-z0-9_.-]", "_", as.character(row$GSM))))
  write_csv_utf8(row, fn)
}

get_counts_obj <- function(obj, assay=NULL) {
  a <- if (is.null(assay)) DefaultAssay(obj) else assay
  out <- tryCatch(SeuratObject::GetAssayData(obj, assay=a, layer="counts"), error=function(e) NULL)
  if (is.null(out)) out <- tryCatch(SeuratObject::GetAssayData(obj, assay=a, slot="counts"), error=function(e) NULL)
  if (is.null(out)) stop("無法取得 counts"); out
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
  md_ <- function(x) if (length(x)) median(x, na.rm=TRUE) else NA_real_
  data.frame(state=label, n_cells=ncol(obj),
    n_genes_detected=sum(Matrix::rowSums(get_counts_obj(obj) > 0) > 0),
    mean_nCount=sm(md$nCount_RNA), median_nCount=md_(md$nCount_RNA),
    mean_nFeature=sm(md$nFeature_RNA), median_nFeature=md_(md$nFeature_RNA),
    mean_percent_mt=sm(md$percent.mt), median_percent_mt=md_(md$percent.mt),
    mean_percent_ribo=if(!is.null(md$percent.ribo)) sm(md$percent.ribo) else NA_real_,
    stringsAsFactors=FALSE)
}
mad_bounds <- function(x, nmads=3, log=FALSE) {
  v <- if (log) log1p(x) else x
  med <- median(v, na.rm=TRUE); m <- mad(v, na.rm=TRUE)
  lo <- med-nmads*m; hi <- med+nmads*m
  if (log) { lo <- expm1(lo); hi <- expm1(hi) }; c(lower=lo, upper=hi)
}

## --- 核心:讀一個 .h5 -> 回傳 gene×cell sparse matrix ---
## 自動處理 Read10X_h5 回傳單矩陣 or list(多 feature type)的情況
read_h5_matrix <- function(path) {
  raw <- Read10X_h5(path)
  note <- "single matrix"
  if (is.list(raw) && !inherits(raw, c("dgCMatrix", "Matrix", "matrix"))) {
    ## 多 feature type / 多 genome -> 取 Gene Expression
    if (GEX_FEATURE_TYPE %in% names(raw)) {
      mat <- raw[[GEX_FEATURE_TYPE]]
      note <- sprintf("list[%s] -> 取'%s'", paste(names(raw), collapse=","), GEX_FEATURE_TYPE)
    } else {
      mat <- raw[[1]]
      note <- sprintf("list[%s] -> 無'%s',取第一個'%s'",
                      paste(names(raw), collapse=","), GEX_FEATURE_TYPE, names(raw)[1])
    }
  } else {
    mat <- raw
  }
  ## 基因/barcode 去重,轉 sparse
  rownames(mat) <- make.unique(as.character(rownames(mat)))
  colnames(mat) <- make.unique(as.character(colnames(mat)))
  mat <- as(mat, "CsparseMatrix")
  list(mat = mat, note = note)
}

## 找一個 GSE 資料夾下所有 .h5
find_h5_files <- function(gse_dir) {
  fs <- list.files(gse_dir, full.names=TRUE, recursive=FALSE)
  fs <- fs[!dir.exists(fs)]
  fs[grepl("\\.h5$", tolower(fs))]
}

## ----------------------------- 3. 掃描來源 -------------------------------- ##
log_line("==== .h5 per-sample QC 開始 ====")
log_line(sprintf("來源: %s", SRC_ROOT)); log_line(sprintf("輸出: %s", OUT_ROOT))
ensure_dir(OUT_ROOT); ensure_dir(LOG_DIR)
if (!dir.exists(SRC_ROOT)) stop(sprintf("找不到來源資料夾: %s", SRC_ROOT))

structure_rows <- list()
add_struct <- function(gse_raw, gse_clean, gsm, status, note, has_hash=FALSE, src_path="") {
  structure_rows[[length(structure_rows)+1]] <<- data.frame(
    GSE_folder_raw=gse_raw, GSE=gse_clean, GSM=gsm, structure_status=status,
    structure_note=note, gse_has_hash_prefix=has_hash, source_path=src_path, stringsAsFactors=FALSE)
}

top_items <- list.files(SRC_ROOT, full.names=TRUE)
gse_dirs <- character(0)
for (it in top_items) {
  base <- basename(it)
  if (base %in% IGNORE_FOLDERS) {
    add_struct(base, sub("^#+","",base), NA, "IGNORED_BY_USER", sprintf("使用者指定忽略(%s)", base), grepl("^#", base), it)
    log_line(sprintf("[ignore] 跳過資料夾: %s", base)); next
  }
  if (!dir.exists(it)) next
  gse_dirs <- c(gse_dirs, it)
}
log_line(sprintf("偵測到 %d 個 GSE 資料夾", length(gse_dirs)))

## ----------------------------- 4. 樣本處理函式 (QC 與 01 相同) ------------ ##
process_one_sample <- function(mat, gse_clean, gse_raw, gse_hash,
                               sample_name, cond, src_path, load_note) {
  out_dir   <- file.path(OUT_ROOT, gse_clean, sample_name)
  plots_dir <- file.path(out_dir, "plots")
  pre_rds   <- file.path(out_dir, paste0(sample_name, "_pre_QC.rds"))
  post_rds  <- file.path(out_dir, paste0(sample_name, "_post_QC.rds"))
  stats_csv <- file.path(out_dir, paste0(sample_name, "_QC_stats.csv"))
  expected_plots <- c(
    file.path(plots_dir, paste0(sample_name, "_violin_preQC.pdf")),
    file.path(plots_dir, paste0(sample_name, "_violin_postQC.pdf")),
    file.path(plots_dir, paste0(sample_name, "_scatter_preQC.pdf")),
    file.path(plots_dir, paste0(sample_name, "_beforeafter.pdf")))

  if (file.exists(pre_rds) && file.exists(post_rds) && file.exists(stats_csv) &&
      all(file.exists(expected_plots))) {
    log_line(sprintf("  [skip] %s/%s 已完整,跳過", gse_clean, sample_name)); return("SKIPPED")
  }
  ensure_dir(out_dir); ensure_dir(plots_dir)

  tryCatch({
    obj <- CreateSeuratObject(counts=mat, project=sample_name,
                              min.cells=MIN_CELLS_GENE, min.features=MIN_FEATURES_CELL)
    obj$orig.ident1 <- gse_clean; obj$orig.ident2 <- sample_name
    if (!is.na(cond)) obj$condition <- cond

    det <- detect_mt_pattern(rownames(obj))
    obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern=det$pattern)
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern=det$ribo)
    log_line(sprintf("    [%s] %s | 物種=%s MT=%s(%d) | 細胞=%d 基因=%d",
                     sample_name, load_note, det$species, det$pattern, det$n_mt, ncol(obj), nrow(obj)))

    pre_sum  <- summarise_obj(obj, "pre_QC")
    mad_feat <- mad_bounds(obj$nFeature_RNA, 3, TRUE)
    mad_cnt  <- mad_bounds(obj$nCount_RNA,   3, TRUE)
    mad_mt   <- mad_bounds(obj$percent.mt,   3, FALSE)
    saveRDS(obj, pre_rds)

    if (HAVE_GGPLOT) {
      feats <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), colnames(obj@meta.data))
      pdf(expected_plots[1], width=4*length(feats), height=5); print(VlnPlot(obj, features=feats, ncol=length(feats), pt.size=0.1)); dev.off()
      s1 <- FeatureScatter(obj,"nCount_RNA","nFeature_RNA")+ggtitle("nCount vs nFeature")
      s2 <- FeatureScatter(obj,"nCount_RNA","percent.mt")+ggtitle("nCount vs percent.mt")
      pdf(expected_plots[3], width=12, height=5); if (HAVE_PATCH) print(s1+s2) else { print(s1); print(s2) }; dev.off()
    }

    filter_expr <- sprintf("nFeature_RNA > %d & nFeature_RNA < %d & percent.mt < %d", MIN_FEATURE, MAX_FEATURE, MAX_MT_PCT)
    n_before <- ncol(obj)
    obj_post <- subset(obj, subset = nFeature_RNA > MIN_FEATURE & nFeature_RNA < MAX_FEATURE & percent.mt < MAX_MT_PCT)
    n_after <- ncol(obj_post)
    post_sum <- summarise_obj(obj_post, "post_QC")

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

    stats_df <- rbind(pre_sum, post_sum)
    stats_df$GSE <- gse_clean; stats_df$GSM <- sample_name; stats_df$condition <- cond
    stats_df$species <- det$species; stats_df$mt_pattern <- det$pattern; stats_df$n_mt_genes <- det$n_mt
    stats_df$filter_used <- filter_expr
    stats_df$MAD_nFeature_lower <- round(mad_feat["lower"],1); stats_df$MAD_nFeature_upper <- round(mad_feat["upper"],1)
    stats_df$MAD_nCount_lower <- round(mad_cnt["lower"],1);    stats_df$MAD_nCount_upper <- round(mad_cnt["upper"],1)
    stats_df$MAD_percent_mt_upper <- round(mad_mt["upper"],2)
    stats_df$cells_removed <- n_before-n_after
    stats_df$pct_cells_removed <- if (n_before>0) round(100*(n_before-n_after)/n_before,2) else NA
    stats_df$n_doublets <- NA_integer_; stats_df$doublet_method <- if (!RUN_DOUBLETFINDER) "skipped(user)" else "none"
    stats_df$load_note <- load_note; stats_df$processed_time <- ts()
    write_csv_utf8(stats_df, stats_csv)

    write_metrics_row(data.frame(
      GSE=gse_clean, GSE_folder_raw=gse_raw, GSM=sample_name, condition=cond,
      status="OK", error_note="", gse_has_hash_prefix=gse_hash,
      species=det$species, mt_pattern=det$pattern,
      pre_n_cells=pre_sum$n_cells, pre_n_genes=pre_sum$n_genes_detected,
      pre_mean_nCount=round(pre_sum$mean_nCount,1), pre_mean_nFeature=round(pre_sum$mean_nFeature,1),
      pre_mean_pct_mt=round(pre_sum$mean_percent_mt,3), pre_median_pct_mt=round(pre_sum$median_percent_mt,3),
      post_n_cells=post_sum$n_cells, post_n_genes=post_sum$n_genes_detected,
      post_mean_nCount=round(post_sum$mean_nCount,1), post_mean_nFeature=round(post_sum$mean_nFeature,1),
      post_mean_pct_mt=round(post_sum$mean_percent_mt,3), post_median_pct_mt=round(post_sum$median_percent_mt,3),
      cells_removed=n_before-n_after,
      pct_cells_removed=if (n_before>0) round(100*(n_before-n_after)/n_before,2) else NA,
      n_doublets=NA, doublet_method=if (!RUN_DOUBLETFINDER) "skipped(user)" else "none",
      filter_used=filter_expr, source_path=src_path, output_path=out_dir,
      processed_time=ts(), stringsAsFactors=FALSE))

    log_line(sprintf("    完成 %s: 細胞 %d -> %d (移除 %.1f%%)",
                     sample_name, n_before, n_after, if (n_before>0) 100*(n_before-n_after)/n_before else 0))
    "OK"
  }, error=function(e) {
    msg <- conditionMessage(e)
    log_line(sprintf("    [FAILED] %s/%s: %s", gse_clean, sample_name, msg))
    write_metrics_row(data.frame(
      GSE=gse_clean, GSE_folder_raw=gse_raw, GSM=sample_name, condition=cond,
      status="FAILED", error_note=msg, gse_has_hash_prefix=gse_hash, species=NA, mt_pattern=NA,
      pre_n_cells=NA,pre_n_genes=NA,pre_mean_nCount=NA,pre_mean_nFeature=NA,pre_mean_pct_mt=NA,pre_median_pct_mt=NA,
      post_n_cells=NA,post_n_genes=NA,post_mean_nCount=NA,post_mean_nFeature=NA,post_mean_pct_mt=NA,post_median_pct_mt=NA,
      cells_removed=NA,pct_cells_removed=NA,n_doublets=NA,doublet_method=NA,filter_used=NA,
      source_path=src_path, output_path=out_dir, processed_time=ts(), stringsAsFactors=FALSE))
    "FAILED"
  })
}

## ----------------------------- 5. 主迴圈 ---------------------------------- ##
for (gse_dir in gse_dirs) {
  gse_raw <- basename(gse_dir); gse_hash <- grepl("^#", gse_raw); gse_clean <- sub("^#+", "", gse_raw)
  if (gse_hash) log_line(sprintf("[!] GSE '%s' 有 # 前綴,仍嘗試處理", gse_raw))

  h5s <- find_h5_files(gse_dir)
  if (length(h5s) == 0) {
    add_struct(gse_raw, gse_clean, NA, "NO_H5_FILE", "GSE 內找不到 .h5 檔", gse_hash, gse_dir)
    log_line(sprintf("[!] GSE '%s' 無 .h5 檔,跳過", gse_raw)); next
  }

  for (h5 in h5s) {
    fbase <- basename(h5)
    ## 樣本名:檔名去 .h5 與常見後綴
    sample_name <- sub("\\.h5$", "", fbase, ignore.case=TRUE)
    sample_name <- sub("(_filtered_feature_bc_matrix|_raw_feature_bc_matrix|_filtered_gene_bc_matrices|_feature_bc_matrix)$",
                       "", sample_name, ignore.case=TRUE)
    ## condition 由檔名推測
    cond <- if (grepl("keloid|_KL|_K[0-9]", sample_name, ignore.case=TRUE)) "Keloid"
            else if (grepl("normal|_NS|nskin|scar", sample_name, ignore.case=TRUE)) "NormalSkin" else NA

    log_line(sprintf("==> 讀取 %s / %s", gse_clean, fbase))
    loaded <- tryCatch(read_h5_matrix(h5), error=function(e) {
      log_line(sprintf("  [FAILED] 讀取 .h5 失敗 %s: %s", fbase, conditionMessage(e)))
      add_struct(gse_raw, gse_clean, fbase, "READ_FAILED", conditionMessage(e), gse_hash, h5); NULL
    })
    if (is.null(loaded)) next
    add_struct(gse_raw, gse_clean, fbase, "OK", loaded$note, gse_hash, h5)

    process_one_sample(loaded$mat, gse_clean, gse_raw, gse_hash,
                       sample_name=sample_name, cond=cond, src_path=h5, load_note=loaded$note)
    suppressWarnings(rm(loaded)); gc(verbose=FALSE)
  }
}

## ----------------------------- 6. 結構掃描表 ------------------------------ ##
if (length(structure_rows)) {
  struct_df <- do.call(rbind, structure_rows)
  write_csv_utf8(struct_df, file.path(LOG_DIR, "structure_scan_h5.csv"))
  log_line(sprintf("結構掃描表已寫出 (%d 列)", nrow(struct_df)))
}
log_line("==== .h5 per-sample QC 全部完成。可接 02_Aggregate_Summary.R ====")
