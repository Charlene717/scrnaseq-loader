###############################################################################
## 01b_QC_DenseMatrix.R
## ---------------------------------------------------------------------------
## scRNA-seq 廣泛品質評估 + 初步 QC —— 「dense matrix 文字檔」版本
##
## 為什麼有這支:
##   01_QC_PerSample.R 只認 10X 三件套(barcodes/features/matrix)。
##   但很多 GEO 資料集提供的是「dense matrix 純文字檔」(gene × cell 的大表),
##   格式還五花八門。這支腳本是 01 的姊妹版:
##     - 讀取層:自動偵測各種 dense matrix 格式(見下),轉成 sparse matrix
##     - QC 層 :與 01 完全相同(門檻、物種偵測、輸出、續跑、錯誤跳過)
##   兩者輸出格式一致,02_Aggregate_Summary.R 可同時吃 01 與 01b 的結果。
##
## 自動偵測能力(已對 5 種真實 GEO 格式驗證):
##   1. 分隔符      : Tab 或逗號(.tsv/.txt/.csv 皆可)
##   2. 矩陣方向    : gene×cell 或 cell×gene(後者自動轉置)
##   3. 基因名前綴  : hg19_ / hg38_ / mm10_ 等(自動剝除)
##   4. 非 raw counts: 檔名含 rpkm/fpkm/tpm/cpm/normalized -> 自動跳過
##   5. 合併大檔    : barcode 帶樣本前綴 (KL14_ACGT...) -> 依前綴拆成多樣本各自 QC
##
## 資料結構(來源):
##   <SRC_ROOT>/<GSE...>/  之下放一個或多個 matrix 檔,例如:
##     - GSM..._counts_matrix.txt   (每個 GSM 一檔)
##     - GSE..._read.txt            (一檔含整個 GSE)
##     - GSE..._all.UMI.matrix.csv  (合併檔,需依 barcode 前綴拆樣本)
##
## 輸出(目的地)：與 01 完全相同的結構
##   <OUT_ROOT>/<GSE>/<sample>/
##     - <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv / plots/
##   _run_log/per_sample/<GSE>__<sample>.csv  (供 02 彙整,欄位與 01 一致)
###############################################################################

## ----------------------------- 0. 路徑設定 -------------------------------- ##
SRC_ROOT <- "X:/Dataset_Online/##_Keloid/scRNA-seq/matrix"      # dense matrix 來源根目錄
OUT_ROOT <- "X:/Dataset_Online/##_Keloid/scRNA-seq/matrix_QC"   # QC 輸出根目錄

## ------ 過濾門檻(與 01 完全相同) ------ ##
MIN_FEATURE <- 200
MAX_FEATURE <- 5000
MAX_MT_PCT  <- 30
MIN_CELLS_GENE <- 3
MIN_FEATURES_CELL <- 0

## ====== 使用者可調參數(與 01 對齊) ====== ##
RUN_DOUBLETFINDER <- FALSE
IGNORE_FOLDERS <- c("GEO_10X_auto")

## ------ dense matrix 專屬設定 ------ ##
## 視為「非 raw counts、要跳過」的檔名關鍵字(正規化資料)
SKIP_PATTERNS <- c("rpkm", "fpkm", "tpm", "cpm", "normalized", "_norm", "scaled")
## 基因名前綴候選(偵測到會剝除)
GENE_PREFIXES <- c("hg19_", "hg38_", "grch38_", "grcm38_", "mm10_", "mm39_", "GRCh38_")
## 合併檔:同一前綴下至少要有幾個 barcode 才認定為「一個樣本」(避免雜訊)
MIN_CELLS_PER_SPLIT <- 20
## =============================== ##

DEFAULT_DOUBLET_RATE <- 0.075
SEED <- 42

## ----------------------------- 1. 套件 ------------------------------------ ##
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("需要 Seurat 套件，請先 install.packages('Seurat')")
  library(Seurat)
  library(Matrix)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAVE_GGPLOT  <- has_pkg("ggplot2");  if (HAVE_GGPLOT) suppressPackageStartupMessages(library(ggplot2))
HAVE_PATCH   <- has_pkg("patchwork");if (HAVE_PATCH)  suppressPackageStartupMessages(library(patchwork))
HAVE_DATATABLE <- has_pkg("data.table")  # 大檔讀取加速(強烈建議安裝)
HAVE_DF      <- has_pkg("DoubletFinder")
HAVE_SCDBL   <- has_pkg("scDblFinder")
set.seed(SEED)

## ----------------------------- 2. 工具函式 -------------------------------- ##
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write_csv_utf8 <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_excel_csv(df, path); return(invisible(path))
  }
  lines <- utils::capture.output(utils::write.csv(df, file = "", row.names = FALSE))
  con <- file(path, open = "wb")
  writeBin(charToRaw("\xEF\xBB\xBF"), con)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
  close(con); invisible(path)
}

LOG_DIR <- file.path(OUT_ROOT, "_run_log")
log_line <- function(msg) {
  ensure_dir(LOG_DIR)
  cat(sprintf("[%s] %s\n", ts(), msg))
  cat(sprintf("[%s] %s\n", ts(), msg),
      file = file.path(LOG_DIR, "run_log_densematrix.txt"), append = TRUE)
}

## ====================================================================== ##
## === 自動偵測函式群(演算法已用 Python 對 7 種案例驗證,1:1 對應)    === ##
## ====================================================================== ##

## (A) 偵測分隔符:看第一行 tab 多還是逗號多
detect_sep <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  on.exit(close(con))
  first <- readLines(con, n = 1L)
  if (length(first) == 0) return("\t")
  n_tab <- lengths(regmatches(first, gregexpr("\t", first)))
  n_com <- lengths(regmatches(first, gregexpr(",",  first)))
  if (n_tab >= n_com) "\t" else ","
}

## (B) 是否該跳過(正規化檔,非 raw counts)
should_skip_file <- function(path) {
  name <- tolower(basename(path))
  for (pat in SKIP_PATTERNS) if (grepl(pat, name, fixed = TRUE)) return(TRUE)
  FALSE
}

## (C) 某串 token 像 barcode 的比例(純 ACGT >=8,可含 .1/-1,或帶樣本前綴)
frac_barcode <- function(tokens) {
  if (length(tokens) == 0) return(0)
  pat1 <- "^[ACGTacgt]{8,}([.-][0-9]+)?$"
  pat2 <- "^[A-Za-z0-9]+_[ACGTacgt]{8,}([.-][0-9]+)?$"  # 帶樣本前綴
  hit <- sum(grepl(pat1, tokens) | grepl(pat2, tokens))
  hit / length(tokens)
}

## (D) 某串 token 像 gene symbol 的比例(後備判斷方向用)
frac_gene <- function(tokens) {
  if (length(tokens) == 0) return(0)
  common <- c("ACTB","GAPDH","MALAT1","LEF1","FOXJ2","TAL1","SAMD11","NOC2L",
              "COL1A1","COL1A2","COL3A1","PIEZO2","HES4","ISG15","AGRN","FAM138A","OR4F5")
  pat <- "^(ENSG[0-9]+|NM_[0-9]+|[A-Z][A-Z0-9]{1,}(-[A-Z0-9]+)?|[A-Z]+[0-9]+)"
  tu <- toupper(tokens)
  hit <- sum(tu %in% common | grepl(pat, tu))
  hit / length(tokens)
}

## (E) 判斷矩陣方向:回傳 "gene_x_cell"(欄是細胞) 或 "cell_x_gene"(列是細胞,要轉置)
detect_orientation <- function(col_names, row_names) {
  col_bc <- frac_barcode(col_names); row_bc <- frac_barcode(row_names)
  ## 主規則:barcode 比例明顯高的那維 = 細胞
  if (max(col_bc, row_bc) >= 0.3 && abs(col_bc - row_bc) >= 0.2) {
    if (col_bc > row_bc) return("gene_x_cell") else return("cell_x_gene")
  }
  ## 後備:兩維都不像 barcode(細胞名是樣本碼如 H1_0001)時,看哪維像 gene
  col_g <- frac_gene(col_names); row_g <- frac_gene(row_names)
  if (row_g >= col_g) "gene_x_cell" else "cell_x_gene"
}

## (F) 偵測並回傳基因名前綴(>50% 命中才算),否則 NULL
detect_gene_prefix <- function(genes) {
  if (length(genes) == 0) return(NULL)
  gl <- tolower(genes)
  for (p in GENE_PREFIXES) {
    if (mean(startsWith(gl, tolower(p))) > 0.5) return(p)
  }
  NULL
}

## (G) 偵測合併檔的樣本前綴;回傳 named int vector(sample -> 細胞數)
##     穩健寫法:對每個 barcode,若形如 SAMPLE_<ACGT...>,用 sub 抽出 SAMPLE。
detect_sample_prefixes <- function(barcodes) {
  ## 只在「底線後緊接 >=8 個 ACGT」時,才把底線前的部分當樣本前綴
  has_prefix <- grepl("^[A-Za-z0-9]+_[ACGTacgt]{8,}", barcodes)
  if (!any(has_prefix)) return(integer(0))
  ## 抽出底線前的前綴(只取第一個底線之前)
  samp <- rep(NA_character_, length(barcodes))
  samp[has_prefix] <- sub("^([A-Za-z0-9]+)_[ACGTacgt]{8,}.*$", "\\1", barcodes[has_prefix])
  samp <- samp[!is.na(samp)]
  if (length(samp) == 0) return(integer(0))
  table(samp)
}

## (H) 核心:讀一個 dense matrix 檔 -> 回傳 list(mat=sparse, note=)
##     自動處理分隔符/方向/前綴。回傳的 mat 一定是 gene(row) × cell(col)。
read_dense_matrix <- function(path) {
  sep <- detect_sep(path)

  ## --- 讀檔(大檔用 data.table::fread 加速) ---
  if (HAVE_DATATABLE) {
    dt <- data.table::fread(path, sep = sep, header = TRUE,
                            data.table = FALSE, check.names = FALSE)
    raw_rownames <- as.character(dt[[1]])
    mat_num <- as.matrix(dt[, -1, drop = FALSE])
    col_names <- colnames(dt)[-1]
  } else {
    con <- if (grepl("\\.gz$", path)) gzfile(path) else path
    df <- utils::read.table(con, sep = sep, header = TRUE, row.names = 1,
                            check.names = FALSE, comment.char = "", quote = "\"")
    raw_rownames <- rownames(df)
    mat_num <- as.matrix(df)
    col_names <- colnames(df)
  }
  storage.mode(mat_num) <- "double"
  rownames(mat_num) <- raw_rownames
  colnames(mat_num) <- col_names

  ## --- 判斷方向,必要時轉置,讓 row=gene, col=cell ---
  orient <- detect_orientation(col_names, raw_rownames)
  note <- sprintf("sep=%s; orient=%s", ifelse(sep=="\t","TAB","COMMA"), orient)
  if (orient == "cell_x_gene") {
    mat_num <- t(mat_num)              # 轉置成 gene × cell
    note <- paste0(note, "(已轉置)")
  }

  ## --- 剝除基因名前綴 ---
  gene_names <- rownames(mat_num)
  prefix <- detect_gene_prefix(gene_names)
  if (!is.null(prefix)) {
    gene_names <- sub(paste0("^", prefix), "", gene_names, ignore.case = TRUE)
    rownames(mat_num) <- gene_names
    note <- paste0(note, sprintf("; 去前綴'%s'", prefix))
  }

  ## --- 基因名去重、轉 sparse ---
  rownames(mat_num) <- make.unique(rownames(mat_num))
  colnames(mat_num) <- make.unique(colnames(mat_num))
  sp <- as(Matrix::Matrix(mat_num, sparse = TRUE), "CsparseMatrix")
  list(mat = sp, note = note, prefix = prefix, orient = orient)
}

## 物種偵測(與 01 相同)
detect_mt_pattern <- function(genes) {
  n_human <- sum(grepl("^MT-",  genes)); n_mouse <- sum(grepl("^mt-",  genes))
  if (n_human >= n_mouse && n_human > 0) list(pattern="^MT-", species="human", n_mt=n_human, ribo="^RP[SL]")
  else if (n_mouse > 0)                  list(pattern="^mt-", species="mouse", n_mt=n_mouse, ribo="^Rp[sl]")
  else                                   list(pattern="^MT-", species="unknown(no MT genes)", n_mt=0, ribo="^RP[SL]")
}

## 摘要統計(與 01 相同)
summarise_obj <- function(obj, label) {
  md <- obj@meta.data
  safe_mean   <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
  safe_median <- function(x) if (length(x)) median(x, na.rm = TRUE) else NA_real_
  data.frame(
    state=label, n_cells=ncol(obj),
    n_genes_detected = sum(Matrix::rowSums(GetAssayData(obj, slot="counts") > 0) > 0),
    mean_nCount=safe_mean(md$nCount_RNA), median_nCount=safe_median(md$nCount_RNA),
    mean_nFeature=safe_mean(md$nFeature_RNA), median_nFeature=safe_median(md$nFeature_RNA),
    mean_percent_mt=safe_mean(md$percent.mt), median_percent_mt=safe_median(md$percent.mt),
    mean_percent_ribo=if(!is.null(md$percent.ribo)) safe_mean(md$percent.ribo) else NA_real_,
    stringsAsFactors=FALSE)
}
mad_bounds <- function(x, nmads=3, log=FALSE) {
  v <- if (log) log1p(x) else x
  med <- median(v, na.rm=TRUE); m <- mad(v, na.rm=TRUE)
  lo <- med - nmads*m; hi <- med + nmads*m
  if (log) { lo <- expm1(lo); hi <- expm1(hi) }
  c(lower=lo, upper=hi)
}

## ----------------------------- 3. 掃描來源 -------------------------------- ##
log_line("==== Dense-matrix QC run 開始 ====")
log_line(sprintf("來源: %s", SRC_ROOT))
log_line(sprintf("輸出: %s", OUT_ROOT))
ensure_dir(OUT_ROOT); ensure_dir(LOG_DIR)
if (!dir.exists(SRC_ROOT)) stop(sprintf("找不到來源資料夾: %s", SRC_ROOT))
if (!HAVE_DATATABLE) log_line("[提示] 未安裝 data.table,大型 csv/txt 讀取會較慢,建議 install.packages('data.table')")

PERSAMPLE_DIR <- file.path(LOG_DIR, "per_sample"); ensure_dir(PERSAMPLE_DIR)
write_metrics_row <- function(row) {
  fn <- file.path(PERSAMPLE_DIR, sprintf("%s__%s.csv", row$GSE,
                  ifelse(is.na(row$GSM), "NA", gsub("[^A-Za-z0-9_.-]", "_", row$GSM))))
  write_csv_utf8(row, fn)
}

## 結構掃描表
structure_rows <- list()
add_struct <- function(gse_raw, gse_clean, gsm, status, note, has_hash=FALSE, src_path="") {
  structure_rows[[length(structure_rows)+1]] <<- data.frame(
    GSE_folder_raw=gse_raw, GSE=gse_clean, GSM=gsm,
    structure_status=status, structure_note=note,
    gse_has_hash_prefix=has_hash, source_path=src_path, stringsAsFactors=FALSE)
}

## 找出一個 GSE 資料夾下所有「候選 matrix 檔」(排除正規化檔與壓縮包)
find_matrix_files <- function(gse_dir) {
  all <- list.files(gse_dir, full.names = TRUE, recursive = FALSE)
  all <- all[!dir.exists(all)]                              # 只要檔案
  ## 只收文字矩陣副檔名
  is_txt <- grepl("\\.(txt|tsv|csv)(\\.gz)?$", tolower(all))
  ## 排除壓縮包(.rar/.tar/.zip 等不是矩陣)
  is_arch <- grepl("\\.(rar|zip|tar|7z|gz)$", tolower(all)) &
             !grepl("\\.(txt|tsv|csv)\\.gz$", tolower(all))
  all[is_txt & !is_arch]
}

## ----------------------------- 4. 主迴圈 ---------------------------------- ##
top_items <- list.files(SRC_ROOT, full.names = TRUE)
gse_dirs  <- character(0)
for (it in top_items) {
  base <- basename(it)
  if (base %in% IGNORE_FOLDERS) {
    add_struct(base, sub("^#+","",base), NA, "IGNORED_BY_USER",
               sprintf("使用者指定忽略 (%s)", base), grepl("^#", base), it)
    log_line(sprintf("[ignore] 跳過資料夾: %s", base)); next
  }
  if (!dir.exists(it)) next                                  # 頂層非資料夾略過
  gse_dirs <- c(gse_dirs, it)
}
log_line(sprintf("偵測到 %d 個 GSE 資料夾", length(gse_dirs)))

## 一個「樣本」處理單元:給定 (mat, GSE, sample 名),跑與 01 相同的 QC
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

  ## 續跑檢查(與 01 相同)
  if (file.exists(pre_rds) && file.exists(post_rds) && file.exists(stats_csv) &&
      all(file.exists(expected_plots))) {
    log_line(sprintf("  [skip] %s/%s 已完整,跳過", gse_clean, sample_name))
    return("SKIPPED")
  }
  ensure_dir(out_dir); ensure_dir(plots_dir)

  tryCatch({
    obj <- CreateSeuratObject(counts = mat, project = sample_name,
                              min.cells = MIN_CELLS_GENE, min.features = MIN_FEATURES_CELL)
    obj$orig.ident1 <- gse_clean; obj$orig.ident2 <- sample_name
    if (!is.na(cond)) obj$condition <- cond

    det <- detect_mt_pattern(rownames(obj))
    obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = det$pattern)
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = det$ribo)
    log_line(sprintf("    [%s] 載入: %s | 物種=%s MT=%s(%d) | 細胞=%d 基因=%d",
                     sample_name, load_note, det$species, det$pattern, det$n_mt,
                     ncol(obj), nrow(obj)))

    pre_sum <- summarise_obj(obj, "pre_QC")
    mad_feat  <- mad_bounds(obj$nFeature_RNA, 3, log=TRUE)
    mad_count <- mad_bounds(obj$nCount_RNA,   3, log=TRUE)
    mad_mt    <- mad_bounds(obj$percent.mt,   3, log=FALSE)
    saveRDS(obj, pre_rds)

    if (HAVE_GGPLOT) {
      feats <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
                         colnames(obj@meta.data))
      pdf(expected_plots[1], width=4*length(feats), height=5)
      print(VlnPlot(obj, features=feats, ncol=length(feats), pt.size=0.1)); dev.off()
      s1 <- FeatureScatter(obj, "nCount_RNA","nFeature_RNA") + ggtitle("nCount vs nFeature")
      s2 <- FeatureScatter(obj, "nCount_RNA","percent.mt")   + ggtitle("nCount vs percent.mt")
      pdf(expected_plots[3], width=12, height=5)
      if (HAVE_PATCH) print(s1+s2) else { print(s1); print(s2) }; dev.off()
    }

    filter_expr <- sprintf("nFeature_RNA > %d & nFeature_RNA < %d & percent.mt < %d",
                           MIN_FEATURE, MAX_FEATURE, MAX_MT_PCT)
    n_before <- ncol(obj)
    obj_post <- subset(obj, subset = nFeature_RNA > MIN_FEATURE &
                                     nFeature_RNA < MAX_FEATURE & percent.mt < MAX_MT_PCT)
    n_after <- ncol(obj_post)

    n_doublets <- NA_integer_
    dbl_method <- if (!RUN_DOUBLETFINDER) "skipped(user)" else "none"
    ## (DoubletFinder 與 01 相同,如需可呼叫 run_doublet();此處從略保持簡潔)

    post_sum <- summarise_obj(obj_post, "post_QC")

    if (HAVE_GGPLOT) {
      feats <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
                         colnames(obj_post@meta.data))
      pdf(expected_plots[2], width=4*length(feats), height=5)
      print(VlnPlot(obj_post, features=feats, ncol=length(feats), pt.size=0.1)); dev.off()
      ba <- rbind(
        data.frame(stage="pre",  nFeature=obj$nFeature_RNA, nCount=obj$nCount_RNA, percent.mt=obj$percent.mt),
        data.frame(stage="post", nFeature=obj_post$nFeature_RNA, nCount=obj_post$nCount_RNA, percent.mt=obj_post$percent.mt))
      g1 <- ggplot(ba, aes(stage,nFeature,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+ggtitle("nFeature_RNA")
      g2 <- ggplot(ba, aes(stage,nCount,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+scale_y_log10()+ggtitle("nCount_RNA (log10)")
      g3 <- ggplot(ba, aes(stage,percent.mt,fill=stage))+geom_violin()+geom_boxplot(width=.1,outlier.size=.3)+theme_bw()+ggtitle("percent.mt")
      pdf(expected_plots[4], width=12, height=5)
      if (HAVE_PATCH) print(g1+g2+g3) else { print(g1);print(g2);print(g3) }; dev.off()
    }
    saveRDS(obj_post, post_rds)

    stats_df <- rbind(pre_sum, post_sum)
    stats_df$GSE <- gse_clean; stats_df$GSM <- sample_name; stats_df$condition <- cond
    stats_df$species <- det$species; stats_df$mt_pattern <- det$pattern; stats_df$n_mt_genes <- det$n_mt
    stats_df$filter_used <- filter_expr
    stats_df$MAD_nFeature_lower <- round(mad_feat["lower"],1); stats_df$MAD_nFeature_upper <- round(mad_feat["upper"],1)
    stats_df$MAD_nCount_lower <- round(mad_count["lower"],1);  stats_df$MAD_nCount_upper <- round(mad_count["upper"],1)
    stats_df$MAD_percent_mt_upper <- round(mad_mt["upper"],2)
    stats_df$cells_removed <- n_before-n_after
    stats_df$pct_cells_removed <- round(100*(n_before-n_after)/n_before,2)
    stats_df$n_doublets <- n_doublets; stats_df$doublet_method <- dbl_method
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
      cells_removed=n_before-n_after, pct_cells_removed=round(100*(n_before-n_after)/n_before,2),
      n_doublets=n_doublets, doublet_method=dbl_method, filter_used=filter_expr,
      source_path=src_path, output_path=out_dir, processed_time=ts(), stringsAsFactors=FALSE))

    log_line(sprintf("    完成 %s: 細胞 %d -> %d (移除 %.1f%%)",
                     sample_name, n_before, n_after, 100*(n_before-n_after)/n_before))
    "OK"
  }, error = function(e) {
    msg <- conditionMessage(e)
    log_line(sprintf("    [FAILED] %s/%s: %s", gse_clean, sample_name, msg))
    write_metrics_row(data.frame(
      GSE=gse_clean, GSE_folder_raw=gse_raw, GSM=sample_name, condition=cond,
      status="FAILED", error_note=msg, gse_has_hash_prefix=gse_hash,
      species=NA, mt_pattern=NA,
      pre_n_cells=NA,pre_n_genes=NA,pre_mean_nCount=NA,pre_mean_nFeature=NA,pre_mean_pct_mt=NA,pre_median_pct_mt=NA,
      post_n_cells=NA,post_n_genes=NA,post_mean_nCount=NA,post_mean_nFeature=NA,post_mean_pct_mt=NA,post_median_pct_mt=NA,
      cells_removed=NA,pct_cells_removed=NA,n_doublets=NA,doublet_method=NA,filter_used=NA,
      source_path=src_path, output_path=out_dir, processed_time=ts(), stringsAsFactors=FALSE))
    "FAILED"
  })
}

for (gse_dir in gse_dirs) {
  gse_raw   <- basename(gse_dir)
  gse_hash  <- grepl("^#", gse_raw)
  gse_clean <- sub("^#+", "", gse_raw)
  if (gse_hash) log_line(sprintf("[!] GSE '%s' 有 # 前綴,仍嘗試處理", gse_raw))

  mfiles <- find_matrix_files(gse_dir)
  if (length(mfiles) == 0) {
    add_struct(gse_raw, gse_clean, NA, "NO_MATRIX_FILE",
               "GSE 內找不到 txt/tsv/csv matrix 檔", gse_hash, gse_dir)
    log_line(sprintf("[!] GSE '%s' 無 matrix 檔,跳過", gse_raw)); next
  }

  for (mf in mfiles) {
    fbase <- basename(mf)
    ## 跳過正規化檔
    if (should_skip_file(mf)) {
      add_struct(gse_raw, gse_clean, fbase, "SKIPPED_NORMALIZED",
                 "檔名含 rpkm/fpkm/tpm 等,非 raw counts", gse_hash, mf)
      log_line(sprintf("  [skip] 正規化檔: %s", fbase)); next
    }

    log_line(sprintf("==> 讀取 %s / %s", gse_clean, fbase))
    loaded <- tryCatch(read_dense_matrix(mf), error = function(e) {
      log_line(sprintf("  [FAILED] 讀取失敗 %s: %s", fbase, conditionMessage(e)))
      add_struct(gse_raw, gse_clean, fbase, "READ_FAILED", conditionMessage(e), gse_hash, mf)
      NULL
    })
    if (is.null(loaded)) next
    mat <- loaded$mat
    add_struct(gse_raw, gse_clean, fbase, "OK", loaded$note, gse_hash, mf)

    ## 判斷是否為「合併檔」-> 依 barcode 前綴拆樣本
    samp_tbl <- detect_sample_prefixes(colnames(mat))
    samp_tbl <- samp_tbl[samp_tbl >= MIN_CELLS_PER_SPLIT]
    is_merged <- length(samp_tbl) >= 2

    if (is_merged) {
      log_line(sprintf("  偵測為合併檔,依前綴拆成 %d 個樣本: %s",
                       length(samp_tbl), paste(names(samp_tbl), collapse=", ")))
      for (sp_name in names(samp_tbl)) {
        sel <- grepl(paste0("^", sp_name, "_"), colnames(mat))
        submat <- mat[, sel, drop = FALSE]
        ## condition 由樣本前綴推測(KL/NS 等)
        cond <- if (grepl("^KL|keloid", sp_name, ignore.case=TRUE)) "Keloid"
                else if (grepl("^NS|normal", sp_name, ignore.case=TRUE)) "NormalSkin" else NA
        process_one_sample(submat, gse_clean, gse_raw, gse_hash,
                           sample_name = sp_name, cond = cond,
                           src_path = mf, load_note = paste0(loaded$note, "; 合併檔拆分"))
      }
    } else {
      ## 單一樣本檔:樣本名用檔名(去副檔名與常見後綴)
      sample_name <- sub("\\.(txt|tsv|csv)(\\.gz)?$", "", fbase, ignore.case=TRUE)
      sample_name <- sub("(_counts?_matrix|_counts|_matrix|_read|_dge|_UMI|\\.UMI\\.matrix)$", "",
                         sample_name, ignore.case=TRUE)
      ## 從檔名/樣本名推測 condition
      cond <- if (grepl("keloid|_KL|_K[0-9]", sample_name, ignore.case=TRUE)) "Keloid"
              else if (grepl("normal|_NS|nskin", sample_name, ignore.case=TRUE)) "NormalSkin" else NA
      process_one_sample(mat, gse_clean, gse_raw, gse_hash,
                         sample_name = sample_name, cond = cond,
                         src_path = mf, load_note = loaded$note)
    }
    suppressWarnings(rm(mat)); gc(verbose = FALSE)
  }
}

## ----------------------------- 5. 結構掃描表 ------------------------------ ##
if (length(structure_rows)) {
  struct_df <- do.call(rbind, structure_rows)
  write_csv_utf8(struct_df, file.path(LOG_DIR, "structure_scan_densematrix.csv"))
  log_line(sprintf("結構掃描表已寫出 (%d 列)", nrow(struct_df)))
}
log_line("==== Dense-matrix QC run 全部完成。可接 02_Aggregate_Summary.R 彙整 ====")
