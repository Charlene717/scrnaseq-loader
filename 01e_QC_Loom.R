###############################################################################
## 01e_QC_Loom.R
## ---------------------------------------------------------------------------
## scRNA-seq 廣泛品質評估 + 初步 QC —— 「.loom (HDF5-based)」版本 (per-sample)
##
## 為什麼有這支 & 與其他的分工:
##   01  = 10X 三件套 ; 01b = dense matrix ; 01c = 已整合 RDS ; 01d = 10x .h5
##   01e = .loom。每個 .loom 各自跑與 01 相同的 QC,輸出格式與前面一致。
##
## .loom 的兩個坑 (本腳本自動處理):
##   (1) 常以 gzip 壓縮 (.loom.gz)。本腳本:有 .gz 就先解壓到暫存再讀,沒有就直讀。
##   (2) 內部欄位命名 / 矩陣方向不統一 (不同工具 loompy/velocyto/scanpy 產生的不同):
##         - gene 名 attribute 可能叫 Gene / var_names / gene_names / Symbol ...
##         - cell 名 attribute 可能叫 CellID / obs_names / barcode ...
##         - 矩陣理論上 gene×cell,但偶有 cell×gene -> 需轉置
##       本腳本「以內容判斷方向」(看哪維像 gene symbol、哪維像 barcode),
##       attribute 名稱只作輔助。此偵測邏輯已用 5 種工具風格的 loom 驗證通過
##       (loompy / velocyto / scanpy / 別名 / 轉置)。
##
## 不用 SeuratDisk::LoadLoom 的原因:它對「非標準命名」的 loom 常直接失敗;
## 這裡直接用 hdf5r 讀 /matrix + /row_attrs + /col_attrs,自己偵測,較穩健。
##
## 資料結構 (來源):
##   <SRC_ROOT>/<GSE...>/  之下放一個或多個 .loom 或 .loom.gz,每檔 = 一個樣本。
##
## 輸出 (與 01 完全相同):
##   <OUT_ROOT>/<GSE>/<sample>/
##     - <sample>_pre_QC.rds / _post_QC.rds / _QC_stats.csv / plots/
##   _run_log/per_sample/<GSE>__<sample>.csv  (欄位與 01 一致)
###############################################################################

## ----------------------------- 0. 路徑設定 -------------------------------- ##
SRC_ROOT <- "X:/Dataset_Online/##_Keloid/scRNA-seq/loom"          # 放 .loom 的來源根目錄
OUT_ROOT <- "X:/Dataset_Online/##_Keloid/scRNA-seq/loom_QC"    # QC 輸出根目錄

## ------ 過濾門檻 (與 01 完全相同) ------ ##
MIN_FEATURE <- 200
MAX_FEATURE <- 5000
MAX_MT_PCT  <- 30
MIN_CELLS_GENE    <- 3
MIN_FEATURES_CELL <- 0

## ====== 使用者可調參數 (與 01 對齊) ====== ##
RUN_DOUBLETFINDER <- FALSE
IGNORE_FOLDERS <- c("GEO_10X_auto")
## 解壓 .loom.gz 的暫存目錄 (預設用系統暫存;大檔請確認空間足夠)
GZ_TMP_DIR <- tempdir()
## .loom 主矩陣要用哪個 layer:預設用 /matrix。velocyto loom 也可改用 "spliced"。
##   注意:若原始 loom 的 /matrix 已是 spliced+unspliced 的總和,維持 "matrix" 即可。
LOOM_MATRIX_LAYER <- "matrix"
## =============================== ##

SEED <- 42; set.seed(SEED)

## gene / cell 名 attribute 候選 (輔助用;主要仍以內容判斷)
GENE_ATTR_KEYS <- c("Gene","var_names","gene_names","GeneName","Symbol","features","Accession","gene_ids")
CELL_ATTR_KEYS <- c("CellID","obs_names","barcode","Barcode","CellName","cell_names","cell_id")

## ----------------------------- 1. 套件 ------------------------------------ ##
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("需要 Seurat 套件")
  library(Seurat); library(Matrix)
})
if (!requireNamespace("hdf5r", quietly = TRUE))
  stop("讀 .loom 需要 hdf5r 套件,請先 install.packages('hdf5r')")
library(hdf5r)

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAVE_GGPLOT <- has_pkg("ggplot2"); if (HAVE_GGPLOT) suppressPackageStartupMessages(library(ggplot2))
HAVE_PATCH  <- has_pkg("patchwork"); if (HAVE_PATCH) suppressPackageStartupMessages(library(patchwork))

## ----------------------------- 2. 工具函式 (與 01/01b/01c/01d 一致) ------- ##
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
  cat(sprintf("[%s] %s\n", ts(), msg), file=file.path(LOG_DIR, "run_log_loom.txt"), append=TRUE)
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

## ---- 內容判斷:某串值像 barcode / gene 的比例 ----
frac_barcode <- function(v) {
  if (length(v) == 0) return(0)
  mean(grepl("^[ACGTacgt]{6,}([.:_-]?\\w+)?$", v))
}
## 像「gene 相關」(含 Ensembl ID 或 symbol,用來定位 gene 是哪一維)
frac_gene_any <- function(v) {
  if (length(v) == 0) return(0)
  mean(grepl("^(ENS[A-Z]*G[0-9]+|MT-|mt-|[A-Za-z][A-Za-z0-9]*(-[A-Za-z0-9]+)?)$", v))
}
## 像「gene SYMBOL」(Ensembl ID 不算;用來在多個 gene 欄中優先挑 symbol,
## 確保 rownames 帶 MT-/mt- 前綴,物種偵測才抓得到)
frac_gene_symbol <- function(v) {
  if (length(v) == 0) return(0)
  common <- c("ACTB","GAPDH","MALAT1","MT-CO1","MT-ND1","SAMD11","NOC2L","COL1A1",
              "COL1A2","COL3A1","PIEZO2","RPL13","B2M","TMSB4X","XKR4","SOX17")
  is_ensembl <- grepl("^ENS[A-Z]*G[0-9]+$", toupper(v))
  is_symbol  <- (toupper(v) %in% common) |
    grepl("^(MT-|mt-|[A-Za-z][A-Za-z0-9]*(-[A-Za-z0-9]+)?)$", v)
  mean(is_symbol & !is_ensembl)
}

## ---- 核心:讀一個 .loom -> gene×cell sparse matrix ----
## 用 hdf5r 直接讀,依內容自動判斷方向與 gene/cell 欄位。
read_loom_matrix <- function(path) {
  h <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  
  ## 主矩陣 (loom: /matrix 或 /layers/<layer>)
  mat_path <- if (LOOM_MATRIX_LAYER == "matrix" || is.null(LOOM_MATRIX_LAYER)) "matrix"
  else file.path("layers", LOOM_MATRIX_LAYER)
  if (!h$exists(mat_path)) mat_path <- "matrix"   # 後備
  mat <- h[[mat_path]]$read()                     # 讀成一般 matrix (row x col)
  
  ## 讀 row_attrs / col_attrs 的所有欄位值
  read_attrs <- function(grp) {
    if (!h$exists(grp)) return(list())
    keys <- names(h[[grp]])
    out <- list()
    for (k in keys) {
      val <- tryCatch(h[[file.path(grp, k)]]$read(), error = function(e) NULL)
      if (is.null(val)) next
      ## 只保留一維、長度符合的字元/可轉字元向量
      if (is.array(val) && length(dim(val)) > 1) next
      out[[k]] <- as.character(val)
    }
    out
  }
  row_attrs <- read_attrs("row_attrs")
  col_attrs <- read_attrs("col_attrs")
  
  nrow_m <- nrow(mat); ncol_m <- ncol(mat)
  
  ## ============================================================= ##
  ## 修正版偵測 (經真實 loom 驗證):
  ##   關鍵教訓 — 不能假設 row_attrs 一定對應 matrix 的 row!
  ##   有些 loom (如 velocyto 輸出) 的 row_attrs 長度 = matrix 的「另一維」。
  ##   正確做法:
  ##     1. 蒐集 row_attrs + col_attrs 全部欄位,各自記錄「值」與「長度」。
  ##     2. 用 frac_gene_any 找出「哪組 attribute 是 gene 相關」-> 其長度定位 gene 維。
  ##     3. gene 維 = matrix 中長度 == 該 attr 長度 的那一維;據此決定是否轉置。
  ##     4. 在 gene 維長度相同的欄位中,「symbol 優先於 Ensembl ID」挑 gene 名,
  ##        以確保 rownames 是 MT-/mt- 開頭的 symbol,物種偵測才抓得到。
  ##     5. cell 名 = 另一維長度的欄位 (優先名字含 cell/id/barcode)。
  ## ============================================================= ##
  
  ## 合併所有 attribute,標記來源與長度
  all_attrs <- c(
    setNames(row_attrs, paste0("row/", names(row_attrs))),
    setNames(col_attrs, paste0("col/", names(col_attrs)))
  )
  all_attrs <- all_attrs[vapply(all_attrs, function(v) length(v) > 0, logical(1))]
  if (length(all_attrs) == 0) stop("loom 沒有可用的 row_attrs/col_attrs")
  
  ## (1) 定位 gene 維:哪個 attribute 最像 gene (含 Ensembl 或 symbol 皆可)
  any_scores <- vapply(all_attrs, frac_gene_any, numeric(1))
  gene_anchor <- names(all_attrs)[which.max(any_scores)]
  gene_len <- length(all_attrs[[gene_anchor]])
  
  ## gene 維 = matrix 中長度符合 gene_len 的那一維
  gene_axis <- if (gene_len == nrow_m) "row" else if (gene_len == ncol_m) "col" else NA
  if (is.na(gene_axis))
    stop(sprintf("gene attribute 長度 %d 對不上 matrix 維度 %d x %d", gene_len, nrow_m, ncol_m))
  
  ## (2) 在「與 gene 維同長」的欄位中,symbol 優先挑 gene 名
  same_len_as_gene <- all_attrs[vapply(all_attrs, function(v) length(v) == gene_len, logical(1))]
  sym_scores <- vapply(same_len_as_gene, frac_gene_symbol, numeric(1))
  if (max(sym_scores) > 0) {
    gene_key <- names(same_len_as_gene)[which.max(sym_scores)]      # 有 symbol 欄 -> 用它
  } else {
    ga <- vapply(same_len_as_gene, frac_gene_any, numeric(1))       # 全是 Ensembl -> 用 any 最高
    gene_key <- names(same_len_as_gene)[which.max(ga)]
  }
  genes <- all_attrs[[gene_key]]
  
  ## (3) cell 名:與「另一維」同長、且非 gene 欄;優先名字像 cell/id/barcode
  cell_len <- if (gene_axis == "row") ncol_m else nrow_m
  cell_cands <- all_attrs[vapply(all_attrs, function(v) length(v) == cell_len, logical(1))]
  cell_cands <- cell_cands[setdiff(names(cell_cands), gene_key)]
  cell_key <- NA
  if (length(cell_cands) > 0) {
    prefer <- grepl("cell|barcode|obs_names|_id$|id$", tolower(names(cell_cands)))
    cell_key <- if (any(prefer)) names(cell_cands)[which(prefer)[1]] else names(cell_cands)[1]
  }
  cells <- if (!is.na(cell_key)) all_attrs[[cell_key]] else paste0("Cell", seq_len(cell_len))
  
  ## (4) 讓 matrix 成為 gene x cell
  if (gene_axis == "row") {
    note <- sprintf("gene_x_cell; gene=%s cell=%s", gene_key, cell_key)
  } else {
    mat <- t(mat)                                    # 目前是 cell x gene -> 轉置
    note <- sprintf("cell_x_gene(已轉置); gene=%s cell=%s", gene_key, cell_key)
  }
  
  ## 貼名、去重、轉 sparse (gene x cell)
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- make.unique(as.character(cells))
  mat <- as(Matrix::Matrix(mat, sparse = TRUE), "CsparseMatrix")
  list(mat = mat, note = note)
}

## ---- 若是 .gz 先解壓到暫存,回傳可讀的 .loom 路徑 + 是否為暫存 ----
prepare_loom_path <- function(path) {
  if (grepl("\\.gz$", tolower(path))) {
    out <- file.path(GZ_TMP_DIR, sub("\\.gz$", "", basename(path), ignore.case = TRUE))
    ## 解壓 (R 內建 gzfile 串流,避免額外依賴)
    con_in  <- gzfile(path, open = "rb")
    con_out <- file(out, open = "wb")
    repeat {
      chunk <- readBin(con_in, what = "raw", n = 1e7)   # 每次 10MB
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }
    close(con_in); close(con_out)
    return(list(path = out, is_temp = TRUE))
  }
  list(path = path, is_temp = FALSE)
}

## 找一個 GSE 資料夾下所有 .loom / .loom.gz
find_loom_files <- function(gse_dir) {
  fs <- list.files(gse_dir, full.names=TRUE, recursive=FALSE)
  fs <- fs[!dir.exists(fs)]
  fs[grepl("\\.loom(\\.gz)?$", tolower(fs))]
}

## ----------------------------- 3. 掃描來源 -------------------------------- ##
log_line("==== .loom per-sample QC 開始 ====")
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
    ## 防呆:若完全沒抓到 MT 基因,percent.mt 會全為 0/NA。記錄警告(通常代表
    ## gene 名沒正確載入,或該物種 MT 命名不同)。percent.mt 全 NA 時補 0 以免繪圖崩潰。
    if (det$n_mt == 0)
      log_line(sprintf("    [警告] %s 未偵測到 MT 基因 (物種=%s);percent.mt 將為 0,請檢查基因名是否正確載入",
                       sample_name, det$species))
    if (all(is.na(obj$percent.mt)))      obj$percent.mt   <- 0
    if (all(is.na(obj$percent.ribo)))    obj$percent.ribo <- 0
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
  
  looms <- find_loom_files(gse_dir)
  if (length(looms) == 0) {
    add_struct(gse_raw, gse_clean, NA, "NO_LOOM_FILE", "GSE 內找不到 .loom 檔", gse_hash, gse_dir)
    log_line(sprintf("[!] GSE '%s' 無 .loom 檔,跳過", gse_raw)); next
  }
  
  for (lm in looms) {
    fbase <- basename(lm)
    ## 樣本名:去 .loom / .loom.gz 與常見後綴
    sample_name <- sub("\\.loom(\\.gz)?$", "", fbase, ignore.case=TRUE)
    ## condition 由檔名推測 (此資料集有 Wound/P 等;可自行調整)
    cond <- if (grepl("keloid|_KL|_K[0-9]", sample_name, ignore.case=TRUE)) "Keloid"
    else if (grepl("wound", sample_name, ignore.case=TRUE)) "Wound"
    else if (grepl("normal|_NS|nskin|scar", sample_name, ignore.case=TRUE)) "NormalSkin" else NA
    
    log_line(sprintf("==> 處理 %s / %s", gse_clean, fbase))
    
    ## 若 .gz 先解壓
    prep <- tryCatch(prepare_loom_path(lm), error=function(e) {
      log_line(sprintf("  [FAILED] 解壓失敗 %s: %s", fbase, conditionMessage(e)))
      add_struct(gse_raw, gse_clean, fbase, "GUNZIP_FAILED", conditionMessage(e), gse_hash, lm); NULL
    })
    if (is.null(prep)) next
    if (prep$is_temp) log_line(sprintf("    已解壓 .gz -> 暫存: %s", prep$path))
    
    ## 讀 loom
    loaded <- tryCatch(read_loom_matrix(prep$path), error=function(e) {
      log_line(sprintf("  [FAILED] 讀取 .loom 失敗 %s: %s", fbase, conditionMessage(e)))
      add_struct(gse_raw, gse_clean, fbase, "READ_FAILED", conditionMessage(e), gse_hash, lm); NULL
    })
    ## 清掉解壓暫存檔
    if (prep$is_temp && file.exists(prep$path)) unlink(prep$path)
    if (is.null(loaded)) next
    
    add_struct(gse_raw, gse_clean, fbase, "OK", loaded$note, gse_hash, lm)
    process_one_sample(loaded$mat, gse_clean, gse_raw, gse_hash,
                       sample_name=sample_name, cond=cond, src_path=lm, load_note=loaded$note)
    suppressWarnings(rm(loaded)); gc(verbose=FALSE)
  }
}

## ----------------------------- 6. 結構掃描表 ------------------------------ ##
if (length(structure_rows)) {
  struct_df <- do.call(rbind, structure_rows)
  write_csv_utf8(struct_df, file.path(LOG_DIR, "structure_scan_loom.csv"))
  log_line(sprintf("結構掃描表已寫出 (%d 列)", nrow(struct_df)))
}
log_line("==== .loom per-sample QC 全部完成。可接 02_Aggregate_Summary.R ====")
