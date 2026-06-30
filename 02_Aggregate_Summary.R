###############################################################################
## 02_Aggregate_Summary.R
## ---------------------------------------------------------------------------
## 在 01_QC_PerSample.R 全部跑完後執行。
## 把 _run_log/per_sample/*.csv (每個樣本一列) 合併成「大表格」，
## 並併入結構掃描資訊 (不符規定者的註記)。
##
## 輸出:
##   Z:\Dataset_Online\scRNA-seq\#Keloid_QC\QC_master_summary.csv
##   Z:\Dataset_Online\scRNA-seq\#Keloid_QC\QC_master_summary.xlsx (若有 openxlsx)
##
## 大表格欄位 (必要):
##   GSE (GSE 開頭, 大資料夾編號)、GSM (GSM 開頭, 子資料夾編號)、condition、status
##   QC 前後: Gene 數、細胞數、平均 Count 數、平均/中位 粒線體比例
##   cells_removed / pct_cells_removed、doublets、物種、過濾門檻、錯誤註記...
##   跑不出東西的樣本 (status=FAILED / 結構問題) 也會列出並註記。
###############################################################################

OUT_ROOT <- "Z:/Dataset_Online/scRNA-seq/#Keloid_QC"
LOG_DIR        <- file.path(OUT_ROOT, "_run_log")
PERSAMPLE_DIR  <- file.path(LOG_DIR, "per_sample")
STRUCT_CSV     <- file.path(LOG_DIR, "structure_scan.csv")

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
message(sprintf("[%s] 開始彙整大表", ts()))

## 跨平台穩健的 UTF-8 (帶 BOM) CSV 寫出，避免 locale 轉換破壞中文，
## 且 Windows Excel 直接開啟即正確顯示中文。
## 優先用 readr::write_excel_csv (原生 UTF-8 BOM)，否則用 base R 自行加 BOM。
write_csv_utf8 <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_excel_csv(df, path)            # 自帶 UTF-8 BOM
    return(invisible(path))
  }
  ## base R fallback: 用 capture.output 取得 CSV 文字，binary 寫出 + BOM
  lines <- utils::capture.output(utils::write.csv(df, file = "", row.names = FALSE))
  con <- file(path, open = "wb")                 # binary，不做 locale 轉換
  writeBin(charToRaw("\xEF\xBB\xBF"), con)        # UTF-8 BOM
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
  close(con)
  invisible(path)
}

if (!dir.exists(PERSAMPLE_DIR))
  stop(sprintf("找不到 per-sample 目錄: %s (請先跑 01_QC_PerSample.R)", PERSAMPLE_DIR))

## ---- 1. 讀取所有 per-sample row ----
files <- list.files(PERSAMPLE_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("per_sample 目錄內沒有任何 CSV，請確認 01 腳本有跑出結果")

read_safe <- function(f) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                          fileEncoding = "UTF-8-BOM"),
                 error = function(e)
                   tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                            error = function(e2) NULL))
  ## 保險: 若第一欄名殘留 BOM，清掉
  if (!is.null(df) && length(names(df)))
    names(df)[1] <- sub("^\ufeff", "", names(df)[1])
  df
}
rows <- lapply(files, read_safe)
rows <- rows[!vapply(rows, is.null, logical(1))]

## 對齊欄位 (不同 row 可能欄位數略異 -> 取聯集，缺的補 NA)
all_cols <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(df) {
  miss <- setdiff(all_cols, names(df))
  for (m in miss) df[[m]] <- NA
  df[all_cols]
})
master <- do.call(rbind, rows)

## ---- 2. 併入結構掃描的「不符規定」註記 ----
if (file.exists(STRUCT_CSV)) {
  struct <- tryCatch(
    read.csv(STRUCT_CSV, stringsAsFactors = FALSE, check.names = FALSE,
             fileEncoding = "UTF-8-BOM"),
    error = function(e)
      read.csv(STRUCT_CSV, stringsAsFactors = FALSE, check.names = FALSE))
  if (length(names(struct))) names(struct)[1] <- sub("^\ufeff", "", names(struct)[1])
  ## 標記非標準結構 (含 # 前綴、非資料夾、缺檔...) 給每個 GSM
  struct_key <- struct
  struct_key$GSM <- ifelse(is.na(struct_key$GSM), "", struct_key$GSM)
  master$GSM_join <- ifelse(is.na(master$GSM), "", master$GSM)

  master <- merge(
    master,
    struct_key[, c("GSE", "GSM", "structure_status", "structure_note")],
    by.x = c("GSE", "GSM_join"), by.y = c("GSE", "GSM"),
    all.x = TRUE, sort = FALSE
  )
  master$GSM_join <- NULL
} else {
  master$structure_status <- NA
  master$structure_note   <- NA
}

## ---- 3. 整理欄位順序 (必要欄位置前) ----
preferred <- c(
  "GSE", "GSM", "condition", "status",
  "gse_has_hash_prefix", "structure_status", "structure_note", "error_note",
  "species", "mt_pattern",
  ## QC 前
  "pre_n_cells", "pre_n_genes", "pre_mean_nCount", "pre_mean_nFeature",
  "pre_mean_pct_mt", "pre_median_pct_mt",
  ## QC 後
  "post_n_cells", "post_n_genes", "post_mean_nCount", "post_mean_nFeature",
  "post_mean_pct_mt", "post_median_pct_mt",
  ## 變化
  "cells_removed", "pct_cells_removed",
  "n_doublets", "doublet_method",
  "filter_used", "source_path", "output_path", "processed_time",
  "GSE_folder_raw"
)
ordered_cols <- c(intersect(preferred, names(master)),
                  setdiff(names(master), preferred))
master <- master[, ordered_cols]

## ---- 4. 排序: 先 GSE 再 GSM；FAILED / 結構問題排在各 GSE 內後面以便檢視 ----
master <- master[order(master$GSE,
                       master$status != "OK",   # OK 在前
                       master$GSM), ]

## ---- 5. 寫出 ----
out_csv <- file.path(OUT_ROOT, "QC_master_summary.csv")
write_csv_utf8(master, out_csv)
message(sprintf("[%s] 大表已寫出 (%d 列, %d 欄): %s",
                ts(), nrow(master), ncol(master), out_csv))

## (選用) 也寫一份 xlsx，方便直接開
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(master, file.path(OUT_ROOT, "QC_master_summary.xlsx"),
                       overwrite = TRUE)
  message(sprintf("[%s] 同時寫出 xlsx 版本", ts()))
}

## ---- 6. 簡短統計摘要印到 console ----
n_total  <- nrow(master)
n_ok     <- sum(master$status == "OK", na.rm = TRUE)
n_failed <- sum(master$status == "FAILED", na.rm = TRUE)
n_other  <- n_total - n_ok - n_failed
message("---------------- 彙整摘要 ----------------")
message(sprintf("總列數: %d", n_total))
message(sprintf("  成功 (OK):        %d", n_ok))
message(sprintf("  失敗 (FAILED):    %d", n_failed))
message(sprintf("  其他/結構問題:    %d", n_other))
if (n_ok > 0) {
  ok_rows <- master[master$status == "OK", ]
  message(sprintf("成功樣本 QC 後總細胞數: %s",
                  format(sum(as.numeric(ok_rows$post_n_cells), na.rm = TRUE),
                         big.mark = ",")))
}
message("------------------------------------------")
