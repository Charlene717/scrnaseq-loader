###############################################################################
## diagnose_loom.R  —  攤開一個 .loom 的真實內部結構 (除錯用)
## ---------------------------------------------------------------------------
## 用途:01e 對真實 loom 偵測失敗 (gene=NA cell=NA)。這支只做「檢查」,
##       印出 loom 的實際結構,讓我們看清楚 gene/cell 名到底存在哪、什麼型別。
##       不做 QC、不寫任何檔、不依賴 01e。
##
## 用法:改下面 LOOM_GZ 成你任一個 .loom.gz 的完整路徑,然後 source 這支。
##       (P21_2 那種「只有 1866 基因」的、和 P211Wound 那種「13127 基因」的
##        最好各跑一次貼給我,因為兩者差異很大。)
###############################################################################

## >>>>>>>>>> 改這裡:指向一個 .loom.gz <<<<<<<<<<
LOOM_GZ <- "X:/Dataset_Online/##_Keloid/scRNA-seq/loom/#GSE153596/GSM4647789_P21_2.loom.gz"
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if (!requireNamespace("hdf5r", quietly = TRUE)) stop("需要 hdf5r")
library(hdf5r)

cat("========================================================\n")
cat("診斷檔案:", LOOM_GZ, "\n")
cat("========================================================\n\n")

## --- 先解壓 .gz 到暫存 ---
if (grepl("\\.gz$", tolower(LOOM_GZ))) {
  tmp <- file.path(tempdir(), sub("\\.gz$", "", basename(LOOM_GZ), ignore.case = TRUE))
  ci <- gzfile(LOOM_GZ, "rb"); co <- file(tmp, "wb")
  repeat { ch <- readBin(ci, "raw", n = 1e7); if (!length(ch)) break; writeBin(ch, co) }
  close(ci); close(co)
  loom_path <- tmp
  cat("已解壓到暫存:", tmp, "\n\n")
} else loom_path <- LOOM_GZ

h <- H5File$new(loom_path, mode = "r")

## --- 1. 頂層有哪些群組/資料集 ---
cat("── 1. 頂層結構 ──────────────────────\n")
cat("頂層名稱:", paste(names(h), collapse = ", "), "\n\n")

## --- 2. /matrix 的維度與型別 ---
cat("── 2. /matrix ───────────────────────\n")
if (h$exists("matrix")) {
  d <- h[["matrix"]]
  cat("維度 (dims):", paste(d$dims, collapse = " x "), "\n")
  cat("(loom 規範: 第一維=genes, 第二維=cells;若相反則需轉置)\n")
  cat("dtype:", d$get_type()$to_text(), "\n\n")
} else cat("!! 沒有 /matrix\n\n")

## --- 3. row_attrs:每個欄位的型別 + 前 5 個值 ---
peek_group <- function(grp) {
  if (!h$exists(grp)) { cat("!! 沒有", grp, "\n\n"); return(invisible()) }
  keys <- names(h[[grp]])
  cat("欄位 (", length(keys), "個):", paste(keys, collapse = ", "), "\n\n")
  for (k in keys) {
    obj <- h[[file.path(grp, k)]]
    cat("  ▸", k, "\n")
    cat("      dims :", paste(obj$dims, collapse = " x "),
        "| dtype:", tryCatch(obj$get_type()$to_text(), error = function(e) "?"), "\n")
    val <- tryCatch(obj$read(), error = function(e) paste("讀取失敗:", conditionMessage(e)))
    ## 若是多維 (matrix),標示出來
    if (is.array(val) && length(dim(val)) > 1) {
      cat("      (多維陣列, dim =", paste(dim(val), collapse = "x"), ")\n")
      cat("      前幾個值:", paste(utils::head(as.vector(val), 5), collapse = " | "), "\n")
    } else {
      cat("      class:", class(val)[1], "\n")
      cat("      前5值:", paste(utils::head(as.character(val), 5), collapse = " | "), "\n")
    }
    cat("\n")
  }
}

cat("── 3. /row_attrs (通常放 gene 資訊) ─────────────\n")
peek_group("row_attrs")

cat("── 4. /col_attrs (通常放 cell 資訊) ─────────────\n")
peek_group("col_attrs")

## --- 5. layers (velocyto 會有 spliced/unspliced) ---
cat("── 5. /layers ───────────────────────\n")
if (h$exists("layers")) {
  cat("layers:", paste(names(h[["layers"]]), collapse = ", "), "\n\n")
} else cat("(無 layers)\n\n")

## --- 6. 關鍵自動判斷:哪個欄位像 gene、哪個像 barcode ---
cat("── 6. 內容判斷測試 (看分數為什麼會是 0) ──────────\n")
frac_barcode <- function(v) if (!length(v)) 0 else mean(grepl("^[ACGTacgt]{6,}([.:_-]?\\w+)?$", v))
frac_gene <- function(v) {
  if (!length(v)) return(0)
  common <- c("ACTB","GAPDH","MALAT1","MT-CO1","SAMD11","COL1A1","PIEZO2","B2M")
  mean(toupper(v) %in% common |
       grepl("^(ENSG[0-9]+|ENSMUSG[0-9]+|[A-Z][A-Z0-9]{1,}(-[A-Z0-9]+)?)$", toupper(v)))
}
test_group <- function(grp) {
  if (!h$exists(grp)) return()
  for (k in names(h[[grp]])) {
    val <- tryCatch(as.character(h[[file.path(grp,k)]]$read()), error = function(e) character(0))
    if (!length(val)) next
    cat(sprintf("  %s/%-15s  frac_gene=%.2f  frac_barcode=%.2f  (前2: %s)\n",
                grp, k, frac_gene(val), frac_barcode(val),
                paste(utils::head(val,2), collapse=",")))
  }
}
test_group("row_attrs")
test_group("col_attrs")

h$close_all()
if (exists("tmp") && file.exists(tmp)) unlink(tmp)
cat("\n========== 診斷完成,請把以上整段輸出貼給我 ==========\n")
