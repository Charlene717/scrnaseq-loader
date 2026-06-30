#!/usr/bin/env bash
#
# download_sra.sh — 下載 SRA run 並轉成壓縮的 FASTQ
#
# 用法:
#   ./download_sra.sh SRR16475068
#   ./download_sra.sh SRR16475068 SRR16475069 SRR16475070
#   ./download_sra.sh -o /path/to/outdir -t 8 SRR16475068
#
set -euo pipefail

# ---- 預設參數 ----
OUTDIR="."
THREADS=4

usage() {
    cat <<EOF
用法: $0 [選項] <SRR_ID> [SRR_ID ...]

選項:
  -o DIR    輸出目錄 (預設: 目前目錄)
  -t N      fasterq-dump 使用的執行緒數 (預設: 4)
  -h        顯示此說明

範例:
  $0 SRR16475068
  $0 -o ./fastq -t 8 SRR16475068 SRR16475069
EOF
}

while getopts "o:t:h" opt; do
    case "$opt" in
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ "$#" -eq 0 ]; then
    echo "錯誤: 請至少提供一個 SRR ID" >&2
    usage
    exit 1
fi

# ---- 檢查工具是否安裝 ----
for tool in prefetch fasterq-dump; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "錯誤: 找不到 '$tool',請先安裝 SRA Toolkit" >&2
        echo "  conda install -c bioconda sra-tools" >&2
        echo "  或 sudo apt install sra-toolkit" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

# ---- 主迴圈 ----
for ACC in "$@"; do
    echo "==> [$ACC] 開始下載..."

    # 1. prefetch 抓 .sra
    prefetch "$ACC" --output-directory "$OUTDIR"

    # 2. fasterq-dump 轉 FASTQ (--split-files 處理雙端)
    echo "==> [$ACC] 轉換為 FASTQ..."
    fasterq-dump "$OUTDIR/$ACC/$ACC.sra" \
        --threads "$THREADS" \
        --split-files \
        --outdir "$OUTDIR"

    # 3. gzip 壓縮
    echo "==> [$ACC] 壓縮..."
    gzip -f "$OUTDIR/$ACC"*.fastq

    # 4. (選用) 清掉中間的 .sra 以節省空間,如不需要可註解掉下一行
    rm -rf "${OUTDIR:?}/$ACC"

    echo "==> [$ACC] 完成: $OUTDIR/$ACC*.fastq.gz"
    echo
done

echo "全部完成。"
