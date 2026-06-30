#!/usr/bin/env bash
###############################################################################
## inspect_sra.sh
## ---------------------------------------------------------------------------
## 在「真正開跑 Cell Ranger 之前」,先用一小撮 reads 偵測一筆 SRA run 的:
##   1. 是不是 10x Genomics  (看 R1 長度結構)
##   2. 化學版本 v2 / v3      (R1 = 26bp -> v2 ; 28bp -> v3/v3.1)
##   3. 物種 human / mouse    (抽樣 reads 比對 human/mouse 看誰 mapping 高)
##   4. read 配置 (哪個是 barcode+UMI、哪個是 cDNA、哪個是 index)
##
## 設計成「唯讀偵測」:只轉前 N 條 reads,不動原始 .sra,不產生大檔。
## 物種偵測為「選用」:需要 minimap2 + 兩個物種的 transcriptome(可只放一個)。
## 沒有 minimap2 / 參考檔時,仍會輸出 1~2 與 4(這幾項不需比對)。
##
## 用法:
##   ./inspect_sra.sh <SRA_PATH_或_SRR_ID> [輸出目錄]
##
## 範例:
##   ./inspect_sra.sh /path/PRJNA772373/SRR16475068
##   ./inspect_sra.sh SRR16475068 ./inspect_out
##
## 物種偵測(選用)— 用環境變數指定參考(transcriptome FASTA,可 .gz):
##   HUMAN_REF=/ref/human_cdna.fa.gz MOUSE_REF=/ref/mouse_cdna.fa.gz \
##     ./inspect_sra.sh SRR16475068
###############################################################################
set -euo pipefail

# ---- 參數 ----
SRA_INPUT="${1:-}"
OUT_DIR="${2:-./inspect_$(basename "${SRA_INPUT%.*}")}"
N_READS="${N_READS:-100000}"     # 抽樣 reads 數(偵測用,夠判斷即可)
THREADS="${THREADS:-4}"

# 物種參考(選用):transcriptome / cDNA FASTA
HUMAN_REF="${HUMAN_REF:-}"
MOUSE_REF="${MOUSE_REF:-}"

if [ -z "$SRA_INPUT" ]; then
    echo "用法: $0 <SRA_PATH 或 SRR_ID> [輸出目錄]" >&2
    exit 1
fi

# ---- 工具檢查 ----
need() { command -v "$1" >/dev/null 2>&1; }
if ! need fasterq-dump; then
    echo "錯誤: 找不到 fasterq-dump,請先安裝 SRA Toolkit" >&2
    echo "  conda install -c bioconda sra-tools" >&2
    exit 1
fi
HAVE_MINIMAP=false; need minimap2 && HAVE_MINIMAP=true
HAVE_SEQKIT=false;  need seqkit   && HAVE_SEQKIT=true

mkdir -p "$OUT_DIR"
PEEK_DIR="$OUT_DIR/peek_fastq"
mkdir -p "$PEEK_DIR"
REPORT="$OUT_DIR/inspect_report.txt"
: > "$REPORT"

log() { echo "$@" | tee -a "$REPORT"; }

SRR_ID="$(basename "${SRA_INPUT%.*}")"
# 若傳進來是不帶副檔名的 .sra 檔(像你截圖那樣),basename 還是會對
[ -e "$SRA_INPUT" ] && SRR_ID="$(basename "$SRA_INPUT")"
SRR_ID="${SRR_ID%.sra}"

log "================================================================"
log "  SRA Inspect Report : $SRR_ID"
log "  時間: $(date '+%Y-%m-%d %H:%M:%S')"
log "================================================================"

# =========================================================================
# 1) 抽樣轉出前 N 條 reads(含 technical reads,才看得到 barcode read)
# =========================================================================
log ""
log "[1/4] 抽樣前 $N_READS 條 reads ..."
# -X N : 只取前 N 個 spot ; --include-technical : 保留 barcode/index read
fasterq-dump "$SRA_INPUT" \
    --split-files \
    --include-technical \
    -X "$N_READS" \
    --threads "$THREADS" \
    --outdir "$PEEK_DIR" \
    --force >/dev/null 2>&1 || {
        log "  ! fasterq-dump 抽樣失敗,請確認路徑/檔案正確"
        exit 1
    }

mapfile -t FQS < <(ls "$PEEK_DIR"/*.fastq 2>/dev/null | sort)
if [ "${#FQS[@]}" -eq 0 ]; then
    log "  ! 沒有轉出任何 fastq,無法偵測"
    exit 1
fi
log "  轉出 ${#FQS[@]} 個 fastq 檔:"
for f in "${FQS[@]}"; do log "    - $(basename "$f")"; done

# =========================================================================
# 2) 量每個 fastq 的 read 長度(取眾數),判斷 10x 與化學版本
# =========================================================================
log ""
log "[2/4] 量測 read 長度 ..."

# 回傳某 fastq 最常見的 read 長度
modal_len() {
    awk 'NR%4==2 { print length($0) }' "$1" \
        | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}'
}

declare -A LEN
for f in "${FQS[@]}"; do
    base="$(basename "$f")"
    L="$(modal_len "$f")"
    LEN["$base"]="$L"
    log "    $base : 眾數長度 = ${L} bp"
done

# 判斷各檔角色。
# 關鍵:index read(~8bp)比 barcode read(26/28bp)更短,不能用「最短=barcode」。
#   cDNA(R2)    = 最長
#   index(I1)   = 8~10bp
#   barcode(R1) = 20~32bp(排除 cDNA/index 後最短)
longest_len=0; longest_file=""
for f in "${FQS[@]}"; do
    base="$(basename "$f")"; L_="${LEN[$base]}"
    [ "$L_" -gt "$longest_len" ] && { longest_len="$L_"; longest_file="$base"; }
done
# index
index_file=""
for f in "${FQS[@]}"; do
    base="$(basename "$f")"; L_="${LEN[$base]}"
    [ "$base" = "$longest_file" ] && continue
    if [ "$L_" -ge 8 ] && [ "$L_" -le 10 ]; then index_file="$base"; fi
done
# barcode = 20~32bp,排除 cDNA/index 後最短
shortest_len=999999; shortest_file=""
for f in "${FQS[@]}"; do
    base="$(basename "$f")"; L_="${LEN[$base]}"
    [ "$base" = "$longest_file" ] && continue
    [ "$base" = "$index_file" ] && continue
    if [ "$L_" -ge 20 ] && [ "$L_" -le 32 ] && [ "$L_" -lt "$shortest_len" ]; then
        shortest_len="$L_"; shortest_file="$base"
    fi
done
# 後備:barcode 不在典型區間時,取排除 cDNA/index 後最短
if [ -z "$shortest_file" ]; then
    shortest_len=999999
    for f in "${FQS[@]}"; do
        base="$(basename "$f")"; L_="${LEN[$base]}"
        [ "$base" = "$longest_file" ] && continue
        [ "$base" = "$index_file" ] && continue
        if [ "$L_" -lt "$shortest_len" ]; then shortest_len="$L_"; shortest_file="$base"; fi
    done
fi

# ---- 判斷化學版本 ----
CHEM="unknown"
IS_10X="否(不確定)"
case "$shortest_len" in
    24|26) CHEM="10x 3' v2 (16bp barcode + 10bp UMI)"; IS_10X="是" ;;
    27|28) CHEM="10x 3' v3/v3.1 (16bp barcode + 12bp UMI)"; IS_10X="是" ;;
    *)
        # 5' 或其他情況:barcode read 也是 16+10/12,但有時長度略不同
        if [ "$shortest_len" -ge 24 ] && [ "$shortest_len" -le 30 ]; then
            CHEM="疑似 10x(barcode read ${shortest_len}bp,非典型,建議用 --chemistry=auto)"
            IS_10X="很可能"
        else
            CHEM="非典型(barcode read ${shortest_len}bp,可能不是標準 10x)"
            IS_10X="否(不確定)"
        fi
        ;;
esac

log ""
log "  ── 判讀 ──────────────────────────────────"
log "  最短 read : $shortest_file (${shortest_len}bp)  -> 候選 barcode+UMI (R1)"
log "  最長 read : $longest_file (${longest_len}bp)  -> 候選 cDNA (R2)"
log "  是否 10x  : $IS_10X"
log "  化學版本  : $CHEM"

# 8bp 的通常是 sample index (I1)
if [ -n "$index_file" ]; then
    log "  $index_file (${LEN[$index_file]}bp) -> 候選 sample index (I1)"
fi

# =========================================================================
# 3) 物種偵測(選用):抽 cDNA read 比對 human / mouse
# =========================================================================
log ""
log "[3/4] 物種偵測 ..."

if ! $HAVE_MINIMAP; then
    log "  (略過) 未安裝 minimap2。如需自動判物種:conda install -c bioconda minimap2"
elif [ -z "$HUMAN_REF" ] && [ -z "$MOUSE_REF" ]; then
    log "  (略過) 未提供參考。設定 HUMAN_REF / MOUSE_REF 環境變數即可啟用。"
else
    CDNA="$PEEK_DIR/$longest_file"
    # 只抽 2000 條 cDNA read 來比對,夠看 mapping 比例
    SUB="$OUT_DIR/_cdna_subset.fq"
    awk 'NR%4==1{i++} i<=2000{print}' "$CDNA" > "$SUB"
    n_sub=$(( $(wc -l < "$SUB") / 4 ))

    map_rate() {  # $1 = ref fasta
        [ -z "$1" ] && { echo "NA"; return; }
        [ -e "$1" ] || { echo "NA(ref不存在)"; return; }
        local mapped
        mapped=$(minimap2 -ax map-ont -t "$THREADS" "$1" "$SUB" 2>/dev/null \
                  | awk '$1!~/^@/ && $3!="*"' | wc -l)
        awk -v m="$mapped" -v n="$n_sub" 'BEGIN{ if(n>0) printf "%.1f", 100*m/n; else print "NA" }'
    }

    log "  抽樣 cDNA reads: $n_sub 條"
    H_RATE="NA"; M_RATE="NA"
    [ -n "$HUMAN_REF" ] && { H_RATE="$(map_rate "$HUMAN_REF")"; log "    human mapping rate ≈ ${H_RATE}%"; }
    [ -n "$MOUSE_REF" ] && { M_RATE="$(map_rate "$MOUSE_REF")"; log "    mouse mapping rate ≈ ${M_RATE}%"; }

    SPECIES="unknown"
    if [ "$H_RATE" != "NA" ] && [ "$M_RATE" != "NA" ]; then
        SPECIES=$(awk -v h="$H_RATE" -v m="$M_RATE" \
            'BEGIN{ if(h>m && h>20) print "human"; else if(m>h && m>20) print "mouse"; else print "unknown(兩者皆低,檢查參考)" }')
    elif [ "$H_RATE" != "NA" ]; then
        SPECIES=$(awk -v h="$H_RATE" 'BEGIN{ print (h>20)?"human(僅測human)":"unknown" }')
    elif [ "$M_RATE" != "NA" ]; then
        SPECIES=$(awk -v m="$M_RATE" 'BEGIN{ print (m>20)?"mouse(僅測mouse)":"unknown" }')
    fi
    log "  ── 物種判定: $SPECIES"
    rm -f "$SUB"
fi

# =========================================================================
# 4) 給 Cell Ranger 的建議
# =========================================================================
log ""
log "[4/4] 建議 ──────────────────────────────────"
log "  • Cell Ranger reads 對應:"
log "      R1 (barcode+UMI) = 最短檔 $shortest_file"
log "      R2 (cDNA)        = 最長檔 $longest_file"
log "  • 跑 count 時用 --chemistry=auto 最保險(讓 Cell Ranger 自己確認 v2/v3)"
log "  • 參考基因組:human -> refdata-gex-GRCh38-2020-A"
log "                mouse -> refdata-gex-GRCm39-2024-A"
log ""
log "報告已寫出: $REPORT"

# 清掉抽樣 fastq(偵測完不需要),保留報告
if [ "${KEEP_PEEK:-false}" != "true" ]; then
    rm -rf "$PEEK_DIR"
fi
