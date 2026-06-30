#!/usr/bin/env bash
###############################################################################
## sra_to_cellranger.sh
## ---------------------------------------------------------------------------
## 一條龍:SRA  ->  10x FASTQ(改成 Cell Ranger 命名)  ->  cellranger count
##         -> 產出 barcodes/features/matrix(可直接接 01_QC_PerSample.R)
##
## 對應你的資料結構:
##   <SRC_ROOT>/<PRJNA...>/<SRR...>            (prefetch 下載的 .sra,可能無副檔名)
##   或 <SRC_ROOT>/<PRJNA...>/<SRR...>.sra
##
## 主要處理:
##   1. fasterq-dump 把 .sra 拆成 R1/R2(/I1),--include-technical 保留 barcode read
##   2. 自動判斷哪個檔是 barcode(短)/cDNA(長)/index,改成 Cell Ranger 要求的命名:
##        <Sample>_S1_L001_R1_001.fastq.gz   (barcode+UMI)
##        <Sample>_S1_L001_R2_001.fastq.gz   (cDNA)
##        <Sample>_S1_L001_I1_001.fastq.gz   (sample index,有才放)
##   3. cellranger count --chemistry=auto 跑比對定量
##   4. 把 filtered_feature_bc_matrix 整理到輸出目錄,接 QC 腳本
##
## 特性:
##   * 批次:可一次處理多個 SRR(或整個 PRJNA 資料夾)
##   * 可續跑:已有 cellranger outs/ 的樣本自動跳過
##   * 物種:可自動偵測(需 inspect_sra.sh 同目錄 + minimap2 + 參考),
##           或用 --species human/mouse 手動指定
##   * 失敗的樣本記錄後跳過,不中斷整批
##
## 用法:
##   # 單筆,手動指定物種與參考
##   ./sra_to_cellranger.sh \
##       --species human \
##       --ref /ref/refdata-gex-GRCh38-2020-A \
##       /data/PRJNA772373/SRR16475068
##
##   # 整個 PRJNA 資料夾(每個 SRR 一個樣本),自動偵測物種
##   ./sra_to_cellranger.sh \
##       --human-ref /ref/refdata-gex-GRCh38-2020-A \
##       --mouse-ref /ref/refdata-gex-GRCm39-2024-A \
##       --auto-species \
##       /data/PRJNA772373/
##
## 依賴:sra-tools(fasterq-dump)、cellranger(在 PATH)、pigz 或 gzip
##       自動物種偵測另需 minimap2 + transcriptome FASTA(給 inspect_sra.sh 用)
###############################################################################
set -uo pipefail   # 注意:不開 -e,讓單一樣本失敗不會中斷整批

# ----------------------------- 預設參數 ---------------------------------- #
SPECIES=""                 # human / mouse;留空且開 --auto-species 則自動測
REF=""                     # 直接指定 Cell Ranger transcriptome 路徑(優先)
HUMAN_REF=""               # human Cell Ranger reference
MOUSE_REF=""               # mouse Cell Ranger reference
AUTO_SPECIES=false
OUT_ROOT="./cellranger_out"
THREADS="${THREADS:-8}"
LOCALMEM="${LOCALMEM:-32}" # Cell Ranger localmem (GB)
CHEMISTRY="auto"
KEEP_FASTQ=false           # 跑完是否保留中間 fastq(預設刪,省空間)
EXPECT_CELLS=""            # 選用:--expect-cells N

# 物種偵測用的 transcriptome(給 inspect_sra.sh),選用
INSPECT_HUMAN_REF="${INSPECT_HUMAN_REF:-}"
INSPECT_MOUSE_REF="${INSPECT_MOUSE_REF:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    sed -n '2,60p' "$0" | sed 's/^## \{0,1\}//'
    exit "${1:-0}"
}

# ----------------------------- 解析參數 ---------------------------------- #
INPUTS=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --species)      SPECIES="$2"; shift 2 ;;
        --ref)          REF="$2"; shift 2 ;;
        --human-ref)    HUMAN_REF="$2"; shift 2 ;;
        --mouse-ref)    MOUSE_REF="$2"; shift 2 ;;
        --auto-species) AUTO_SPECIES=true; shift ;;
        --out)          OUT_ROOT="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --localmem)     LOCALMEM="$2"; shift 2 ;;
        --chemistry)    CHEMISTRY="$2"; shift 2 ;;
        --expect-cells) EXPECT_CELLS="$2"; shift 2 ;;
        --keep-fastq)   KEEP_FASTQ=true; shift ;;
        -h|--help)      usage 0 ;;
        -*)             echo "未知參數: $1" >&2; usage 1 ;;
        *)              INPUTS+=("$1"); shift ;;
    esac
done

if [ "${#INPUTS[@]}" -eq 0 ]; then
    echo "錯誤: 請提供至少一個 SRR 路徑或一個 PRJNA 資料夾" >&2
    usage 1
fi

# ----------------------------- 工具檢查 ---------------------------------- #
need() { command -v "$1" >/dev/null 2>&1; }
for t in fasterq-dump cellranger; do
    if ! need "$t"; then
        echo "錯誤: 找不到 '$t'。" >&2
        [ "$t" = "fasterq-dump" ] && echo "  conda install -c bioconda sra-tools" >&2
        [ "$t" = "cellranger" ]   && echo "  請依 10x 官網下載 Cell Ranger 並加入 PATH" >&2
        exit 1
    fi
done
GZIP="gzip"; need pigz && GZIP="pigz -p $THREADS"

mkdir -p "$OUT_ROOT"
LOG_DIR="$OUT_ROOT/_logs"; mkdir -p "$LOG_DIR"
SUMMARY="$LOG_DIR/run_summary.csv"
if [ ! -f "$SUMMARY" ]; then
    echo "srr,sample,species,ref_used,chemistry,status,n_cells,note,time" > "$SUMMARY"
fi
ts() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(ts)] $*"; }

# ----------------------------- 展開輸入清單 ------------------------------ #
## 把資料夾展開成裡面的每個 SRR;單一檔/資料夾(本身是 SRR)直接收進來
SRR_PATHS=()
for inp in "${INPUTS[@]}"; do
    inp="${inp%/}"
    base="$(basename "$inp")"
    if [[ "$base" =~ ^(SRR|ERR|DRR)[0-9]+ ]]; then
        # 本身就是一個 SRR(檔案或資料夾)
        SRR_PATHS+=("$inp")
    elif [ -d "$inp" ]; then
        # 視為 PRJNA 資料夾,撈裡面的 SRR/ERR/DRR
        while IFS= read -r p; do SRR_PATHS+=("$p"); done < <(
            find "$inp" -maxdepth 1 \( -type f -o -type d \) \
                 -regextype posix-extended \
                 -regex '.*/(SRR|ERR|DRR)[0-9]+(\.sra)?$' | sort -u
        )
    else
        log "! 跳過無法辨識的輸入: $inp"
    fi
done

if [ "${#SRR_PATHS[@]}" -eq 0 ]; then
    log "! 找不到任何 SRR,結束"
    exit 1
fi
log "共 ${#SRR_PATHS[@]} 筆要處理"

# ----------------------------- 物種 -> 參考 對應 ------------------------- #
resolve_ref() {  # $1 = species ; echo ref path
    case "$1" in
        human) [ -n "$HUMAN_REF" ] && echo "$HUMAN_REF" || echo "" ;;
        mouse) [ -n "$MOUSE_REF" ] && echo "$MOUSE_REF" || echo "" ;;
        *)     echo "" ;;
    esac
}

# 自動偵測物種(呼叫 inspect_sra.sh,解析它報告裡的「物種判定」)
detect_species() {  # $1 = sra path ; echo human/mouse/unknown
    local sra="$1"
    if [ ! -x "$SCRIPT_DIR/inspect_sra.sh" ] && [ ! -f "$SCRIPT_DIR/inspect_sra.sh" ]; then
        echo "unknown"; return
    fi
    local tmp; tmp="$(mktemp -d)"
    HUMAN_REF="$INSPECT_HUMAN_REF" MOUSE_REF="$INSPECT_MOUSE_REF" \
        bash "$SCRIPT_DIR/inspect_sra.sh" "$sra" "$tmp" >/dev/null 2>&1 || true
    local sp="unknown"
    if [ -f "$tmp/inspect_report.txt" ]; then
        sp="$(grep -oE '物種判定: ?(human|mouse)' "$tmp/inspect_report.txt" \
               | grep -oE '(human|mouse)' | head -1 || true)"
        [ -z "$sp" ] && sp="unknown"
    fi
    rm -rf "$tmp"
    echo "$sp"
}

# =========================================================================
# 主迴圈:逐筆處理
# =========================================================================
for SRA in "${SRR_PATHS[@]}"; do
    SRR="$(basename "$SRA")"; SRR="${SRR%.sra}"
    SAMPLE="$SRR"                          # Cell Ranger sample 名就用 SRR
    WORK="$OUT_ROOT/$SRR"
    FASTQ_DIR="$WORK/fastq"
    CR_ID="${SRR}_cellranger"              # cellranger --id
    OUTS="$WORK/$CR_ID/outs"
    MATRIX_OUT="$WORK/filtered_feature_bc_matrix"

    log "========================================================"
    log "處理 $SRR"

    # ---- 續跑檢查 ----
    if [ -d "$OUTS/filtered_feature_bc_matrix" ] || [ -d "$MATRIX_OUT" ]; then
        log "  [skip] 已有 Cell Ranger 輸出,跳過"
        echo "$SRR,$SAMPLE,-,-,-,SKIPPED,-,已有輸出,$(ts)" >> "$SUMMARY"
        continue
    fi

    mkdir -p "$FASTQ_DIR"

    # ---- 決定物種與參考 ----
    sp="$SPECIES"; ref="$REF"
    if [ -z "$ref" ]; then
        if [ -z "$sp" ] && $AUTO_SPECIES; then
            log "  偵測物種中 ..."
            sp="$(detect_species "$SRA")"
            log "  偵測結果: $sp"
        fi
        ref="$(resolve_ref "$sp")"
    fi
    if [ -z "$ref" ] || [ ! -d "$ref" ]; then
        log "  ! 無可用參考基因組 (species=$sp, ref=$ref),跳過此樣本"
        echo "$SRR,$SAMPLE,$sp,NONE,$CHEMISTRY,FAILED,-,無參考基因組,$(ts)" >> "$SUMMARY"
        continue
    fi
    log "  物種=$sp  參考=$ref"

    # ---- (1) fasterq-dump 拆檔 ----
    log "  [1] fasterq-dump 拆檔 ..."
    if ! fasterq-dump "$SRA" \
            --split-files --include-technical \
            --threads "$THREADS" --outdir "$FASTQ_DIR" --force \
            > "$LOG_DIR/${SRR}_fasterqdump.log" 2>&1; then
        log "  ! fasterq-dump 失敗,見 ${SRR}_fasterqdump.log"
        echo "$SRR,$SAMPLE,$sp,$ref,$CHEMISTRY,FAILED,-,fasterq-dump失敗,$(ts)" >> "$SUMMARY"
        $KEEP_FASTQ || rm -rf "$FASTQ_DIR"
        continue
    fi

    # ---- (2) 判斷 R1/R2/I1 並改名 ----
    log "  [2] 判斷 read 角色並改成 Cell Ranger 命名 ..."
    mapfile -t RAW < <(ls "$FASTQ_DIR"/*.fastq 2>/dev/null | sort)
    if [ "${#RAW[@]}" -lt 2 ]; then
        log "  ! 拆出的 fastq 少於 2 個,可能不是雙端 10x,跳過"
        echo "$SRR,$SAMPLE,$sp,$ref,$CHEMISTRY,FAILED,-,fastq數<2,$(ts)" >> "$SUMMARY"
        $KEEP_FASTQ || rm -rf "$FASTQ_DIR"
        continue
    fi

    modal_len() { awk 'NR%4==2{print length($0)}' "$1" | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}'; }

    # 算每個檔的長度,判斷角色。
    # 注意:index read(~8bp)比 barcode read(26/28bp)更短,
    #       所以「最短=barcode」是錯的!正確順序:
    #         cDNA(R2)    = 最長
    #         index(I1)   = 8~10bp
    #         barcode(R1) = 20~32bp(10x v2=26 / v3=28),排除 cDNA 與 index 後最短
    declare -A L
    for f in "${RAW[@]}"; do L["$f"]="$(modal_len "$f")"; done
    r1=""; r2=""; i1=""; smax=0
    # (a) cDNA = 最長
    for f in "${RAW[@]}"; do l="${L[$f]}"; (( l > smax )) && { smax=$l; r2="$f"; }; done
    # (b) index = 8~10bp(且非 cDNA)
    for f in "${RAW[@]}"; do
        l="${L[$f]}"; [ "$f" = "$r2" ] && continue
        if (( l >= 8 && l <= 10 )); then i1="$f"; fi
    done
    # (c) barcode = 20~32bp,排除 cDNA/index 後取最短
    smin=999999
    for f in "${RAW[@]}"; do
        l="${L[$f]}"; [ "$f" = "$r2" ] && continue; [ "$f" = "$i1" ] && continue
        if (( l >= 20 && l <= 32 )) && (( l < smin )); then smin=$l; r1="$f"; fi
    done
    # (d) 後備:若 barcode 長度不在典型區間,改取「排除 cDNA/index 後最短」
    if [ -z "$r1" ]; then
        smin=999999
        for f in "${RAW[@]}"; do
            l="${L[$f]}"; [ "$f" = "$r2" ] && continue; [ "$f" = "$i1" ] && continue
            (( l < smin )) && { smin=$l; r1="$f"; }
        done
    fi

    if [ -z "$r1" ] || [ -z "$r2" ]; then
        log "  ! 無法判定 R1/R2 角色,跳過(建議先用 inspect_sra.sh 檢查)"
        echo "$SRR,$SAMPLE,$sp,$ref,$CHEMISTRY,FAILED,-,無法判定R1R2,$(ts)" >> "$SUMMARY"
        $KEEP_FASTQ || rm -rf "$FASTQ_DIR"
        continue
    fi

    log "      R1(barcode ${smin}bp) = $(basename "$r1")"
    log "      R2(cDNA ${smax}bp)    = $(basename "$r2")"
    [ -n "$i1" ] && log "      I1(index ${L[$i1]}bp)  = $(basename "$i1")"

    # sanity:barcode 長度應在 24~30 之間(10x v2/v3)
    if (( smin < 20 || smin > 32 )); then
        log "  ! 警告:barcode read=${smin}bp,不像 10x v2/v3。仍嘗試,但建議先用 inspect_sra.sh 檢查"
    fi

    CR_FASTQ="$WORK/cr_fastq"; mkdir -p "$CR_FASTQ"
    log "      壓縮並改名 ..."
    $GZIP -c "$r1" > "$CR_FASTQ/${SAMPLE}_S1_L001_R1_001.fastq.gz"
    $GZIP -c "$r2" > "$CR_FASTQ/${SAMPLE}_S1_L001_R2_001.fastq.gz"
    [ -n "$i1" ] && $GZIP -c "$i1" > "$CR_FASTQ/${SAMPLE}_S1_L001_I1_001.fastq.gz"

    # 拆檔用的原始 fastq 不再需要
    rm -rf "$FASTQ_DIR"

    # ---- (3) cellranger count ----
    log "  [3] cellranger count(--chemistry=$CHEMISTRY) ..."
    EC_ARG=""; [ -n "$EXPECT_CELLS" ] && EC_ARG="--expect-cells=$EXPECT_CELLS"
    (
        cd "$WORK" || exit 1
        cellranger count \
            --id="$CR_ID" \
            --transcriptome="$ref" \
            --fastqs="$CR_FASTQ" \
            --sample="$SAMPLE" \
            --chemistry="$CHEMISTRY" \
            --localcores="$THREADS" \
            --localmem="$LOCALMEM" \
            $EC_ARG \
            --create-bam=true \
            > "$LOG_DIR/${SRR}_cellranger.log" 2>&1
    )
    CR_RC=$?

    if [ "$CR_RC" -ne 0 ] || [ ! -d "$OUTS/filtered_feature_bc_matrix" ]; then
        log "  ! cellranger 失敗 (rc=$CR_RC),見 ${SRR}_cellranger.log"
        echo "$SRR,$SAMPLE,$sp,$ref,$CHEMISTRY,FAILED,-,cellranger失敗,$(ts)" >> "$SUMMARY"
        $KEEP_FASTQ || rm -rf "$CR_FASTQ"
        continue
    fi

    # ---- (4) 整理輸出:把 matrix 複製到好找的位置 ----
    log "  [4] 整理輸出 ..."
    cp -r "$OUTS/filtered_feature_bc_matrix" "$MATRIX_OUT"
    # 細胞數(從 barcodes.tsv.gz 行數)
    bc="$MATRIX_OUT/barcodes.tsv.gz"
    n_cells="NA"
    [ -f "$bc" ] && n_cells="$(zcat "$bc" | wc -l)"

    $KEEP_FASTQ || rm -rf "$CR_FASTQ"

    log "  ✔ 完成 $SRR : 細胞數=$n_cells"
    log "    matrix -> $MATRIX_OUT"
    echo "$SRR,$SAMPLE,$sp,$ref,$CHEMISTRY,OK,$n_cells,,$(ts)" >> "$SUMMARY"
done

log "========================================================"
log "全部處理完畢。彙整表: $SUMMARY"
log ""
log "下一步:把各樣本的 filtered_feature_bc_matrix(barcodes/features/matrix)"
log "        放進你的 #Keloid/<GSE>/<GSM>/ 結構,接著跑 01_QC_PerSample.R"
