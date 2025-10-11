#!/usr/bin/env bash
# --------------------------------------------------------------
# Author: Md. Jubayer Hossain
# Organization: DeepBio Limited / CHIRAL Bangladesh
# Project: RNA-seq Quantification with Salmon
# Description: SRA ? FASTQ ? Salmon ? Quantification
# Date: 2025-10-06
# Version: 1.2 (repo-aware paths)

# --------------------------------------------------------------
# Data 
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
# SRA link: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=229998

# --------------------------------------------------------------
# RNA-seq Quantification with Salmon (repo-root aware)
# Layout expected:
#   repo/
#    +- input/SRR_Acc_List.txt
#    +- scripts/00_salmon_pipeline.sh  <-- this file
#    +- (generated at repo root) sra/, fastq/, QC_fastqc/, QC_multiqc/, Salmon.out/, logs/


# --------------------------------------------------------------
# Expected layout:
#   repo/
#    +- input/srr_acc_list.txt
#    +- scripts/00_salmon_pipeline.sh
#    +- (generated at repo root)
#       +- sra/
#       +- fastq/
#       +- qc_fastqc/
#       +- qc_multiqc/
#       +- salmon.out/
#       +- logs/


# ------- Locate repo root (parent of this script) -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(cd "$SCRIPT_DIR/.." && pwd)"       # repo root regardless of CWD
cd "$WORKDIR"

# ------- User config (defaults; override via env if needed) -------------------
INPUT_DIR="${INPUT_DIR:-$WORKDIR/input}"
SAMPLES_FILE="${SAMPLES_FILE:-$INPUT_DIR/SRR_Acc_List.txt}"

SRA_DIR="${SRA_DIR:-$WORKDIR/sra}"
FASTQ_DIR="${FASTQ_DIR:-$WORKDIR/fastq}"
QC_DIR="${QC_DIR:-$WORKDIR/qc_fastqc}"
QC_SUMMARY_DIR="${QC_SUMMARY_DIR:-$WORKDIR/qc_multiqc}"
QUANT_DIR="${QUANT_DIR:-$WORKDIR/salmon.out}"
LOG_DIR="${LOG_DIR:-$WORKDIR/logs}"

# Salmon index (adjust to your installed index)
SALMON_INDEX="${SALMON_INDEX:-$HOME/hg38/salmon_partial_sa_index/default}"

THREADS="${THREADS:-12}"
DO_QC="${DO_QC:-1}"        # 1 = run FastQC + MultiQC, 0 = skip QC
FORCE="${FORCE:-0}"        # 1 = re-run even if quant.sf exists
PIGZ="${PIGZ:-1}"          # 1 = compress FASTQs (.gz)
SALMON_OPTS="${SALMON_OPTS:---validateMappings --seqBias --gcBias --posBias}"
NCBI_PREFETCH_CART="${NCBI_PREFETCH_CART:-}"

# ------- Setup directories & logging -----------------------------------------
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$QC_DIR" "$QC_SUMMARY_DIR" "$QUANT_DIR" "$LOG_DIR"
RUN_ID="$(date -u +'%Y%m%dT%H%M%SZ')"
MASTER_LOG="$LOG_DIR/run_$RUN_ID.log"
exec > >(tee -a "$MASTER_LOG") 2>&1

echo ">>> Start: $RUN_ID"
echo "Repo root : $WORKDIR"
echo "SRR list  : $SAMPLES_FILE"
echo "Threads   : $THREADS | QC: $DO_QC | FORCE: $FORCE | pigz: $PIGZ"
echo "Salmon idx: $SALMON_INDEX"

# ------- Tool checks + version capture ---------------------------------------
REQ=(prefetch fasterq-dump salmon)
[[ "$DO_QC" -eq 1 ]] && REQ+=(fastqc multiqc)
[[ "$PIGZ" -eq 1 ]] && REQ+=(pigz)

for c in "${REQ[@]}"; do
  command -v "$c" >/dev/null || { echo "ERROR: '$c' not in PATH"; exit 1; }
done

VERS_FILE="$LOG_DIR/tool_versions_$RUN_ID.txt"
{
  echo "UTC: $(date -u +'%F %T')"
  prefetch --version 2>&1 | head -n1 || true
  fasterq-dump --version 2>&1 | head -n1 || true
  salmon --version 2>&1 || true
  if [[ "$DO_QC" -eq 1 ]]; then fastqc --version 2>&1 || true; multiqc --version 2>&1 || true; fi
  if [[ "$PIGZ" -eq 1 ]]; then pigz --version 2>&1 || true; fi
} | tee "$VERS_FILE"

# ------- Read SRR list (ignore blanks/comments) ------------------------------
if [[ ! -s "$SAMPLES_FILE" ]]; then
  echo "ERROR: SRR list not found or empty: $SAMPLES_FILE"; exit 1;
fi
mapfile -t SAMPLES < <(grep -v -E '^\s*(#|$)' "$SAMPLES_FILE" | tr -d '\r')
[[ "${#SAMPLES[@]}" -gt 0 ]] || { echo "ERROR: No SRR IDs in $SAMPLES_FILE"; exit 1; }

# ------- Helpers --------------------------------------------------------------
finish(){ local c=$?; [[ $c -eq 0 ]] && echo "SUCCESS $(date -u +'%F %T')" || echo "FAILED (code $c)"; exit $c; }
trap finish EXIT

download_sra(){ # echo path to .sra
  local srr="$1" out="$SRA_DIR/$srr"
  mkdir -p "$out"
  if [[ -n "$NCBI_PREFETCH_CART" ]]; then
    prefetch --cart "$NCBI_PREFETCH_CART" --output-directory "$out" "$srr"
  else
    prefetch --output-directory "$out" "$srr"
  fi
  echo "$out/$srr.sra"
}

sra_to_fastq(){ # echo "R1 R2" or "SE -"
  local srr="$1" sra="$2" fqdir="$FASTQ_DIR/$srr"
  mkdir -p "$fqdir"
  # --split-3: paired ? _1/_2; single ? single file
  fasterq-dump --split-3 --threads "$THREADS" --temp "$fqdir" --outdir "$fqdir" "$sra"
  if [[ "$PIGZ" -eq 1 ]]; then
    find "$fqdir" -maxdepth 1 -type f -name "*.fastq" -print0 | xargs -0 -r pigz -p "$THREADS" -f
    local R1="$fqdir/${srr}_1.fastq.gz"
    local R2="$fqdir/${srr}_2.fastq.gz"
    local SE="$fqdir/${srr}.fastq.gz"
  else
    local R1="$fqdir/${srr}_1.fastq"
    local R2="$fqdir/${srr}_2.fastq"
    local SE="$fqdir/${srr}.fastq"
  fi
  if [[ -s "$R1" && -s "$R2" ]]; then
    echo "$R1 $R2"
  elif [[ -s "$SE" ]]; then
    echo "$SE -"
  else
    echo "ERROR: No FASTQ files found for $srr in $fqdir" >&2
    return 1
  fi
}

run_fastqc(){
  local f1="$1" f2="$2"
  [[ "$DO_QC" -eq 1 ]] || return 0
  if [[ "$f2" == "-" ]]; then fastqc -t "$THREADS" -o "$QC_DIR" "$f1"
  else fastqc -t "$THREADS" -o "$QC_DIR" "$f1" "$f2"
  fi
}

run_multiqc(){ [[ "$DO_QC" -eq 1 ]] && multiqc "$QC_DIR" -o "$QC_SUMMARY_DIR"; }

run_salmon(){
  local srr="$1" f1="$2" f2="$3" out="$QUANT_DIR/$srr"
  mkdir -p "$out"
  if [[ "$f2" == "-" ]]; then
    salmon quant -l A -r "$f1" -i "$SALMON_INDEX" -p "$THREADS" $SALMON_OPTS -o "$out"
  else
    salmon quant -l A -1 "$f1" -2 "$f2" -i "$SALMON_INDEX" -p "$THREADS" $SALMON_OPTS -o "$out"
  fi
}

# ------- Main loop ------------------------------------------------------------
for srr in "${SAMPLES[@]}"; do
  echo "------------------------------------------------------------"
  echo ">>> $srr  ($(date -u +'%F %T'))"
  SAMPLE_LOG="$LOG_DIR/${srr}_$RUN_ID.log"
  {
    set -x
    if [[ "$FORCE" -eq 0 && -s "$QUANT_DIR/$srr/quant.sf" ]]; then
      echo "quant.sf exists; skipping $srr (FORCE=1 to re-run)"; set +x; continue
    fi
    sra_path="$(download_sra "$srr")"
    read -r fq1 fq2 < <(sra_to_fastq "$srr" "$sra_path")
    run_fastqc "$fq1" "$fq2"
    run_salmon "$srr" "$fq1" "$fq2"
    set +x
  } 2>&1 | tee -a "$SAMPLE_LOG"
done

# ------- One-time MultiQC -----------------------------------------------------
if [[ "$DO_QC" -eq 1 && -n "$(ls -A "$QC_DIR" 2>/dev/null || true)" ]]; then
  echo "Generating MultiQC summary..."
  run_multiqc
fi

echo "Outputs at repo root:"
echo "  sra/          (downloaded .sra)"
echo "  fastq/<SRR>/  (FASTQs)"
echo "  salmon.out/<SRR>/quant.sf"
echo "  logs/"
[[ "$DO_QC" -eq 1 ]] && { echo "  qc_fastqc/"; echo "  qc_multiqc/"; }
