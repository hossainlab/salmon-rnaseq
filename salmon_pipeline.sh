#!/bin/bash
# --------------------------------------------------------------
# Author: Md. Jubayer Hossain
# Organization: DeepBio Limited / CHIRAL Bangladesh
# Project: RNA-seq Quantification with Salmon
# Description: Download SRA ? FASTQ ? Salmon quant with
#              per-step timing and CSV logs.
# Date: 2025-10-06
# --------------------------------------------------------------

set -euo pipefail

# ====== CONFIG ======
BASE_DIR="$HOME/RNAseq_Project"
RAW_DIR="$BASE_DIR/raw_data"
FASTQ_DIR="$BASE_DIR/fastq"
SALMON_DIR="$BASE_DIR/Salmon_out"
LOG_DIR="$BASE_DIR/logs"
INDEX_DIR="$HOME/hg38/salmon_partial_sa_index/default"

# SRA Samples
samples_to_quant=(SRR1039508 SRR1039509)

# ====== PREP ======
mkdir -p "$RAW_DIR" "$FASTQ_DIR" "$SALMON_DIR" "$LOG_DIR"

MASTER_CSV="$LOG_DIR/timings_master.csv"
if [[ ! -f "$MASTER_CSV" ]]; then
  echo "sample,step,start_iso,end_iso,duration_seconds" > "$MASTER_CSV"
fi

# Command presence checks (nice to have)
need_cmds=(prefetch fasterq-dump salmon)
for c in "${need_cmds[@]}"; do
  command -v "$c" >/dev/null 2>&1 || { echo "ERROR: '$c' not found in PATH"; exit 127; }
done

# ====== HELPERS ======
now_epoch() { date +%s; }
now_iso()   { date -u +"%Y-%m-%dT%H:%M:%SZ"; }   # UTC ISO-8601

log_step_csv() {
  local sample="$1" step="$2" start_iso="$3" end_iso="$4" dur="$5" sample_csv="$6"
  echo "$sample,$step,$start_iso,$end_iso,$dur" | tee -a "$sample_csv" >> "$MASTER_CSV"
}

# ====== MAIN LOOP ======
for sample in "${samples_to_quant[@]}"; do
  echo "========================================"
  echo "Processing sample: $sample"
  echo "========================================"

  # Paths per sample
  SAMPLE_RAW="$RAW_DIR/$sample"
  SAMPLE_FASTQ="$FASTQ_DIR/$sample"
  SAMPLE_OUT="$SALMON_DIR/$sample"
  SAMPLE_LOG_DIR="$LOG_DIR/$sample"
  SAMPLE_CSV="$SAMPLE_LOG_DIR/${sample}_timing.csv"
  SAMPLE_STD_LOG="$SAMPLE_LOG_DIR/${sample}.log"

  mkdir -p "$SAMPLE_RAW" "$SAMPLE_FASTQ" "$SAMPLE_OUT" "$SAMPLE_LOG_DIR"

  # Fresh CSV with header for this sample
  echo "sample,step,start_iso,end_iso,duration_seconds" > "$SAMPLE_CSV"

  sample_start_epoch=$(now_epoch)
  sample_start_iso=$(now_iso)

  # ---------- Step 1: Prefetch ----------
  echo "[1/3] Prefetch: $sample"
  s1_start_iso=$(now_iso); s1_start_epoch=$(now_epoch)
  {
    prefetch -O "$SAMPLE_RAW/" "$sample"
  } >>"$SAMPLE_STD_LOG" 2>&1
  s1_end_iso=$(now_iso); s1_end_epoch=$(now_epoch)
  s1_dur=$(( s1_end_epoch - s1_start_epoch ))
  log_step_csv "$sample" "prefetch" "$s1_start_iso" "$s1_end_iso" "$s1_dur" "$SAMPLE_CSV"

  # ---------- Step 2: fasterq-dump ----------
  echo "[2/3] fasterq-dump: $sample"
  s2_start_iso=$(now_iso); s2_start_epoch=$(now_epoch)
  {
    fasterq-dump -e 14 -p -O "$SAMPLE_FASTQ/" "$SAMPLE_RAW/$sample/$sample.sra"
  } >>"$SAMPLE_STD_LOG" 2>&1
  s2_end_iso=$(now_iso); s2_end_epoch=$(now_epoch)
  s2_dur=$(( s2_end_epoch - s2_start_epoch ))
  log_step_csv "$sample" "fasterq-dump" "$s2_start_iso" "$s2_end_iso" "$s2_dur" "$SAMPLE_CSV"

  # FASTQ paths (NCBI naming style: ${sample}_1.fastq / ${sample}_2.fastq)
  R1="$SAMPLE_FASTQ/${sample}_1.fastq"
  R2="$SAMPLE_FASTQ/${sample}_2.fastq"

  # ---------- Step 3: Salmon quant ----------
  echo "[3/3] Salmon quant: $sample"
  s3_start_iso=$(now_iso); s3_start_epoch=$(now_epoch)
  {
    salmon quant \
      -l A \
      -1 "$R1" \
      -2 "$R2" \
      --validateMappings \
      -i "$INDEX_DIR" \
      -o "$SAMPLE_OUT" \
      -p 15
  } >>"$SAMPLE_STD_LOG" 2>&1
  s3_end_iso=$(now_iso); s3_end_epoch=$(now_epoch)
  s3_dur=$(( s3_end_epoch - s3_start_epoch ))
  log_step_csv "$sample" "salmon_quant" "$s3_start_iso" "$s3_end_iso" "$s3_dur" "$SAMPLE_CSV"

  # ---------- Totals ----------
  sample_end_epoch=$(now_epoch)
  sample_end_iso=$(now_iso)
  total_dur=$(( sample_end_epoch - sample_start_epoch ))
  log_step_csv "$sample" "total" "$sample_start_iso" "$sample_end_iso" "$total_dur" "$SAMPLE_CSV"

  echo "? Completed $sample | Total seconds: $total_dur"
  echo "Logs: $SAMPLE_STD_LOG"
  echo "Timing CSV: $SAMPLE_CSV"
done

echo "----------------------------------------"
echo "?? All samples processed."
echo "Master timing CSV: $MASTER_CSV"
echo "Per-sample logs & CSVs: $LOG_DIR/<sample>/"
echo "----------------------------------------"
