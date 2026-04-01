#!/usr/bin/env bash
#SBATCH --account=massey04083
#SBATCH --job-name=1.QC_with_nanoQC
#SBATCH --time=08:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=6
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err


set -euo pipefail
shopt -s nullglob

# -------------------------
# Config 
# -------------------------
WD=$(pwd)
INPUT_DIR="${WD}/raw_merged_reads_all"
QC_DIR="${WD}/Quality_control_raw"    

#   Nanoplot/, NanoQC/, MultiQC/, logs/
NANOPLOT_DIR="${QC_DIR}/Nanoplot"
NANOQC_DIR="${QC_DIR}/NanoQC"
MULTIQC_DIR="${QC_DIR}/MultiQC"
LOG_DIR="${QC_DIR}/logs"

THREADS=${SLURM_CPUS_PER_TASK:-6}
MAX_JOBS=${MAX_JOBS:-1}   # number of concurrent NanoPlot / nanoQC jobs

# create main directories
mkdir -p "$INPUT_DIR" "$QC_DIR" "$NANOPLOT_DIR" "$NANOQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# -------------------------
# Helpers
# -------------------------
log_step() { printf '\033[1;32m[%s] %s\033[0m\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }

# Capture version output once
log_module_version() {
  local name="$1"; shift
  local outfile="$LOG_DIR/pipeline_modules_$(date +%Y%m%d_%H%M%S).log"
  if out=$("$@" 2>&1); then
    echo "[$name] $(echo "$out" | head -n 1)" >> "$outfile"
  else
    echo "[$name] ERROR running version command" >> "$outfile"
  fi
}

# -------------------------
# QC Suite
# -------------------------
run_qc_suite() {
  local stage="$1"
  local input_dir="$2"
  local qc_dir="$3"

  log_step "🚀 Starting QC suite for $stage reads"
  # use pre-defined tool directories
  mkdir -p "$NANOPLOT_DIR" "$NANOQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

  #======================================================
  # NanoPlot
  #======================================================
  module purge
  module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
  log_module_version "NanoPlot ($stage)" NanoPlot --version

  job_count=0
  for fq in "$input_dir"/*.fastq.gz; do
    [[ -e "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    sample_out_dir="$NANOPLOT_DIR/$sample"
    [[ -f "$sample_out_dir/NanoPlot-report.html" ]] && { log_step "Skipping NanoPlot $sample (exists)"; continue; }
    mkdir -p "$sample_out_dir"

    log_step "NanoPlot: $sample -> $sample_out_dir"
    NanoPlot \
      --fastq "$fq" \
      --threads "$THREADS" \
      --plots hex dot \
      --maxlength 1000000 \
      --outdir "$sample_out_dir" \
      >"$LOG_DIR/${sample}_${stage}_NanoPlot.log" 2>&1 &

    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
      wait -n
      job_count=$((job_count-1))
    fi
  done
  wait

  #======================================================
  # nanoQC
  #======================================================
  module purge
  module load nanoQC/0.9.4-gimkl-2022a-Python-3.10.5
  log_module_version "nanoQC ($stage)" nanoQC --version || true

  job_count=0
  for fq in "$input_dir"/*.fastq.gz; do
    [[ -e "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    sample_out_dir="$NANOQC_DIR/$sample"
    mkdir -p "$sample_out_dir"

    # Skip if per-sample nanoQC log exists and output dir contains any file (simple up-to-date check)
    if [[ -f "$LOG_DIR/${sample}_${stage}_nanoQC.log" && -n "$(ls -A "$sample_out_dir" 2>/dev/null || true)" ]]; then
      log_step "Skipping nanoQC $sample (already run)"
      continue
    fi

    log_step "nanoQC: $sample -> $sample_out_dir"
    # nanoQC may write outputs into current working dir; run it from the sample dir to keep outputs local
    (
      cd "$sample_out_dir"
      nanoQC "$fq" >"$LOG_DIR/${sample}_${stage}_nanoQC.log" 2>&1 || {
        echo "nanoQC failed for $sample (see $LOG_DIR/${sample}_${stage}_nanoQC.log)" >&2
      }
    ) &

    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
      wait -n
      job_count=$((job_count-1))
    fi
  done
  wait

  #======================================================
  # MultiQC (per-sample and global)
  #======================================================
  module purge
  module load MultiQC/1.24.1-foss-2023a-Python-3.11.6
  log_module_version "MultiQC ($stage)" multiqc --version

  # Per-sample MultiQC pages (scan each sample's NanoPlot dir)
  for fq in "$input_dir"/*.fastq.gz; do
    [[ -e "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    sample_nanoplot_dir="$NANOPLOT_DIR/$sample"
    sample_multi_out="$MULTIQC_DIR/$sample"
    mkdir -p "$sample_multi_out"

    if [[ -d "$sample_nanoplot_dir" && -n "$(ls -A "$sample_nanoplot_dir" 2>/dev/null || true)" ]]; then
      log_step "MultiQC (per-sample): $sample -> $sample_multi_out"
      multiqc "$sample_nanoplot_dir" \
        --force --filename "multiqc_${sample}.html" \
        -o "$sample_multi_out" \
        >"$LOG_DIR/${sample}_${stage}_multiqc.log" 2>&1 || {
          echo "⚠️ MultiQC failed for $sample (see $LOG_DIR/${sample}_${stage}_multiqc.log)"
        }
    else
      log_step "Skipping MultiQC for $sample (no NanoPlot output found)"
    fi
  done

  # Global MultiQC across all NanoPlot outputs
  log_step "MultiQC (global) -> $MULTIQC_DIR"
  multiqc "$NANOPLOT_DIR" \
    --force --filename "multiqc_all_barcodes.html" \
    -o "$MULTIQC_DIR" \
    >"$LOG_DIR/multiqc_${stage}_all.log" 2>&1 || {
      echo "⚠️ MultiQC (all) failed — see $LOG_DIR/multiqc_${stage}_all.log"
    }

  log_step "✅ $stage QC complete"
}

# Run QC suite
run_qc_suite "raw" "$INPUT_DIR" "$QC_DIR"
