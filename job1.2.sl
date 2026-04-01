#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=1.QC_Raw
#SBATCH --time=08:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=6
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
shopt -s nullglob

#------------------------------------------------------
# 🗂️ Setup
#------------------------------------------------------
WD=$(pwd)
INPUT_DIR="$WD/raw_merged_reads_all"
QC_DIR="$WD/Quality_control_raw_4"
LOG_DIR="$QC_DIR/Quality_control_logs_4"

mkdir -p "$QC_DIR" "$LOG_DIR" "$INPUT_DIR"

THREADS=${SLURM_CPUS_PER_TASK:-6}
MAX_JOBS=1

# Prefer pigz for parallel compression if available
if command -v pigz >/dev/null 2>&1; then
  COMPRESS_CMD="pigz -p $THREADS -c"
else
  COMPRESS_CMD="gzip -c"
fi

#------------------------------------------------------
# 🎨 Logging helper
#------------------------------------------------------
log_step() { echo -e "\033[1;32m[$(date '+%Y-%m-%d %H:%M:%S')] $1\033[0m"; }

#------------------------------------------------------
# 🧾 Version log
#------------------------------------------------------
VER_FILE="$LOG_DIR/pipeline_versions_$(date '+%Y%m%d_%H%M%S').log"
echo "📜 QC Pipeline Log - $(date)" > "$VER_FILE"

# Improved version logger: capture stdout/stderr once and record result
log_module_version() {
  local name="$1"; shift
  local out
  if out=$("$@" 2>&1); then
    echo "[$name] $(echo "$out" | head -n 1)" >> "$VER_FILE"
  else
    echo "[$name] ERROR: $(echo "$out" | head -n 1)" >> "$VER_FILE"
  fi
}

# ====================================================================================
# Step 0: Merge barcode fastq files
# ====================================================================================
log_step "🔄 Step 0: Merging barcode fastq.gz files"

for folder in barcode*; do
    if [[ -d "$folder" ]]; then
        barcode_number=$(echo "$folder" | sed -E 's/barcode([0-9]+)/\1/')
        formatted_barcode=$(printf "%02d" $((10#$barcode_number)))

        out_fq="$INPUT_DIR/all_barcode${formatted_barcode}.fastq.gz"
        echo "Merging $folder/*.fastq.gz → $out_fq"

        # Collect matching files safely (nullglob set, so pattern becomes empty array if none)
        files=( "$folder"/*.fastq.gz )
        if (( ${#files[@]} == 0 )); then
            echo "⚠️  No FASTQ files found in $folder — skipping"
            continue
        fi

        # Use zcat on the list of files and compress with pigz/gzip
        zcat "${files[@]}" | $COMPRESS_CMD > "$out_fq"
    fi
done

log_step "✅ Step 0 complete."

#------------------------------------------------------
# 🧪 Input check
#------------------------------------------------------
fastq_files=( "$INPUT_DIR"/*.fastq.gz )
if (( ${#fastq_files[@]} == 0 )); then
  echo "❌ No FASTQ files found in $INPUT_DIR"
  exit 1
fi

#======================================================
# Function: Run QC suite
#======================================================
run_qc_suite() {
  local stage="$1"
  local input_dir="$2"
  local output_dir="$3"

  log_step "🚀 Starting QC suite for $stage reads"

  local NANOPLOT_DIR="$output_dir/nanoplot"
  local MULTIQC_DIR="$output_dir/multiqc"

  mkdir -p "$NANOPLOT_DIR" "$MULTIQC_DIR" "$LOG_DIR"

  #======================================================
  # NanoPlot
  #======================================================
  module purge
  module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
  log_module_version "NanoPlot ($stage)" NanoPlot --version

  job_count=0

  for fq in "$input_dir"/*.fastq.gz; do
    # If nullglob removed the pattern then loop won't run; this guard is mostly defensive
    [[ -e "$fq" ]] || continue

    sample=$(basename "$fq" .fastq.gz)
    out_dir="$NANOPLOT_DIR/$sample"
    [[ -f "$out_dir/NanoPlot-report.html" ]] && continue
    mkdir -p "$out_dir"

    NanoPlot \
      --fastq "$fq" \
      --threads "$THREADS" \
      --plots hex dot \
      --maxlength 1000000 \
      --outdir "$out_dir" \
      >"$LOG_DIR/${sample}_${stage}_NanoPlot.log" 2>&1 &

    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
      # Wait for at least one background job to finish before starting more
      wait -n
      job_count=$((job_count-1))
    fi
  done

  # Wait for any remaining background NanoPlot jobs to complete
  wait

  #======================================================
  # MultiQC
  #======================================================
  module purge
  module load MultiQC/1.15-gimkl-2022a-Python-3.10.5
  log_module_version "MultiQC ($stage)" multiqc --version

  INDEX="$MULTIQC_DIR/index.html"

  echo "<html><head><title>${stage} QC Dashboard</title>
  <style>
      body{font-family:sans-serif;}
      a{margin:8px;padding:6px;background:#eef;border-radius:5px;text-decoration:none;}
      iframe{width:100%;height:90vh;border:none;}
  </style>
  </head><body><h1>📊 ${stage} QC Dashboard</h1>" > "$INDEX"

  for fq in "$input_dir"/*.fastq.gz; do
    [[ -e "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    out_sub="$MULTIQC_DIR/$sample"

    mkdir -p "$out_sub"

    # Run multiqc on the NanoPlot output for the sample. Capture logs.
    multiqc "$NANOPLOT_DIR/$sample" \
      --force --filename "multiqc_${sample}.html" \
      -o "$out_sub" \
      >"$LOG_DIR/${sample}_${stage}_multiqc.log" 2>&1 || {
        echo "⚠️ MultiQC failed for $sample (see $LOG_DIR/${sample}_${stage}_multiqc.log)"
      }

    echo "<a href=\"${sample}/multiqc_${sample}.html\" target=\"content\">${sample}</a>" >> "$INDEX"
  done

  # Global multiqc across all samples
  multiqc "$NANOPLOT_DIR" \
    --force --filename "multiqc_all_barcodes.html" \
    -o "$MULTIQC_DIR" \
    >"$LOG_DIR/multiqc_${stage}_all.log" 2>&1 || {
      echo "⚠️ MultiQC (all) failed — see $LOG_DIR/multiqc_${stage}_all.log"
    }

  echo "<hr><iframe name=\"content\" src=\"multiqc_all_barcodes.html\"></iframe></body></html>" >> "$INDEX"

  log_step "✅ $stage QC complete"
}

# Run QC suite
run_qc_suite "raw" "$INPUT_DIR" "$QC_DIR"