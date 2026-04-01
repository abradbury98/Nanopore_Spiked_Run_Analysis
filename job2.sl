#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=2.Filtering-hostremove
#SBATCH --time=24:00:00
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=8
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
shopt -s nullglob

# -------------------------
# Config 
# -------------------------
WD=$(pwd)

INPUT_DIR="${WD}/raw_merged_reads_all"

# Intermediate dirs
PORECHOP_DIR="${WD}/porechop_2"
NANOFILT_DIR="${WD}/nanofilt_2"
FILTER_DIR="$NANOFILT_DIR"   # inputs for host removal

# Host removal outputs
OUT_DIR="${WD}/host_removed_reads_2"

# QC outputs
QC_DIR="${WD}/Quality_control_cleaned"
NANOPLOT_DIR="${QC_DIR}/Nanoplot"   # per-sample nanoplot/<sample>
NANOQC_DIR="${QC_DIR}/NanoQC"       # per-sample nanoQC/<sample>
MULTIQC_DIR="${QC_DIR}/MultiQC"     # multiqc outputs

# Single central log directory for all steps
LOG_DIR="${WD}/logs"

# References (edit)
CHICKEN_REF="${WD}/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna"
CTRL_REF="${WD}/DCS.fasta"
HUMAN_REF="${WD}/GCF_000001405.40_GRCh38.p14_genomic.fna"
COMBINED_REF="${WD}/combined_host_refs.fna"

# Resources
THREADS=${SLURM_CPUS_PER_TASK:-8}
MAX_JOBS=${MAX_JOBS:-4}

mkdir -p "$INPUT_DIR" "$PORECHOP_DIR" "$NANOFILT_DIR" "$OUT_DIR" "$NANOPLOT_DIR" "$NANOQC_DIR" "$MULTIQC_DIR" "$LOG_DIR"

# -------------------------
# Helpers
# -------------------------
log_step() { printf '\033[1;32m[%s] %s\033[0m\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }
err_exit() { echo "ERROR: $1" >&2; exit 1; }

count_reads_gz() {
  local f="$1"
  [[ -f "$f" ]] || { echo 0; return; }
  local lines
  lines=$(zcat -f "$f" | wc -l)
  echo $((lines/4))
}

log_module_version() {
  local name="$1"; shift
  if "$@" &>/dev/null; then
    echo "[$name] $("$@" 2>&1 | head -n 1)" >> "$LOG_DIR/module_versions.log"
  else
    echo "[$name] ERROR running version command" >> "$LOG_DIR/module_versions.log"
  fi
}

# -------------------------
# Validate refs (for hostremoval)
# -------------------------
check_refs() {
  [[ -f "$CHICKEN_REF" ]] || err_exit "Missing chicken reference: $CHICKEN_REF"
  [[ -f "$CTRL_REF" ]]   || err_exit "Missing control FASTA: $CTRL_REF"
  [[ -f "$HUMAN_REF" ]]  || err_exit "Missing human FASTA: $HUMAN_REF"
  if [[ ! -s "$COMBINED_REF" ]]; then
    log_step "Combining references into $COMBINED_REF"
    { sed -e '$a\' "$CHICKEN_REF"; sed -e '$a\' "$CTRL_REF"; sed -e '$a\' "$HUMAN_REF"; } > "$COMBINED_REF"
    log_step "Combined reference created."
  fi
}

# -------------------------
# Step functions
# -------------------------
step_porechop() {
  log_step "Step: porechop (adapter trimming)"
  module purge
  module load Porechop/0.2.4-gimkl-2022a-Python-3.11.3
  log_module_version "Porechop" porechop --version

  job_count=0
  for fq in "$INPUT_DIR"/*.fastq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    out="$PORECHOP_DIR/${sample}_porechop.fastq.gz"
    if [[ -f "$out" ]]; then
      log_step "Skipping porechop $sample (exists)"
      continue
    fi
    log_step "Porechop -> $out"
    porechop -i "$fq" -o "$out" > "$LOG_DIR/${sample}.porechop.log" 2>&1 &
    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then wait -n; job_count=$((job_count-1)); fi
  done
  wait
  log_step "porechop done"
}

step_nanofilt() {
  log_step "Step: nanofilt (quality filtering of porechop outputs)"
  module purge
  module load nanofilt/2.6.0-gimkl-2020a-Python-3.8.2
  log_module_version "NanoFilt" NanoFilt --version

  job_count=0
  for fq in "$PORECHOP_DIR"/*_porechop.fastq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    out="$NANOFILT_DIR/${sample}_nanofilt.fastq.gz"
    if [[ -f "$out" ]]; then
      log_step "Skipping NanoFilt $sample (exists)"
      continue
    fi
    log_step "NanoFilt $sample -> $out"
    gunzip -c "$fq" | NanoFilt -q 10 -l 250 | gzip > "$out" 2> "$LOG_DIR/${sample}.nanofilt.log" &
    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then wait -n; job_count=$((job_count-1)); fi
  done
  wait
  log_step "nanofilt done"
}

step_hostremoval() {
  log_step "Step: hostremoval (minimap2 -> remove mapped reads)"
  check_refs
  module purge
  module load minimap2/2.28-GCC-12.3.0 SAMtools/1.22-GCC-12.3.0 SeqKit/2.4.0
  log_module_version "minimap2" minimap2 --version
  log_module_version "samtools" samtools --version
  log_module_version "seqkit" seqkit version

  job_count=0
  for fq in "$FILTER_DIR"/*_nanofilt.fastq.gz "$FILTER_DIR"/*_porechop_nanofilt.fastq.gz; do
    [[ -f "$fq" ]] || continue
    base=$(basename "$fq" .fastq.gz)
    sample="${base%_porechop}"
    sample="${sample%_nanofilt}"
    out_fastq="$OUT_DIR/${sample}_host_removed.fastq.gz"
    stats_log="$LOG_DIR/${sample}.hostremoval.log"
    temp_sam="$LOG_DIR/${sample}.mapped.sam"
    mapped_ids="$LOG_DIR/${sample}.mapped_ids.txt"

    if [[ -f "$out_fastq" ]]; then
      log_step "Skipping host removal $sample (exists)"
      continue
    fi

    reads_before=$(count_reads_gz "$fq")
    log_step "Mapping $sample (reads: $reads_before)"
    minimap2 -ax map-ont --secondary=yes -t "$THREADS" "$COMBINED_REF" "$fq" > "$temp_sam" 2> "$LOG_DIR/${sample}.minimap2.log"

    samtools view -F 4 "$temp_sam" | cut -f1 | sort -u > "$mapped_ids"
    mapped_count=$(wc -l < "$mapped_ids" 2>/dev/null || echo 0)

    log_step "Removing $mapped_count mapped reads for $sample"
    if [[ "$mapped_count" -gt 0 ]]; then
      seqkit grep -v -f "$mapped_ids" "$fq" | gzip -c > "$out_fastq"
    else
      gzip -c "$fq" > "$out_fastq"
    fi

    reads_after=$(count_reads_gz "$out_fastq")
    {
      echo "Sample: $sample"
      echo "Reads before: $reads_before"
      echo "Mapped (removed): $mapped_count"
      echo "Reads after: $reads_after"
    } > "$stats_log"

    rm -f "$temp_sam"
    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then wait -n; job_count=$((job_count-1)); fi
  done
  wait
  log_step "hostremoval done"
}

step_nanoplot() {
  log_step "Step: nanoplot (run on Porechop-trimmed reads)"
  module purge
  module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
  log_module_version "NanoPlot" NanoPlot --version

  job_count=0
  for fq in "$NANOFILT_DIR"/*_porechop_nanofilt.fastq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    np_out_dir="$NANOPLOT_DIR/$sample"
    if [[ -f "$np_out_dir/NanoPlot-report.html" ]]; then
      log_step "Skipping NanoPlot $sample (exists)"
      continue
    fi
    mkdir -p "$np_out_dir"
    log_step "NanoPlot -> $np_out_dir"
    NanoPlot --fastq "$fq" --threads "$THREADS" --plots hex dot --maxlength 1000000 --outdir "$np_out_dir" > "$LOG_DIR/${sample}.nanoplot.log" 2>&1 &
    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then wait -n; job_count=$((job_count-1)); fi
  done
  wait
  log_step "nanoplot done"
}

step_nanoqc() {
  log_step "Step: nanoqc (run on Porechop-trimmed reads)"
  module purge
  module load nanoQC/0.9.4-gimkl-2022a-Python-3.10.5
  log_module_version "nanoQC" nanoQC --version || true

  job_count=0
  for fq in "$NANOFILT_DIR"/*_porechop_nanofilt.fastq.gz; do
    [[ -f "$fq" ]] || continue
    sample=$(basename "$fq" .fastq.gz)
    nanoqc_sample_dir="$NANOQC_DIR/$sample"
    mkdir -p "$nanoqc_sample_dir"
    # up-to-date check
    if [[ -f "$LOG_DIR/${sample}.nanoQC.log" && -n "$(ls -A "$nanoqc_sample_dir" 2>/dev/null || true)" ]]; then
      log_step "Skipping nanoQC $sample (already run)"
      continue
    fi
    log_step "nanoQC -> $nanoqc_sample_dir"
    (
      cd "$nanoqc_sample_dir" || exit 1
      nanoQC "$fq" > "$LOG_DIR/${sample}.nanoQC.log" 2>&1 || echo "nanoQC failed for $sample (see $LOG_DIR/${sample}.nanoQC.log)" >&2
    ) &
    job_count=$((job_count+1))
    if [ "$job_count" -ge "$MAX_JOBS" ]; then wait -n; job_count=$((job_count-1)); fi
  done
  wait
  log_step "nanoqc done"
}

step_multiqc() {
  log_step "Step: multiqc (aggregate NanoPlot + NanoQC)"
  module purge
  module load MultiQC/1.15-gimkl-2022a-Python-3.10.5
  log_module_version "MultiQC" multiqc --version

  # per-sample MultiQC
  for sample_dir in "$NANOPLOT_DIR"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")
    out_sub="$MULTIQC_DIR/$sample"
    mkdir -p "$out_sub"
    multiqc "$sample_dir" "$NANOQC_DIR/$sample" --force --filename "multiqc_${sample}.html" -o "$out_sub" > "$LOG_DIR/${sample}.multiqc.log" 2>&1 || true
  done

  # global MultiQC
  multiqc "$NANOPLOT_DIR" "$NANOQC_DIR" --force --filename "multiqc_all_samples.html" -o "$MULTIQC_DIR" > "$LOG_DIR/multiqc_all.log" 2>&1 || true

  log_step "multiqc done"
}

# -------------------------
# Argument parsing / dispatch
# -------------------------
usage() {
  cat <<EOF
Usage: $0 [--step NAME]
Where NAME is one of: porechop, nanofilt, hostremoval, nanoplot, nanoqc, multiqc, all
Default: all
Examples:
  sbatch $0 --step porechop
  sbatch $0 --step multiqc
EOF
  exit 1
}

STEP="all"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --step) STEP="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1"; usage ;;
  esac
done

case "$STEP" in
  porechop) step_porechop ;;
  nanofilt) step_nanofilt ;;
  hostremoval) step_hostremoval ;;
  nanoplot) step_nanoplot ;;
  nanoqc) step_nanoqc ;;
  multiqc) step_multiqc ;;
  all)
    step_porechop
    step_nanofilt
    step_hostremoval
    step_nanoplot
    step_nanoqc
    step_multiqc
    ;;
  *)
    echo "Unknown step: $STEP"
    usage
    ;;
esac

log_step "Requested step ($STEP) complete. Logs in $LOG_DIR"