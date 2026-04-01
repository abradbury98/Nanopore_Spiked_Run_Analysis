#!/bin/bash -e
#SBATCH --account=massey04083
#SBATCH --job-name=3-assembly-contig_analysis
#SBATCH --time=08:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=6
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
shopt -s nullglob

# -------------------------
# Minimal helpers (standalone)
# -------------------------
log_step() { printf '\033[1;32m[%s] %s\033[0m\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }
log_module_version() {
  local name="$1"; shift
  if "$@" &>/dev/null; then
    echo "[$name] $("$@" 2>&1 | head -n 1)"
  else
    echo "[$name] ERROR running version command"
  fi
}

# -------------------------
# Required env vars (caller must set these)
# -------------------------
WD=$(pwd)

NANOFILT_DIR="${WD}/host_removed_reads_2"
FLYE_DIR="${WD}/flye_assembly"
THREADS=${SLURM_CPUS_PER_TASK:-8}

# central log dir for this script
LOG_DIR="${FLYE_DIR}/logs"
mkdir -p "$LOG_DIR"

# ============================
# Step 4: Flye assembly
# ============================
log_step "🧩 Step 4: Assembly with Flye"

module purge
module load Flye/2.9.5-foss-2023a-Python-3.11.6
log_module_version "Flye" flye --version | tee -a "$LOG_DIR/module_versions.txt"

mkdir -p "$FLYE_DIR"

# NOTE: running Flye for many samples in one job can be slow; consider per-sample jobs if assemblies are large.
for fq in "$NANOFILT_DIR"/*_porechop_host_removed.fastq.gz; do
    [[ -f "$fq" ]] || continue
    base=$(basename "$fq" porechop_host_removed.fastq.gz)
    outdir="$FLYE_DIR/flye_output_${base}"
    mkdir -p "$outdir"
    log_step "Running Flye for sample: $base -> $outdir"
    flye --nano-raw "$fq" --out-dir "$outdir" --meta --threads "$THREADS" > "$LOG_DIR/${base}.flye.log" 2>&1 || {
      log_step "⚠️ Flye failed for $base (see $LOG_DIR/${base}.flye.log)"
      continue
    }
done

log_step "✅ Step 4 complete: Assemblies saved in $FLYE_DIR."
unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

# ============================================
# Step 5: Seqmagick filtering and BBMap stats
# ============================================
log_step "🧹 Step 5.1: Seqmagick filtering"

module purge
module load seqmagick/0.8.4-gimkl-2020a-Python-3.8.2
log_module_version "Seqmagick" seqmagick --version | tee -a "$LOG_DIR/module_versions.txt"

export JAVA_TOOL_OPTIONS="-Xmx120g"

for dir in "$FLYE_DIR"/flye_output_*; do
    [[ -d "$dir" ]] || continue
    input_fasta="$dir/assembly.fasta"
    sample_dirname=$(basename "$dir")           # e.g. flye_output_<sample>
    sample_short="${sample_dirname#flye_output_}"
    output_fasta="$dir/${sample_short}_assembly.m1000.fasta"

    if [[ -f "$input_fasta" ]]; then
        log_step "Filtering contigs <1000bp for sample $sample_short"
        seqmagick convert --min-length 1000 "$input_fasta" "$output_fasta" > "$LOG_DIR/${sample_short}.seqmagick.log" 2>&1 || {
          log_step "⚠️ Seqmagick failed for $sample_short (see log)"
          continue
        }
        log_step "Original contigs: $(grep -c '^>' "$input_fasta" || true)"
        log_step "Filtered contigs: $(grep -c '^>' "$output_fasta" || true)"
    else
        log_step "⚠️ No assembly.fasta in $dir, skipping seqmagick."
    fi
done

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

# Start BBMap stats
log_step "📊 Step 5.2: BBMap stats (stats.sh)"

module purge
module load minimap2/2.28-GCC-12.3.0
module load BBMap/39.01-GCC-11.3.0

# stats.sh may not support --version on all installs; guard it
if command -v stats.sh >/dev/null 2>&1; then
  log_module_version "BBMap.stats.sh" stats.sh --version | tee -a "$LOG_DIR/module_versions.txt"
fi

export JAVA_TOOL_OPTIONS="-Xmx120g"

# Run BBMap stats on both original and filtered assemblies
for dir in "$FLYE_DIR"/flye_output_*/; do
    [[ -d "$dir" ]] || continue
    sample_dirname=$(basename "$dir")
    sample_short="${sample_dirname#flye_output_}"
    input_orig="${dir}assembly.fasta"
    input_filt="${dir}/${sample_short}_assembly.m1000.fasta"

    if [[ -f "$input_orig" ]]; then
        stats.sh in="$input_orig" > "${dir}${sample_short}_stats.txt" 2> "$LOG_DIR/${sample_short}.bbmap_stats.log" || log_step "⚠️ stats.sh failed for ${sample_short} original"
    fi
    if [[ -f "$input_filt" ]]; then
        stats.sh in="$input_filt" > "${dir}${sample_short}_m1000_stats.txt" 2> "$LOG_DIR/${sample_short}.bbmap_stats_m1000.log" || log_step "⚠️ stats.sh failed for ${sample_short} filtered"
    fi
done

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

log_step "✅ Step 5 complete."

# ================================================
# Step 6: Minimap2 mapping of reads to assembly
# ================================================
log_step "🧮 Step 6: Minimap2 mapping"

module purge
module load minimap2/2.28-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0

log_module_version "minimap2" minimap2 --version | tee -a "$LOG_DIR/module_versions.txt"
log_module_version "SAMtools" samtools --version | tee -a "$LOG_DIR/module_versions.txt"

export JAVA_TOOL_OPTIONS="-Xmx24g"

mkdir -p "${FLYE_DIR}/mapping_bams"

# Map reads to their corresponding assembly.
for file in "$NANOFILT_DIR"/*_porechop_nanofilt.fastq.gz "$NANOFILT_DIR"/*_nanofilt.fastq.gz; do
    [[ -e "$file" ]] || continue
    sample=$(basename "$file" | sed -E 's/_porechop_nanofilt\.fastq\.gz$//; s/_nanofilt\.fastq\.gz$//')
    assembly="$FLYE_DIR/flye_output_${sample}/assembly.fasta"

    if [[ -f "$assembly" ]]; then
        log_step "Mapping reads for sample: $sample"
        bam_out="${FLYE_DIR}/mapping_bams/${sample}_aligned.bam"
        # create index
        minimap2 -d "${assembly}.mmi" "$assembly" >/dev/null 2>&1 || true

        minimap2 -ax map-ont -t "$THREADS" "${assembly}.mmi" "$file" | \
            samtools view -bS - | \
            samtools sort -@ "$THREADS" -o "$bam_out" - || { log_step "⚠️ Error mapping/sorting for $sample"; continue; }

        samtools index "$bam_out" || { log_step "⚠️ Error indexing BAM for $sample"; continue; }

        log_step "BAM for $sample written: $bam_out"
    else
        log_step "Assembly not found for $sample ($assembly), skipping."
    fi
done

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

log_step "✅ Step 6.1: Mapping complete."

# ============================================
# Step 7: jgi_summarize_bam_contig_depths
# ============================================
log_step "\033Step 6.2: Generating depth file\033"

# ============================================
# Step 7: jgi_summarize_bam_contig_depths
# ============================================
# Set up Packages needed
module purge
module load SAMtools/1.21-GCC-12.3.0
module load MetaBAT/2.17-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

# Generate depth files and MaxBin2 abundance files
for bam in *.sorted.bam; do
    # Get the sample name (remove .sorted.bam)
    sample=$(basename "$bam" .sorted.bam)

    echo "Processing sample: $sample"

    # Index the BAM if not already present
    if [[ ! -f "${bam}.bai" ]]; then
        echo "Indexing BAM for $sample..."
        samtools index "$bam"
    fi

    # Run jgi_summarize_bam_contig_depths
    echo "Generating depth file for $sample..."
    jgi_summarize_bam_contig_depths --outputDepth "${sample}.depth.txt" "$bam"

    if [ $? -eq 0 ]; then
        echo "Successfully generated depth for $sample"
        # Generate MaxBin2-compatible abundance file
        awk 'NR>1 {sum[$1]+=$4; n[$1]++} END{for (c in sum) print c, sum[c]/n[c]}' "${sample}.depth.txt" > "${sample}.maxbin_abund.txt"
        echo "MaxBin2 abundance file generated: ${sample}.maxbin_abund.txt"
    else
        echo "⚠️ERROR generating depth for $sample"
    fi
done

unset JAVA_TOOL_OPTIONS

# ================================
# Step 8: Kraken2 classification
# ================================
log_step "🧫 Step 7: Kraken2 classification"

module purge
module load Kraken2/2.1.3-GCC-11.3.0

log_module_version "Kraken2" kraken2 --version | tee -a "$LOG_DIR/module_versions.txt"

export JAVA_TOOL_OPTIONS="-Xmx100g"

mkdir -p "${FLYE_DIR}/kraken_output"

# NOTE: set KRAKEN2_DB environment variable or edit --db path below if required.
for file in "$FLYE_DIR"/flye_output_*/assembly.fasta; do
    [[ -f "$file" ]] || continue
    sample_dir=$(dirname "$file")
    sample=$(basename "$sample_dir" | sed 's/^flye_output_//')
    out_prefix="${FLYE_DIR}/kraken_output/${sample}"
    log_step "Running Kraken2 for sample $sample"
    kraken2 --use-names --output "${out_prefix}_kraken2.out" --report "${out_prefix}_kraken2.report" "$file" 2> "$LOG_DIR/${sample}.kraken2.log" || log_step "⚠️ kraken2 failed for $sample"
done

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

log_step "✅ Step 7: Kraken2 complete."

# =========================================
# Step 9: ABRicate on Contigs
# =========================================
log_step "🧪 Step 8: Running ABRicate on contigs"

module purge
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

log_module_version "ABRicate" abricate --version | tee -a "$LOG_DIR/module_versions.txt"

export JAVA_TOOL_OPTIONS="-Xmx100g"

# prepare database (runs once). If you already ran setupdb, you can skip this command.
abricate --setupdb > "$LOG_DIR/abricate_setup.log" 2>&1 || log_step "⚠️ abricate --setupdb may have failed or been run earlier"

mkdir -p "${FLYE_DIR}/abricate_results_contigs"

for file in "$FLYE_DIR"/flye_output_*/assembly.fasta; do
    [[ -f "$file" ]] || continue
    sample_dir=$(dirname "$file")
    sample=$(basename "$sample_dir" | sed 's/^flye_output_//')
    out_tsv="${FLYE_DIR}/abricate_results_contigs/${sample}_abricate.tsv"
    log_step "Running ABRicate for $sample"
    abricate "$file" > "$out_tsv" 2> "$LOG_DIR/${sample}.abricate.log" || log_step "⚠️ abricate failed for $sample"
done

unset JAVA_TOOL_OPTIONS
unset _JAVA_OPTIONS || true

log_step "✅ Step 8: ABRicate complete."

# =========================================
# Step 10: Link ABRicate AMR genes to Kraken2 species per contig
# =========================================
mkdir -p "${FLYE_DIR}/combined_results_2"
mkdir -p "${FLYE_DIR}/kraken_output"
mkdir -p "${FLYE_DIR}/abricate_results_contigs"
LOG="${FLYE_DIR}/combined_results_2/join_abr_vs_kraken.log"

log_step "🧬 Step 10 (fixed): Linking ABRicate AMR genes to Kraken2 species per contig"

for d in "$FLYE_DIR"/flye_output_*; do
    [[ -d "$d" ]] || continue
    sample="${d##*/}"
    sample="${sample#flye_output_}"
    log_step "Processing sample: $sample" | tee -a "$LOG"

    kraken_out="${FLYE_DIR}/kraken_output/${sample}_kraken2.out"
    abricate_raw="${FLYE_DIR}/abricate_results_contigs/${sample}_abricate.tsv"
    kraken_tax="${FLYE_DIR}/kraken_output/${sample}_contig_taxonomy.tsv"
    abricate_clean="${FLYE_DIR}/abricate_results_contigs/${sample}_abricate_clean.tsv"
    output="${FLYE_DIR}/combined_results_2/${sample}_AMR_with_species.tsv"

    if [[ -f "$kraken_out" ]]; then
        echo -e "contig\tspecies" > "$kraken_tax"
        awk -F'\t' '
            function norm(s) {
                gsub(/^ +| +$/, "", s);
                gsub(/^>/, "", s);
                sub(/ .*/, "", s);
                return s
            }
            $1 == "C" {
                contig = norm($2)
                species = $3
                # remove trailing " (taxid NNNNN)"
                sub(/[[:space:]]*\(taxid[[:space:]]*[0-9]+\)[[:space:]]*$/, "", species)
                gsub(/^ +| +$/, "", species)
                if (species == "") species = "unclassified"
                print contig "\t" species
            }
        ' "$kraken_out" >> "$kraken_tax"
    else
        log_step "⚠️ No Kraken2 .out for $sample at $kraken_out — skipping sample" | tee -a "$LOG"
        continue
    fi

    if [[ -f "$abricate_raw" ]]; then
        # detect header if present and produce normalized 5-column file: contig,gene,cov,id,db
        awk '
            BEGIN{
                FS=OFS="\t"
                seq_col=gene_col=cov_col=id_col=db_col=0
            }
            NR==1 && /^#/ {
                # header line
                for(i=1;i<=NF;i++){
                    h=$i
                    gsub(/^#/,"",h)
                    H=toupper(h)
                    if(H=="SEQUENCE" || H=="SEQID" || H=="SEQUENCEID" || H=="SEQUENCE_ID" || H=="CONTIG") seq_col=i
                    if(H=="GENE" || H=="GENE_NAME" || H=="GENEID") gene_col=i
                    if(H=="COVERAGE" || H=="COV" || H=="%COVERAGE" || H=="COVERAGE%") cov_col=i
                    if(H=="%IDENTITY" || H=="IDENTITY" || H=="PIDENT" || H=="%ID") id_col=i
                    if(H=="DATABASE" || H=="DB") db_col=i
                }
                next
            }
            !/^#/ {
                if(seq_col==0) seq_col=1
                if(gene_col==0) gene_col=2
                if(cov_col==0) cov_col=4
                if(id_col==0) id_col=5
                if(db_col==0) db_col=NF
                seq = $seq_col
                gene = $gene_col
                cov = $cov_col
                id = $id_col
                db = $db_col
                gsub(/^>/,"",seq)
                sub(/ .*/,"",seq)
                n = split(seq, parts, /[\/:|]/)
                seq = parts[n]
                gsub(/^ +| +$/,"",seq)
                print seq, gene, cov, id, db
            }
        ' "$abricate_raw" > "$abricate_clean"
    else
        log_step "⚠️ No ABRicate output for $sample at $abricate_raw — skipping sample" | tee -a "$LOG"
        continue
    fi

    echo -e "contig\tspecies\tgene\tcov%\tidentity%\tdatabase" > "$output"

    matched_file="${FLYE_DIR}/combined_results_2/${sample}_matched.txt"
    unmatched_file="${FLYE_DIR}/combined_results_2/${sample}_unmatched.txt"
    : > "$matched_file"
    : > "$unmatched_file"

    awk -v out="$output" -v unmatched="$unmatched_file" -v logf="$LOG" '
        BEGIN{FS=OFS="\t"}
        NR==FNR {
            k = $1
            tax[k] = $2
            next
        }
        {
            cont = $1
            contn = cont
            gsub(/^>/,"",contn)
            sub(/ .*/,"",contn)
            if (contn in tax) {
                sp = tax[contn]
                print cont, sp, $2, $3, $4, $5 >> out
                print cont >> "'$matched_file'"
            } else {
                # try relaxed match: look for tax key equal to contn suffix/prefix
                found = 0
                for (k in tax) {
                    if (k == contn || index(contn, k) == 1 || index(k, contn) == 1) {
                        sp = tax[k]
                        print cont, sp, $2, $3, $4, $5 >> out
                        print cont >> "'$matched_file'"
                        found = 1
                        break
                    }
                }
                if (!found) {
                    print cont > unmatched
                    print cont, "unclassified", $2, $3, $4, $5 >> out
                }
            }
        }
    ' "$kraken_tax" "$abricate_clean"

    total=$(awk 'END{print NR}' "$output" 2>/dev/null || echo 0)
    unmatched_count=$(wc -l < "$unmatched_file" 2>/dev/null || echo 0)
    log_step "Sample $sample: lines in final output: $total; unmatched contigs: $unmatched_count" | tee -a "$LOG"

    if [[ -s "$unmatched_file" ]]; then
        log_step "First 20 unmatched contigs (sample $sample):" | tee -a "$LOG"
        head -n 20 "$unmatched_file" | sed 's/^/  /' | tee -a "$LOG"
    fi
done

log_step "✅ Step 10 (fixed) complete. See $LOG for details."