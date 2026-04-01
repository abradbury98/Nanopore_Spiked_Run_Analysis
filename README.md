# 🧬 ONT Metagenomics Pipeline

> Oxford Nanopore long-read metagenomic analysis pipeline for NeSI HPC — from raw basecalled reads to taxonomically-annotated AMR gene tables.

---

## Overview

This pipeline processes Oxford Nanopore Technology (ONT) sequencing data through a series of SLURM batch jobs on the [NeSI](https://www.nesi.org.nz/) HPC platform. It covers the full metagenomic workflow: raw read QC → adapter trimming → quality filtering → host read removal → assembly → contig analysis → taxonomic classification → AMR gene detection.

```
Raw ONT reads
    │
    ▼
[Script 1] Raw QC ──────────────────── NanoPlot · nanoQC · MultiQC
    │
    ▼
[Script 2] Filtering & Host Removal ── Porechop · NanoFilt · minimap2 · seqkit
    │
    ▼
[Script 3] Assembly & Analysis ──────── Flye · Kraken2 · ABRicate · MetaBAT2
```

---

## Pipeline Scripts

| # | File | Job Name | Time | Mem | CPUs | Description |
|---|------|----------|------|-----|------|-------------|
| 1a | `job1-checked-05-12-25.sl` | `1.QC_with_nanoQC` | 8h | 64 GB | 6 | Raw QC — NanoPlot + nanoQC + MultiQC (reads pre-merged) |
| 1b | `job1_2.sl` | `1.QC_Raw` | 8h | 64 GB | 6 | Merge barcodes then run NanoPlot + MultiQC with HTML dashboard |
| 2 | `job2.sl` | `2.Filtering-hostremove` | 24h | 200 GB | 8 | Adapter trim, quality filter, host removal, post-filter QC |
| 3 | `job3.sl` | `3-assembly-contig_analysis` | 8h | 64 GB | 6 | Assembly, contig stats, mapping, taxonomy, AMR, integration |

> **Script 1a vs 1b:** These are alternative approaches to the same QC step — run one, not both. Script 1b includes barcode merging and generates an HTML navigation dashboard; Script 1a additionally runs nanoQC but assumes reads are pre-merged.

---

## Prerequisites

### NeSI Account

All scripts use `module load` via the NeSI Lmod system. No manual software installation is needed. Set your account in each script:

```bash
#SBATCH --account=<your_account>
```

### Required Input Files

**For all scripts:**
- Raw ONT reads in `barcode01/`, `barcode02/`, ... directories (Script 1b will merge these), **or** pre-merged FASTQ files in `raw_merged_reads_all/`

**For Script 2 (host removal):**

Download and place these three reference files in your working directory before submitting Script 2. The script will auto-combine them into `combined_host_refs.fna` on first run.

#### 🐔 Chicken genome — *Gallus gallus* GRCg7b

| Field | Value |
|-------|-------|
| Assembly | bGalGal1.mat.broiler.GRCg7b |
| Accession | [GCF_016699485.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/) |
| Submitted by | Vertebrate Genomes Project, January 2021 |
| Description | Chromosome-level broiler chicken assembly; 41 chromosomes, contig N50 ~18.8 Mb |

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
```

#### 🧬 Human genome — *Homo sapiens* GRCh38.p14

| Field | Value |
|-------|-------|
| Assembly | GRCh38.p14 |
| Accession | [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) |
| Submitted by | Genome Reference Consortium, patch release 14 (February 2022) |
| Description | The standard human reference genome (hg38); 69 patch scaffolds added in p14 |

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

> ⚠️ The human genome is ~3 GB uncompressed. Download time and disk space should be planned accordingly.

#### 🧪 ONT DNA Control Strand (DCS)

| Field | Value |
|-------|-------|
| File | `DCS.fasta` |
| Source | [Oxford Nanopore Technologies Community](https://nanoporetech.com/support/library-prep/kit-contents-and-composition/what-is-dna-cs-dcs) |
| Description | A 3.6 kb amplicon mapping to the 3′ end of the Lambda phage genome, supplied as a positive control in ONT ligation sequencing kits (e.g. SQK-LSK114). Reads derived from DCS are removed to prevent them being assembled as genuine metagenomic signal. |

The DCS FASTA sequence is available from the ONT Community portal (login required). Rename the downloaded file to `DCS.fasta` and place it in your working directory.

Alternatively, the sequence can be retrieved from the ONT assets CDN (no login required):

```bash
wget -O DCS.fasta "https://assets.ctfassets.net/hkzaxo8a05x5/2IX56YmF5ug0kAQYoAg2Uk/159523e3a26d31f5b1f21d5cc94c5b61/DNA_CS.fasta"
```

**For Script 3 (Kraken2):**

Set the Kraken2 database path before submitting:
```bash
export KRAKEN2_DB=/path/to/kraken2_db
sbatch job3.sl
```

### Directory Setup

Create the SLURM log directory before submitting any job:

```bash
mkdir -p log
```

---

## Usage

### Recommended Submission Order

Run scripts sequentially, waiting for each job to complete before submitting the next.

```bash
# 1. Create the SLURM log directory
mkdir -p log

# 2. Raw QC (choose one)
sbatch job1-checked-05-12-25.sl   # includes nanoQC
# OR
sbatch job1_2.sl                  # includes barcode merging + HTML dashboard

# 3. Monitor job completion
squeue --me

# 4. Filtering and host removal
sbatch job2.sl

# 5. Assembly and downstream analysis
export KRAKEN2_DB=/path/to/kraken2_db
sbatch job3.sl
```

### Running Individual Steps (Script 2 only)

Script 2 supports running specific steps via `--step`, useful for re-running after partial failures:

```bash
sbatch job2.sl --step porechop      # adapter trimming only
sbatch job2.sl --step nanofilt      # quality filtering only
sbatch job2.sl --step hostremoval   # host removal only
sbatch job2.sl --step nanoplot      # QC plots only
sbatch job2.sl --step nanoqc        # nanoQC only
sbatch job2.sl --step multiqc       # MultiQC only
sbatch job2.sl --step all           # full pipeline (default)
```

### Monitoring Jobs

```bash
squeue --me                  # view running/queued jobs
sacct -j <JOBID>             # view completed job stats
tail -f log/<jobname>_<jobid>.out   # follow job output live
```

---

## Script Details

### Script 1a — Raw QC with nanoQC (`job1-checked-05-12-25.sl`)

Runs a three-tool QC suite on pre-merged raw reads.

**Steps:**
1. **NanoPlot** — per-sample read length and quality plots (hex, dot)
2. **nanoQC** — read-end quality assessment
3. **MultiQC** — per-sample and global aggregated reports

**Key settings:**
```bash
INPUT_DIR="raw_merged_reads_all"
QC_DIR="Quality_control_raw"
MAX_JOBS=1   # concurrent NanoPlot/nanoQC jobs
```

**NeSI modules:**

| Module | Version |
|--------|---------|
| NanoPlot | `1.43.0-foss-2023a-Python-3.11.6` |
| nanoQC | `0.9.4-gimkl-2022a-Python-3.10.5` |
| MultiQC | `1.24.1-foss-2023a-Python-3.11.6` |

**Outputs:**
```
Quality_control_raw/
├── Nanoplot/<sample>/NanoPlot-report.html
├── NanoQC/<sample>/
├── MultiQC/<sample>/multiqc_<sample>.html
├── MultiQC/multiqc_all_barcodes.html
└── logs/
```

---

### Script 1b — Raw QC with Barcode Merging (`job1_2.sl`)

Merges per-barcode FASTQ files, then runs NanoPlot and MultiQC. Generates an HTML navigation dashboard linking all per-sample reports.

**Steps:**
1. **Merge** — concatenates `barcode*/` FASTQ files; uses `pigz` for parallel compression when available
2. **NanoPlot** — per-sample quality plots
3. **MultiQC** — per-sample and global reports with HTML index dashboard

**Key settings:**
```bash
INPUT_DIR="raw_merged_reads_all"
QC_DIR="Quality_control_raw_4"
MAX_JOBS=1
```

**NeSI modules:**

| Module | Version |
|--------|---------|
| NanoPlot | `1.43.0-foss-2023a-Python-3.11.6` |
| MultiQC | `1.15-gimkl-2022a-Python-3.10.5` |

**Outputs:**
```
raw_merged_reads_all/all_barcode<NN>.fastq.gz
Quality_control_raw_4/
├── nanoplot/<sample>/
├── multiqc/index.html          ← HTML navigation dashboard
├── multiqc/<sample>/multiqc_<sample>.html
└── multiqc/multiqc_all_barcodes.html
```

---

### Script 2 — Filtering & Host Removal (`job2.sl`)

The most resource-intensive step. Trims adapters, filters by quality and length, and removes host-derived reads by alignment to a combined reference.

**Steps:**
1. **Porechop** — adapter trimming
2. **NanoFilt** — quality (≥ Q10) and length (≥ 250 bp) filtering
3. **Host removal** — `minimap2 -ax map-ont` → `samtools` (extract mapped IDs) → `seqkit grep -v` (remove)
4. **NanoPlot + nanoQC + MultiQC** — QC on cleaned reads

**Key settings:**
```bash
# NanoFilt thresholds
NanoFilt -q 10 -l 250

# Concurrent jobs
MAX_JOBS=4
```

**NeSI modules:**

| Module | Version |
|--------|---------|
| Porechop | `0.2.4-gimkl-2022a-Python-3.11.3` |
| nanofilt | `2.6.0-gimkl-2020a-Python-3.8.2` |
| minimap2 | `2.28-GCC-12.3.0` |
| SAMtools | `1.22-GCC-12.3.0` |
| SeqKit | `2.4.0` |
| NanoPlot | `1.43.0-foss-2023a-Python-3.11.6` |
| nanoQC | `0.9.4-gimkl-2022a-Python-3.10.5` |
| MultiQC | `1.15-gimkl-2022a-Python-3.10.5` |

**Outputs:**
```
porechop_2/<sample>_porechop.fastq.gz
nanofilt_2/<sample>_porechop_nanofilt.fastq.gz
host_removed_reads_2/<sample>_host_removed.fastq.gz
logs/<sample>.hostremoval.log      ← read counts before/after
Quality_control_cleaned/
├── Nanoplot/
├── NanoQC/
└── MultiQC/
```

---

### Script 3 — Assembly & Contig Analysis (`job3.sl`)

Full downstream pipeline from assembly through to a taxonomy-annotated AMR gene table.

**Steps:**

| Step | Tool | Description |
|------|------|-------------|
| 4 | Flye | Metagenomic assembly (`--meta --nano-raw`) |
| 5.1 | Seqmagick | Filter contigs < 1,000 bp |
| 5.2 | BBMap (`stats.sh`) | Assembly statistics for raw and filtered contigs |
| 6 | minimap2 + SAMtools | Map reads back to assembly; sorted + indexed BAM |
| 6.2 | MetaBAT2 (`jgi_summarize_bam_contig_depths`) | Per-contig depth and MaxBin2 abundance files |
| 7 | Kraken2 | Taxonomic classification of contigs |
| 8 | ABRicate | AMR gene screening across standard databases |
| 10 | awk | Join ABRicate + Kraken2 per contig → annotated AMR table |

**Key settings:**
```bash
NANOFILT_DIR="host_removed_reads_2"   # input from Script 2
FLYE_DIR="flye_assembly"
# Kraken2 DB: set via environment variable
export KRAKEN2_DB=/path/to/db
```

**NeSI modules:**

| Module | Version |
|--------|---------|
| Flye | `2.9.5-foss-2023a-Python-3.11.6` |
| seqmagick | `0.8.4-gimkl-2020a-Python-3.8.2` |
| BBMap | `39.01-GCC-11.3.0` |
| minimap2 | `2.28-GCC-12.3.0` |
| SAMtools | `1.21-GCC-12.3.0` |
| MetaBAT | `2.17-GCC-12.3.0` |
| Kraken2 | `2.1.3-GCC-11.3.0` |
| ABRicate | `1.0.0-GCC-11.3.0-Perl-5.34.1` |

**Outputs:**
```
flye_assembly/
├── flye_output_<sample>/
│   ├── assembly.fasta
│   ├── <sample>_assembly.m1000.fasta    ← contigs ≥ 1000 bp
│   ├── <sample>_stats.txt               ← BBMap stats (all contigs)
│   └── <sample>_m1000_stats.txt         ← BBMap stats (filtered)
├── mapping_bams/<sample>_aligned.bam
├── kraken_output/
│   ├── <sample>_kraken2.out
│   └── <sample>_kraken2.report
├── abricate_results_contigs/
│   └── <sample>_abricate.tsv
└── combined_results_2/
    └── <sample>_AMR_with_species.tsv    ← final integrated output
```

---

## Directory Structure

Full expected layout after running the complete pipeline:

```
working_directory/
├── barcode01/                             # raw input (if using Script 1b)
├── barcode02/
├── ...
│
├── raw_merged_reads_all/                  # merged per-barcode FASTQs
│
├── Quality_control_raw/                   # Script 1a outputs
├── Quality_control_raw_4/                 # Script 1b outputs
│
├── porechop_2/                            # adapter-trimmed reads
├── nanofilt_2/                            # quality/length filtered reads
├── host_removed_reads_2/                  # host-depleted reads
├── Quality_control_cleaned/               # post-filter QC
│
├── flye_assembly/                         # all Script 3 outputs
│   ├── flye_output_<sample>/
│   ├── mapping_bams/
│   ├── kraken_output/
│   ├── abricate_results_contigs/
│   ├── combined_results_2/
│   └── logs/
│
├── combined_host_refs.fna                 # auto-generated host reference
├── GCF_016699485.2_*.fna                 # chicken reference (user-provided)
├── DCS.fasta                             # control reference (user-provided)
├── GCF_000001405.40_*.fna               # human reference (user-provided)
│
└── log/                                   # SLURM stdout/stderr logs
```

---

## Notes

- **`set -euo pipefail`** is active in all scripts — any unhandled error will abort the job. Check `log/<jobname>_<jobid>.err` if a job exits unexpectedly.
- **Skip logic** is implemented throughout: if expected output files already exist for a sample, that sample is skipped. This makes it safe to re-submit after a partial failure.
- **Script 3 filename dependency:** Step 4 expects reads named `*_porechop_host_removed.fastq.gz`. Verify the filenames produced by Script 2 match before submitting Script 3.
- **Step 6.2 working directory:** The `jgi_summarize_bam_contig_depths` loop in Script 3 processes `*.sorted.bam` files in the current working directory — ensure you submit from the correct location.
- **Memory / time limits:** The values set are conservative estimates. Adjust `#SBATCH --mem` and `#SBATCH --time` based on dataset size and observed usage from `sacct`.
- **Module versions are pinned** — if a module is not available on NeSI, check `module spider <tool>` for the closest available version and update the script accordingly.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Job exits immediately | Missing `log/` directory | `mkdir -p log` before submitting |
| Host removal finds no reads | Combined reference not built | Delete `combined_host_refs.fna` to force rebuild |
| Kraken2 fails | `KRAKEN2_DB` not set | `export KRAKEN2_DB=/path/to/db` before `sbatch` |
| Script 3 skips all samples | Filename mismatch from Script 2 | Check filenames in `host_removed_reads_2/` match `*_porechop_host_removed.fastq.gz` |
| nanoQC output missing | nanoQC writes to current dir | Script runs nanoQC in a subshell from the output dir — check per-sample log |

---

## SLURM Account

All scripts are configured for account `massey04083`. Update `#SBATCH --account` in each script if running under a different NeSI project.

---

## References

### Reference Genomes

**Chicken (*Gallus gallus*) — GRCg7b**
> Vertebrate Genomes Project. *Gallus gallus* genome assembly bGalGal1.mat.broiler.GRCg7b. NCBI Assembly accession [GCF_016699485.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016699485.2/). Submitted January 2021.

**Human (*Homo sapiens*) — GRCh38.p14**
> Genome Reference Consortium. *Homo sapiens* genome assembly GRCh38.p14. NCBI Assembly accession [GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). Patch release 14, February 2022.

**ONT DNA Control Strand (DCS)**
> Oxford Nanopore Technologies. *DNA Control Strand (DCS)*. A 3.6 kb standard amplicon mapping to the 3′ end of the Lambda phage (*Enterobacteria phage lambda*) genome. Supplied in ONT ligation sequencing kits as a positive library preparation control. Sequence available at: https://nanoporetech.com/support/library-prep/kit-contents-and-composition/what-is-dna-cs-dcs

### Key Software

| Tool | Citation |
|------|----------|
| NanoPlot / NanoFilt / nanoQC | De Coster W, et al. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15):2666–2669. [doi:10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149) |
| MultiQC | Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19):3047–3048. [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |
| Porechop | Wick RR. Porechop. GitHub: https://github.com/rrwick/Porechop |
| minimap2 | Li H. (2021). New strategies to improve minimap2 alignment accuracy. *Bioinformatics*, 37(23):4572–4574. [doi:10.1093/bioinformatics/btab705](https://doi.org/10.1093/bioinformatics/btab705) |
| SAMtools | Danecek P, et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2):giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| SeqKit | Shen W, et al. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLOS ONE*, 11(10):e0163962. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962) |
| Flye | Kolmogorov M, et al. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37:540–546. [doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8) |
| Kraken2 | Wood DE, et al. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20:257. [doi:10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0) |
| ABRicate | Seemann T. ABRicate. GitHub: https://github.com/tseemann/abricate |
| MetaBAT2 | Kang DD, et al. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*, 7:e7359. [doi:10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359) |
| BBMap | Bushnell B. BBMap. https://sourceforge.net/projects/bbmap/ |
