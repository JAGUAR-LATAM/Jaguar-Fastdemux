#!/bin/bash
# =============================================================================
# fastdemux_general_pipeline.sh  —  Genotype-based demultiplexing pipeline for scATAC-seq
# =============================================================================
# Project      : JAGUAR multimodal atlas
# Affiliation  : LIIGH-UNAM / Wellcome Sanger Institute
#
# Description
# -----------
# Runs the popscle demuxlet pipeline on CellRanger ATAC data (no UMI).
# For each sample listed in CONFIG_TSV the script:
#   Step -1  Verify BAM exists and is indexed
#   Step  0  Extract chromosome order from BAM header (@SQ lines)
#   Step  1  Subset joint germline VCF to the donors of interest
#   Step  1.5 Filter to biallelic SNPs only (exclude indels)
#   Step  2  Restrict VCF to ATAC peak regions
#   Step  3  Reorder VCF ##contig header to match BAM chromosome order
#   Step  4  Run fastdemux
#
# All steps are idempotent: if the expected output already exists the step
# is skipped, making it safe to re-run after a failure.
#
# Input config (CONFIG_TSV)
# -------------------------
# Tab-separated file with a header row and the following columns:
#   SAMPLE_ID       Unique identifier for the sample
#   BAM_FILE        Absolute path to the position-sorted BAM
#   BARCODE_FILE    Absolute path to the valid-barcode list (one per line)
#   DONOR_LIST_FILE Absolute path to the donor-ID list (one per line, must
#                   match sample names in JOINT_DONOR_VCF)
#   NUM_DONORS      Number of donors in the pool (informational; not yet used
#                   by popscle but useful for downstream vireo)
#
# Modules loaded (versions as of 2026-03-04)
# -------------------------------------------
#   cellgen/bcftools/1.23     — bcftools 1.23 + samtools 1.23 + bgzip + tabix
#   HGI/softpack/groups/jaguar_analysis/jaguar-fastdemux/2 — fastdemux
#
# Usage
# -----
#   # Submit to LSF:
#   bsub < fastdemux_general_pipeline.sh
#
#   # Dry-run locally (no LSF):
#   bash fastdemux_general_pipeline.sh
#
# Notes
# -----
# - Genome reference: hg38
# - Cell barcode tag in BAM: CB (CellRanger ATAC default)
# - bcftools sort --threads is supported since bcftools 1.12
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# LSF job parameters
# ─────────────────────────────────────────────────────────────────────────────
#BSUB -G jaguar_analysis
#BSUB -J fastdemux_scATAC
#BSUB -q normal
#BSUB -n 16
#BSUB -R "select[mem>16000] rusage[mem=64000] span[hosts=1]"
#BSUB -M 16000
#BSUB -q long
#BSUB -o logs/fastdemux_%J.out
#BSUB -e logs/fastdemux_%J.err

# ─────────────────────────────────────────────────────────────────────────────
# Shell options
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
export TERM=xterm
export LSB_DEFAULT_USERGROUP=jaguar_analysis

# ─────────────────────────────────────────────────────────────────────────────
# Load software modules
# ─────────────────────────────────────────────────────────────────────────────
# cellgen/bcftools/1.23 provides: bcftools 1.23, samtools 1.23, bgzip, tabix
# (replaces the former HGI/softpack/users/bh18/bcftools/1.0  and
#  HGI/softpack/groups/delscreen/Align/1 which were used for these tools)
echo "--- Loading software modules ---"
module load HGI/softpack/groups/jaguar_analysis/jaguar-fastdemux/2

# Print tool versions for reproducibility
echo "bcftools : $(bcftools  --version | head -1)"
echo "samtools : $(samtools  --version | head -1)"
echo "bgzip    : $(bgzip     --version | head -1)"
echo "Modules loaded."

# ─────────────────────────────────────────────────────────────────────────────
# Thread count — respect LSF job allocation
# ─────────────────────────────────────────────────────────────────────────────
THREADS="${LSB_DJOB_NUMPROC:-16}"
echo "--- Using THREADS=$THREADS ---"

# ─────────────────────────────────────────────────────────────────────────────
# Global configuration — edit these paths to match your run
# ─────────────────────────────────────────────────────────────────────────────
CONFIG_TSV="config_pipeline.tsv"

# Directory that will hold per-sample output sub-directories
GLOBAL_ANALYSIS_DIR="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jv8/fastdemux-scATAC-07-14--17-20"

# Joint germline VCF (multi-sample, bgzipped + tabix-indexed)
JOINT_DONOR_VCF="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/data_working/bge/1_variant_calling/BGE05/snp_calling_batches_1-5/joint_variant_calling/joint_germline_recalibrated.vcf.gz"

# Directory containing per-sample peaks BED files
# Expected layout: ${BED_FOLDER}/${SAMPLE_ID}/peaks.bed.gz
BED_FOLDER="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jv8/beds"

# BAM tag that stores the cell barcode (CellRanger ATAC default: CB)
CELL_TAG="CB"

# ─────────────────────────────────────────────────────────────────────────────
# Ensure log directory exists (required by #BSUB -o/-e directives above)
# ─────────────────────────────────────────────────────────────────────────────
mkdir -p logs

# ─────────────────────────────────────────────────────────────────────────────
# Helper: print a timestamped message
# ─────────────────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ─────────────────────────────────────────────────────────────────────────────
# Main loop — one iteration per sample in the config TSV
# ─────────────────────────────────────────────────────────────────────────────
log "Starting pipeline from config: $CONFIG_TSV"

tail -n +2 "$CONFIG_TSV" | while IFS=$'\t' read -r SAMPLE_ID BAM_FILE BARCODE_FILE DONOR_LIST_FILE NUM_DONORS
do
  # Skip incomplete rows
  if [ -z "${SAMPLE_ID:-}" ] || [ -z "${BAM_FILE:-}" ] || [ -z "${DONOR_LIST_FILE:-}" ]; then
    log "WARNING: incomplete row — skipping."
    continue
  fi

  log "============================================================"
  log "Processing sample : $SAMPLE_ID"
  log "  BAM              : $BAM_FILE"
  log "  Donors           : $DONOR_LIST_FILE"
  log "  Barcodes         : ${BARCODE_FILE:-(none)}"
  log "  Num donors       : ${NUM_DONORS:-(unset)}"
  log "============================================================"

  BASE_DIR="${GLOBAL_ANALYSIS_DIR}/${SAMPLE_ID}"
  FASTDEMUX_DIR="${BASE_DIR}/fastdemux_out"
  DONOR_VCF_DIR="${BASE_DIR}/donor_vcf_prep"
  mkdir -p "$FASTDEMUX_DIR" "$DONOR_VCF_DIR"

  # -------------------------------------------------------------------------
  # STEP -1: Verify BAM exists and is indexed
  # -------------------------------------------------------------------------
                                                                                                                                                                       146,12        37%
  if [ ! -f "$BAM_FILE" ]; then
    log "ERROR: BAM not found: $BAM_FILE" >&2
    exit 1
  fi
  if [ ! -f "${BAM_FILE}.bai" ]; then
    log "BAM index missing — creating it..."
    samtools index -@ "$THREADS" "$BAM_FILE"
  fi

  # -------------------------------------------------------------------------
  # STEP 0: Extract chromosome order from BAM header (@SQ lines)
  # Used later to reorder VCF ##contig entries to match the BAM.
  # -------------------------------------------------------------------------
  BAM_ORDER_FILE="${DONOR_VCF_DIR}/bam.chrom.order.txt"
  if [ ! -f "$BAM_ORDER_FILE" ]; then
    log "Extracting chromosome order from BAM header..."
    samtools view -H "$BAM_FILE" \
      | grep '^@SQ' \
      | cut -f2 \
      | sed 's/SN://' \
      > "$BAM_ORDER_FILE"
  else
    log "Chromosome order file already exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # STEP 1: Subset joint germline VCF to the donors of interest
  # -------------------------------------------------------------------------
  SUBSET_VCF="${DONOR_VCF_DIR}/donors.subset.vcf.gz"
  if [ ! -f "$SUBSET_VCF" ]; then
    log "Subsetting VCF to donors in $DONOR_LIST_FILE ..."
    bcftools view \
      --threads "$THREADS" \
      -S "$DONOR_LIST_FILE" \
      "$JOINT_DONOR_VCF" \
      -Oz -o "${DONOR_VCF_DIR}/tmp.subset.vcf.gz"
    bcftools sort \
      "${DONOR_VCF_DIR}/tmp.subset.vcf.gz" \
      -Oz -o "$SUBSET_VCF"
    tabix -p vcf "$SUBSET_VCF"
    rm -f "${DONOR_VCF_DIR}/tmp.subset.vcf.gz"
  else
    log "Donor-subset VCF already exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # STEP 1.5: Filter to biallelic SNPs only (exclude indels and multi-allelic)
  # -v snps  : keep SNPs only
  # -m2 -M2  : exactly 2 alleles (strictly biallelic)
  # -------------------------------------------------------------------------
  BIALLELIC_SNPS_VCF="${DONOR_VCF_DIR}/donors.subset.biallelic_snps.vcf.gz"
  if [ ! -f "$BIALLELIC_SNPS_VCF" ]; then
    log "Filtering to biallelic SNPs (removing indels and multi-allelic sites)..."
    bcftools view \
      --threads "$THREADS" \
      -v snps -m2 -M2 \
      "$SUBSET_VCF" \
      -Oz -o "${DONOR_VCF_DIR}/tmp.biallelic_snps.vcf.gz"
    bcftools sort \
      "${DONOR_VCF_DIR}/tmp.biallelic_snps.vcf.gz" \
      -Oz -o "$BIALLELIC_SNPS_VCF"
    tabix -p vcf "$BIALLELIC_SNPS_VCF"
    rm -f "${DONOR_VCF_DIR}/tmp.biallelic_snps.vcf.gz"
  else
    log "Biallelic SNP VCF already exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # STEP 2: Restrict VCF to ATAC peak regions
  # Only SNPs overlapping accessible regions are informative for scATAC-seq
  # demultiplexing; restricting here greatly reduces the pileup runtime.
  # -------------------------------------------------------------------------
  PEAKS_BED="${BED_FOLDER}/${SAMPLE_ID}/peaks.bed.gz"
  if [ ! -f "$PEAKS_BED" ]; then
    log "ERROR: peaks BED not found: $PEAKS_BED" >&2
    exit 1
  fi
  # Index if missing (required by bcftools -R for random access)
  if [ ! -f "${PEAKS_BED}.tbi" ] && [ ! -f "${PEAKS_BED}.csi" ]; then
    log "Peaks BED has no index — attempting tabix indexing..."
    tabix -p bed "$PEAKS_BED" \
      || log "WARNING: could not index peaks.bed.gz; bcftools -R may be slower."
  fi

  FILTERED_VCF="${DONOR_VCF_DIR}/donors.atac_filtered.biallelic_snps.vcf.gz"
  if [ ! -f "$FILTERED_VCF" ]; then
    log "Restricting VCF to ATAC peak regions..."
    bcftools view \
      --threads "$THREADS" \
      -R "$PEAKS_BED" \
      "$BIALLELIC_SNPS_VCF" \
      -Oz -o "${DONOR_VCF_DIR}/tmp.atac.vcf.gz"
    bcftools sort \
      "${DONOR_VCF_DIR}/tmp.atac.vcf.gz" \
      -Oz -o "$FILTERED_VCF"
    tabix -p vcf "$FILTERED_VCF"
    rm -f "${DONOR_VCF_DIR}/tmp.atac.vcf.gz"
  else
    log "Peak-filtered VCF already exists — skipping."
  fi

  # -------------------------------------------------------------------------
  # STEP 3: Reorder VCF ##contig header entries to match the BAM chromosome
  # order. 
  # Note: this step reorders the ##contig header lines only; the variant
  # records themselves are already sorted by bcftools sort above.
  # -------------------------------------------------------------------------
  REORDERED_VCF="${DONOR_VCF_DIR}/donors.atac_filtered.biallelic_snps.reordered.vcf.gz"
  if [ ! -f "$REORDERED_VCF" ]; then
    log "Reordering VCF ##contig header entries to match BAM chromosome order..."

    TMP_HEADER="${DONOR_VCF_DIR}/vcf.header.meta.txt"
    TMP_CONTIGS_ALL="${DONOR_VCF_DIR}/vcf.contigs.all.txt"
    TMP_CONTIGS_ORDERED="${DONOR_VCF_DIR}/vcf.contigs.ordered.txt"
    TMP_CHROMLINE="${DONOR_VCF_DIR}/vcf.chromline.txt"
    TMP_BODY="${DONOR_VCF_DIR}/vcf.body.tmp"

    # Split the VCF header into its components
    bcftools view -h "$FILTERED_VCF" | grep '^##contig'                  > "$TMP_CONTIGS_ALL"
    bcftools view -h "$FILTERED_VCF" | grep '^##' | grep -v '^##contig'  > "$TMP_HEADER"
    bcftools view -h "$FILTERED_VCF" | grep '^#CHROM'                    > "$TMP_CHROMLINE"

    # Re-build ##contig section in BAM chromosome order
    : > "$TMP_CONTIGS_ORDERED"
    while read -r CHR; do
      if grep -m1 -E "ID=${CHR}([,>])" "$TMP_CONTIGS_ALL" >> "$TMP_CONTIGS_ORDERED"; then
        true
      else
        log "WARNING: chromosome '$CHR' not found in VCF header — omitting ##contig entry."
      fi
    done < "$BAM_ORDER_FILE"

    # Extract variant records (non-header lines)
    bgzip -cd "$FILTERED_VCF" | awk '!/^#/' > "$TMP_BODY"

    # Reassemble: meta headers + reordered contigs + #CHROM line + records
    cat "$TMP_HEADER" "$TMP_CONTIGS_ORDERED" "$TMP_CHROMLINE" "$TMP_BODY" \
      | bgzip -c > "$REORDERED_VCF"
    tabix -p vcf "$REORDERED_VCF"

    rm -f "$TMP_HEADER" "$TMP_CONTIGS_ALL" "$TMP_CONTIGS_ORDERED" "$TMP_CHROMLINE" "$TMP_BODY"
    log "Reordered VCF written to: $REORDERED_VCF"
  else
    log "Reordered VCF already exists — skipping."
  fi


  # -------------------------------------------------------------------------
  # STEP 4: Assign cells to donors (FASTDEMUX)
  # Uses genotype likelihoods (GT field) and the pileup counts to assign
  # each cell to a donor or flag it as a doublet.
  # -------------------------------------------------------------------------

  if [ ! -f "${FASTDEMUX_DIR}/fastdemux.fdout.raw.corr.txt.gz" ]; then
    log "Running fastdemux..."
    fastdemux -t 2 \
        "$BAM_FILE" \
        "$REORDERED_VCF" \
        "$BARCODE_FILE" \
        "${FASTDEMUX_DIR}/fastdemux.fdout.raw"
    log "fastdemux finished for sample: $SAMPLE_ID"
  else
    log "fastdemux output already exists — skipping."
  fi

done

log "=== scATAC fastdemux pipeline finished ==="
