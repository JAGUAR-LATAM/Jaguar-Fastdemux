# Jaguar-Fastdemux 🐆

**A streamlined workflow for donor assignment in pooled scATAC-seq experiments.**

This pipeline optimizes computational efficiency by intersecting SNPs with accessible chromatin regions (peaks) and standardizing chromosome ordering between BAM and VCF files, ensuring seamless processing with demultiplexing tools such as **FASTDEMUX**.

## Key Features

-   **Memory Optimization:** Only processes SNPs overlapping with accessible regions (`peaks.bed.gz`).

-   **Automatic Synchronization:** Corrects chromosome ordering discrepancies between the BAM header and the VCF file.

-   **Batch Processing:** Uses a `.tsv` configuration structure to handle multiple samples simultaneously.

-   **Quality Filters:** Strict selection of biallelic SNPs to ensure the accuracy of the additive model.

## Required Data Structure

### 1. BED Files Directory (`/beds`)

Each sample must have its peaks in a dedicated folder. Files must be **sorted by genomic coordinates** and indexed.

```         
|- beds
   |- ATLAS00/
      |- peaks.bed.gz      # Sorted peaks (bgzip)
      |- peaks.bed.gz.tbi  # Index (tabix)
   |- ATLASXX/
      |- peaks.bed.gz      
      |- peaks.bed.gz.tbi
```

### 2. Configuration File (`Config_pipeline.tsv`)

A tab-separated values file defining execution paths and parameters.

| SampleID | BAM_File | Barcode_File | Donor_List_File | Num_Donors |
|----|----|----|----|----|
| **ATLAS00** | `/path/to/possorted_bam.bam` | `/path/to/barcodes.tsv` | `/path/to/donors.txt` | `N` |
| **ATLASXX** | `/path/to/possorted_bam.bam` | `/path/to/barcodes.tsv` | `/path/to/donors.txt` | `N` |

-   **SampleID**: Unique identifier for the sample.

-   **BAM_File**: CellRanger output BAM file (position-sorted).

-   **Barcode_File**: List of barcodes (typically `barcodes.tsv`).

-   **Donor_List_File**: List of donor IDs present in the joint VCF file.

-   **Num_Donors**: Total number of expected donors.

### 3. Joint VCF

A population-level VCF file containing the genotypic information for all possible donors is required. The pipeline handles the extraction of relevant subsets.

## Requirements and Dependencies

To ensure stability during large-scale genomic data processing, the following is recommended:

-   **Hardware:** Minimum 16GB RAM.

-   **Bioinformatics Tools:**

    -   `samtools` & `bcftools` (BAM/VCF handling).

    -   `tabix` (BED indexing).

-   **Programming Environment:**

    -   **R/Python** for orchestration and filtering logic.

## Pipeline Logic

The workflow follows a sequence of steps to maximize the demultiplexing signal:

1.  **Order Extraction:** Reads the BAM header to determine the exact chromosome order (e.g., chr1..22, X, Y).

2.  **Donor Subsetting:** Filters the joint VCF to include only the donors specified for the current sample.

3.  **Variant Filtering:** - Removes Indels and multi-allelic variants.

    -   Retains only biallelic SNPs (additive model).

4.  **Peak Intersection:** Restricts SNPs to accessible regions defined in the BED file, reducing computational noise.

5.  **Reordering:** Adjusts the VCF header to match the BAM chromosome order exactly.

6.  **Execution:** Prepared files are passed to the assignment tool to generate the `fastdemux.fdout.raw` report.

## References
*Ranjbaran, A., Luca, F., & Pique-Regi, R. (2026). fastdemux: Robust SNP-based demultiplexing of single-cell population genomics data. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.64898/2026.02.10.705082*

