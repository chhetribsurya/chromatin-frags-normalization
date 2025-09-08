# Chromatin Fragment Count and Normalization

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.14+-green.svg)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A chromatin fragment counter and normalization pipeline for counting chromatin fragments at user-defined genomic sites with flexible normalization options.

## Overview

`chromatin_count_norm_v2.R` is a streamlined tool designed for quantifying chromatin accessibility or histone modification signals at specific genomic regions. It supports both batch processing of multiple samples and single-sample analysis, generating three types of quantitative matrices:

1. **Raw Counts**: Direct fragment overlap counts with target sites
2. **CPM (Counts Per Million)**: Library size-normalized counts  
3. **Reference-Normalized**: Counts normalized to user-defined reference regions

## Key Features

- **Dual Mode Operation**: Batch processing or single sample analysis
- **Flexible Normalization**: Optional reference-based normalization using housekeeping genes, DHS sites, or custom regions
- **Chromosome Filtering**: Choose between autosomes only or all chromosomes
- **Simple Input Format**: Standard BED files and minimal samplesheet requirements
- **Multiple Output Formats**: Both RDS (R native) and tab-separated text files

## Installation

### Prerequisites

- **R (version ≥ 4.0.0)**
- **Bioconductor** (for genomic analysis packages)

### Required R Packages

The script automatically installs missing packages, but you can install them manually if preferred:

#### CRAN Packages
```r
# Install CRAN packages
install.packages(c("dplyr", "readr", "stringr"))
```

#### Bioconductor Packages
```r
# Install BiocManager if not available
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("rtracklayer", "GenomicRanges"))
```

#### Complete Manual Installation
```r
# Complete installation script
cran_packages <- c("dplyr", "readr", "stringr")
bioc_packages <- c("rtracklayer", "GenomicRanges")

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(bioc_packages, dependencies = TRUE)
```

### Setup

1. **Download the script:**
```bash
wget https://github.com/chhetribsurya/chromatin-frags-normalization/blob/main/chromatin_count_norm_v2.R
chmod +x chromatin_count_norm_v2.R
```

2. **Verify installation (automatic dependency installation):**
```bash
Rscript chromatin_count_norm_v2.R --help
```

3. **Test with sample data:**
```bash
# The script will install missing packages on first run
Rscript chromatin_count_norm_v2.R --help
```

### Alternative Installation Methods

#### Using conda/mamba
```bash
# Create conda environment with R and dependencies
conda create -n chromatin_analysis r-base r-dplyr r-readr r-stringr
conda activate chromatin_analysis

# Install Bioconductor packages
R -e "BiocManager::install(c('rtracklayer', 'GenomicRanges'))"
```

#### Using renv (Reproducible Environments)
```r
# Initialize renv project
renv::init()

# Install packages
install.packages(c("dplyr", "readr", "stringr"))
BiocManager::install(c("rtracklayer", "GenomicRanges"))

# Save package state
renv::snapshot()
```

## Usage

### Command Line Syntax

**Batch Mode (Multiple Samples):**
```bash
Rscript chromatin_count_norm_v2.R --samplesheet SAMPLES.tsv --target-sites PEAKS.bed [OPTIONS]
```

**Single Sample Mode:**
```bash
Rscript chromatin_count_norm_v2.R --sample-name NAME --fragment-file FRAGS.bed --target-sites PEAKS.bed [OPTIONS]
```

### Command Line Options

| Option | Description | Default | Required |
|--------|-------------|---------|----------|
| `--samplesheet PATH` | Sample list file (batch mode) | - | Batch mode |
| `--sample-name NAME` | Sample identifier (single mode) | - | Single mode |
| `--fragment-file PATH` | Fragment BED file (single mode) | - | Single mode |
| `--target-sites PATH` | Target regions BED file | - | Yes |
| `--reference-sites PATH` | Reference regions for normalization | - | No |
| `--frags-dir PATH` | Fragment files directory (batch mode) | `./frags` | No |
| `--output-dir PATH` | Output directory | `./output` | No |
| `--include-all-chr` | Use all chromosomes (not just autosomes) | False | No |
| `--regenerate-counts` | Force regeneration of count matrices | False | No |
| `--save-intermediate` | Save intermediate files as text | False | No |
| `--verbose` | Detailed logging output | False | No |
| `--help` | Show detailed help message | - | No |

## Input File Formats

### 1. Samplesheet (Batch Mode Only)

**Format:** Tab-separated file with header

**Required Column:**
- `sample_name`: Unique identifier for each sample

**Example (`samples.tsv`):**
```
sample_name
sample001
sample002  
sample003
sample004
```

**Alternative with additional columns (ignored):**
```
sample_name	condition	batch
sample001	control	1
sample002	control	2
sample003	treated	1  
sample004	treated	2
```

### 2. Target Sites BED File

**Format:** Standard BED format (minimum 3 columns)

**Columns:**
1. Chromosome (chr1, chr2, etc.)
2. Start position (0-based)
3. End position (1-based)  
4. Peak name/ID (optional, auto-generated if missing)

**Example (`target_peaks.bed`):**
```
chr1	1000	2000	peak_1
chr1	5000	6000	peak_2
chr2	10000	11000	peak_3
chr3	15000	16000	peak_4
chr4	20000	21000	peak_5
```

### 3. Reference Sites BED File (Optional)

**Format:** Same as target sites

**Use Cases:**
- Housekeeping gene promoters
- DNase I hypersensitive sites (DHS)
- Constitutively accessible regions
- User-defined control regions

**Example (`housekeeping_genes.bed`):**
```
chr1	50000	51000	ACTB_promoter
chr2	100000	101000	GAPDH_promoter  
chr3	150000	151000	TUBB_promoter
chr4	200000	201000	RPL13_promoter
```

#### Comprehensive Housekeeping Genes Reference (`housekeeping_genes_comprehensive.bed`)

**Note: This is an example reference file only. Example coordinates are based on hg19/GRCh37 assembly, but should be verified for your specific analysis.**

```
chr7	5527148	5530601	ACTB
chr12	6643556	6647525	GAPDH
chr6	30720181	30723570	TUBB
chr16	89572621	89576362	RPL13
chrX	133466730	133778278	HPRT1
chr6	170863787	170896729	TBP
chr15	42451504	42457019	B2M
chr3	197394478	197412710	RPLP0
chr7	65240563	65263311	GUSB
chr7	44806747	44813876	PPIA
chr2	206639015	206644247	YWHAZ
chr19	49495789	49500052	TFRC
chr12	124911604	124954850	UBC
chr1	203409813	203418583	NONO
chrX	77348022	77385349	PGK1
```

### 4. Fragment Files

**Location:** `--frags-dir` (batch mode) or `--fragment-file` (single mode)

**Format:** Standard BED format with fragment coordinates

**Naming Convention (Batch Mode):** `{sample_name}*.bed`

**Example (`sample001_fragments.bed`):**
```
chr1	1500	1650
chr1	2100	2250
chr1	5200	5350
chr2	10500	10650
chr3	15200	15400
```

## Output Files

All matrices are saved in both RDS (R binary) and tab-separated text formats.

### 1. Raw Counts Matrix (`raw_counts_matrix.txt`)

Direct fragment overlap counts with target sites.

```
peak_coordinates	sample001	sample002	sample003	sample004
chr1_1000_2000	125	134	89	102
chr1_5000_6000	67	78	45	58
chr2_10000_11000	89	95	112	98
chr3_15000_16000	45	52	38	41
```

### 2. CPM Matrix (`cpm_matrix.txt`)

Counts normalized per million mapped fragments.

`CPM = (raw_counts / total_fragments_per_sample) × 1,000,000`

```
peak_coordinates	sample001	sample002	sample003	sample004
chr1_1000_2000	12.5	13.4	8.9	10.2
chr1_5000_6000	6.7	7.8	4.5	5.8
chr2_10000_11000	8.9	9.5	11.2	9.8
chr3_15000_16000	4.5	5.2	3.8	4.1
```

### 3. Reference-Normalized Matrix (`normalized_matrix.txt`)

Raw counts divided by total reference site counts per sample.

`Normalized = raw_counts / total_reference_counts_per_sample`

```
peak_coordinates	sample001	sample002	sample003	sample004
chr1_1000_2000	0.125	0.134	0.089	0.102
chr1_5000_6000	0.067	0.078	0.045	0.058
chr2_10000_11000	0.089	0.095	0.112	0.098
chr3_15000_16000	0.045	0.052	0.038	0.041
```

## Usage Examples

### Example 1: Basic Batch Analysis

Generate raw counts and CPM for multiple samples:

```bash
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites enhancer_peaks.bed \
  --verbose
```

### Example 2: Batch Analysis with Reference Normalization

Include housekeeping gene normalization:

```bash
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites enhancer_peaks.bed \
  --reference-sites housekeeping_promoters.bed \
  --output-dir ./results \
  --save-intermediate
```

### Example 3: Single Sample Analysis

Process one sample without samplesheet:

```bash
Rscript chromatin_count_norm_v2.R \
  --sample-name sample005 \
  --fragment-file ./data/sample005_fragments.bed \
  --target-sites disease_associated_peaks.bed \
  --reference-sites dhs_sites.bed
```

### Example 4: Include All Chromosomes

Process all chromosomes including sex chromosomes and contigs:

```bash
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --include-all-chr \
  --verbose
```

### Example 5: Re-run Analysis with Different Parameters

Force regeneration with new parameters:

```bash
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --reference-sites new_reference.bed \
  --regenerate-counts \
  --include-all-chr
```

## Directory Structure

### Input Directory Layout
```
project/
├── samples.tsv
├── target_peaks.bed
├── housekeeping_genes.bed (optional)
└── frags/
    ├── sample001_fragments.bed
    ├── sample002_fragments.bed
    ├── sample003_fragments.bed
    └── sample004_fragments.bed
```

### Output Directory Layout
```
output/
├── matrices/
│   ├── raw_counts_matrix.RDS
│   ├── raw_counts_matrix.txt
│   ├── cpm_matrix.RDS
│   ├── cpm_matrix.txt
│   ├── normalized_matrix.RDS (if reference provided)
│   ├── normalized_matrix.txt (if reference provided)
│   └── reference_counts_matrix.RDS (if reference provided)
└── logs/
    └── analysis_summary.txt
```

## Use Cases

### ChIP-seq Signal at Candidate Regions
```bash
# Measure H3K27ac signal at enhancer candidates  
Rscript chromatin_count_norm_v2.R \
  --samplesheet chip_samples.tsv \
  --target-sites enhancer_candidates.bed \
  --reference-sites active_promoters.bed
```

## Advanced Usage

### Working with Large Datasets

For datasets with many samples or large peak sets:

```bash
# Use parallel processing and save intermediate files
Rscript chromatin_count_norm_v2.R \
  --samplesheet large_cohort.tsv \
  --target-sites genome_wide_peaks.bed \
  --save-intermediate \
  --verbose \
  --output-dir ./large_analysis
```

### Quality Control Checks

Monitor analysis progress and check for issues:

```bash
# Enable detailed logging for troubleshooting
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --verbose \
  2>&1 | tee analysis.log
```

### Reprocessing Subsets

Rerun analysis with different normalization:

```bash
# First run: basic analysis
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed

# Second run: add reference normalization (reuses raw counts)
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --reference-sites new_reference.bed
```

## Data Integration

### Loading Results in R

```r
# Load count matrices
raw_counts <- readRDS("output/matrices/raw_counts_matrix.RDS")
cpm_matrix <- readRDS("output/matrices/cpm_matrix.RDS")
normalized_matrix <- readRDS("output/matrices/normalized_matrix.RDS")

# Basic exploration
dim(raw_counts)
summary(colSums(raw_counts))
cor(cpm_matrix)
```

### Loading Results in Python

```python
import pandas as pd
import numpy as np

# Load text format matrices
raw_counts = pd.read_csv("output/matrices/raw_counts_matrix.txt", 
                        sep="\t", index_col=0)
cpm_matrix = pd.read_csv("output/matrices/cpm_matrix.txt", 
                        sep="\t", index_col=0)

# Basic analysis
print(raw_counts.shape)
print(raw_counts.sum(axis=0))
```

## Troubleshooting

### Common Issues

**1. Fragment files not found:**
```
Error: Missing fragment files for samples: sample_001
```
- Check file naming convention: `{sample_name}*.bed`
- Verify `--frags-dir` path
- Ensure files have `.bed` extension

**2. Empty count matrices:**
```
Warning: No overlaps found for sample: sample_001
```
- Verify chromosome naming consistency (chr1 vs 1)
- Check coordinate systems (0-based vs 1-based)
- Try `--include-all-chr` flag

**3. Memory issues with large datasets:**
```
Error: Cannot allocate vector of size X
```
- Process samples in batches
- Use single sample mode for very large datasets
- Increase system memory or use computing cluster

**4. Reference normalization fails:**
```
Error: Zero reference counts found for samples: sample_001
```
- Check reference sites overlap with fragments
- Verify reference BED file format
- Consider different reference regions

### Debugging Tips

1. **Enable verbose output:** Use `--verbose` for detailed logging
2. **Check intermediate files:** Use `--save-intermediate` to inspect matrices
3. **Test with small datasets:** Verify workflow with subset of data
4. **Validate input formats:** Ensure BED files follow standard format

### Performance Optimization

- **Parallel processing:** R automatically uses multiple cores for large matrices
- **Disk space:** RDS files are more space-efficient than text files
- **Memory usage:** Process samples individually for very large datasets

## Technical Details

### Fragment Counting Algorithm

1. Import fragment BED files using `rtracklayer::import()`
2. Filter chromosomes based on `--include-all-chr` flag  
3. Count overlaps using `GenomicRanges::countOverlaps()`
4. No fragment collapsing (preserves full fragment overlaps)

### Normalization Methods

- **CPM:** Standard library size normalization
- **Reference-based:** Custom normalization using user-defined regions
- **No filtering:** Raw counts preserved without quality filters

### Chromosome Filtering

- **Default:** Autosomes only (chr1-chr22)
- **All chromosomes:** Includes sex chromosomes (chrX, chrY) and contigs

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Version History

- **v2.0:** Added single sample mode, flexible chromosome filtering, simplified input requirements
- **v1.0:** Initial batch processing implementation

---

**Last updated:** September 2025
