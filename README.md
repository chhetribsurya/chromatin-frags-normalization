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

### Option 1: Docker (Recommended)

The easiest way to run the analysis is using Docker, which handles all dependencies automatically.

#### Prerequisites
- **Docker** (version 20.10 or higher)
- **Docker Compose** (optional, for easier management)

#### Quick Start with Docker

1. **Clone the repository:**
```bash
git clone https://github.com/chhetribsurya/chromatin-frags-normalization.git
cd chromatin-frags-normalization
```

2. **Set up directory structure:**
```bash
./run_analysis.sh setup
```

3. **Build Docker image:**
```bash
./run_analysis.sh build
```

4. **Run analysis:**
```bash
# Batch mode
./run_analysis.sh run-batch

# Single sample mode
./run_analysis.sh run-single --sample-name sample1 --fragment-file /workspace/frags/sample1.bed --target-sites /workspace/input/target_peaks.bed
```

#### Manual Docker Usage

```bash
# Build the image
docker build -t chromatin-frags-normalization .

# Run batch analysis
docker run --rm \
  -v $(pwd)/input:/workspace/input:ro \
  -v $(pwd)/output:/workspace/output \
  -v $(pwd)/frags:/workspace/frags:ro \
  chromatin-frags-normalization \
  --samplesheet /workspace/input/samples.tsv \
  --target-sites /workspace/input/target_peaks.bed \
  --frags-dir /workspace/frags \
  --output-dir /workspace/output \
  --verbose

# Run batch analysis with reference site normalization
docker run --rm \
  -v $(pwd)/input:/workspace/input:ro \
  -v $(pwd)/output:/workspace/output \
  -v $(pwd)/frags:/workspace/frags:ro \
  chromatin-frags-normalization \
  --samplesheet /workspace/input/samples.tsv \
  --target-sites /workspace/input/target_peaks.bed \
  --reference-sites /workspace/input/housekeeping_genes.bed \
  --frags-dir /workspace/frags \
  --output-dir /workspace/output \
  --verbose

# Run single sample analysis
docker run --rm \
  -v $(pwd)/input:/workspace/input:ro \
  -v $(pwd)/output:/workspace/output \
  -v $(pwd)/frags:/workspace/frags:ro \
  chromatin-frags-normalization \
  --sample-name sample1 \
  --fragment-file /workspace/frags/sample1.bed \
  --target-sites /workspace/input/target_peaks.bed \
  --output-dir /workspace/output \
  --verbose

# Run single sample analysis with reference site normalization
docker run --rm \
  -v $(pwd)/input:/workspace/input:ro \
  -v $(pwd)/output:/workspace/output \
  -v $(pwd)/frags:/workspace/frags:ro \
  chromatin-frags-normalization \
  --sample-name sample1 \
  --fragment-file /workspace/frags/sample1.bed \
  --target-sites /workspace/input/target_peaks.bed \
  --reference-sites /workspace/input/housekeeping_genes.bed \
  --output-dir /workspace/output \
  --verbose

# Run analysis with all chromosomes
docker run --rm \
  -v $(pwd)/input:/workspace/input:ro \
  -v $(pwd)/output:/workspace/output \
  -v $(pwd)/frags:/workspace/frags:ro \
  chromatin-frags-normalization \
  --samplesheet /workspace/input/samples.tsv \
  --target-sites /workspace/input/target_peaks.bed \
  --reference-sites /workspace/input/housekeeping_genes.bed \
  --frags-dir /workspace/frags \
  --output-dir /workspace/output \
  --include-all-chr \
  --verbose
```

#### Using Helper Script for Docker Execution

The helper script provides convenient Docker commands with full parameter specification:

```bash
# Set script path and variables
SCRIPT_BASH="./run_analysis.sh"
SAMPLESHEET="input/samples.tsv"
FRAGS_DIR="frags"
TARGET_SITES="input/target_peaks.bed"
REFERENCE_SITES="input/housekeeping_genes.bed"
OUT_DIR_BATCH_AUTO="output/batch_auto"
OUT_DIR_BATCH_AUTO_WITHOUTREF="output/batch_auto_withoutref"
OUT_DIR_SINGLE_AUTO="output/single_auto"
OUT_DIR_SINGLE_AUTO_WITHOUTREF="output/single_auto_withoutref"
SAMPLE_NAME="sample001"
FRAG_FILE="frags/sample001.bed"

# Batch Analysis Commands (Docker)

# Basic batch analysis (Docker)
bash $SCRIPT_BASH \
  run-batch \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Batch analysis with reference site normalization (Docker)
bash $SCRIPT_BASH \
  run-batch-with-ref \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Batch analysis without reference site normalization (Docker)
bash $SCRIPT_BASH \
  run-batch-without-ref \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO_WITHOUTREF" \
  --verbose

# Batch analysis including all chromosomes (Docker)
bash $SCRIPT_BASH \
  run-batch-all-chr \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Single Sample Analysis Commands (Docker)

# Basic single sample analysis (Docker)
bash $SCRIPT_BASH \
  run-single \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose

# Single sample with reference site normalization (Docker)
bash $SCRIPT_BASH \
  run-single-with-ref \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose

# Single sample without reference site normalization (Docker)
bash $SCRIPT_BASH \
  run-single-without-ref \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO_WITHOUTREF" \
  --verbose

# Single sample including all chromosomes (Docker)
bash $SCRIPT_BASH \
  run-single-all-chr \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose
```

#### Using Docker Compose

```bash
# Edit docker-compose.yml to set your analysis parameters
# Then run:
docker-compose up

# Or run with custom command:
docker-compose run chromatin-counter --samplesheet /workspace/input/samples.tsv --target-sites /workspace/input/target_peaks.bed --verbose
```

### Option 2: Local Installation (No Docker Required)

Run the analysis directly on your system without Docker containers.

#### Prerequisites

- **R (version ≥ 4.0.0)**
- **Bioconductor** (for genomic analysis packages)

#### Quick Start with Local Installation

1. **Clone the repository:**
```bash
git clone https://github.com/chhetribsurya/chromatin-frags-normalization.git
cd chromatin-frags-normalization
```

2. **Set up directory structure:**
```bash
./run_analysis.sh setup
```

3. **Verify R installation:**
```bash
./run_analysis.sh validate
```

4. **Run analysis:**
```bash
# Batch mode (local)
./run_analysis.sh run-batch-local

# Single sample mode (local)
./run_analysis.sh run-single-local --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed
```

#### Manual Local Usage

Run the R script directly without the helper script:

```bash
# Run batch analysis
Rscript chromatin_count_norm_v2.R \
  --samplesheet input/samples.tsv \
  --target-sites input/target_peaks.bed \
  --frags-dir frags \
  --output-dir output \
  --verbose

# Run batch analysis with reference site normalization
Rscript chromatin_count_norm_v2.R \
  --samplesheet input/samples.tsv \
  --target-sites input/target_peaks.bed \
  --reference-sites input/housekeeping_genes.bed \
  --frags-dir frags \
  --output-dir output \
  --verbose

# Run single sample analysis
Rscript chromatin_count_norm_v2.R \
  --sample-name sample1 \
  --fragment-file frags/sample1.bed \
  --target-sites input/target_peaks.bed \
  --output-dir output \
  --verbose

# Run single sample analysis with reference site normalization
Rscript chromatin_count_norm_v2.R \
  --sample-name sample1 \
  --fragment-file frags/sample1.bed \
  --target-sites input/target_peaks.bed \
  --reference-sites input/housekeeping_genes.bed \
  --output-dir output \
  --verbose

# Run analysis with all chromosomes
Rscript chromatin_count_norm_v2.R \
  --samplesheet input/samples.tsv \
  --target-sites input/target_peaks.bed \
  --reference-sites input/housekeeping_genes.bed \
  --frags-dir frags \
  --output-dir output \
  --include-all-chr \
  --verbose
```

#### Using Helper Script for Local Execution

The helper script provides convenient local commands with full parameter specification:

```bash
# Set script path and variables
SCRIPT_BASH="./run_analysis.sh"
SAMPLESHEET="input/samples.tsv"
FRAGS_DIR="frags"
TARGET_SITES="input/target_peaks.bed"
REFERENCE_SITES="input/housekeeping_genes.bed"
OUT_DIR_BATCH_AUTO="output/batch_auto"
OUT_DIR_BATCH_AUTO_WITHOUTREF="output/batch_auto_withoutref"
OUT_DIR_SINGLE_AUTO="output/single_auto"
OUT_DIR_SINGLE_AUTO_WITHOUTREF="output/single_auto_withoutref"
SAMPLE_NAME="sample001"
FRAG_FILE="frags/sample001.bed"

# Batch Analysis Commands (Local)

# Basic batch analysis (local)
bash $SCRIPT_BASH \
  run-batch-local \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Batch analysis with reference site normalization (local)
bash $SCRIPT_BASH \
  run-batch-local-with-ref \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Batch analysis without reference site normalization (local)
bash $SCRIPT_BASH \
  run-batch-local-without-ref \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO_WITHOUTREF" \
  --verbose

# Batch analysis including all chromosomes (local)
bash $SCRIPT_BASH \
  run-batch-local-all-chr \
  --samplesheet "$SAMPLESHEET" \
  --frags-dir "$FRAGS_DIR" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_BATCH_AUTO" \
  --verbose

# Single Sample Analysis Commands (Local)

# Basic single sample analysis (local)
bash $SCRIPT_BASH \
  run-single-local \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose

# Single sample with reference site normalization (local)
bash $SCRIPT_BASH \
  run-single-local-with-ref \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose

# Single sample without reference site normalization (local)
bash $SCRIPT_BASH \
  run-single-local-without-ref \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO_WITHOUTREF" \
  --verbose

# Single sample including all chromosomes (local)
bash $SCRIPT_BASH \
  run-single-local-all-chr \
  --sample-name "$SAMPLE_NAME" \
  --fragment-file "$FRAG_FILE" \
  --target-sites "$TARGET_SITES" \
  --reference-sites "$REFERENCE_SITES" \
  --output-dir "$OUT_DIR_SINGLE_AUTO" \
  --verbose
```

### Option 3: Local Installation (Detailed Setup)

#### Prerequisites

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

### Example 1: Quick Start with Docker

1. **Set up the project:**
```bash
git clone https://github.com/chhetribsurya/chromatin-frags-normalization.git
cd chromatin-frags-normalization
./run_analysis.sh setup
```

2. **Copy example data to input directory:**
```bash
cp example_data/* input/
cp -r example_data/frags/* frags/
```

3. **Run analysis:**
```bash
./run_analysis.sh build
./run_analysis.sh run-batch
```

### Example 2: Basic Batch Analysis (Local)

Generate raw counts and CPM for multiple samples:

```bash
# Using helper script (local)
./run_analysis.sh run-batch-local --verbose

# Using R directly
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites enhancer_peaks.bed \
  --verbose
```

### Example 3: Batch Analysis with Reference Normalization

Include housekeeping gene normalization:

```bash
# Using Docker
./run_analysis.sh run-batch-with-ref --save-intermediate

# Using local helper script
./run_analysis.sh run-batch-local-with-ref --save-intermediate

# Using local R
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites enhancer_peaks.bed \
  --reference-sites housekeeping_promoters.bed \
  --output-dir ./results \
  --save-intermediate
```

### Example 4: Single Sample Analysis

Process one sample without samplesheet:

```bash
# Using Docker
./run_analysis.sh run-single-with-ref \
  --sample-name sample005 \
  --fragment-file frags/sample005_fragments.bed \
  --target-sites input/disease_associated_peaks.bed \
  --reference-sites input/dhs_sites.bed

# Using local helper script
./run_analysis.sh run-single-local-with-ref \
  --sample-name sample005 \
  --fragment-file frags/sample005_fragments.bed \
  --target-sites input/disease_associated_peaks.bed \
  --reference-sites input/dhs_sites.bed

# Using local R
Rscript chromatin_count_norm_v2.R \
  --sample-name sample005 \
  --fragment-file ./data/sample005_fragments.bed \
  --target-sites disease_associated_peaks.bed \
  --reference-sites dhs_sites.bed
```

### Example 5: Include All Chromosomes

Process all chromosomes including sex chromosomes and contigs:

```bash
# Using Docker
./run_analysis.sh run-batch-all-chr --verbose

# Using local helper script
./run_analysis.sh run-batch-local-all-chr --verbose

# Using local R
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --include-all-chr \
  --verbose
```

### Example 6: Re-run Analysis with Different Parameters

Force regeneration with new parameters:

```bash
# Using Docker
./run_analysis.sh run-batch-with-ref --regenerate-counts --include-all-chr

# Using local helper script
./run_analysis.sh run-batch-local-with-ref --regenerate-counts --include-all-chr

# Using local R
Rscript chromatin_count_norm_v2.R \
  --samplesheet samples.tsv \
  --target-sites peaks.bed \
  --reference-sites new_reference.bed \
  --regenerate-counts \
  --include-all-chr
```

### Example 7: Environment Variable Usage

Simplify commands using environment variables:

```bash
# Set environment variables
export SAMPLESHEET="my_samples.tsv"
export TARGET_SITES="my_peaks.bed"
export REFERENCE_SITES="my_ref_sites.bed"
export OUTPUT_DIR="my_results"

# Run with simplified commands
./run_analysis.sh run-batch-local-with-ref --verbose
./run_analysis.sh run-single-local-with-ref --sample-name sample001 --fragment-file frags/sample001.bed
```

## Repository Structure

```
chromatin-frags-normalization/
├── chromatin_count_norm_v2.R    # Main analysis script
├── run_analysis.sh              # Convenience script for running analyses
├── Dockerfile                   # Docker configuration
├── docker-compose.yml          # Docker Compose configuration
├── requirements.txt            # R package dependencies
├── .gitignore                  # Git ignore rules
├── .dockerignore              # Docker ignore rules
├── README.md                  # This file
├── example_data/              # Example input files
│   ├── samples.tsv
│   ├── target_peaks.bed
│   ├── housekeeping_genes.bed
│   └── frags/
│       ├── sample001_fragments.bed
│       ├── sample002_fragments.bed
│       ├── sample003_fragments.bed
│       └── sample004_fragments.bed
├── input/                     # Your input files (created by setup)
├── output/                    # Analysis results (created by setup)
├── frags/                     # Fragment files (created by setup)
└── logs/                      # Log files (created by setup)
```

### Input Directory Layout
```
input/
├── samples.tsv                # Sample list (batch mode)
├── target_peaks.bed           # Target regions
└── housekeeping_genes.bed     # Reference regions (optional)

frags/
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

## Helper Script Usage

The `run_analysis.sh` script provides comprehensive commands for both Docker and local execution:

### Available Commands

```bash
# Show all available commands
./run_analysis.sh help

# Show comprehensive examples
./run_analysis.sh examples

# Set up directory structure
./run_analysis.sh setup

# Build Docker image
./run_analysis.sh build

# Validate input files and configuration
./run_analysis.sh validate

# Clean up Docker resources
./run_analysis.sh clean
```

### Batch Analysis Commands

#### Docker Commands (Requires Docker)
```bash
# Basic batch analysis
./run_analysis.sh run-batch

# With reference site normalization
./run_analysis.sh run-batch-with-ref

# Without reference site normalization
./run_analysis.sh run-batch-without-ref

# Including all chromosomes
./run_analysis.sh run-batch-all-chr
```

#### Local Commands (No Docker Required)
```bash
# Basic batch analysis (local)
./run_analysis.sh run-batch-local

# With reference site normalization (local)
./run_analysis.sh run-batch-local-with-ref

# Without reference site normalization (local)
./run_analysis.sh run-batch-local-without-ref

# Including all chromosomes (local)
./run_analysis.sh run-batch-local-all-chr
```

### Single Sample Analysis Commands

#### Docker Commands (Requires Docker)
```bash
# Basic single sample analysis
./run_analysis.sh run-single --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# With reference site normalization
./run_analysis.sh run-single-with-ref --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# Without reference site normalization
./run_analysis.sh run-single-without-ref --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# Including all chromosomes
./run_analysis.sh run-single-all-chr --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed
```

#### Local Commands (No Docker Required)
```bash
# Basic single sample analysis (local)
./run_analysis.sh run-single-local --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# With reference site normalization (local)
./run_analysis.sh run-single-local-with-ref --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# Without reference site normalization (local)
./run_analysis.sh run-single-local-without-ref --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed

# Including all chromosomes (local)
./run_analysis.sh run-single-local-all-chr --sample-name sample1 --fragment-file frags/sample1.bed --target-sites input/peaks.bed
```

### Advanced Commands

```bash
# Run with custom parameters (pass through to R script)
./run_analysis.sh run-custom --samplesheet my_samples.tsv --target-sites my_peaks.bed --reference-sites my_ref.bed --verbose

# Run analysis locally with custom parameters
./run_analysis.sh run-local --samplesheet samples.tsv --target-sites peaks.bed --verbose
```

### Environment Variable Support

You can set environment variables to simplify command usage:

```bash
# Set environment variables
export SAMPLESHEET="my_samples.tsv"
export TARGET_SITES="my_peaks.bed"
export REFERENCE_SITES="my_ref_sites.bed"
export OUTPUT_DIR="my_results"
export SAMPLE_NAME="sample001"
export FRAGMENT_FILE="frags/sample001.bed"

# Then run with simplified commands
./run_analysis.sh run-batch-local-with-ref --verbose
./run_analysis.sh run-single-local-with-ref --verbose
```

## Use Cases

### ChIP-seq Signal at Candidate Regions
```bash
# Using Docker
./run_analysis.sh run-batch --reference-sites /workspace/input/active_promoters.bed

# Using local R
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

**1. Docker command not found:**
```
Error: docker: command not found
```
- Install Docker or use local commands instead
- Use `./run_analysis.sh run-batch-local` instead of `./run_analysis.sh run-batch`
- Use `./run_analysis.sh run-single-local` instead of `./run_analysis.sh run-single`

**2. Fragment files not found:**
```
Error: Missing fragment files for samples: sample_001
```
- Check file naming convention: `{sample_name}*.bed`
- Verify `--frags-dir` path
- Ensure files have `.bed` extension

**3. Empty count matrices:**
```
Warning: No overlaps found for sample: sample_001
```
- Verify chromosome naming consistency (chr1 vs 1)
- Check coordinate systems (0-based vs 1-based)
- Try `--include-all-chr` flag

**4. Memory issues with large datasets:**
```
Error: Cannot allocate vector of size X
```
- Process samples in batches
- Use single sample mode for very large datasets
- Increase system memory or use computing cluster

**5. Reference normalization fails:**
```
Error: Zero reference counts found for samples: sample_001
```
- Check reference sites overlap with fragments
- Verify reference BED file format
- Consider different reference regions

**6. R or Rscript not found (local execution):**
```
Error: Rscript not found
```
- Install R (version 4.0.0 or higher)
- Use Docker commands instead: `./run_analysis.sh run-batch`
- Add R to your system PATH

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

- **v2.1:** Added comprehensive local execution support, environment variable support, enhanced helper script with Docker and local commands, improved error handling and validation
- **v2.0:** Added single sample mode, flexible chromosome filtering, simplified input requirements
- **v1.0:** Initial batch processing implementation

---

**Last updated:** September 2025
