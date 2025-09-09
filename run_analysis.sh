#!/bin/bash

# =========================================================
# Chromatin Fragments Normalization Pipeline Script
# =========================================================
# 
# This script provides comprehensive commands for running the chromatin
# fragment normalization pipeline with support for:
# - Batch mode (multiple samples via samplesheet)
# - Single sample mode
# - Reference site normalization
# - Docker and local execution
# - All chromosome inclusion
# - Comprehensive logging and error handling
#
# =========================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/chromatin_count_norm_v2.R"
DOCKER_IMAGE="chromatin-frags-normalization"

# Default paths (can be overridden via environment variables)
DEFAULT_SAMPLESHEET="input/samples.tsv"
DEFAULT_FRAGS_DIR="frags"
DEFAULT_TARGET_SITES="input/target_peaks.bed"
DEFAULT_REFERENCE_SITES="input/housekeeping_genes.bed"
DEFAULT_OUTPUT_DIR="output"

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${PURPLE}[HEADER]${NC} $1"
}

print_command() {
    echo -e "${CYAN}[COMMAND]${NC} $1"
}

# Function to show comprehensive help
show_help() {
    cat << 'EOF'
=============================================================================
Chromatin Fragment Normalization - Comprehensive Analysis Runner
=============================================================================

DESCRIPTION:
  This script provides convenient commands for running the chromatin fragment
  normalization pipeline with support for both batch and single sample modes,
  reference site normalization, and comprehensive analysis options.

USAGE:
  ./run_analysis.sh [COMMAND] [OPTIONS]

COMMANDS:
  help                    Show this comprehensive help message
  setup                   Set up directory structure and example files
  build                   Build Docker image
  clean                   Clean up Docker containers and images
  
  # Batch Analysis Commands (Docker)
  run-batch               Run batch analysis with Docker
  run-batch-with-ref      Run batch analysis with reference site normalization (Docker)
  run-batch-without-ref   Run batch analysis without reference site normalization (Docker)
  run-batch-all-chr       Run batch analysis including all chromosomes (Docker)
  
  # Batch Analysis Commands (Local)
  run-batch-local         Run batch analysis locally (without Docker)
  run-batch-local-with-ref      Run batch analysis with reference site normalization (Local)
  run-batch-local-without-ref   Run batch analysis without reference site normalization (Local)
  run-batch-local-all-chr       Run batch analysis including all chromosomes (Local)
  
  # Single Sample Analysis Commands (Docker)
  run-single              Run single sample analysis with Docker
  run-single-with-ref     Run single sample analysis with reference normalization (Docker)
  run-single-without-ref  Run single sample analysis without reference normalization (Docker)
  run-single-all-chr      Run single sample analysis including all chromosomes (Docker)
  
  # Single Sample Analysis Commands (Local)
  run-single-local        Run single sample analysis locally (without Docker)
  run-single-local-with-ref     Run single sample analysis with reference normalization (Local)
  run-single-local-without-ref  Run single sample analysis without reference normalization (Local)
  run-single-local-all-chr      Run single sample analysis including all chromosomes (Local)
  
  # Advanced Commands
  run-custom              Run with custom parameters (pass through to R script)
  validate                Validate input files and configuration
  examples                Show comprehensive usage examples

ENVIRONMENT VARIABLES:
  SAMPLESHEET             Path to samplesheet file (default: input/samples.tsv)
  FRAGS_DIR               Path to fragments directory (default: frags)
  TARGET_SITES            Path to target sites BED file (default: input/target_peaks.bed)
  REFERENCE_SITES         Path to reference sites BED file (default: input/housekeeping_genes.bed)
  OUTPUT_DIR              Path to output directory (default: output)
  SAMPLE_NAME             Name for single sample analysis
  FRAGMENT_FILE           Path to fragment file for single sample analysis

EXAMPLES:
  # Basic setup
  ./run_analysis.sh setup
  ./run_analysis.sh build
  
  # Batch analysis examples
  ./run_analysis.sh run-batch
  ./run_analysis.sh run-batch-with-ref
  ./run_analysis.sh run-batch-without-ref
  ./run_analysis.sh run-batch-all-chr
  
  # Single sample analysis examples
  ./run_analysis.sh run-single --sample-name sample001 --fragment-file frags/sample001.bed
  ./run_analysis.sh run-single-with-ref --sample-name sample001 --fragment-file frags/sample001.bed
  ./run_analysis.sh run-single-without-ref --sample-name sample001 --fragment-file frags/sample001.bed
  
  # Custom analysis with specific parameters
  ./run_analysis.sh run-custom --samplesheet my_samples.tsv --target-sites my_peaks.bed --reference-sites my_ref.bed --verbose
  
  # Local execution (without Docker)
  ./run_analysis.sh run-batch-local
  ./run_analysis.sh run-single-local --sample-name sample001 --fragment-file frags/sample001.bed

=============================================================================
EOF
}

# Function to show comprehensive examples
show_examples() {
    cat << 'EOF'
=============================================================================
COMPREHENSIVE USAGE EXAMPLES
=============================================================================

1. BASIC SETUP AND PREPARATION:
   # Set up directory structure
   ./run_analysis.sh setup
   
   # Build Docker image
   ./run_analysis.sh build
   
   # Validate your input files
   ./run_analysis.sh validate

2. BATCH ANALYSIS (Multiple Samples):
   # Basic batch analysis (raw counts + CPM)
   ./run_analysis.sh run-batch
   
   # Batch analysis with reference site normalization
   ./run_analysis.sh run-batch-with-ref
   
   # Batch analysis without reference site normalization
   ./run_analysis.sh run-batch-without-ref
   
   # Batch analysis including all chromosomes (including sex chromosomes)
   ./run_analysis.sh run-batch-all-chr
   
   # Batch analysis with custom parameters
   ./run_analysis.sh run-custom \
     --samplesheet my_samples.tsv \
     --frags-dir my_frags \
     --target-sites my_peaks.bed \
     --reference-sites my_ref_sites.bed \
     --output-dir my_results \
     --include-all-chr \
     --verbose

3. SINGLE SAMPLE ANALYSIS:
   # Basic single sample analysis
   ./run_analysis.sh run-single \
     --sample-name sample001 \
     --fragment-file frags/sample001.bed
   
   # Single sample with reference normalization
   ./run_analysis.sh run-single-with-ref \
     --sample-name sample001 \
     --fragment-file frags/sample001.bed
   
   # Single sample without reference normalization
   ./run_analysis.sh run-single-without-ref \
     --sample-name sample001 \
     --fragment-file frags/sample001.bed
   
   # Single sample including all chromosomes
   ./run_analysis.sh run-single-all-chr \
     --sample-name sample001 \
     --fragment-file frags/sample001.bed

4. LOCAL EXECUTION (Without Docker):
   # Batch analysis locally
   ./run_analysis.sh run-batch-local
   
   # Single sample analysis locally
   ./run_analysis.sh run-single-local \
     --sample-name sample001 \
     --fragment-file frags/sample001.bed

5. ADVANCED USAGE:
   # Custom analysis with all options
   ./run_analysis.sh run-custom \
     --samplesheet samples.tsv \
     --frags-dir frags \
     --target-sites peaks.bed \
     --reference-sites housekeeping.bed \
     --output-dir results \
     --include-all-chr \
     --regenerate-counts \
     --save-intermediate \
     --verbose

6. ENVIRONMENT VARIABLE USAGE:
   # Set environment variables for repeated use
   export SAMPLESHEET="my_samples.tsv"
   export TARGET_SITES="my_peaks.bed"
   export REFERENCE_SITES="my_ref.bed"
   export OUTPUT_DIR="my_results"
   
   # Then run with defaults
   ./run_analysis.sh run-batch-with-ref

=============================================================================
EOF
}

# Function to set up directory structure
setup_directories() {
    print_header "Setting up directory structure and example files..."
    
    # Create directories
    mkdir -p input output frags logs
    
    # Create example samplesheet if it doesn't exist
    if [ ! -f "input/samples.tsv" ]; then
        cat > input/samples.tsv << 'EOF'
sample_name
sample001
sample002
sample003
sample004
EOF
        print_success "Created example samplesheet: input/samples.tsv"
    fi
    
    # Create example target sites if it doesn't exist
    if [ ! -f "input/target_peaks.bed" ]; then
        cat > input/target_peaks.bed << 'EOF'
chr1	1000	2000	peak1
chr1	5000	6000	peak2
chr2	1000	2000	peak3
chr2	5000	6000	peak4
EOF
        print_success "Created example target sites: input/target_peaks.bed"
    fi
    
    # Create example reference sites if it doesn't exist
    if [ ! -f "input/housekeeping_genes.bed" ]; then
        cat > input/housekeeping_genes.bed << 'EOF'
chr1	10000	11000	housekeeping1
chr1	20000	21000	housekeeping2
chr2	10000	11000	housekeeping3
chr2	20000	21000	housekeeping4
EOF
        print_success "Created example reference sites: input/housekeeping_genes.bed"
    fi
    
    print_success "Directory structure created:"
    echo "  input/     - Place your input files here (samplesheets, BED files)"
    echo "  output/    - Analysis results will be saved here"
    echo "  frags/     - Place your fragment BED files here"
    echo "  logs/      - Log files will be saved here"
}

# Function to build Docker image
build_docker() {
    print_header "Building Docker image..."
    docker build -t "$DOCKER_IMAGE" .
    print_success "Docker image built successfully"
}

# Function to clean up Docker resources
clean_docker() {
    print_header "Cleaning up Docker resources..."
    
    # Remove containers
    docker ps -a --filter "name=chromatin-frags-normalization" --format "table {{.Names}}" | grep -v NAMES | xargs -r docker rm
    
    # Remove images
    docker images --filter "reference=chromatin-frags-normalization" --format "table {{.Repository}}" | grep -v REPOSITORY | xargs -r docker rmi
    
    print_success "Docker cleanup completed"
}

# Function to validate input files
validate_inputs() {
    print_header "Validating input files and configuration..."
    
    local errors=0
    
    # Check if R script exists
    if [ ! -f "$SCRIPT" ]; then
        print_error "R script not found: $SCRIPT"
        ((errors++))
    else
        print_success "R script found: $SCRIPT"
    fi
    
    # Check if Docker is available
    if command -v docker &> /dev/null; then
        print_success "Docker is available"
    else
        print_warning "Docker not found - local execution only"
    fi
    
    # Check if R is available for local execution
    if command -v Rscript &> /dev/null; then
        print_success "Rscript is available for local execution"
    else
        print_warning "Rscript not found - Docker execution only"
    fi
    
    # Check input files if they exist
    if [ -f "${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}" ]; then
        print_success "Samplesheet found: ${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    else
        print_warning "Samplesheet not found: ${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    fi
    
    if [ -f "${TARGET_SITES:-$DEFAULT_TARGET_SITES}" ]; then
        print_success "Target sites found: ${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    else
        print_warning "Target sites not found: ${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    fi
    
    if [ -f "${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}" ]; then
        print_success "Reference sites found: ${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    else
        print_warning "Reference sites not found: ${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    fi
    
    if [ -d "${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}" ]; then
        local frag_count=$(find "${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}" -name "*.bed" | wc -l)
        print_success "Fragments directory found: ${FRAGS_DIR:-$DEFAULT_FRAGS_DIR} ($frag_count .bed files)"
    else
        print_warning "Fragments directory not found: ${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    fi
    
    if [ $errors -eq 0 ]; then
        print_success "Validation completed successfully"
    else
        print_error "Validation completed with $errors errors"
        exit 1
    fi
}

# Function to run Docker command
run_docker() {
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    print_command "Running Docker command: $*"
    
    docker run --rm \
        -v "$(pwd)/input:/workspace/input:ro" \
        -v "$(pwd)/output:/workspace/output" \
        -v "$(pwd)/frags:/workspace/frags:ro" \
        "$DOCKER_IMAGE" \
        "$@"
}

# Function to run local command
run_local() {
    print_command "Running local command: Rscript $SCRIPT $*"
    
    if ! command -v Rscript &> /dev/null; then
        print_error "Rscript not found. Please install R (version 4.0.0 or higher)"
        exit 1
    fi
    
    Rscript "$SCRIPT" "$@"
}

# Function to run batch analysis with Docker
run_batch_docker() {
    print_header "Running batch analysis with Docker..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --samplesheet "/workspace/$samplesheet" \
        --frags-dir "/workspace/$frags_dir" \
        --target-sites "/workspace/$target_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis completed"
}

# Function to run batch analysis locally
run_batch_local() {
    print_header "Running batch analysis locally..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --samplesheet "$samplesheet" \
        --frags-dir "$frags_dir" \
        --target-sites "$target_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis completed"
}

# Function to run batch analysis locally with reference sites
run_batch_local_with_ref() {
    print_header "Running batch analysis locally with reference site normalization..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    if [ ! -f "$reference_sites" ]; then
        print_error "Reference sites not found: $reference_sites"
        exit 1
    fi
    
    run_local \
        --samplesheet "$samplesheet" \
        --frags-dir "$frags_dir" \
        --target-sites "$target_sites" \
        --reference-sites "$reference_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis with reference normalization completed"
}

# Function to run batch analysis locally without reference sites
run_batch_local_without_ref() {
    print_header "Running batch analysis locally without reference site normalization..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --samplesheet "$samplesheet" \
        --frags-dir "$frags_dir" \
        --target-sites "$target_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis without reference normalization completed"
}

# Function to run batch analysis locally with all chromosomes
run_batch_local_all_chr() {
    print_header "Running batch analysis locally including all chromosomes..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --samplesheet "$samplesheet" \
        --frags-dir "$frags_dir" \
        --target-sites "$target_sites" \
        --reference-sites "$reference_sites" \
        --output-dir "$output_dir" \
        --include-all-chr \
        --verbose \
        "$@"
    
    print_success "Batch analysis with all chromosomes completed"
}

# Function to run batch analysis with reference sites
run_batch_with_ref() {
    print_header "Running batch analysis with reference site normalization..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    if [ ! -f "$reference_sites" ]; then
        print_error "Reference sites not found: $reference_sites"
        exit 1
    fi
    
    run_docker \
        --samplesheet "/workspace/$samplesheet" \
        --frags-dir "/workspace/$frags_dir" \
        --target-sites "/workspace/$target_sites" \
        --reference-sites "/workspace/$reference_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis with reference normalization completed"
}

# Function to run batch analysis without reference sites
run_batch_without_ref() {
    print_header "Running batch analysis without reference site normalization..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --samplesheet "/workspace/$samplesheet" \
        --frags-dir "/workspace/$frags_dir" \
        --target-sites "/workspace/$target_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Batch analysis without reference normalization completed"
}

# Function to run batch analysis with all chromosomes
run_batch_all_chr() {
    print_header "Running batch analysis including all chromosomes..."
    
    local samplesheet="${SAMPLESHEET:-$DEFAULT_SAMPLESHEET}"
    local frags_dir="${FRAGS_DIR:-$DEFAULT_FRAGS_DIR}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required files
    if [ ! -f "$samplesheet" ]; then
        print_error "Samplesheet not found: $samplesheet"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --samplesheet "/workspace/$samplesheet" \
        --frags-dir "/workspace/$frags_dir" \
        --target-sites "/workspace/$target_sites" \
        --reference-sites "/workspace/$reference_sites" \
        --output-dir "/workspace/$output_dir" \
        --include-all-chr \
        --verbose \
        "$@"
    
    print_success "Batch analysis with all chromosomes completed"
}

# Function to run single sample analysis with Docker
run_single_docker() {
    print_header "Running single sample analysis with Docker..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --sample-name "$sample_name" \
        --fragment-file "/workspace/$fragment_file" \
        --target-sites "/workspace/$target_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis completed"
}

# Function to run single sample analysis locally
run_single_local() {
    print_header "Running single sample analysis locally..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --sample-name "$sample_name" \
        --fragment-file "$fragment_file" \
        --target-sites "$target_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis completed"
}

# Function to run single sample analysis locally with reference sites
run_single_local_with_ref() {
    print_header "Running single sample analysis locally with reference site normalization..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    if [ ! -f "$reference_sites" ]; then
        print_error "Reference sites not found: $reference_sites"
        exit 1
    fi
    
    run_local \
        --sample-name "$sample_name" \
        --fragment-file "$fragment_file" \
        --target-sites "$target_sites" \
        --reference-sites "$reference_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis with reference normalization completed"
}

# Function to run single sample analysis locally without reference sites
run_single_local_without_ref() {
    print_header "Running single sample analysis locally without reference site normalization..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --sample-name "$sample_name" \
        --fragment-file "$fragment_file" \
        --target-sites "$target_sites" \
        --output-dir "$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis without reference normalization completed"
}

# Function to run single sample analysis locally with all chromosomes
run_single_local_all_chr() {
    print_header "Running single sample analysis locally including all chromosomes..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_local \
        --sample-name "$sample_name" \
        --fragment-file "$fragment_file" \
        --target-sites "$target_sites" \
        --reference-sites "$reference_sites" \
        --output-dir "$output_dir" \
        --include-all-chr \
        --verbose \
        "$@"
    
    print_success "Single sample analysis with all chromosomes completed"
}

# Function to run single sample analysis with reference sites
run_single_with_ref() {
    print_header "Running single sample analysis with reference site normalization..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    if [ ! -f "$reference_sites" ]; then
        print_error "Reference sites not found: $reference_sites"
        exit 1
    fi
    
    run_docker \
        --sample-name "$sample_name" \
        --fragment-file "/workspace/$fragment_file" \
        --target-sites "/workspace/$target_sites" \
        --reference-sites "/workspace/$reference_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis with reference normalization completed"
}

# Function to run single sample analysis without reference sites
run_single_without_ref() {
    print_header "Running single sample analysis without reference site normalization..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --sample-name "$sample_name" \
        --fragment-file "/workspace/$fragment_file" \
        --target-sites "/workspace/$target_sites" \
        --output-dir "/workspace/$output_dir" \
        --verbose \
        "$@"
    
    print_success "Single sample analysis without reference normalization completed"
}

# Function to run single sample analysis with all chromosomes
run_single_all_chr() {
    print_header "Running single sample analysis including all chromosomes..."
    
    local sample_name="${SAMPLE_NAME}"
    local fragment_file="${FRAGMENT_FILE}"
    local target_sites="${TARGET_SITES:-$DEFAULT_TARGET_SITES}"
    local reference_sites="${REFERENCE_SITES:-$DEFAULT_REFERENCE_SITES}"
    local output_dir="${OUTPUT_DIR:-$DEFAULT_OUTPUT_DIR}"
    
    # Check required parameters
    if [ -z "$sample_name" ]; then
        print_error "Sample name not provided. Use --sample-name or set SAMPLE_NAME environment variable"
        exit 1
    fi
    
    if [ -z "$fragment_file" ]; then
        print_error "Fragment file not provided. Use --fragment-file or set FRAGMENT_FILE environment variable"
        exit 1
    fi
    
    if [ ! -f "$target_sites" ]; then
        print_error "Target sites not found: $target_sites"
        exit 1
    fi
    
    run_docker \
        --sample-name "$sample_name" \
        --fragment-file "/workspace/$fragment_file" \
        --target-sites "/workspace/$target_sites" \
        --reference-sites "/workspace/$reference_sites" \
        --output-dir "/workspace/$output_dir" \
        --include-all-chr \
        --verbose \
        "$@"
    
    print_success "Single sample analysis with all chromosomes completed"
}

# Function to run custom analysis
run_custom() {
    print_header "Running custom analysis with provided parameters..."
    
    # Check if R script exists
    if [ ! -f "$SCRIPT" ]; then
        print_error "R script not found: $SCRIPT"
        exit 1
    fi
    
    # Run with all provided arguments
    run_local "$@"
    
    print_success "Custom analysis completed"
}

# Parse command line arguments for environment variables
parse_env_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --sample-name)
                export SAMPLE_NAME="$2"
                shift 2
                ;;
            --fragment-file)
                export FRAGMENT_FILE="$2"
                shift 2
                ;;
            --samplesheet)
                export SAMPLESHEET="$2"
                shift 2
                ;;
            --frags-dir)
                export FRAGS_DIR="$2"
                shift 2
                ;;
            --target-sites)
                export TARGET_SITES="$2"
                shift 2
                ;;
            --reference-sites)
                export REFERENCE_SITES="$2"
                shift 2
                ;;
            --output-dir)
                export OUTPUT_DIR="$2"
                shift 2
                ;;
            *)
                # Unknown option, pass through
                break
                ;;
        esac
    done
}

# Main script logic
case "${1:-help}" in
    help)
        show_help
        ;;
    examples)
        show_examples
        ;;
    setup)
        setup_directories
        ;;
    build)
        build_docker
        ;;
    clean)
        clean_docker
        ;;
    validate)
        validate_inputs
        ;;
    run-batch)
        shift
        parse_env_args "$@"
        run_batch_docker
        ;;
    run-batch-local)
        shift
        parse_env_args "$@"
        run_batch_local
        ;;
    run-batch-local-with-ref)
        shift
        parse_env_args "$@"
        run_batch_local_with_ref
        ;;
    run-batch-local-without-ref)
        shift
        parse_env_args "$@"
        run_batch_local_without_ref
        ;;
    run-batch-local-all-chr)
        shift
        parse_env_args "$@"
        run_batch_local_all_chr
        ;;
    run-batch-with-ref)
        shift
        parse_env_args "$@"
        run_batch_with_ref
        ;;
    run-batch-without-ref)
        shift
        parse_env_args "$@"
        run_batch_without_ref
        ;;
    run-batch-all-chr)
        shift
        parse_env_args "$@"
        run_batch_all_chr
        ;;
    run-single)
        shift
        parse_env_args "$@"
        run_single_docker
        ;;
    run-single-local)
        shift
        parse_env_args "$@"
        run_single_local
        ;;
    run-single-local-with-ref)
        shift
        parse_env_args "$@"
        run_single_local_with_ref
        ;;
    run-single-local-without-ref)
        shift
        parse_env_args "$@"
        run_single_local_without_ref
        ;;
    run-single-local-all-chr)
        shift
        parse_env_args "$@"
        run_single_local_all_chr
        ;;
    run-single-with-ref)
        shift
        parse_env_args "$@"
        run_single_with_ref
        ;;
    run-single-without-ref)
        shift
        parse_env_args "$@"
        run_single_without_ref
        ;;
    run-single-all-chr)
        shift
        parse_env_args "$@"
        run_single_all_chr
        ;;
    run-custom)
        shift
        run_custom "$@"
        ;;
    *)
        print_error "Unknown command: $1"
        echo ""
        show_help
        exit 1
        ;;
esac