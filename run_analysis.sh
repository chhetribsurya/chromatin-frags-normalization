#!/bin/bash

# Chromatin Fragment Normalization - Analysis Runner Script
# This script provides convenient commands for running the analysis

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

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

# Function to show help
show_help() {
    echo "Chromatin Fragment Normalization - Analysis Runner"
    echo "=================================================="
    echo ""
    echo "Usage: $0 [COMMAND] [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  help                    Show this help message"
    echo "  build                   Build Docker image"
    echo "  run-batch              Run batch analysis with Docker"
    echo "  run-single             Run single sample analysis with Docker"
    echo "  run-local              Run analysis locally (without Docker)"
    echo "  clean                  Clean up Docker containers and images"
    echo "  setup                  Set up directory structure"
    echo ""
    echo "Examples:"
    echo "  $0 setup"
    echo "  $0 build"
    echo "  $0 run-batch --samplesheet samples.tsv --target-sites peaks.bed"
    echo "  $0 run-single --sample-name sample1 --fragment-file frags/sample1.bed --target-sites peaks.bed"
    echo "  $0 run-local --samplesheet samples.tsv --target-sites peaks.bed --verbose"
    echo ""
}

# Function to set up directory structure
setup_directories() {
    print_status "Setting up directory structure..."
    
    mkdir -p input
    mkdir -p output
    mkdir -p frags
    mkdir -p logs
    
    print_success "Directory structure created:"
    echo "  input/     - Place your input files here (samplesheets, BED files)"
    echo "  output/    - Analysis results will be saved here"
    echo "  frags/     - Place your fragment BED files here"
    echo "  logs/      - Log files will be saved here"
}

# Function to build Docker image
build_docker() {
    print_status "Building Docker image..."
    docker build -t chromatin-frags-normalization .
    print_success "Docker image built successfully"
}

# Function to run batch analysis with Docker
run_batch_docker() {
    print_status "Running batch analysis with Docker..."
    
    # Check if required files exist
    if [ ! -f "input/samples.tsv" ]; then
        print_error "Samplesheet not found: input/samples.tsv"
        print_status "Please create input/samples.tsv with your sample names"
        exit 1
    fi
    
    if [ ! -f "input/target_peaks.bed" ]; then
        print_error "Target sites BED file not found: input/target_peaks.bed"
        print_status "Please create input/target_peaks.bed with your target regions"
        exit 1
    fi
    
    # Run the analysis
    docker run --rm \
        -v "$(pwd)/input:/workspace/input:ro" \
        -v "$(pwd)/output:/workspace/output" \
        -v "$(pwd)/frags:/workspace/frags:ro" \
        chromatin-frags-normalization \
        --samplesheet /workspace/input/samples.tsv \
        --target-sites /workspace/input/target_peaks.bed \
        --frags-dir /workspace/frags \
        --output-dir /workspace/output \
        --verbose \
        "$@"
    
    print_success "Batch analysis completed"
}

# Function to run single sample analysis with Docker
run_single_docker() {
    print_status "Running single sample analysis with Docker..."
    
    # Check if required files exist
    if [ ! -f "input/target_peaks.bed" ]; then
        print_error "Target sites BED file not found: input/target_peaks.bed"
        print_status "Please create input/target_peaks.bed with your target regions"
        exit 1
    fi
    
    # Run the analysis
    docker run --rm \
        -v "$(pwd)/input:/workspace/input:ro" \
        -v "$(pwd)/output:/workspace/output" \
        -v "$(pwd)/frags:/workspace/frags:ro" \
        chromatin-frags-normalization \
        --target-sites /workspace/input/target_peaks.bed \
        --output-dir /workspace/output \
        --verbose \
        "$@"
    
    print_success "Single sample analysis completed"
}

# Function to run analysis locally
run_local() {
    print_status "Running analysis locally..."
    
    # Check if R is available
    if ! command -v Rscript &> /dev/null; then
        print_error "Rscript not found. Please install R (version 4.0.0 or higher)"
        exit 1
    fi
    
    # Run the analysis
    Rscript chromatin_count_norm_v2.R "$@"
    
    print_success "Local analysis completed"
}

# Function to clean up Docker resources
clean_docker() {
    print_status "Cleaning up Docker resources..."
    
    # Remove containers
    docker ps -a --filter "name=chromatin-frags-normalization" --format "table {{.Names}}" | grep -v NAMES | xargs -r docker rm
    
    # Remove images
    docker images --filter "reference=chromatin-frags-normalization" --format "table {{.Repository}}" | grep -v REPOSITORY | xargs -r docker rmi
    
    print_success "Docker cleanup completed"
}

# Main script logic
case "${1:-help}" in
    help)
        show_help
        ;;
    setup)
        setup_directories
        ;;
    build)
        build_docker
        ;;
    run-batch)
        shift
        run_batch_docker "$@"
        ;;
    run-single)
        shift
        run_single_docker "$@"
        ;;
    run-local)
        shift
        run_local "$@"
        ;;
    clean)
        clean_docker
        ;;
    *)
        print_error "Unknown command: $1"
        echo ""
        show_help
        exit 1
        ;;
esac
