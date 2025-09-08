# Use official R base image with Bioconductor
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('dplyr', 'readr', 'stringr'), dependencies=TRUE, repos='https://cran.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('rtracklayer', 'GenomicRanges'), dependencies=TRUE, ask=FALSE)"

# Copy the main script
COPY chromatin_count_norm_v2.R /workspace/

# Make the script executable
RUN chmod +x /workspace/chromatin_count_norm_v2.R

# Create directories for input and output
RUN mkdir -p /workspace/input /workspace/output /workspace/frags

# Set the default command
ENTRYPOINT ["Rscript", "/workspace/chromatin_count_norm_v2.R"]
