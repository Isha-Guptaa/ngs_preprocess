#!/bin/bash

# Setup directory structure
echo "Setting up directory structure..."
mkdir -p data results/fastqc results/fastp results/feature_extraction logs/fastqc logs/feature_extraction logs/fastp

echo "Directory structure created successfully!"
echo "Please place your FASTQ files in the 'data' directory."

# Check if data directory has FASTQ files
if [ ! "$(ls -A data 2>/dev/null)" ]; then
    echo "Warning: No files found in the data directory. Please add your FASTQ files."
fi

# Run the pipeline with increased latency wait
echo "Running the pipeline with increased latency wait..."
snakemake --cores 4 --use-conda --latency-wait 60 -p --force
