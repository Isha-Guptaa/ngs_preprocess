#!/bin/bash

# Create directory structure (if not already exists)
mkdir -p data
mkdir -p results/fastqc
mkdir -p results/fastp
mkdir -p results/feature_extraction
mkdir -p logs/fastqc
mkdir -p logs/feature_extraction
mkdir -p scripts

# Check if scripts exist and copy them if they don't
if [ ! -f "scripts/extract_features.py" ]; then
    echo "Copying extract_features.py to scripts directory..."
    # Add the code for extract_features.py here if needed
fi

if [ ! -f "scripts/summarize_features.py" ]; then
    echo "Copying summarize_features.py to scripts directory..."
    # Add the code for summarize_features.py here if needed
fi

# Make Python scripts executable
if [ -d "scripts" ]; then
    chmod +x scripts/*.py 2>/dev/null || true
fi

echo "Directory structure created successfully!"
echo "Please place your FASTQ files in the 'data' directory."