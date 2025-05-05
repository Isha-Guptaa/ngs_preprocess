# FastQC Quality Control and Feature Extraction Pipeline

A Snakemake-based bioinformatics workflow for analyzing the quality of next-generation sequencing (NGS) data. This pipeline runs FastQC on sequencing files and extracts important quality metrics into structured CSV files for further analysis.

## Overview

This pipeline performs two main tasks:

1. **Quality Control**: Runs FastQC on FASTQ files to generate quality reports
2. **Feature Extraction**: Extracts key metrics from FastQC reports for downstream analysis

## Requirements

- Python 3.6+
- Snakemake
- Conda (for environment management)
- Required Python packages:
  - pandas
  - zipfile (standard library)

## Directory Structure

```
.
├── data/              # Input FASTQ files (.fastq.gz, .fq.gz)
├── results/           # Output directory
│   ├── fastqc/        # FastQC HTML and ZIP reports
│   └── feature_extraction/ # Extracted quality metrics (CSV)
├── logs/              # Log files
│   ├── fastqc/        # FastQC logs
│   └── feature_extraction/ # Feature extraction logs
├── scripts/
│   └── extract_features.py # Script for extracting metrics from FastQC ZIP files
├── envs/
│   └── fastqc.yaml    # Conda environment for FastQC
└── Snakefile          # Main Snakemake workflow file
```

## Installation

1. Clone this repository:
   ```bash
   git clone <https://github.com/Isha-Guptaa/ngs_preprocess.git>
   cd fastqc-pipeline
   ```

2. Create the conda environment:
   ```bash
   conda env create -f envs/fastqc.yaml
   ```

## Usage

1. Place your FASTQ files in the `data/` directory (supported formats: `.fastq.gz`, `.fq.gz`)

2. Run the pipeline:
   ```bash
   snakemake --use-conda -j <cores>
   ```
   Replace `<cores>` with the number of CPU cores to use.

3. For a dry run to check the workflow:
   ```bash
   snakemake -n
   ```

## Input

The pipeline auto-detects paired-end or single-end sequencing data based on filenames. Common naming conventions are supported:
- Paired-end: `*_R1.fastq.gz` and `*_R2.fastq.gz` or `*_1.fastq.gz` and `*_2.fastq.gz`
- Single-end: Any other FASTQ file that doesn't match paired-end patterns

If no input files are found, the pipeline will create dummy files for testing purposes.

## Output

### FastQC Reports

FastQC HTML reports and ZIP files will be generated in the `results/fastqc/` directory:
- `{sample}_R1_fastqc.html` and `{sample}_R1_fastqc.zip` for forward reads
- `{sample}_R2_fastqc.html` and `{sample}_R2_fastqc.zip` for reverse reads (paired-end data only)

### Feature Extraction

For each sample, a CSV file containing extracted quality metrics is created in the `results/feature_extraction/` directory:
- `{sample}_features.csv`

These CSV files contain metrics such as:
- Sequence quality scores (mean, min, max)
- GC content
- Sequence length distribution
- Adapter content
- Overrepresented sequences

## Feature Extraction Details

The `extract_features.py` script extracts the following metrics from FastQC reports:

- **Basic Statistics**: Total sequences, sequence length, %GC
- **Sequence Quality**: Mean, minimum, and maximum quality scores
- **GC Content**: Mean GC content across all sequences
- **Sequence Length**: Average sequence length
- **Overrepresented Sequences**: Count and presence of overrepresented sequences
- **Adapter Content**: Maximum and mean adapter content

For paired-end data, metrics are extracted separately for forward (R1) and reverse (R2) reads.

