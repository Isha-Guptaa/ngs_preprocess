# QC Pipeline

A Snakemake-based bioinformatics workflow for analyzing the quality of next-generation sequencing (NGS) data. This pipeline runs FastQC on sequencing files and extracts important quality metrics into structured CSV files for further analysis.

⚙️ Step 1 - Create the Conda Environment

Create a conda environment named qc_pipeline using the following command:

conda create --name qc_pipeline

📄 Step 2 - Create the Environment YAML File

qc_pipeline.yaml defines all dependencies:

name: qc_pipeline
channels:
 - conda-forge
 - bioconda
 - defaults
dependencies:
 - python=3.9
 - pandas=1.5.*
 - fastqc=0.12.*
 - biopython=1.81
 - numpy=1.24.*
 - scipy=1.10.*
 - matplotlib=3.7.*
 - pip=23.1.*
 - pip:
   - seaborn==0.12.*

🔁 Step 3 - Update Conda Environment

Once the YAML file is created, update the environment:

conda env update --name qc_pipeline --file qc_pipeline.yaml

Activate the environment:

conda activate qc_pipeline

🐍 Step 4 - Run the Pipeline

To execute the Snakemake pipeline:

snakemake --cores <num_cores> --use-conda

Example:

snakemake --cores 4 --use-conda

🧠 Scripts Overview

    extract_features.py
    Extracts basic statistics and per-base sequence quality from FastQC .zip files.

    quality_dispersion.py
    Computes variability in sequence quality.

    summarize_features.py
    Aggregates extracted features into a summarized CSV for downstream analysis.
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
   snakemake --cores 4 --use-conda
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

