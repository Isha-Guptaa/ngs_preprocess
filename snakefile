import os
import glob
from pathlib import Path

# Define output directories
FASTQC_DIR = "results/fastqc"
FASTP_DIR = "results/fastp"
FEATURE_DIR = "results/feature_extraction"

# Function to detect if samples are single-end or paired-end
def detect_paired_end():
    fastq_files = glob.glob("data/*.fastq.gz") + glob.glob("data/*.fq.gz")
    
    # If no files found, create dummy files for testing
    if not fastq_files:
        print("No FASTQ files found. Creating dummy files for testing...")
        os.makedirs("data", exist_ok=True)
        
        # Create small dummy FASTQ files
        with open("data/SRR12345678_R1.fastq.gz", "wb") as f:
            f.write(b"dummy fastq content")
        with open("data/SRR12345678_R2.fastq.gz", "wb") as f:
            f.write(b"dummy fastq content")
            
        # Update the list of files
        fastq_files = glob.glob("data/*.fastq.gz") + glob.glob("data/*.fq.gz")
    
    # Dictionary to store sample information
    samples_info = {}
    
    for file_path in fastq_files:
        file_name = os.path.basename(file_path)
        
        # Look for common paired-end patterns in the filename
        if "_R1" in file_name or "_1" in file_name:
            # This is a forward read
            base_name = file_name.replace("_R1", "").replace("_1", "").split(".")[0]
            if base_name not in samples_info:
                samples_info[base_name] = {"R1": file_path, "is_paired": False}
            else:
                samples_info[base_name]["R1"] = file_path
        
        elif "_R2" in file_name or "_2" in file_name:
            # This is a reverse read
            base_name = file_name.replace("_R2", "").replace("_2", "").split(".")[0]
            if base_name not in samples_info:
                samples_info[base_name] = {"R2": file_path, "is_paired": False}
            else:
                samples_info[base_name]["R2"] = file_path
                samples_info[base_name]["is_paired"] = True
        
        else:
            # This is likely a single-end read
            base_name = file_name.split(".")[0]
            samples_info[base_name] = {"R1": file_path, "is_paired": False}
    
    # Determine if any sample is paired-end
    is_paired = any(info["is_paired"] for info in samples_info.values())
    return is_paired, list(samples_info.keys()), samples_info

# Get sample information
IS_PAIRED, SAMPLES, SAMPLES_INFO = detect_paired_end()

# Ensure output directories exist
os.makedirs(FASTQC_DIR, exist_ok=True)
os.makedirs(FEATURE_DIR, exist_ok=True)
os.makedirs("logs/fastqc", exist_ok=True)
os.makedirs("logs/feature_extraction", exist_ok=True)

# Rule all to define the final targets
rule all:
    input:
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.zip"), sample=SAMPLES),
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.html"), sample=SAMPLES) if IS_PAIRED else [],
        expand(os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.zip"), sample=SAMPLES) if IS_PAIRED else [],
        expand(os.path.join(FEATURE_DIR, "{sample}_features.csv"), sample=SAMPLES)

# Function to get the correct input file for a given sample and read
def get_input_fastq(wildcards):
    sample = wildcards.sample
    read = wildcards.read
    
    if read in SAMPLES_INFO[sample]:
        return SAMPLES_INFO[sample][read]
    else:
        raise ValueError(f"No {read} file found for sample {sample}")

# Function to get FastQC output name from input file
def get_fastqc_name(input_file):
    base = os.path.basename(input_file)
    return f"{os.path.splitext(os.path.splitext(base)[0])[0]}_fastqc"

# Rule to run FastQC
rule fastqc:
    input:
        get_input_fastq
    output:
        html=os.path.join(FASTQC_DIR, "{sample}_{read}_fastqc.html"),
        zip=os.path.join(FASTQC_DIR, "{sample}_{read}_fastqc.zip")
    params:
        outdir=FASTQC_DIR
    log:
        "logs/fastqc/{sample}_{read}.log"
    threads: 4
    conda: "envs/fastqc.yaml"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p {params.outdir}
        
        # Run FastQC
        fastqc -o {params.outdir} -t {threads} {input} &> {log}
        
        # Get the name of the output files that FastQC actually created
        fastqc_basename=$(basename {input} | sed 's/.fastq.gz/_fastqc/g' | sed 's/.fq.gz/_fastqc/g')
        
        # If the files don't match the expected output names, move them
        if [ -f "{params.outdir}/${{fastqc_basename}}.html" ] && [ "{params.outdir}/${{fastqc_basename}}.html" != "{output.html}" ]; then
            mv "{params.outdir}/${{fastqc_basename}}.html" "{output.html}"
        fi
        
        if [ -f "{params.outdir}/${{fastqc_basename}}.zip" ] && [ "{params.outdir}/${{fastqc_basename}}.zip" != "{output.zip}" ]; then
            mv "{params.outdir}/${{fastqc_basename}}.zip" "{output.zip}"
        fi
        """

# Rule to extract features from FastQC data
rule extract_features:
    input:
        zip=lambda wildcards: [os.path.join(FASTQC_DIR, f"{wildcards.sample}_R1_fastqc.zip")] + 
            ([os.path.join(FASTQC_DIR, f"{wildcards.sample}_R2_fastqc.zip")] if IS_PAIRED and "R2" in SAMPLES_INFO[wildcards.sample] else [])
    output:
        features=os.path.join(FEATURE_DIR, "{sample}_features.csv")
    log:
        "logs/feature_extraction/{sample}.log"
    conda: "envs/fastqc.yaml"
    script:
        "scripts/extract_features.py"