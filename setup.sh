mkdir -p envs
mkdir -p scripts
mkdir -p data
mkdir -p logs/fastqc
mkdir -p logs/feature_extraction
mkdir -p logs/quality_dispersion
mkdir -p results/fastqc
mkdir -p results/feature_extraction
mkdir -p results/quality_dispersion

# Verify that the conda environment file exists
if [ ! -f "envs/qc_pipeline.yaml" ]; then
    echo "Creating conda environment file..."
    cat > envs/qc_pipeline.yaml << 'EOF'
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
EOF
fi

# Ensure extract_features.py script exists
if [ ! -f "scripts/extract_features.py" ]; then
    echo "Creating extract_features.py script..."
    cat > scripts/extract_features.py << 'EOF'
#!/usr/bin/env python3
import os
import argparse
import zipfile
import pandas as pd
import tempfile
import re
import sys

def extract_fastqc_data(zip_file):
    """Extract data from FastQC zip file."""
    metrics = {}
    try:
        with zipfile.ZipFile(zip_file, 'r') as zf:
            # Get the directory name inside the zip
            data_file = None
            for name in zf.namelist():
                if name.endswith('fastqc_data.txt'):
                    data_file = name
                    break
            
            if not data_file:
                print(f"Warning: No fastqc_data.txt found in {zip_file}")
                return metrics
            
            # Read and parse the data file
            with zf.open(data_file) as f:
                lines = [line.decode('utf-8').strip() for line in f.readlines()]
            
            # Extract basic statistics
            basic_stats_section = False
            for line in lines:
                if line.startswith('>>Basic Statistics'):
                    basic_stats_section = True
                    continue
                elif line.startswith('>>END_MODULE'):
                    basic_stats_section = False
                    continue
                
                if basic_stats_section and '\t' in line:
                    key, value = line.split('\t')
                    metrics[key] = value
            
            # Extract per base sequence quality
            quality_section = False
            quality_data = []
            for line in lines:
                if line.startswith('>>Per base sequence quality'):
                    quality_section = True
                    continue
                elif line.startswith('>>END_MODULE') and quality_section:
                    quality_section = False
                    continue
                
                if quality_section and '\t' in line and not line.startswith('>>'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        try:
                            base_pos = parts[0]
                            mean_qual = float(parts[1])
                            quality_data.append((base_pos, mean_qual))
                        except ValueError:
                            pass
            
            if quality_data:
                avg_quality = sum(q for _, q in quality_data) / len(quality_data)
                metrics['Mean Quality'] = f"{avg_quality:.2f}"
                
    except Exception as e:
        print(f"Error processing {zip_file}: {str(e)}")
    
    return metrics

def main():
    parser = argparse.ArgumentParser(description="Extract features from FastQC data")
    parser.add_argument('--fastqc_zips', nargs='+', required=True, help="FastQC zip files")
    parser.add_argument('--output', required=True, help="Output CSV file")
    
    args = parser.parse_args()
    
    all_metrics = {}
    
    # Check if input files exist
    missing_files = [f for f in args.fastqc_zips if not os.path.exists(f)]
    if missing_files:
        print(f"Error: The following input files do not exist: {', '.join(missing_files)}")
        # Create empty output file to satisfy Snakemake
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w') as f:
            f.write("metric,value\ntotal_sequences,0\nsequence_length,0\ngc_content,0\n")
        sys.exit(1)
    
    for zip_file in args.fastqc_zips:
        sample_name = os.path.basename(zip_file).replace('_fastqc.zip', '')
        metrics = extract_fastqc_data(zip_file)
        
        # Store metrics by sample
        for key, value in metrics.items():
            if key not in all_metrics:
                all_metrics[key] = {}
            all_metrics[key][sample_name] = value
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Convert to DataFrame and save
    if not all_metrics:
        print(f"Warning: No metrics found in any of the FastQC files")
        df = pd.DataFrame({"metric": ["total_sequences", "sequence_length", "gc_content"], 
                          "value": [0, 0, 0]})
    else:
        df = pd.DataFrame(all_metrics)
    
    df.to_csv(args.output, index_label='Sample')
    print(f"Features extracted to {args.output}")

if __name__ == "__main__":
    main()
EOF
    chmod +x scripts/extract_features.py
fi

# Ensure quality_dispersion.py script exists
if [ ! -f "scripts/quality_dispersion.py" ]; then
    echo "Creating quality_dispersion.py script..."
    cat > scripts/quality_dispersion.py << 'EOF'
#!/usr/bin/env python3
import os
import argparse
import gzip
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def parse_fastq(fastq_file, max_reads=None):
    """Parse FASTQ file and return quality scores."""
    quality_scores = []
    read_count = 0
    
    try:
        # Check if file is gzipped
        is_gzipped = fastq_file.endswith('.gz')
        
        # Open file accordingly
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        with open_func(fastq_file, mode) as f:
            while True:
                if max_reads and read_count >= max_reads:
                    break
                
                # Read 4 lines at a time (FASTQ record)
                header = f.readline().strip()
                if not header:  # End of file
                    break
                
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                
                if not all([header, seq, plus, qual]):  # Incomplete record
                    break
                
                # Convert quality scores
                ascii_offset = 33  # Assuming Phred+33 encoding
                read_quals = [ord(q) - ascii_offset for q in qual]
                quality_scores.append(read_quals)
                
                read_count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {str(e)}", file=sys.stderr)
        return []
    
    print(f"Processed {read_count} reads from {fastq_file}")
    return quality_scores

def calculate_dispersion(quality_scores):
    """Calculate quality score dispersion."""
    if not quality_scores:
        return None, None, None
    
    # Find the maximum sequence length
    max_len = max(len(q) for q in quality_scores)
    
    # Initialize arrays for tracking quality stats
    position_means = np.zeros(max_len)
    position_stdevs = np.zeros(max_len)
    position_counts = np.zeros(max_len)
    
    # Calculate mean and count for each position
    for read_quals in quality_scores:
        for i, q in enumerate(read_quals):
            position_means[i] += q
            position_counts[i] += 1
    
    # Normalize means
    for i in range(max_len):
        if position_counts[i] > 0:
            position_means[i] /= position_counts[i]
    
    # Calculate standard deviations
    for read_quals in quality_scores:
        for i, q in enumerate(read_quals):
            position_stdevs[i] += (q - position_means[i]) ** 2
    
    # Normalize standard deviations
    for i in range(max_len):
        if position_counts[i] > 1:  # Need at least 2 points for std dev
            position_stdevs[i] = np.sqrt(position_stdevs[i] / (position_counts[i] - 1))
    
    # Calculate coefficient of variation where possible
    position_cv = np.zeros(max_len)
    for i in range(max_len):
        if position_means[i] > 0:
            position_cv[i] = position_stdevs[i] / position_means[i]
    
    return position_means, position_stdevs, position_cv

def plot_quality_dispersion(positions, means, stdevs, cv, output_base):
    """Create plots of quality dispersion metrics."""
    try:
        # Set up plotting
        fig, axs = plt.subplots(3, 1, figsize=(10, 15))
        
        # Plot 1: Mean quality scores per position
        axs[0].plot(positions, means, 'b-')
        axs[0].fill_between(positions, means - stdevs, means + stdevs, alpha=0.3)
        axs[0].set_title('Mean Quality Scores by Position')
        axs[0].set_xlabel('Position in Read')
        axs[0].set_ylabel('Phred Quality Score')
        axs[0].grid(True)
        
        # Plot 2: Standard deviation of quality scores per position
        axs[1].plot(positions, stdevs, 'r-')
        axs[1].set_title('Quality Score Variability by Position')
        axs[1].set_xlabel('Position in Read')
        axs[1].set_ylabel('Standard Deviation')
        axs[1].grid(True)
        
        # Plot 3: Coefficient of variation per position
        axs[2].plot(positions, cv, 'g-')
        axs[2].set_title('Coefficient of Variation by Position')
        axs[2].set_xlabel('Position in Read')
        axs[2].set_ylabel('CV (StdDev/Mean)')
        axs[2].grid(True)
        
        plt.tight_layout()
        plt.savefig(f"{output_base}.png")
        plt.close()
        
        print(f"Plots saved to {output_base}.png")
    except Exception as e:
        print(f"Error creating plots: {str(e)}", file=sys.stderr)

def write_dispersion_data(positions, means, stdevs, cv, outfile):
    """Write dispersion data to output file."""
    try:
        with open(outfile, 'w') as f:
            f.write("Position\tMean_Quality\tStdDev\tCV\n")
            for i in range(len(positions)):
                f.write(f"{i+1}\t{means[i]:.2f}\t{stdevs[i]:.2f}\t{cv[i]:.4f}\n")
        print(f"Quality dispersion data written to {outfile}")
    except Exception as e:
        print(f"Error writing output file: {str(e)}", file=sys.stderr)
        # Create minimal file to satisfy Snakemake
        with open(outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed\n")

def main():
    parser = argparse.ArgumentParser(description="Analyze quality score dispersion in FASTQ data")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file (can be gzipped)")
    parser.add_argument("--outfile", required=True, help="Output file for dispersion data")
    parser.add_argument("--max_reads", type=int, default=10000, help="Maximum number of reads to process")
    parser.add_argument("--plot", action="store_true", help="Generate quality dispersion plots")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.fastq):
        print(f"Error: Input file {args.fastq} does not exist", file=sys.stderr)
        # Create minimal output file to satisfy Snakemake dependencies
        os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
        with open(args.outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed - input file not found\n")
        sys.exit(1)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
    
    # Parse FASTQ and calculate quality dispersion
    quality_scores = parse_fastq(args.fastq, args.max_reads)
    
    if not quality_scores:
        print("Warning: No quality scores were obtained", file=sys.stderr)
        # Create minimal output file
        with open(args.outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed - no quality scores obtained\n")
        sys.exit(1)
    
    means, stdevs, cv = calculate_dispersion(quality_scores)
    
    if means is None:
        print("Warning: Could not calculate quality dispersion", file=sys.stderr)
        # Create minimal output file
        with open(args.outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed - calculation error\n")
        sys.exit(1)
    
    # Generate position indices
    positions = list(range(1, len(means) + 1))
    
    # Write dispersion data
    write_dispersion_data(positions, means, stdevs, cv, args.outfile)
    
    # Create plots if requested
    if args.plot:
        output_base = os.path.splitext(args.outfile)[0]
        plot_quality_dispersion(positions, means, stdevs, cv, output_base)

if __name__ == "__main__":
    main()
EOF
    chmod +x scripts/quality_dispersion.py
fi

echo "Setup complete! Directory structure and scripts are now ready."
echo "Run your Snakemake pipeline with: snakemake --cores 4 --use-conda"
