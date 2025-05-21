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
