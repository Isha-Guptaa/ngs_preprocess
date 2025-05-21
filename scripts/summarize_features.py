#!/usr/bin/env python3
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import re
from pathlib import Path

def parse_quality_dispersion_files(qual_disp_dir):
    """Parse quality dispersion files and combine into a single DataFrame."""
    quality_data = {}
    
    # Find all QUAL_DISP files
    qual_disp_files = glob.glob(os.path.join(qual_disp_dir, "*.QUAL_DISP"))
    
    for file_path in qual_disp_files:
        filename = os.path.basename(file_path)
        match = re.match(r'(.+)_(R[12])\.QUAL_DISP', filename)
        if match:
            sample_name = match.group(1)
            read = match.group(2)
            
            # Read the data
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Skip files with no data or headers only
            if len(lines) <= 1:
                continue
            
            # Extract position and quality data
            positions = []
            qualities = []
            
            for line in lines:
                if line.startswith('#') or not line.strip():
                    continue
                    
                try:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        position = int(parts[0])
                        quality = float(parts[1])
                        positions.append(position)
                        qualities.append(quality)
                except (ValueError, IndexError):
                    continue
            
            if positions and qualities:
                if sample_name not in quality_data:
                    quality_data[sample_name] = {}
                quality_data[sample_name][read] = {'positions': positions, 'qualities': qualities}
    
    return quality_data

def plot_mean_quality_by_position(quality_data, output_file):
    """Plot mean quality by position across all samples."""
    plt.figure(figsize=(12, 8))
    
    for sample_name, reads in quality_data.items():
        for read, data in reads.items():
            positions = data['positions']
            qualities = data['qualities']
            plt.plot(positions, qualities, label=f"{sample_name}_{read}")
    
    plt.xlabel('Position in read (bp)')
    plt.ylabel('Mean Quality Score')
    plt.title('Mean Quality Score by Position')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def generate_summary_plots(features_df, output_prefix):
    """Generate summary plots from feature data."""
    # Extract numeric values for certain metrics
    try:
        # Total sequences
        if 'Total Sequences' in features_df.columns:
            total_seq = features_df['Total Sequences'].astype(int)
            plt.figure(figsize=(10, 6))
            plt.bar(total_seq.index, total_seq)
            plt.title('Total Sequences per Sample')
            plt.ylabel('Number of Sequences')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_total_sequences.png")
            plt.close()
        
        # GC content
        if '%GC' in features_df.columns:
            gc_content = features_df['%GC'].str.rstrip('%').astype(float)
            plt.figure(figsize=(10, 6))
            plt.bar(gc_content.index, gc_content)
            plt.title('GC Content per Sample')
            plt.ylabel('GC Content (%)')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_gc_content.png")
            plt.close()
        
        # Mean Quality
        if 'Mean Quality' in features_df.columns:
            mean_quality = features_df['Mean Quality'].astype(float)
            plt.figure(figsize=(10, 6))
            plt.bar(mean_quality.index, mean_quality)
            plt.title('Mean Quality per Sample')
            plt.ylabel('Mean Quality Score')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"{output_prefix}_mean_quality.png")
            plt.close()
    except Exception as e:
        print(f"Error generating plots: {str(e)}")

def generate_summary_report(features_df, output_file):
    """Generate a summary report in Markdown format."""
    with open(output_file, 'w') as f:
        f.write("# NGS Quality Control Summary Report\n\n")
        
        # Basic statistics
        f.write("## Basic Statistics\n\n")
        f.write("| Sample | Total Sequences | Sequence Length | GC Content | Mean Quality |\n")
        f.write("|--------|----------------|----------------|------------|-------------|\n")
        
        for sample in features_df.index:
            total_seq = features_df.loc[sample, 'Total Sequences'] if 'Total Sequences' in features_df.columns else 'N/A'
            seq_len = features_df.loc[sample, 'Sequence length'] if 'Sequence length' in features_df.columns else 'N/A'
            gc = features_df.loc[sample, '%GC'] if '%GC' in features_df.columns else 'N/A'
            mean_qual = features_df.loc[sample, 'Mean Quality'] if 'Mean Quality' in features_df.columns else 'N/A'
            
            f.write(f"| {sample} | {total_seq} | {seq_len} | {gc} | {mean_qual} |\n")
        
        # Image references
        f.write("\n## Quality Control Plots\n\n")
        
        # Use correct filenames for the plots
        f.write("### Total Sequences\n\n")
        f.write(f"![Total Sequences](qc_summary_total_sequences.png)\n\n")
        
        f.write("### GC Content\n\n")
        f.write(f"![GC Content](qc_summary_gc_content.png)\n\n")
        
        f.write("### Mean Quality\n\n")
        f.write(f"![Mean Quality](qc_summary_mean_quality.png)\n\n")
        
        f.write("### Mean Quality by Position\n\n")
        f.write(f"![Mean Quality by Position](qc_summary_mean_quality.png)\n\n")

def main():
    parser = argparse.ArgumentParser(description="Summarize QC features and generate reports")
    parser.add_argument('--feature_dir', required=True, help="Directory containing feature files")
    parser.add_argument('--qual_disp_dir', required=True, help="Directory containing quality dispersion files")
    parser.add_argument('--output_prefix', required=True, help="Prefix for output files")
    
    args = parser.parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)
        
        # Read feature files
        feature_files = glob.glob(os.path.join(args.feature_dir, "*_features.csv"))
        if not feature_files:
            print(f"No feature files found in {args.feature_dir}")
            return
        
        # Combine all feature data
        all_features = pd.concat([pd.read_csv(file, index_col=0) for file in feature_files])
        
        # Parse quality dispersion files
        quality_data = parse_quality_dispersion_files(args.qual_disp_dir)
        
        # Generate summary plots
        generate_summary_plots(all_features, args.output_prefix)
        
        # Plot mean quality by position
        plot_mean_quality_by_position(quality_data, f"{args.output_prefix}_mean_quality.png")
        
        # Generate summary report
        generate_summary_report(all_features, f"{args.output_prefix}_summary_report.md")
        
        print(f"Summary report and plots generated with prefix {args.output_prefix}")
        
    except Exception as e:
        print(f"Error in summarize_features.py: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
