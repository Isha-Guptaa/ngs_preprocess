#!/usr/bin/env python3
import os
import argparse
import gzip
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def parse_fastq(fastq_file, max_reads=None, verbose=False):
    """Parse FASTQ file and return quality scores."""
    quality_scores = []
    read_count = 0
    error_count = 0
    
    try:
        if verbose:
            print(f"Attempting to open file: {fastq_file}")
        
        # Check if file is gzipped
        is_gzipped = fastq_file.endswith('.gz')
        
        # Open file accordingly
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        with open_func(fastq_file, mode) as f:
            line_count = 0
            while True:
                if max_reads and read_count >= max_reads:
                    if verbose:
                        print(f"Reached maximum read count of {max_reads}")
                    break
                
                # Read 4 lines at a time (FASTQ record)
                try:
                    header = f.readline().strip()
                    line_count += 1
                    
                    if not header:  # End of file
                        break
                    
                    if not header.startswith('@'):
                        if verbose:
                            print(f"Warning: Invalid header at line {line_count}: {header[:50]}...")
                        error_count += 1
                        # Try to resync by reading until we find another header
                        continue
                    
                    seq = f.readline().strip()
                    line_count += 1
                    
                    plus = f.readline().strip()
                    line_count += 1
                    
                    if not plus.startswith('+'):
                        if verbose:
                            print(f"Warning: Invalid '+' line at line {line_count}: {plus[:50]}...")
                        error_count += 1
                        continue
                    
                    qual = f.readline().strip()
                    line_count += 1
                    
                    if len(seq) != len(qual):
                        if verbose:
                            print(f"Warning: Sequence and quality length mismatch at read {read_count+1}: {len(seq)} vs {len(qual)}")
                        error_count += 1
                        if error_count > 10:  # Give up after too many errors
                            raise ValueError(f"Too many format errors in {fastq_file}")
                        continue
                    
                    # Convert quality scores
                    ascii_offset = 33  # Assuming Phred+33 encoding
                    try:
                        read_quals = [ord(q) - ascii_offset for q in qual]
                        quality_scores.append(read_quals)
                        read_count += 1
                        
                        if read_count % 1000000 == 0 and verbose:
                            print(f"Processed {read_count} reads...")
                    
                    except Exception as e:
                        if verbose:
                            print(f"Error processing quality string: {str(e)}")
                        error_count += 1
                
                except Exception as e:
                    if verbose:
                        print(f"Error reading FASTQ entry: {str(e)}")
                    error_count += 1
                    if error_count > 100:  # Give up after too many errors
                        raise ValueError(f"Too many reading errors in {fastq_file}")
    
    except Exception as e:
        print(f"Error reading {fastq_file}: {str(e)}", file=sys.stderr)
        return []
    
    if verbose or read_count > 0:
        print(f"Processed {read_count} reads from {fastq_file}, encountered {error_count} errors")
    
    return quality_scores

def calculate_dispersion(quality_scores, verbose=False):
    """Calculate quality score dispersion."""
    if not quality_scores:
        if verbose:
            print("No quality scores provided")
        return None, None, None
    
    try:
        # Find the maximum sequence length
        max_len = max(len(q) for q in quality_scores)
        
        if verbose:
            print(f"Max read length: {max_len}")
            print(f"Number of reads: {len(quality_scores)}")
        
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
    
    except Exception as e:
        print(f"Error calculating dispersion: {str(e)}", file=sys.stderr)
        if verbose:
            import traceback
            traceback.print_exc()
        return None, None, None

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

def check_file_integrity(fastq_file, verbose=False):
    """Quick check for file integrity."""
    try:
        is_gzipped = fastq_file.endswith('.gz')
        open_func = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'
        
        with open_func(fastq_file, mode) as f:
            # Try to read a few lines to check format
            lines = [f.readline().strip() for _ in range(4)]
            
            if not lines[0] or not lines[0].startswith('@'):
                print(f"Warning: First line doesn't look like a FASTQ header: {lines[0][:50]}")
                return False
                
            if not lines[2] or not lines[2].startswith('+'):
                print(f"Warning: Third line doesn't start with '+': {lines[2][:50]}")
                return False
                
            if len(lines[1]) != len(lines[3]):
                print(f"Warning: Sequence and quality length mismatch: {len(lines[1])} vs {len(lines[3])}")
                return False
            
            if verbose:
                print(f"File check passed for {fastq_file}")
                print(f"Header: {lines[0][:50]}...")
                print(f"Sequence length: {len(lines[1])}")
            
            return True
            
    except Exception as e:
        print(f"Error checking file {fastq_file}: {str(e)}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="Analyze quality score dispersion in FASTQ data")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file (can be gzipped)")
    parser.add_argument("--outfile", required=True, help="Output file for dispersion data")
    parser.add_argument("--max_reads", type=int, default=100000, help="Maximum number of reads to process")
    parser.add_argument("--plot", action="store_true", help="Generate quality dispersion plots")
    parser.add_argument("--verbose", action="store_true", help="Print verbose debugging information")
    parser.add_argument("--sample", action="store_true", help="Process only a small sample to check format")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.fastq):
        print(f"Error: Input file {args.fastq} does not exist", file=sys.stderr)
        # Create minimal output file to satisfy Snakemake dependencies
        os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
        with open(args.outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed - input file not found\n")
        sys.exit(1)
    
    # Print file size for debugging
    if args.verbose:
        file_size = os.path.getsize(args.fastq)
        print(f"File size of {args.fastq}: {file_size/1024/1024:.2f} MB")
    
    # Quick file integrity check
    if not check_file_integrity(args.fastq, args.verbose):
        print(f"Warning: File integrity check failed for {args.fastq}", file=sys.stderr)
        if not args.sample:  # Continue anyway unless we're just sampling
            print("Continuing analysis despite integrity check failure...")
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
    
    # If sample mode, just check a few reads
    if args.sample:
        sample_reads = parse_fastq(args.fastq, max_reads=10, verbose=True)
        print(f"Sample analysis complete. Found {len(sample_reads)} reads.")
        if len(sample_reads) > 0:
            print(f"First read length: {len(sample_reads[0])}")
            print(f"First read quality scores: {sample_reads[0][:10]}... (first 10 positions)")
        sys.exit(0)
    
    # Parse FASTQ and calculate quality dispersion
    quality_scores = parse_fastq(args.fastq, args.max_reads, args.verbose)
    
    if not quality_scores:
        print("Warning: No quality scores were obtained", file=sys.stderr)
        # Create minimal output file
        with open(args.outfile, 'w') as f:
            f.write("# Quality dispersion analysis failed - no quality scores obtained\n")
        sys.exit(1)
    
    means, stdevs, cv = calculate_dispersion(quality_scores, args.verbose)
    
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
