#!/usr/bin/env python3
"""
Script to extract quality metrics from FastQC zip files.
"""

import os
import sys
import zipfile
import pandas as pd
import logging

# Set up logging to both console and log file
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler()
    ]
)

def extract_fastqc_data(zip_file_path):
    """
    Extract relevant data from FastQC zip file.
    
    Parameters:
    -----------
    zip_file_path : str
        Path to the FastQC zip file
    
    Returns:
    --------
    dict
        Dictionary of extracted metrics
    """
    metrics = {}
    
    try:
        logging.info(f"Processing {zip_file_path}")
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            # Find the fastqc_data.txt file in the archive
            fastqc_data_files = [f for f in zip_ref.namelist() if f.endswith('fastqc_data.txt')]
            
            if not fastqc_data_files:
                logging.warning(f"No fastqc_data.txt found in {zip_file_path}")
                return {}
                
            # Use the first one found
            fastqc_data_file = fastqc_data_files[0]
            logging.info(f"Found data file: {fastqc_data_file}")
            
            with zip_ref.open(fastqc_data_file) as f:
                content = f.read().decode('utf-8')
                
                # Extract basic statistics
                if ">>Basic Statistics" in content:
                    section = content.split(">>Basic Statistics")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    for line in lines:
                        if '\t' in line:
                            key, value = line.split('\t')
                            metrics[f"basic_{key.strip().lower().replace(' ', '_')}"] = value.strip()
                
                # Extract sequence quality
                if ">>Per base sequence quality" in content:
                    section = content.split(">>Per base sequence quality")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    
                    quality_data = []
                    for line in lines:
                        if '\t' in line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) >= 3:  # Base, Mean, Median, etc.
                                try:
                                    # Skip header row that might contain 'Mean' as text
                                    if parts[1].lower() == 'mean':
                                        continue
                                    quality_data.append(float(parts[1]))  # Mean quality
                                except ValueError:
                                    continue
                    
                    if quality_data:
                        metrics["mean_sequence_quality"] = sum(quality_data) / len(quality_data)
                        metrics["min_sequence_quality"] = min(quality_data)
                        metrics["max_sequence_quality"] = max(quality_data)
                
                # Extract GC content
                if ">>Per sequence GC content" in content:
                    section = content.split(">>Per sequence GC content")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    
                    gc_values = []
                    counts = []
                    for line in lines:
                        if '\t' in line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) >= 2:
                                try:
                                    # Skip header rows
                                    if parts[0].lower() in ['gc', 'gc content']:
                                        continue
                                    gc = float(parts[0])
                                    count = float(parts[1])
                                    gc_values.append(gc)
                                    counts.append(count)
                                except ValueError:
                                    continue
                    
                    if gc_values and counts:
                        weighted_sum = sum(g * c for g, c in zip(gc_values, counts))
                        total_count = sum(counts)
                        if total_count > 0:
                            metrics["mean_gc_content"] = weighted_sum / total_count
                
                # Extract sequence length
                if ">>Sequence Length Distribution" in content:
                    section = content.split(">>Sequence Length Distribution")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    
                    lengths = []
                    counts = []
                    for line in lines:
                        if '\t' in line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) >= 2:
                                try:
                                    # Skip header rows
                                    if parts[0].lower() in ['length', 'sequence length']:
                                        continue
                                        
                                    length_range = parts[0]
                                    count = float(parts[1])
                                    
                                    # Parse the length range (e.g., "40-45" or "75")
                                    if '-' in length_range:
                                        start, end = map(int, length_range.split('-'))
                                        length = (start + end) / 2
                                    else:
                                        length = int(length_range)
                                    
                                    lengths.append(length)
                                    counts.append(count)
                                except (ValueError, TypeError):
                                    continue
                    
                    if lengths and counts:
                        weighted_sum = sum(l * c for l, c in zip(lengths, counts))
                        total_count = sum(counts)
                        if total_count > 0:
                            metrics["mean_sequence_length"] = weighted_sum / total_count
                
                # Extract overrepresented sequences
                if ">>Overrepresented sequences" in content:
                    section = content.split(">>Overrepresented sequences")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    
                    if len(lines) > 1:  # If there are overrepresented sequences (more than just the header)
                        metrics["has_overrepresented_sequences"] = True
                        
                        # Count overrepresented sequences
                        overrep_count = 0
                        for line in lines:
                            if '\t' in line:
                                overrep_count += 1
                        
                        metrics["overrepresented_sequences_count"] = overrep_count
                    else:
                        metrics["has_overrepresented_sequences"] = False
                        metrics["overrepresented_sequences_count"] = 0
                
                # Extract adapter content
                if ">>Adapter Content" in content:
                    section = content.split(">>Adapter Content")[1].split(">>")[0]
                    lines = section.strip().split('\n')[1:]  # Skip the header line
                    
                    adapter_percentages = []
                    for line in lines:
                        if '\t' in line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) > 1:
                                try:
                                    # Skip header rows
                                    if any(h.lower() in ['position', 'adapter'] for h in parts):
                                        continue
                                        
                                    # Sum all adapter percentages for this position
                                    adapter_sum = sum(float(p) for p in parts[1:] if p and p != 'NaN' and p.lower() != 'nan')
                                    adapter_percentages.append(adapter_sum)
                                except ValueError:
                                    continue
                    
                    if adapter_percentages:
                        metrics["max_adapter_content"] = max(adapter_percentages)
                        metrics["mean_adapter_content"] = sum(adapter_percentages) / len(adapter_percentages)
    
    except Exception as e:
        logging.error(f"Error processing {zip_file_path}: {str(e)}")
        return {}
    
    logging.info(f"Successfully extracted {len(metrics)} metrics from {zip_file_path}")
    return metrics

# Get Snakemake input and output variables
try:
    input_files = snakemake.input.zip
    output_file = snakemake.output.features
    
    # Make sure input_files is a list
    if not isinstance(input_files, list):
        input_files = [input_files]
    
    logging.info(f"Processing {len(input_files)} FastQC zip files")
    logging.info(f"Input files: {input_files}")
    logging.info(f"Output file: {output_file}")
    
    # Extract features from all zip files
    all_metrics = {}
    
    for zip_file in input_files:
        # Get base filename without path and extension
        base_name = os.path.basename(zip_file).replace('_fastqc.zip', '')
        logging.info(f"Extracting metrics for {base_name}")
        
        # Extract metrics from this file
        metrics = extract_fastqc_data(zip_file)
        
        # Add file identifier to each metric key
        file_metrics = {f"{base_name}_{key}": value for key, value in metrics.items()}
        
        # Add to combined metrics
        all_metrics.update(file_metrics)
    
    # Create a DataFrame and save it
    if all_metrics:
        df = pd.DataFrame([all_metrics])
        
        # Save the data
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        df.to_csv(output_file, index=False)
        logging.info(f"Features extracted and saved to {output_file}")
    else:
        # Create a simple DataFrame with basic info to avoid empty file errors
        first_sample = os.path.basename(input_files[0]).replace('_fastqc.zip', '') if input_files else "unknown"
        df = pd.DataFrame({
            "sample": [first_sample],
            "status": ["No valid metrics extracted"],
            "timestamp": [pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")]
        })
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        df.to_csv(output_file, index=False)
        logging.info(f"No valid metrics found. Basic file created at {output_file}")

except Exception as e:
    logging.error(f"Critical error in script: {str(e)}")
    # Create an empty output file to avoid Snakemake errors
    os.makedirs(os.path.dirname(snakemake.output.features), exist_ok=True)
    with open(snakemake.output.features, 'w') as f:
        f.write("error,message\n")
        f.write(f"script_error,{str(e)}\n")
    sys.exit(1)