#!/usr/bin/env python3
"""
Script to read JSON files and convert to TSV format.
Reads all *.json files in the current directory and outputs a single TSV file
with columns: foldid, fraction_disordered, has_clash, ptm, iptm, ranking_score
Results are sorted by ranking_score in descending order.
"""

import json
import glob
import os
import sys
from pathlib import Path

def process_json_files(input_dir=".", output_file="results.tsv"):
    """
    Process all JSON files in the specified directory and create a TSV file.
    
    Args:
        input_dir (str): Directory to search for JSON files (default: current directory)
        output_file (str): Output TSV filename (default: results.tsv)
    """
    # Find all JSON files
    json_pattern = os.path.join(input_dir, "*.json")
    json_files = glob.glob(json_pattern)
    
    if not json_files:
        print(f"No JSON files found in {input_dir}")
        return
    
    print(f"Found {len(json_files)} JSON files")
    
    # Store data for all files
    data_rows = []
    
    # Process each JSON file
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Extract basename without extension for foldid
            foldid = Path(json_file).stem
            
            # Extract required fields
            row = {
                'foldid': foldid,
                'fraction_disordered': data.get('fraction_disordered', ''),
                'has_clash': data.get('has_clash', ''),
                'ptm': data.get('ptm', ''),
                'iptm': data.get('iptm', ''),
                'ranking_score': data.get('ranking_score', '')
            }
            
            data_rows.append(row)
            
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error processing {json_file}: {e}")
            continue
    
    if not data_rows:
        print("No valid JSON files processed")
        return
    
    # Sort by ranking_score in descending order
    # Handle cases where ranking_score might be missing or non-numeric
    def safe_sort_key(row):
        score = row['ranking_score']
        if isinstance(score, (int, float)):
            return score
        return -float('inf')  # Put invalid scores at the end
    
    data_rows.sort(key=safe_sort_key, reverse=True)
    
    # Write TSV file
    headers = ['foldid', 'fraction_disordered', 'has_clash', 'ptm', 'iptm', 'ranking_score']
    
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write('\t'.join(headers) + '\n')
            
            # Write data rows
            for row in data_rows:
                values = [str(row[header]) for header in headers]
                f.write('\t'.join(values) + '\n')
        
        print(f"Successfully wrote {len(data_rows)} rows to {output_file}")
        print("Results sorted by ranking_score (highest to lowest)")
        
    except IOError as e:
        print(f"Error writing output file {output_file}: {e}")

def main():
    """Main function with command line argument support"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert JSON files to TSV format')
    parser.add_argument('--input-dir', '-i', default='.', 
                        help='Input directory containing JSON files (default: current directory)')
    parser.add_argument('--output', '-o', default='results.tsv',
                        help='Output TSV filename (default: results.tsv)')
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    
    process_json_files(args.input_dir, args.output)

if __name__ == "__main__":
    main()