#!/usr/bin/env python3
"""
Script to read AF3 summary confidence scores from JSON files and output a single ranked TSV file.
Reads all *.json files in the current directory and outputs a single TSV file
with columns: foldid, fraction_disordered, has_clash, ptm, iptm, ranking_score
Results are sorted by ranking_score in descending order.
"""

import json
import glob
import os
import sys
from pathlib import Path

def process_json_files(input_dir: str, output_file: str, mode: str):
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
            
            if mode == "colabfold":
                # Extract basename without extension and the suffix for foldid
                foldid = Path(json_file).stem.removesuffix("_toprank")

                # Extract required fields
                row = {
                    'foldid': foldid,
                    'iptm': data.get('iptm', '')
                }
                headers = ['foldid', 'iptm']
            elif mode == "alphafold3":
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
                headers = ['foldid', 'fraction_disordered', 'has_clash', 'ptm', 'iptm', 'ranking_score']
            
            data_rows.append(row)
            
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error processing {json_file}: {e}")
            continue
    
    if not data_rows:
        print("No valid JSON files processed")
        return
    
    # Sort by the last entry in the dict
    # Handle cases where sorting key is missing or non-numeric
    def safe_sort_key(row):
        last_key = list(row.keys())[-1]
        score = row[last_key]
        if isinstance(score, (int, float)):
            return score
        return -float('inf')  # Put invalid scores at the end
    
    data_rows.sort(key=safe_sort_key, reverse=True)
    
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
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert JSON files to TSV format')
    parser.add_argument('--input-dir', '-i', default='.', 
                        help='Input directory containing JSON files (default: current directory)')
    parser.add_argument('--output', '-o', default='results.tsv',
                        help='Output TSV filename (default: results.tsv)')
    parser.add_argument('--mode',
                        help='Set the input format. Options: colabfold, alphafold3')
    
    args = parser.parse_args()
    
    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    
    process_json_files(args.input_dir, args.output, args.mode)

if __name__ == "__main__":
    main()