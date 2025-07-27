#!/usr/bin/env python3
"""
AlphaFold3 TSV to JSON Converter
Converts a TSV file with bait/prey entries to separate AlphaFold3 JSON files
for each unique bait-prey combination.
"""

import json
import pandas as pd
import argparse
import sys
import requests
import time
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional, Any


class AF3Converter:
    def __init__(self, workdir: str = ".") -> None:
        self.base_structure = {
            "name": "",
            "modelSeeds": [1],
            "dialect": "alphafold3",
            "version": 1,
            "sequences": []
        }
        self.workdir = Path(workdir)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'AF3Converter/1.0 (Python script for AlphaFold3 conversion)'
        })
    
    def read_tsv(self, tsv_file: str | Path) -> pd.DataFrame:
        """Read the input TSV file."""
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            if 'Entry' not in df.columns or 'bait' not in df.columns:
                raise ValueError("TSV must contain 'Entry' and 'bait' columns")
            return df
        except Exception as e:
            print(f"Error reading TSV file: {e}")
            sys.exit(1)
    
    def fetch_uniprot_sequence(self, uniprot_id: str) -> Optional[str]:
        """Fetch sequence from UniProt."""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        
        try:
            response = self.session.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                sequence = ''.join(lines[1:])  # Skip header line
                return sequence
            else:
                print(f"Warning: Could not fetch UniProt sequence for {uniprot_id}")
                return None
        except Exception as e:
            print(f"Error fetching UniProt sequence for {uniprot_id}: {e}")
            return None
    
    def read_fasta(self, fasta_file: str | Path) -> str:
        """Read a FASTA file and return the sequence. Ensures only 1 entry per file."""
        # Build full path using workdir if fasta_file is relative
        fasta_path = Path(fasta_file)
        if not fasta_path.is_absolute():
            fasta_path = self.workdir / fasta_path
        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()
            
            header_count = 0
            sequence = ""
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    header_count += 1
                    if header_count > 1:
                        raise ValueError(f"FASTA file {fasta_path} contains multiple entries. Only single-entry FASTA files are supported.")
                elif line:  # Non-empty sequence line
                    sequence += line
            
            if header_count == 0:
                raise ValueError(f"FASTA file {fasta_path} contains no header lines.")
            
            if not sequence:
                raise ValueError(f"FASTA file {fasta_path} contains no sequence data.")
            
            return sequence
        except FileNotFoundError:
            raise FileNotFoundError(f"FASTA file {fasta_path} not found.")
        except Exception as e:
            raise RuntimeError(f"Error reading FASTA file {fasta_path}: {e}")
    
    def get_entry_type(self, entry: str) -> str:
        """Determine the entry type."""
        # Check for FASTA files
        if entry.endswith('.fasta') or entry.endswith('.fa'):
            return 'fasta_file'
        
        # Check for CCD codes (3-letter codes, sometimes with numbers)
        if re.match(r'^CCD:', entry):
            return 'ccd_code'
        
        # Check for SMILES (contains chemical notation characters)
        if re.match(r'^SMILES:', entry):
            return 'smiles'
        
        # Default to uniprot
        return 'uniprot'
    
    def process_entry(self, entry: str, sequence_id: str) -> Optional[Dict[str, Any]]:
        """Process a single entry and return the appropriate sequence object."""
        entry_type = self.get_entry_type(entry)
        
        if entry_type == 'ccd_code':
            san_entry = entry.removeprefix('CCD:')            
            return {
                "ligand": {
                    "id": sequence_id,
                    "ccdCodes": [san_entry]
                }
            }
        
        elif entry_type == 'smiles':
            san_entry = entry.removeprefix('SMILES:')
            return {
                "ligand": {
                    "id": sequence_id,
                    "smiles": san_entry
                }
            }
        
        elif entry_type == 'fasta_file':
            sequence = self.read_fasta(entry)
            if sequence:
                if self.is_dna_sequence(sequence):
                    return {
                        "dna": {
                            "id": sequence_id,
                            "sequence": sequence
                        }
                    }
                elif self.is_rna_sequence(sequence):
                    return {
                        "rna": {
                            "id": sequence_id,
                            "sequence": sequence
                        }
                    }
                else:
                    return {
                        "protein": {
                            "id": sequence_id,
                            "sequence": sequence
                        }
                    }
            else:
                print(f"Warning: Could not read FASTA file {entry}")
                return None
        
        elif entry_type == 'uniprot':
            sequence = self.fetch_uniprot_sequence(entry)
            if sequence:
                return {
                    "protein": {
                        "id": sequence_id,
                        "sequence": sequence
                    }
                }
            else:
                print(f"Warning: Could not fetch sequence for UniProt ID {entry}")
                return None
        
        else:
            print(f"Warning: Unknown entry type for {entry}")
            return None
    
    def is_dna_sequence(self, sequence: str) -> bool:
        """Check if sequence is DNA."""
        dna_chars = set('ATCGN')
        sequence_chars = set(sequence.upper())
        return len(sequence_chars - dna_chars) == 0 and 'U' not in sequence.upper()
    
    def is_rna_sequence(self, sequence: str) -> bool:
        """Check if sequence is RNA."""
        rna_chars = set('AUCGN')
        sequence_chars = set(sequence.upper())
        return len(sequence_chars - rna_chars) == 0 and 'T' not in sequence.upper()
    
    def generate_combinations(self, df: pd.DataFrame) -> List[Tuple[str, str]]:
        """Generate all unique bait-prey combinations."""
        baits = df[df['bait'] == 1]['Entry'].tolist()
        preys = df[df['bait'] == 0]['Entry'].tolist()
        
        combinations_list = []
        
        # Generate all bait-prey pairs
        for bait in baits:
            for prey in preys:
                combinations_list.append((bait, prey))
        
        return combinations_list
    
    def create_json_for_combination(self, bait_entry: str, prey_entry: str, output_dir: str | Path) -> Optional[Path]:
        """Create a JSON file for a specific bait-prey combination."""
        sequences = []
        sequence_id = 'A'
        
        # Process bait
        bait_seq = self.process_entry(bait_entry, sequence_id)
        if bait_seq:
            sequences.append(bait_seq)
            sequence_id = chr(ord(sequence_id) + 1)
        
        # Process prey
        prey_seq = self.process_entry(prey_entry, sequence_id)
        if prey_seq:
            sequences.append(prey_seq)
        
        if not sequences:
            print(f"Warning: No valid sequences for combination {bait_entry}-{prey_entry}")
            return None
        
        # Create structure name using basename for FASTA files
        bait_name = Path(bait_entry).stem if bait_entry.endswith(('.fasta', '.fa')) else bait_entry
        prey_name = Path(prey_entry).stem if prey_entry.endswith(('.fasta', '.fa')) else prey_entry
        safe_bait = re.sub(r'[^\w\-_.]', '_', bait_name)
        safe_prey = re.sub(r'[^\w\-_.]', '_', prey_name)
        safe_name = f"{safe_bait}_{safe_prey}"
        
        # Create structure
        structure = self.base_structure.copy()
        structure["name"] = safe_name
        structure["sequences"] = sequences
        
        # Create filename using basename for FASTA files
        filename = f"{safe_name}.json"
        filepath = Path(output_dir) / filename
        
        # Write file
        try:
            with open(filepath, 'w') as f:
                json.dump(structure, f, indent=2)
            print(f"Created: {filepath}")
            return filepath
        except Exception as e:
            print(f"Error writing file {filepath}: {e}")
            return None
    
    def convert(self, tsv_file: str | Path, output_dir: str = "output", workdir: str = ".") -> List[Path] | ValueError:
        """Convert TSV to multiple AlphaFold3 JSON files."""
        df = self.read_tsv(tsv_file)
        
        # Create output directory
        Path(output_dir).mkdir(exist_ok=True)
        
        # Validate data
        if df['bait'].sum() == 0:
            raise ValueError("No bait entries found (bait=1)")
        
        if (df['bait'] == 0).sum() == 0:
            raise ValueError("No prey entries found (bait=0)")
        
        # Generate combinations
        combinations = self.generate_combinations(df)
        
        print(f"Found {len(combinations)} bait-prey combinations")
        
        # Process each combination
        created_files = []
        for i, (bait, prey) in enumerate(combinations, 1):
            print(f"Processing combination {i}/{len(combinations)}: {bait} (bait) + {prey} (prey)")
            
            filepath = self.create_json_for_combination(bait, prey, output_dir)
            if filepath:
                created_files.append(filepath)
            
            # Add small delay to be respectful to UniProt API
            if i % 10 == 0:
                time.sleep(1)
        
        print(f"\nCompleted! Created {len(created_files)} JSON files in '{output_dir}' directory")
        return created_files


def main() -> None:
    parser = argparse.ArgumentParser(description='Convert TSV to AlphaFold3 JSON format')
    parser.add_argument('input_tsv', help='Input TSV file')
    parser.add_argument('-o', '--output-dir', default='output', 
                        help='Output directory for JSON files (default: output)')
    parser.add_argument('--workdir', default='.',
                        help='Work directory for relative paths (default: .)')
    
    args = parser.parse_args()
    
    if not Path(args.input_tsv).exists():
        print(f"Error: Input file {args.input_tsv} does not exist")
        sys.exit(1)
    
    converter = AF3Converter(args.workdir)
    converter.convert(args.input_tsv, args.output_dir)


if __name__ == "__main__":
    main()