#!/usr/bin/env python3
"""
Converts a TSV file with bait/prey entries to separate AlphaFold3 JSON files
for each unique bait-prey combination, or ColabFold FASTA files, or Boltz FASTA files.
"""

import json
import pandas as pd
import argparse
import sys
import requests
import time
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional, Any, Union


class TSV2AFConverter:
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
    
    def read_tsv(self, tsv_file: Union[str, Path]) -> pd.DataFrame:
        """Read the input TSV file."""
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            df.columns = df.columns.str.lower()
            if 'entry' not in df.columns or 'bait' not in df.columns:
                raise ValueError("TSV must contain 'Entry' and 'Bait' columns")
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
    
    def read_fasta(self, fasta_file: Union[str, Path]) -> Union[str, Dict[str, str]]:
        """Read a FASTA file and return the sequence(s). 
        
        Args:
            fasta_file: Path to FASTA file
        """
        # Build full path using workdir if fasta_file is relative
        fasta_path = Path(fasta_file)
        if not fasta_path.is_absolute():
            fasta_path = self.workdir / fasta_path
        
        try:
            with open(fasta_path, 'r') as f:
                lines = f.readlines()
            
            sequences = {}
            current_header = None
            current_sequence = ""
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header is not None:
                        sequences[current_header] = current_sequence
                    
                    # Start new sequence
                    current_header = line.removeprefix('>').split(" ")[0]  # Remove '>' prefix and only take first part of the name
                    current_sequence = ""
                elif line:  # Non-empty sequence line
                    current_sequence += line
            
            # Save last sequence
            if current_header is not None:
                sequences[current_header] = current_sequence
            
            if len(sequences) == 0:
                raise ValueError(f"FASTA file {fasta_path} contains no entries.")
            
            # Check for empty sequences
            for header, seq in sequences.items():
                if not seq:
                    raise ValueError(f"FASTA file {fasta_path} contains empty sequence for entry '{header}'.")
            
            return sequences
                
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
    
    def validate_colabfold_entry(self, entry: str) -> bool:
        """Validate that entry is suitable for ColabFold mode (UniProt ID or FASTA file only)."""
        entry_type = self.get_entry_type(entry)
        return entry_type in ['uniprot', 'fasta_file']
    
    def get_protein_sequence(self, entry: str) -> Union[str, Dict[str, str]]:
        """Get protein sequence(s) from entry.
        
        Args:
            entry: Entry identifier (UniProt ID or FASTA file path)
        
        Returns:
            For colabfold: dict {header: sequence} for FASTA files, string for UniProt
        """
        entry_type = self.get_entry_type(entry)
        
        if entry_type == 'uniprot':
            sequence = self.fetch_uniprot_sequence(entry)
            if sequence is None:
                raise RuntimeError(f"Could not fetch UniProt sequence for {entry}")
            return sequence
        elif entry_type == 'fasta_file':
            sequences = self.read_fasta(entry)
            
            if isinstance(sequences, dict):
                # Multiple sequences - validate all are proteins
                for header, seq in sequences.items():
                    if self.is_dna_sequence(seq):
                        raise ValueError(f"FASTA file {entry} entry '{header}' contains DNA sequence, not protein")
                    if self.is_rna_sequence(seq):
                        raise ValueError(f"FASTA file {entry} entry '{header}' contains RNA sequence, not protein")
                return sequences
            else:
                # Single sequence - validate it's protein
                if self.is_dna_sequence(sequences):
                    raise ValueError(f"FASTA file {entry} contains DNA sequence, not protein")
                if self.is_rna_sequence(sequences):
                    raise ValueError(f"FASTA file {entry} contains RNA sequence, not protein")
                return sequences
        else:
            raise ValueError(f"Entry type {entry_type} not supported in ColabFold mode")
    
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
            try:
                # Allow multiple entries in both modes now
                sequences = self.read_fasta(entry)
                
                if isinstance(sequences, str):
                    # Single sequence
                    if self.is_dna_sequence(sequences):
                        return {
                            "dna": {
                                "id": sequence_id,
                                "sequence": sequences
                            }
                        }
                    elif self.is_rna_sequence(sequences):
                        return {
                            "rna": {
                                "id": sequence_id,
                                "sequence": sequences
                            }
                        }
                    else:
                        return {
                            "protein": {
                                "id": sequence_id,
                                "sequence": sequences
                            }
                        }
                else:
                    # Multiple sequences - return as list for processing
                    sequence_objects = []
                    for i, (header, seq) in enumerate(sequences.items()):
                        current_id = chr(ord(sequence_id) + i)
                        if self.is_dna_sequence(seq):
                            sequence_objects.append({
                                "dna": {
                                    "id": current_id,
                                    "sequence": seq
                                }
                            })
                        elif self.is_rna_sequence(seq):
                            sequence_objects.append({
                                "rna": {
                                    "id": current_id,
                                    "sequence": seq
                                }
                            })
                        else:
                            sequence_objects.append({
                                "protein": {
                                    "id": current_id,
                                    "sequence": seq
                                }
                            })
                    return {"multiple": sequence_objects}  # Special marker for multiple sequences
                    
            except Exception as e:
                raise RuntimeError(f"Could not read FASTA file {entry}: {e}")
        
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
                raise RuntimeError(f"Could not fetch sequence for UniProt ID {entry}")
        
        else:
            raise ValueError(f"Unknown entry type for {entry}")
    
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
        baits = df[df['bait'] == 1]['entry'].tolist()
        preys = df[df['bait'] == 0]['entry'].tolist()
        
        combinations_list = []
        
        # Generate all bait-prey pairs
        for bait in baits:
            for prey in preys:
                combinations_list.append((bait, prey))
        
        return combinations_list
    
    def get_entry_name(self, entry: str) -> str:
        """Get the name to use for the entry in output filenames."""
        entry_type = self.get_entry_type(entry)
        
        if entry_type == 'fasta_file':
            # Use the sequence name from the FASTA file
            try:
                sequences = self.read_fasta(entry)
                if isinstance(sequences, dict):
                    # If multiple sequences, use the first one's name
                    first_header = next(iter(sequences.keys()))
                    return first_header
                else:
                    # This shouldn't happen with but handle it
                    return Path(entry).stem
            except Exception:
                # Fallback to filename if we can't read the FASTA
                return Path(entry).stem
        else:
            # For non-FASTA entries, use the entry itself (e.g., UniProt ID)
            return entry
    def get_entry_name_for_sequence(self, entry: str, sequence_index: int = 0) -> str:
        """Get the name to use for a specific sequence in a FASTA entry."""
        entry_type = self.get_entry_type(entry)

        if entry_type == 'fasta_file':
            try:
                sequences = self.read_fasta(entry)
                if isinstance(sequences, dict):
                    # Get the sequence name at the specified index
                    headers = list(sequences.keys())
                    if sequence_index < len(headers):
                        return headers[sequence_index]
                    else:
                        # Fallback if index is out of range
                        return headers[0] if headers else Path(entry).stem
                else:
                    # Single sequence
                    return Path(entry).stem
            except Exception:
                # Fallback to filename if we can't read the FASTA
                return Path(entry).stem
        else:
            # For non-FASTA entries, use the entry itself (e.g., UniProt ID)
            return entry
    
    def create_json_for_combination(self, bait_entry: str, prey_entry: str, output_dir: Union[str, Path]) -> List[Path]:
        """Create JSON file(s) for a specific bait-prey combination."""
        created_files = []

        # Process bait with temporary sequence ID
        bait_seq = self.process_entry(bait_entry, 'A')
        if bait_seq:
            if "multiple" in bait_seq:
                bait_sequences = bait_seq["multiple"]
            else:
                bait_sequences = [bait_seq]
        else:
            bait_sequences = []

        # Process prey with temporary sequence ID
        prey_seq = self.process_entry(prey_entry, 'B')
        if prey_seq:
            if "multiple" in prey_seq:
                prey_sequences = prey_seq["multiple"]
            else:
                prey_sequences = [prey_seq]
        else:
            prey_sequences = []

        if not bait_sequences or not prey_sequences:
            raise RuntimeError(f"No valid sequences for combination {bait_entry}-{prey_entry}")

        # Generate all combinations for multi-entry FASTA files
        for i, bait_seq_obj in enumerate(bait_sequences):
            for j, prey_seq_obj in enumerate(prey_sequences):
                # Reset sequence IDs for each JSON file: bait = "A", prey = "B"
                bait_seq_corrected = bait_seq_obj.copy()
                prey_seq_corrected = prey_seq_obj.copy()

                # Update the sequence IDs to ensure bait is "A" and prey is "B"
                for seq_type in ['protein', 'dna', 'rna', 'ligand']:
                    if seq_type in bait_seq_corrected:
                        bait_seq_corrected[seq_type]['id'] = 'A'
                    if seq_type in prey_seq_corrected:
                        prey_seq_corrected[seq_type]['id'] = 'B'

                sequences = [bait_seq_corrected, prey_seq_corrected]

                # Create structure name using specific sequence names for FASTA files
                bait_name = self.get_entry_name_for_sequence(bait_entry, i)
                prey_name = self.get_entry_name_for_sequence(prey_entry, j)

                # Create structure
                structure = self.base_structure.copy()
                structure["name"] = f"{bait_name}_{prey_name}"
                structure["sequences"] = sequences

                # Create filename using sequence name for FASTA files
                safe_bait = re.sub(r'[^\w\-_.]', '_', bait_name)
                safe_prey = re.sub(r'[^\w\-_.]', '_', prey_name)
                filename = f"{safe_bait}_{safe_prey}.json"
                filepath = Path(output_dir) / filename

                # Write file
                try:
                    with open(filepath, 'w') as f:
                        json.dump(structure, f, indent=2)
                    print(f"Created: {filepath}")
                    created_files.append(filepath)
                except Exception as e:
                    raise RuntimeError(f"Error writing file {filepath}: {e}")

        return created_files
    
    def create_fasta_for_colab_combination(self, bait_entry: str, prey_entry: str, output_dir: Union[str, Path]) -> List[Path]:
        """Create FASTA file(s) for a specific bait-prey combination (ColabFold mode)."""
        try:
            # Get sequences
            bait_seqs = self.get_protein_sequence(bait_entry)
            prey_seqs = self.get_protein_sequence(prey_entry)
            
            created_files = []
            
            # Convert single sequences to dict format for uniform handling
            if isinstance(bait_seqs, str):
                bait_name = self.get_entry_name(bait_entry)
                bait_seqs = {bait_name: bait_seqs}
            
            if isinstance(prey_seqs, str):
                prey_name = self.get_entry_name(prey_entry)
                prey_seqs = {prey_name: prey_seqs}
            
            # Generate all combinations
            for bait_header, bait_seq in bait_seqs.items():
                for prey_header, prey_seq in prey_seqs.items():
                    # Create concatenated sequence with : separator
                    combined_sequence = f"{bait_seq}:{prey_seq}"
                    
                    # Create filename using sequence names
                    safe_bait = re.sub(r'[^\w\-_.]', '_', bait_header)
                    safe_prey = re.sub(r'[^\w\-_.]', '_', prey_header)
                    filename = f"{safe_bait}_{safe_prey}.fasta"
                    filepath = Path(output_dir) / filename
                    
                    # Write FASTA file
                    with open(filepath, 'w') as f:
                        f.write(f">{bait_header}_{prey_header}\n")
                        # Write sequence with line breaks every 80 characters
                        for i in range(0, len(combined_sequence), 80):
                            f.write(combined_sequence[i:i+80] + '\n')
                    
                    print(f"Created: {filepath}")
                    created_files.append(filepath)
            
            return created_files
            
        except Exception as e:
            raise RuntimeError(f"Error creating FASTA file for combination {bait_entry}-{prey_entry}: {e}")
    
    def create_fasta_for_boltz_combination(self, bait_entry: str, prey_entry: str, output_dir: Union[str, Path]) -> List[Path]:
        """Create FASTA file(s) for a specific bait-prey combination (Boltz mode)."""
        created_files = []

        # Process bait with temporary sequence ID
        bait_seq = self.process_entry(bait_entry, 'A')
        if bait_seq:
            if "multiple" in bait_seq:
                bait_sequences = bait_seq["multiple"]
            else:
                bait_sequences = [bait_seq]
        else:
            bait_sequences = []

        # Process prey with temporary sequence ID
        prey_seq = self.process_entry(prey_entry, 'B')
        if prey_seq:
            if "multiple" in prey_seq:
                prey_sequences = prey_seq["multiple"]
            else:
                prey_sequences = [prey_seq]
        else:
            prey_sequences = []

        if not bait_sequences or not prey_sequences:
            raise RuntimeError(f"No valid sequences for combination {bait_entry}-{prey_entry}")

        # Generate all combinations for multi-entry FASTA files
        for i, bait_seq_obj in enumerate(bait_sequences):
            for j, prey_seq_obj in enumerate(prey_sequences):
                # Start with sequence ID 'A' and increment
                current_id = 'A'
                fasta_content = []

                # Process bait sequence
                bait_seq_corrected = bait_seq_obj.copy()
                for seq_type in ['protein', 'dna', 'rna', 'ligand']:
                    if seq_type in bait_seq_corrected:
                        bait_seq_corrected[seq_type]['id'] = current_id

                        if seq_type == 'ligand':
                            # Handle ligand entries
                            ligand_data = bait_seq_corrected[seq_type]
                            if 'ccdCodes' in ligand_data:
                                fasta_content.append(f">{current_id}|ccd")
                                fasta_content.append(ligand_data['ccdCodes'][0])
                            elif 'smiles' in ligand_data:
                                fasta_content.append(f">{current_id}|smiles")
                                fasta_content.append(ligand_data['smiles'])
                        else:
                            # Handle protein/DNA/RNA sequences
                            sequence = bait_seq_corrected[seq_type]['sequence']
                            fasta_content.append(f">{current_id}|{seq_type}")
                            # Write sequence with line breaks every 80 characters
                            for k in range(0, len(sequence), 80):
                                fasta_content.append(sequence[k:k+80])

                        current_id = chr(ord(current_id) + 1)
                        break  # Only process one sequence type per entry
                    
                # Process prey sequence
                prey_seq_corrected = prey_seq_obj.copy()
                for seq_type in ['protein', 'dna', 'rna', 'ligand']:
                    if seq_type in prey_seq_corrected:
                        prey_seq_corrected[seq_type]['id'] = current_id

                        if seq_type == 'ligand':
                            # Handle ligand entries
                            ligand_data = prey_seq_corrected[seq_type]
                            if 'ccdCodes' in ligand_data:
                                fasta_content.append(f">{current_id}|ccd")
                                fasta_content.append(ligand_data['ccdCodes'][0])
                            elif 'smiles' in ligand_data:
                                fasta_content.append(f">{current_id}|smiles")
                                fasta_content.append(ligand_data['smiles'])
                        else:
                            # Handle protein/DNA/RNA sequences
                            sequence = prey_seq_corrected[seq_type]['sequence']
                            fasta_content.append(f">{current_id}|{seq_type}")
                            # Write sequence with line breaks every 80 characters
                            for k in range(0, len(sequence), 80):
                                fasta_content.append(sequence[k:k+80])

                        break  # Only process one sequence type per entry
                    
                # Create filename using specific sequence names for FASTA files
                bait_name = self.get_entry_name_for_sequence(bait_entry, i)
                prey_name = self.get_entry_name_for_sequence(prey_entry, j)

                safe_bait = re.sub(r'[^\w\-_.]', '_', bait_name)
                safe_prey = re.sub(r'[^\w\-_.]', '_', prey_name)
                filename = f"{safe_bait}_{safe_prey}.fasta"
                filepath = Path(output_dir) / filename

                # Write FASTA file
                try:
                    with open(filepath, 'w') as f:
                        f.write('\n'.join(fasta_content) + '\n')

                    print(f"Created: {filepath}")
                    created_files.append(filepath)
                except Exception as e:
                    raise RuntimeError(f"Error writing file {filepath}: {e}")

        return created_files
    
    def convert(self, tsv_file: Union[str, Path], output_dir: str = "output", mode: str = "alphafold3") -> List[Path]:
        """Convert TSV to multiple AlphaFold3 JSON files or ColabFold FASTA files."""
        df = self.read_tsv(tsv_file)
        
        # Create output directory
        Path(output_dir).mkdir(exist_ok=True)
        
        # Validate data
        if df['bait'].sum() == 0:
            raise ValueError("No bait entries found (bait=1)")
        
        if (df['bait'] == 0).sum() == 0:
            raise ValueError("No prey entries found (bait=0)")
        
        # Validate entries for ColabFold mode
        if mode == "colabfold":
            invalid_entries = []
            for entry in df['entry']:
                if not self.validate_colabfold_entry(entry):
                    invalid_entries.append(entry)
            
            if invalid_entries:
                raise ValueError(f"ColabFold mode only supports UniProt IDs and FASTA files with protein sequences. "
                                f"Invalid entries: {invalid_entries}")
        
        # Generate combinations
        combinations = self.generate_combinations(df)
        
        print(f"Found {len(combinations)} bait-prey combinations")
        print(f"Mode: {mode}")
        
        # Process each combination
        created_files = []
        for i, (bait, prey) in enumerate(combinations, 1):
            print(f"Processing combination {i}/{len(combinations)}: {bait} (bait) + {prey} (prey)")
            
            if mode == "alphafold3":
                filepaths = self.create_json_for_combination(bait, prey, output_dir)
                created_files.extend(filepaths)
            elif mode == "colabfold":
                filepaths = self.create_fasta_for_colab_combination(bait, prey, output_dir)
                created_files.extend(filepaths)
            elif mode == "boltz":
                filepaths = self.create_fasta_for_boltz_combination(bait, prey, output_dir)
                created_files.extend(filepaths)
            else:
                raise ValueError(f"Unknown mode '{mode}'. Supported modes: alphafold3, colabfold")
            
            # Add small delay to be respectful to UniProt API
            if i % 10 == 0:
                time.sleep(1)
        
        file_type = "JSON" if mode == "alphafold3" else "FASTA"
        print(f"\nCompleted! Created {len(created_files)} {file_type} files in '{output_dir}' directory")
        return created_files


def main() -> None:
    parser = argparse.ArgumentParser(description='Convert TSV to AlphaFold3 JSON format or ColabFold FASTA format')
    parser.add_argument('input_tsv', help='Input TSV file')
    parser.add_argument('-o', '--output-dir', default='output', 
                        help='Output directory for JSON/FASTA files (default: output)')
    parser.add_argument('--workdir', default='.',
                        help='Work directory for relative paths (default: .)')
    parser.add_argument('--mode', choices=['alphafold3', 'colabfold', 'boltz'], default='alphafold3',
                        help='Output mode: alphafold3 (JSON files) or colabfold (FASTA files) or boltz (FASTA files) (default: alphafold3)')
    
    args = parser.parse_args()
    
    if not Path(args.input_tsv).exists():
        print(f"Error: Input file {args.input_tsv} does not exist")
        sys.exit(1)
    
    converter = TSV2AFConverter(args.workdir)
    converter.convert(args.input_tsv, args.output_dir, args.mode)


if __name__ == "__main__":
    main()