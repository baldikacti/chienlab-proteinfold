#!/usr/bin/env python3

"""
Script to generate pools of proteins, including the bait protein, to pass into AF3 for folding prediction.
This script generates pools such that each pool has fewer amino acids than the maximum pool depth (specified in config files),
each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in each of its pools.
Written by Jonathan Coombs
"""

import random
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union

# Function to parse FASTA files.
# Input: fasta_file: a FASTA file with one or more sequences.
# Output: a dictionary with protein name as key and sequence as value.
# Protein name is treated as the first word in a name line, represented by a line starting with '>'.
# Raises exception if: FASTA file is not found; FASTA file is empty.
def read_fasta(fasta_file: Union[str, Path]) -> Union[str, Dict[str, str]]:
        """Read a FASTA file and return the sequence(s). 
        
        Args:
            fasta_file: Path to FASTA file
        """
        # Build full path using workdir if fasta_file is relative
        fasta_path = Path(fasta_file)
        if not fasta_path.is_absolute():
            fasta_path = "./" / fasta_path
        
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

# Function to generate a list of protein pools, following the rules specified in this file.
# Inputs:   bait_protein: a dictionary with one entry, the bait protein name and sequence.
#           test_proteins: a dictionary with one or more entries, the test protein names and sequences.
#           max_pool_depth: the maximum number of amino acids in each pool.
# Output: a list of dictionaries, each dictionary representing a pool of proteins, with protein names as keys and sequences as values.
# Raises exception if: bait_protein is not size 1; test_proteins is empty; max_pool_depth is not a positive integer;
# any protein is too long to be in a pool with the bait protein.
def generate_pools(bait_protein: Dict[str, str], test_proteins: Dict[str, str], max_pool_depth: int) -> List[Dict[str, str]]:
        """Generate pools of proteins for AF3 pooled folding prediction.
        
        Args:
            bait_protein: Dictionary with one entry, the bait protein name and sequence.
            test_proteins: Dictionary with one or more entries, the test protein names and sequences.
            max_pool_depth: Maximum number of amino acids in each pool.
        """
        if len(bait_protein) != 1:
            raise ValueError("bait_protein must contain exactly one entry.")
        if not test_proteins:
            raise ValueError("test_proteins must contain at least one entry.")
        if not isinstance(max_pool_depth, int) or max_pool_depth <= 0:
            raise ValueError("max_pool_depth must be a positive integer.")
        if any(len(seq) > max_pool_depth - len(bait_protein[list(bait_protein.keys())[0]]) for seq in test_proteins.values()):
            raise ValueError("One or more test proteins are too long to be in a pool with the bait protein.")

        # Initialize pools
        pools = []
        current_pool = {}
        current_pool_size = 0
        
        # Add bait protein to every pool
        bait_name, bait_sequence = next(iter(bait_protein.items()))

        # Initial approach: greedy algo with shuffled array at each step.

        # Data structures:
        #   - id_lookup: str array lookup table for protein ID
        #   - length_lookup: int array lookup table for protein length
        #   - same_pool: N x N boolean array for protein ID x protein ID, indicating whether a given protein is in the same pool as another protein
        #   - unused_prots: int array containing proteins NOT USED YET
        #   - once_used: int array containing proteins used once

        id_lookup = list(test_proteins.keys())
        length_lookup = [len(test_proteins[prot]) for prot in id_lookup]
        same_pool = [[False for _ in range(len(id_lookup))] for _ in range(len(id_lookup))]
        unused_prots = list(range(len(id_lookup)))
        once_used = []

        # First loop: until all proteins have been used once
        while unused_prots:
            # Start a new pool with the bait protein
            current_pool = {bait_name: bait_sequence}
            current_pool_size = len(bait_sequence)

            # Shuffle the unused_prots list to randomize selection
            random.shuffle(unused_prots)

            # Create temporary lists to hold proteins we've looked at
            tmp = []
            tmp2 = []

            # Step 1: iterate over all unused proteins and try to add them to the current pool
            for _ in range(len(unused_prots)-1, -1, -1):  # Iterate backwards to safely remove elements
                prot_index = unused_prots.pop()
                prot_length = length_lookup[prot_index]

                # Check if adding this protein would exceed the max pool depth
                if current_pool_size + prot_length > max_pool_depth:
                    tmp.append(prot_index)
                    continue  # Skip this protein since it would exceed the max pool depth

                # Check if this protein has been in a pool with any of the proteins already in the current pool
                overlap = False
                for existing_prot in current_pool.keys():
                    # Continue if it's the bait protein, since it is in every pool
                    if existing_prot == bait_name:
                        continue
                    if same_pool[prot_index][id_lookup.index(existing_prot)]:
                        overlap = True
                        break

                if overlap:
                    tmp.append(prot_index)
                    continue  # Skip this protein since it has been in a pool with an existing protein

                # Add the protein to the current pool
                current_pool_size += prot_length
                prot_name = id_lookup[prot_index]
                current_pool[prot_name] = test_proteins[prot_name]

                # Update the same_pool matrix to indicate that these proteins have been in a pool together
                for existing_prot in current_pool.keys():
                    # Continue if it's the bait protein, since it is in every pool
                    if existing_prot == bait_name:
                        continue
                    existing_index = id_lookup.index(existing_prot)
                    same_pool[prot_index][existing_index] = True
                    same_pool[existing_index][prot_index] = True

                # Add this protein to once_used
                once_used.append(prot_index)

            # Step 2: Check if any proteins added once can fit into the current pool
            random.shuffle(once_used)

            # Iterate over all once_used proteins and try to add them to the current pool
            for _ in range(len(once_used)-1, -1, -1):  # Iterate backwards to safely remove elements
                prot_index = once_used.pop()
                prot_length = length_lookup[prot_index]

                # Check if adding this protein would exceed the max pool depth
                if current_pool_size + prot_length > max_pool_depth:
                    tmp2.append(prot_index)
                    continue  # Skip this protein since it would exceed the max pool depth

                # Check if this protein has been in a pool with any of the proteins already in the current pool
                overlap = False
                for existing_prot in current_pool.keys():
                    # Continue if it's the bait protein, since it is in every pool
                    if existing_prot == bait_name:
                        continue
                    if same_pool[prot_index][id_lookup.index(existing_prot)]:
                        overlap = True
                        break

                if overlap:
                    tmp2.append(prot_index)
                    continue  # Skip this protein since it has been in a pool with an existing protein

                # Add the protein to the current pool
                current_pool_size += prot_length
                prot_name = id_lookup[prot_index]
                current_pool[prot_name] = test_proteins[prot_name]

                # Update the same_pool matrix to indicate that these proteins have been in a pool together
                for existing_prot in current_pool.keys():
                    # Continue if it's the bait protein, since it is in every pool
                    if existing_prot == bait_name:
                        continue
                    existing_index = id_lookup.index(existing_prot)
                    same_pool[prot_index][existing_index] = True
                    same_pool[existing_index][prot_index] = True

                # Once a protein is used twice, we remove it from consideration for a third pool.
                # May want to change this later if some of the last few pools have low depth.        

            # Add the completed pool to the list of pools
            pools.append(current_pool)

            # Add back any proteins that were not added to the current pool
            unused_prots = tmp
            once_used = tmp2

        # Second loop: until all proteins have been used twice
        while once_used:
            # Start a new pool with the bait protein
            current_pool = {bait_name: bait_sequence}
            current_pool_size = len(bait_sequence)

            # Shuffle the once_used list to randomize selection
            random.shuffle(once_used)

            # Create a temporary list to hold proteins we've looked at
            tmp = []

            # Iterate over all once_used proteins and try to add them to the current pool
            for _ in range(len(once_used)-1, -1, -1):  # Iterate backwards to safely remove elements
                prot_index = once_used.pop()
                prot_length = length_lookup[prot_index]

                # Check if adding this protein would exceed the max pool depth
                if current_pool_size + prot_length > max_pool_depth:
                    tmp.append(prot_index)
                    continue  # Skip this protein since it would exceed the max pool depth

                # Check if this protein has been in a pool with any of the proteins already in the current pool
                overlap = False
                for existing_prot in current_pool.keys():
                    if existing_prot == bait_name:
                        continue
                    if same_pool[prot_index][id_lookup.index(existing_prot)]:
                        overlap = True
                        break

                if overlap:
                    tmp.append(prot_index)
                    continue  # Skip this protein since it has been in a pool with an existing protein

                # Add the protein to the current pool
                current_pool_size += prot_length
                prot_name = id_lookup[prot_index]
                current_pool[prot_name] = test_proteins[prot_name]

                # Update the same_pool matrix to indicate that these proteins have been in a pool together
                for existing_prot in current_pool.keys():
                    # Continue if it's the bait protein, since it is in every pool
                    if existing_prot == bait_name:
                        continue
                    existing_index = id_lookup.index(existing_prot)
                    same_pool[prot_index][existing_index] = True
                    same_pool[existing_index][prot_index] = True

            # Add the completed pool to the list of pools
            pools.append(current_pool)

            # Add back any proteins that were not added to the current pool
            once_used = tmp

        return pools

# Test run: read in example_fa1.fasta and example_fa2.fasta from "../examples" directory, generate pools with max_pool_depth=1000, and print the resulting pools.
if __name__ == "__main__":
    # Read in the example FASTA files
    test_proteins = read_fasta("examples/simulated_proteins_1000.fasta")

    # Randomly select one protein to be the bait protein and remove it from test_proteins
    bait_protein_name = random.choice(list(test_proteins.keys()))
    bait_protein = {bait_protein_name: test_proteins.pop(bait_protein_name)}

    # # Select the first protein in the dictionary to be the bait protein and remove it from test_proteins
    # bait_protein_name = next(iter(test_proteins))
    # bait_protein = {bait_protein_name: test_proteins.pop(bait_protein_name)}

    # Print the selected bait protein
    print(f"Selected bait protein: {bait_protein_name}")

    # Print proteins with name and length
    for name, seq in test_proteins.items():
        print(f"{name}: {len(seq)} amino acids")
    print(f"Remaining proteins: {len(test_proteins)}")

    # For 10 iterations, generate pools with max_pool_depth=4000. Keep the set of pools with the fewest total number of pools.
    best_pools = None
    for iteration in range(10):
        print(f"\nIteration {iteration + 1}:")
        pools = generate_pools(bait_protein, test_proteins, max_pool_depth=4000)

        # Keep the set of pools with the fewest total number of pools
        if best_pools is None or len(pools) < len(best_pools):
            best_pools = pools

    # Print the best set of pools
    print("\nBest set of pools:")
    for i, pool in enumerate(best_pools):
        print(f"Pool {i+1}:")
        for name, seq in pool.items():
            print(f">{name}")

    # Sanity check: ensure each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in both of its pools.
    protein_pool_count = {name: 0 for name in test_proteins.keys()}
    for pool in best_pools:
        # Check that the bait protein is in every pool
        assert bait_protein_name in pool, f"Bait protein {bait_protein_name} not found in pool."
        for name, seq in pool.items():
            if name != bait_protein_name:
                protein_pool_count[name] += 1
    print(f"Bait protein {bait_protein_name} is in every pool.")

    # Check that each protein is in at least two pools
    for name, count in protein_pool_count.items():
        assert count >= 2, f"Protein {name} is in only {count} pools."
    print("All proteins are in at least two pools.")

    # Check that no protein has any given protein in both of its pools
    for name, seq in pool.items():
        for i, pool in enumerate(best_pools):
            if name != bait_protein_name:
                for other_name, other_seq in pool.items():
                    if other_name != bait_protein_name and other_name != name:
                        assert other_name not in seq, f"Protein {name} and {other_name} are both in pool {i+1}."
    print("No protein has any given protein in both of its pools.")
