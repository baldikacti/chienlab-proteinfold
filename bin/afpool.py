#!/usr/bin/env python3

"""
Script to generate pools of proteins, including the bait protein, to pass into AF3 for folding prediction.
This script generates pools such that each pool has fewer amino acids than the maximum pool depth (specified in config files),
each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in each of its pools.
Written by Jonathan Coombs
"""

import os
import random
from pathlib import Path
import sys
from typing import Dict, List, Union

# Global arrays, defined at the start of the script
depth_lookup = None # Int array from 0 to max_pool_depth, where each entry links a remaining depth to the index of the longest fitting protein

class FenwickTree:
    def __init__(self, size: int, initial_value: int = 1):
        """
        Initializes a 0-based Fenwick Tree of the given size.
        Every index starts as active (1).
        """
        self.size = size
        # Track the active/inactive state of each index directly (1 or 0)
        self.data = [initial_value] * self.size
        # The Fenwick tree structure storing the range sums
        self.tree = [initial_value] * self.size
        
        # O(n) linear-time tree construction
        for i in range(self.size):
            parent = i | (i + 1)
            if parent < self.size:
                self.tree[parent] += self.tree[i]

    def set(self, index: int, value: int) -> None:
        """
        Flips the index to active (1) or inactive (0).
        Propagates the change across the tree if the state changes.
        """
        if not (0 <= index < self.size):
            raise IndexError("Index out of bounds")
        if value not in (0, 1):
            raise ValueError("Value must be either 0 or 1")
            
        delta = value - self.data[index]
        if delta == 0:
            return  # No state change, skip updates
            
        self.data[index] = value
        
        # Propagate up using 0-based bitwise navigation
        i = index
        while i < self.size:
            self.tree[i] += delta
            i = i | (i + 1)

    def prefix_sum(self, index: int) -> int:
        """
        Returns the total number of active indices from 0 up to the given index (inclusive).
        Time Complexity: O(log n)
        """
        if not (0 <= index < self.size):
            return 0
            
        total = 0
        i = index
        while i >= 0:
            total += self.tree[i]
            i = (i & (i + 1)) - 1  # Move down using 0-based bitwise navigation
        return total

    def kth_instance(self, k: int) -> int:
        """
        Finds the 0-based index of the k-th active element (1-indexed 'k').
        Returns -1 if there are fewer than k active instances in the tree.
        Time Complexity: O(log n)
        """
        if k <= 0:
            return -1
        
        index = -1
        power = 1
        while (power << 1) <= self.size:
            power <<= 1
            
        # Binary lift down to find the exact target index
        while power > 0:
            next_idx = index + power
            if next_idx < self.size and self.tree[next_idx] < k:
                k -= self.tree[next_idx]
                index = next_idx
            power >>= 1
            
        final_idx = index + 1
        if final_idx < self.size and k == 1:
            return final_idx
            
        return -1

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

# Helper function to randomly select an available protein index with length <= remaining_depth, or return -1 if no such protein exists.
# Inputs:   fenwick_tree: a FenwickTree object representing the active/inactive state of protein indices.
#           shared_prots: a set of protein indices that are shared and should be avoided.
#           remaining_depth: the maximum number of amino acids that can still be added to the current pool.
# Output: an integer representing the index of a random protein that can be added to the current pool, or -1 if no such protein exists.
def generate_index(fenwick_tree: FenwickTree, shared_prots: set, remaining_depth: int) -> int:
    """Randomly select an available protein index with length <= remaining_depth, or return -1 if no such protein exists.
    
#     Args:
#         fenwick_tree: A FenwickTree object representing the active/inactive state of protein indices.
#         shared_prots: A set of protein indices that are shared and should be avoided.
#         remaining_depth: Maximum number of amino acids that can still be added to the current pool.
#     """
    # Check if remaining_depth is valid
    if remaining_depth < 0 or remaining_depth >= len(depth_lookup):
        raise ValueError(f"remaining_depth must be between 0 and max_pool_depth={len(depth_lookup)-1}, got {remaining_depth}.")

    # Get the index of the longest fitting protein for the given remaining depth
    longest_fitting_index = depth_lookup[remaining_depth]

    # Use the Fenwick tree to find a random active index, repeating if the index is in shared_prots
    num_active = fenwick_tree.prefix_sum(longest_fitting_index)
    hits = [] # Proteins in shared_prots that we have already seen, to avoid infinite loops
    for i in range(num_active):  # Check up to the number of possible proteins, stop if all possible indices are in shared_prots
        index = fenwick_tree.kth_instance(int(random.random() * (num_active + 1))) # +1 since random.random() is 1-exclusive
        if index not in shared_prots:
            return index
        elif index in hits:
            i -= 1  # Decrement i to try again, since we hit a duplicate
        else:
            hits.append(index)  # Add to hits to avoid infinite loops

    return -1

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
        
        # Save bait protein name and sequence
        bait_name, bait_sequence = next(iter(bait_protein.items()))

        # Data structures:
        #   - id_lookup: str array lookup table for protein ID, where the shortest protein is at index 0 and the longest protein is at index len(id_lookup)-1
        #   - depth_lookup (global): int array from 0 to max_pool_depth, where each entry links a remaining depth to the index of the longest fitting protein
        #   - length_lookup: int array lookup table for protein length, in order of sorted_protein_indices
        #   - shared_prots: int list of sets containing protein index for any protein that shares a pool
        #   - used_prots: int array containing proteins that have been used, and can be used again in the next iteration
        #   - extra_prots: int array containing the subset of used_prots that can be added to the current pool, but have been in a pool already

        protein_names = list(test_proteins.keys())
        protein_lengths = [len(seq) for seq in test_proteins.values()]

        sorted_protein_indices = sorted(
            range(len(protein_names)),
            key=protein_lengths.__getitem__
        )

        id_lookup = [protein_names[i] for i in sorted_protein_indices]
        length_lookup = [protein_lengths[i] for i in sorted_protein_indices]

        global depth_lookup
        depth_lookup = [-1] * (max_pool_depth + 1)  # Initialize depth lookup table with -1
        # Fill the depth_lookup table
        index = 0
        for i in range(max_pool_depth + 1):
            while index < len(length_lookup) and length_lookup[index] <= i:
                index += 1
            depth_lookup[i] = index - 1  # Store the index of the longest fitting protein

        # Generate 2 Fenwick trees: one for proteins that have not been used yet, and one for proteins that have been used once.
        unused_fenwick = FenwickTree(len(id_lookup))
        used_fenwick = FenwickTree(len(id_lookup), 0)  # Start with all proteins inactive in the used_fenwick tree
        # saved_fenwick = FenwickTree(len(id_lookup), 1)  # To track proteins that have been used twice, currently not being used

        shared_prots = {i: set() for i in range(len(id_lookup))}

        num_iters = 0 # Counter, up to 2
        num_used_prots = 0 # Counter for the number of proteins that have been used at least once
        num_twice_used_prots = 0 # Counter for the number of proteins that have been used at least twice, currently not being used
        current_pool_id = -1
        print(f"Done with setup, starting iterations")

        while num_iters < 2:
            while num_used_prots < len(id_lookup): # Iterate until all proteins have been used at least once
                # Start a new pool, with the bait protein's length already counted
                current_pool_id += 1
                current_pool = [] # List of indices of proteins in the current pool
                current_pool_size = len(bait_sequence)
                current_shared_prots = set() # Set of indices of proteins that share a pool with any protein in the current pool

                # Step 1: iterate over all unused proteins short enough to fit in the pool
                while (new_index := generate_index(unused_fenwick, current_shared_prots, max_pool_depth - current_pool_size)) != -1:
                    # Add all proteins shared with this protein to the used_prots list
                    current_shared_prots.update(shared_prots[new_index])
                    
                    # Add the protein to the current pool
                    current_pool.append(new_index)
                    current_pool_size += length_lookup[new_index]

                    # Update the shared_prots list to indicate that these proteins have been in a pool together
                    for existing_index in current_pool:
                        shared_prots[new_index].add(existing_index)
                        shared_prots[existing_index].add(new_index)

                    # Update the Fenwick tree to mark this protein as used
                    unused_fenwick.set(new_index, 0)
                    used_fenwick.set(new_index, 1)  # Mark this protein as used in the used_fenwick tree

                    num_used_prots += 1

                # Step 2: If in iteration 0, check if any proteins in extra_prots can fit into the current pool
                while (new_index := generate_index(used_fenwick, current_shared_prots, max_pool_depth - current_pool_size)) != -1:
                    # Add all proteins shared with this protein to the used_prots list
                    current_shared_prots.update(shared_prots[new_index])
                                    
                    # Add the protein to the current pool
                    current_pool.append(new_index)
                    current_pool_size += length_lookup[new_index]
                
                    # Update the shared_prots list to indicate that these proteins have been in a pool together
                    for existing_index in current_pool:
                        shared_prots[new_index].add(existing_index)
                        shared_prots[existing_index].add(new_index)
                
                    # Update the Fenwick tree to mark this protein as used
                    used_fenwick.set(new_index, 0)
                    # saved_fenwick.set(new_index, 0)  # Mark this protein as used in the saved_fenwick tree
                    # num_twice_used_prots += 1

                # Add the completed pool to the list of pools
                tmp_pool = {bait_name: bait_sequence}
                for prot_index in current_pool:
                    prot_name = id_lookup[prot_index]
                    tmp_pool[prot_name] = test_proteins[prot_name]
                pools.append(tmp_pool)

            # Reset for the second iteration
            num_used_prots = len(id_lookup) - used_fenwick.prefix_sum(len(id_lookup) - 1)  # Count the number of proteins that have already been used twice
            unused_fenwick = used_fenwick
            used_fenwick = FenwickTree(len(id_lookup), 0)  # Start with all proteins inactive in the used_fenwick tree
            # used_fenwick = saved_fenwick
            # saved_fenwick = FenwickTree(len(id_lookup), 1)
            num_iters += 1
            # num_twice_used_prots = 0

        return pools

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Read JSON files into a protein dictionary")
    parser.add_argument(
        "--bait-fasta",
        "-b",
        default=".",
        help="FASTA file containing the bait protein",
    )
    parser.add_argument(
        "--pool-fasta",
        "-p",
        default=".",
        help="FASTA file containing the pool of all test proteins",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="./results",
        help="Output path and directory for pool FASTA files (default: ./results)",
    )
    parser.add_argument(
        "--max-pool-depth",
        "-d",
        default=5000,
        help="Set the maximum depth of each pool (default: 5000)",
    )

    args = parser.parse_args()

    # Check if input files exist
    if not os.path.isfile(args.bait_fasta):
        print(f"Error: Bait FASTA file '{args.bait_fasta}' does not exist")
        sys.exit(1)

    if not os.path.isfile(args.pool_fasta):
        print(f"Error: Pool FASTA file '{args.pool_fasta}' does not exist")
        sys.exit(1)

    # Read in the input FASTA files
    bait_protein = read_fasta(args.bait_fasta)
    test_proteins = read_fasta(args.pool_fasta)

    pools = generate_pools(bait_protein, test_proteins, max_pool_depth=int(args.max_pool_depth))

    # Write the resulting pools to FASTA files
    for i, pool in enumerate(pools):
        os.makedirs(args.output, exist_ok=True)
        with open(f"{args.output}/pools_{i}.fasta", "w") as f:
            for name, seq in pool.items():
                f.write(f">{name}\n{seq}\n")

if __name__ == "__main__":
    main()

# Test run: read in a test FASTA file, generate pools with the first protein as bait and max_pool_depth=5000, and print the resulting pools.
# if __name__ == "__main__":
#     # Read in the example FASTA files
#     test_proteins = read_fasta("examples/simulated_proteins_50000.fasta")

#     # # Randomly select one protein to be the bait protein and remove it from test_proteins
#     # bait_protein_name = random.choice(list(test_proteins.keys()))
#     # bait_protein = {bait_protein_name: test_proteins.pop(bait_protein_name)}

#     # Select the first protein in the dictionary to be the bait protein and remove it from test_proteins
#     bait_protein_name = next(iter(test_proteins))
#     bait_protein = {bait_protein_name: test_proteins.pop(bait_protein_name)}

#     # Print the selected bait protein
#     print(f"Selected bait protein: {bait_protein_name}")

#     # Print proteins with name and length
#     # for name, seq in test_proteins.items():
#         # print(f"{name}: {len(seq)} amino acids")
#     print(f"Remaining proteins: {len(test_proteins)}")

#     # For 10 iterations, generate pools with max_pool_depth=4000. Keep the set of pools with the fewest total number of pools.
#     best_pools = None
#     for iteration in range(1):
#         print(f"\nIteration {iteration + 1}:")
#         pools = generate_pools(bait_protein, test_proteins, max_pool_depth=5000)

#         # Keep the set of pools with the fewest total number of pools
#         if best_pools is None or len(pools) < len(best_pools):
#             best_pools = pools

#     # Print the best set of pools
#     print(f"\nBest set of pools has {len(best_pools)} pools.")
#     # print("\nBest set of pools:")
#     # for i, pool in enumerate(best_pools):
#     #     print(f"Pool {i+1}:")
#     #     for name, seq in pool.items():
#     #         print(f">{name}")

#     # Sanity check: ensure each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in both of its pools.
#     protein_pool_count = {name: 0 for name in test_proteins.keys()}
#     for pool in best_pools:
#         # Check that the bait protein is in every pool
#         assert bait_protein_name in pool, f"Bait protein {bait_protein_name} not found in pool."
#         for name, seq in pool.items():
#             if name != bait_protein_name:
#                 protein_pool_count[name] += 1
#     print(f"Bait protein {bait_protein_name} is in every pool.")

#     # Check that each protein is in at least two pools
#     three_count = 0
#     for name, count in protein_pool_count.items():
#         assert count >= 2, f"Protein {name} is in only {count} pools."
#         if count == 3:
#             three_count += 1
#         if count > 3:
#             print(f"Warning: Protein {name} is in {count} pools.")
#     print(f"All proteins are in at least two pools, {three_count} are in exactly three pools.")

#     # Check that no protein has any given protein in both of its pools
#     for name, seq in pool.items():
#         for i, pool in enumerate(best_pools):
#             if name != bait_protein_name:
#                 for other_name, other_seq in pool.items():
#                     if other_name != bait_protein_name and other_name != name:
#                         assert other_name not in seq, f"Protein {name} and {other_name} are both in pool {i+1}."
#     print("No protein has any given protein in both of its pools.")