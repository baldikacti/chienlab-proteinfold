#!/usr/bin/env python3

"""
Script to generate pools of proteins, including the bait protein, to pass into AF3 for folding prediction.
This script generates pools such that each pool has fewer amino acids than the maximum pool depth (specified in config files),
each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in each of its pools.
generate_index, generate_pools written by Jonathan Coombs, read_fasta written by Berent Aldikacti, all_vs_all adapted by Jonathan Coombs from pooled_ppi in afpool paper.
"""

import json
import os
import random
import sys
from pathlib import Path

import yaml

# Global arrays, defined at the start of the script
depth_lookup = None  # Int array from 0 to max_pool_depth, where each entry links a remaining depth to the index of the longest fitting protein


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


def sanitize_header(header: str) -> str:
    """Sanitize a FASTA header to a single unique identifier.

    Args:
        header: FASTA header line, with or without the leading '>'.
    """
    header = header.removeprefix(">").strip()
    if not header:
        raise ValueError("FASTA entry has an empty header.")

    # Everything after the first whitespace is description text
    first_word = header.split()[0]

    # UniProt style: db|UniqueIdentifier|EntryName -> keep the UniqueIdentifier
    fields = first_word.split("|")
    if len(fields) >= 2 and fields[1]:
        return fields[1]

    return first_word


# Function to parse FASTA files.
# Input: fasta_file: a FASTA file with one or more sequences.
# Output: a dictionary with protein name as key and sequence as value.
# Protein name is treated as the first word in a name line, represented by a line starting with '>'.
# Raises exception if: FASTA file is not found; FASTA file is empty.
def read_fasta(fasta_file: str | Path) -> str | dict[str, str]:
    """Read a FASTA file and return the sequence(s).

    Args:
        fasta_file: Path to FASTA file
    """
    # Build full path using workdir if fasta_file is relative
    fasta_path = Path(fasta_file)
    if not fasta_path.is_absolute():
        fasta_path = "./" / fasta_path

    try:
        with open(fasta_path, "r") as f:
            lines = f.readlines()

        sequences = {}
        current_header = None
        current_sequence = ""

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences[current_header] = current_sequence

                # Start new sequence
                current_header = sanitize_header(line)
                if current_header in sequences:
                    raise ValueError(
                        f"FASTA file {fasta_path} contains duplicate entry '{current_header}'."
                    )
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
                raise ValueError(
                    f"FASTA file {fasta_path} contains empty sequence for entry '{header}'."
                )

        return sequences

    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file {fasta_path} not found.")
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Error reading FASTA file {fasta_path}: {e}")


def export_pool(
    proteins: dict[str, str],
    pool_id: int,
    output_dir: Path,
    mode: str,
) -> None:
    """
    Export a protein pool in the specified format.

    Args:
        proteins: Dictionary with protein names as keys and sequences as values.
        pool_id: ID of the pool to export.
        output_dir: Directory to save the output file.
        mode: Export format, one of 'boltz', 'alphafold3', or 'colabfold'.
    """

    mode = mode.lower()

    if mode == "boltz":
        output_path = output_dir / f"pool_{pool_id}.yaml"

        data = {
            "sequences": [
                {
                    "protein": {
                        "id": protein_id,
                        "sequence": sequence,
                    }
                }
                for protein_id, sequence in proteins.items()
            ],
        }

        with output_path.open("w") as f:
            yaml.safe_dump(data, f, sort_keys=False)

    elif mode == "alphafold3":
        output_path = output_dir / f"pool_{pool_id}.json"

        data = {
            "dialect": "alphafold3",
            "version": 1,
            "name": f"pool_{pool_id}",
            "sequences": [
                {
                    "protein": {
                        "id": protein_id,
                        "sequence": sequence,
                    }
                }
                for protein_id, sequence in proteins.items()
            ],
        }

        with output_path.open("w") as f:
            json.dump(data, f, indent=2)

    elif mode == "colabfold":
        output_path = output_dir / f"pool_{pool_id}.fasta"

        protein_ids = "_".join(proteins.keys())
        protein_sequences = ":".join(proteins.values())

        wrap_length = 80
        with output_path.open("w") as f:
            f.write(f">{protein_ids}\n")
            # Split the string every 'wrap_length' characters
            for i in range(0, len(protein_sequences), wrap_length):
                f.write(protein_sequences[i : i + wrap_length] + "\n")

    else:
        raise ValueError(
            f"Unknown mode '{mode}'. Expected 'boltz', 'alphafold3', or 'colabfold'."
        )


# Helper function to randomly select an available protein index with length <= remaining_depth, or return -1 if no such protein exists.
# Inputs:   fenwick_tree: a FenwickTree object representing the active/inactive state of protein indices.
#           shared_prots: a set of protein indices that are shared and should be avoided.
#           remaining_depth: the maximum number of amino acids that can still be added to the current pool.
# Output: an integer representing the index of a random protein that can be added to the current pool, or -1 if no such protein exists.
def generate_index(
    fenwick_tree: FenwickTree, shared_prots: set, remaining_depth: int
) -> int:
    """Randomly select an available protein index with length <= remaining_depth, or return -1 if no such protein exists.

    #     Args:
    #         fenwick_tree: A FenwickTree object representing the active/inactive state of protein indices.
    #         shared_prots: A set of protein indices that are shared and should be avoided.
    #         remaining_depth: Maximum number of amino acids that can still be added to the current pool.
    #"""
    # Check if remaining_depth is valid
    if remaining_depth < 0 or remaining_depth >= len(depth_lookup):
        raise ValueError(
            f"remaining_depth must be between 0 and max_pool_depth={len(depth_lookup) - 1}, got {remaining_depth}."
        )

    # Get the index of the longest fitting protein for the given remaining depth
    longest_fitting_index = depth_lookup[remaining_depth]

    # Use the Fenwick tree to find a random active index that is not in shared_prots.
    # Start from a uniformly random rank, then walk the remaining ranks in order (wrapping
    # around), so each eligible protein is examined at most once and -1 is returned only
    # when every fitting protein is in shared_prots.
    num_active = fenwick_tree.prefix_sum(longest_fitting_index)
    if num_active == 0:
        return -1

    start = random.randrange(num_active)
    for offset in range(num_active):
        index = fenwick_tree.kth_instance(
            (start + offset) % num_active + 1
        )  # kth_instance is 1-indexed, so ranks run from 1 to num_active
        if index not in shared_prots:
            return index

    return -1


# Function to generate a list of protein pools, following the rules specified in this file.
# Inputs:   bait_protein: a dictionary with one entry, the bait protein name and sequence.
#           test_proteins: a dictionary with one or more entries, the test protein names and sequences.
#           max_pool_depth: the maximum number of amino acids in each pool.
# Output: a list of dictionaries, each dictionary representing a pool of proteins, with protein names as keys and sequences as values.
# Raises exception if: bait_protein is not size 1; test_proteins is empty; max_pool_depth is not a positive integer;
# any protein is too long to be in a pool with the bait protein.
def generate_pools(
    bait_protein: dict[str, str], test_proteins: dict[str, str], max_pool_depth: int
) -> list[dict[str, str]]:
    """Generate pools of proteins for AF3 pooled folding prediction.

    Args:
        bait_protein: Dictionary with one entry, the bait protein name and sequence.
        test_proteins: Dictionary with one or more entries, the test protein names and sequences.
        max_pool_depth: Maximum number of amino acids in each pool.
    """
    if len(bait_protein) != 1:
        raise ValueError("bait_protein must contain exactly one entry.")
    if not isinstance(max_pool_depth, int) or max_pool_depth <= 0:
        raise ValueError("max_pool_depth must be a positive integer.")

    # Save bait protein name and sequence
    bait_name, bait_sequence = next(iter(bait_protein.items()))

    # Removes bait from the pool if it exists in the pool
    test_proteins = {
        name: seq
        for name, seq in test_proteins.items()
        if not (name == bait_name and seq == bait_sequence)
    }

    if not test_proteins:
        raise ValueError("test_proteins must contain at least one entry.")
    if any(
        len(seq) > max_pool_depth - len(bait_sequence) for seq in test_proteins.values()
    ):
        raise ValueError(
            "One or more test proteins are too long to be in a pool with the bait protein."
        )

    # Initialize pools
    pools = []

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
        range(len(protein_names)), key=protein_lengths.__getitem__
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
    used_fenwick = FenwickTree(
        len(id_lookup), 0
    )  # Start with all proteins inactive in the used_fenwick tree
    # saved_fenwick = FenwickTree(len(id_lookup), 1)  # To track proteins that have been used twice, currently not being used

    shared_prots = {i: set() for i in range(len(id_lookup))}

    num_iters = 0  # Counter, up to 2
    num_used_prots = (
        0  # Counter for the number of proteins that have been used at least once
    )
    current_pool_id = -1
    print("Done with setup, starting iterations")

    while num_iters < 2:
        while num_used_prots < len(
            id_lookup
        ):  # Iterate until all proteins have been used at least once
            # Start a new pool, with the bait protein's length already counted
            current_pool_id += 1
            current_pool = []  # List of indices of proteins in the current pool
            current_pool_size = len(bait_sequence)
            current_shared_prots = set()  # Set of indices of proteins that share a pool with any protein in the current pool

            # Step 1: iterate over all unused proteins short enough to fit in the pool
            while (
                new_index := generate_index(
                    unused_fenwick,
                    current_shared_prots,
                    max_pool_depth - current_pool_size,
                )
            ) != -1:
                # Exclude this protein and everything it has already shared a pool with
                # from the rest of this pool. The protein itself has to be excluded because
                # it stays active in used_fenwick, which step 2 draws from.
                current_shared_prots.update(shared_prots[new_index])
                current_shared_prots.add(new_index)

                # Update the shared_prots list to indicate that these proteins have been in a pool together
                for existing_index in current_pool:
                    shared_prots[new_index].add(existing_index)
                    shared_prots[existing_index].add(new_index)

                # Add the protein to the current pool
                current_pool.append(new_index)
                current_pool_size += length_lookup[new_index]

                # Update the Fenwick tree to mark this protein as used
                unused_fenwick.set(new_index, 0)
                used_fenwick.set(
                    new_index, 1
                )  # Mark this protein as used in the used_fenwick tree

                num_used_prots += 1

            # Step 2: If in iteration 0, check if any proteins in extra_prots can fit into the current pool
            while (
                new_index := generate_index(
                    used_fenwick,
                    current_shared_prots,
                    max_pool_depth - current_pool_size,
                )
            ) != -1:
                # Exclude this protein and everything it has already shared a pool with
                # from the rest of this pool
                current_shared_prots.update(shared_prots[new_index])
                current_shared_prots.add(new_index)

                # Update the shared_prots list to indicate that these proteins have been in a pool together
                for existing_index in current_pool:
                    shared_prots[new_index].add(existing_index)
                    shared_prots[existing_index].add(new_index)

                # Add the protein to the current pool
                current_pool.append(new_index)
                current_pool_size += length_lookup[new_index]

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
        num_used_prots = len(id_lookup) - used_fenwick.prefix_sum(
            len(id_lookup) - 1
        )  # Count the number of proteins that have already been used twice
        unused_fenwick = used_fenwick
        used_fenwick = FenwickTree(
            len(id_lookup), 0
        )  # Start with all proteins inactive in the used_fenwick tree
        # used_fenwick = saved_fenwick
        # saved_fenwick = FenwickTree(len(id_lookup), 1)
        num_iters += 1
        # num_twice_used_prots = 0

    return pools


"""
Generate pools as in https://doi.org/10.1101/2025.07.01.662654 except that interactions are weighted by the product of the protein sizes
Optimised for maximum performance using incremental updates and Numba JIT.
"""

import argparse  # , pandas as pd
import itertools

import numba
import numpy as np
import tqdm


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


@numba.njit(parallel=True, nogil=True)
def calculate_uncovered_weights(sizes, covered):
    """Calculate the sum of uncovered interaction weights for each protein."""
    n = sizes.shape[0]
    weights = np.zeros(n, dtype=np.float64)
    for i in numba.prange(n):
        w = 0.0
        for j in range(n):
            if i != j and not covered[i, j]:
                w += sizes[i] * sizes[j]
        weights[i] = w
    return weights


@numba.njit(nogil=True)
def get_initial_large_pools_numba(sizes, covered, max_size):
    """Find pairs of proteins that exceed max_size and are not yet covered."""
    n = sizes.shape[0]
    res = []
    for i in range(n):
        for j in range(i + 1, n):
            if sizes[i] + sizes[j] > max_size and not covered[i, j]:
                res.append((i, j))
    return res


@numba.njit(parallel=True, nogil=True)
def find_best_i_parallel(pool_mask, sizes, pool_C, current_pool_size, max_size):
    n = sizes.shape[0]
    scores = np.full(n, -1e18, dtype=np.float64)
    for i in numba.prange(n):
        if not pool_mask[i] and current_pool_size + sizes[i] <= max_size:
            # Score: newly covered interaction weight per unit size of protein i.
            # Newly covered interaction weight for adding i to current pool is 2 * sizes[i] * (current_pool_size - pool_C[i])
            # Dividing by sizes[i] results in (current_pool_size - pool_C[i])
            scores[i] = current_pool_size - pool_C[i]

    best_i = np.argmax(scores)

    if scores[best_i] <= 0:
        return -1
    return best_i


@numba.njit(parallel=True, nogil=True)
def update_pool_C_parallel(pool_C, covered, best_i, size_best_i):
    n = pool_C.shape[0]
    for i in numba.prange(n):
        if covered[i, best_i]:
            pool_C[i] += size_best_i


@numba.njit(nogil=True)
def update_coverage_numba(pool_ix, covered, sizes, uncovered_weight_per_protein):
    """Mark all interactions within the finished pool as covered."""
    newly_covered_count = 0
    for i_idx in range(len(pool_ix)):
        p1 = pool_ix[i_idx]
        for j_idx in range(i_idx + 1, len(pool_ix)):
            p2 = pool_ix[j_idx]
            if not covered[p1, p2]:
                covered[p1, p2] = 1
                covered[p2, p1] = 1
                weight = sizes[p1] * sizes[p2]
                uncovered_weight_per_protein[p1] -= weight
                uncovered_weight_per_protein[p2] -= weight
                newly_covered_count += 1
    return newly_covered_count


def all_vs_all(sizes, max_size=5120, skip_pairs=None):
    if skip_pairs is None:
        skip_pairs = []
    n = len(sizes)
    covered = np.zeros((n, n), dtype=np.uint8)
    np.fill_diagonal(covered, 1)

    for i, j in skip_pairs:
        if i < n and j < n:
            covered[i, j] = 1
            covered[j, i] = 1

    uncovered_weight_per_protein = calculate_uncovered_weights(sizes, covered)

    # Initial large pairs
    large_pairs = get_initial_large_pools_numba(sizes, covered, max_size)
    for i, j in large_pairs:
        yield ({i, j}, sizes[i] + sizes[j])
        covered[i, j] = 1
        covered[j, i] = 1
        weight = sizes[i] * sizes[j]
        uncovered_weight_per_protein[i] -= weight
        uncovered_weight_per_protein[j] -= weight

    total_interactions = n * (n - 1) // 2
    covered_interactions = int(np.sum(covered) // 2 - n // 2)

    pbar = tqdm.tqdm(total=total_interactions)
    pbar.update(covered_interactions)

    pool_mask = np.zeros(n, dtype=np.bool_)
    pool_C = np.zeros(n, dtype=np.float64)

    while covered_interactions < total_interactions:
        avail = np.where(uncovered_weight_per_protein > 1e-6)[0]
        if len(avail) == 0:
            break

        # Greedy: pick protein with most uncovered interaction weight
        avail_choice = avail[np.argmax(uncovered_weight_per_protein[avail])]

        pool_mask.fill(False)
        pool_mask[avail_choice] = True

        # Initialise pool_C for the first protein
        pool_C.fill(0)
        update_pool_C_parallel(pool_C, covered, avail_choice, sizes[avail_choice])

        current_pool_size = sizes[avail_choice]

        while True:
            best_i = find_best_i_parallel(
                pool_mask, sizes, pool_C, current_pool_size, max_size
            )
            if best_i == -1:
                break

            current_pool_size += sizes[best_i]
            pool_mask[best_i] = True
            update_pool_C_parallel(pool_C, covered, best_i, sizes[best_i])

        pool_ix = np.where(pool_mask)[0]
        yield (set(pool_ix.tolist()), float(current_pool_size))

        newly_covered = update_coverage_numba(
            pool_ix, covered, sizes, uncovered_weight_per_protein
        )
        covered_interactions += newly_covered
        pbar.update(newly_covered)

    pbar.close()


def run_bait_vs_all(
    bait_fasta: Path,
    pool_fasta: Path,
    output_dir: Path,
    max_pool_depth: int,
    export_mode: str,
):
    # Read in the input FASTA files
    bait_protein = read_fasta(bait_fasta)
    test_proteins = read_fasta(pool_fasta)

    pools = generate_pools(bait_protein, test_proteins, max_pool_depth)

    print(f"Generated {len(pools)} pools with max pool depth {max_pool_depth}")

    # Sanity check: ensure each protein is in at least two pools, the bait protein is in every pool, and no protein has any given protein in both of its pools.
    bait_protein_name = next(iter(bait_protein))
    protein_pool_count = {name: 0 for name in test_proteins}
    for pool in pools:
        # Check that the bait protein is in every pool
        assert bait_protein_name in pool, (
            f"Bait protein {bait_protein_name} not found in pool."
        )
        for name, seq in pool.items():
            if name != bait_protein_name:
                protein_pool_count[name] += 1
    print(f"Bait protein {bait_protein_name} is in every pool.")

    # Check that each protein is in at least two pools
    three_count = 0
    for name, count in protein_pool_count.items():
        assert count >= 2, f"Protein {name} is in only {count} pools."
        if count == 3:
            three_count += 1
        if count > 3:
            print(f"Warning: Protein {name} is in {count} pools.")
    print(
        f"All proteins are in at least two pools, {three_count} are in exactly three pools."
    )

    # Check that no protein has any given protein in both of its pools
    for name, seq in pool.items():
        for i, pool in enumerate(pools):
            if name != bait_protein_name:
                for other_name, other_seq in pool.items():
                    if other_name != bait_protein_name and other_name != name:
                        assert other_name not in seq, (
                            f"Protein {name} and {other_name} are both in pool {i + 1}."
                        )
    print("No protein shares multiple pools with any other protein.")

    # Export pools to the specified output directory in the chosen format
    # Build a .tsv file with three columns: pool_id, protein_ids (underscore-separated), pool_size (sum of lengths of all proteins in the pool)
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(f"{output_dir}/pools.tsv", "w") as f:
        f.write("pool_id\tprotein_ids\tpool_size\n")
        for i, pool in enumerate(pools):
            protein_ids = "_".join(pool.keys())
            pool_size = sum(len(seq) for seq in pool.values())
            f.write(f"{i}\t{protein_ids}\t{pool_size}\n")
            export_pool(pool, i, output_dir, export_mode)


def run_all_vs_all(
    pool_fasta: Path,
    init_pools: Path | None,
    max_pools: int | None,
    output_dir: Path,
    max_pool_depth: int,
    export_mode: str,
) -> Path:

    proteins = read_fasta(pool_fasta)
    protein_ids = list(proteins.keys())
    protein_id_to_ix = {id_: i for i, id_ in enumerate(protein_ids)}
    sizes = np.fromiter(
        (len(sequence) for sequence in proteins.values()), dtype=np.float64
    )

    eprint(numba.get_num_threads(), "threads available for numba")
    eprint(len(protein_ids), "proteins in input")

    skip_pairs = []
    if init_pools is not None:
        initial_pools = read_fasta(init_pools)
        # Each header of an initial pool FASTA names a pool as underscore-joined protein ids
        for pool_id in initial_pools:
            ids = pool_id.split("_")
            ixs = [protein_id_to_ix[id_] for id_ in ids if id_ in protein_id_to_ix]
            if len(ixs) >= 2:
                skip_pairs.extend(itertools.combinations(sorted(ixs), 2))

    pools_data_id = []
    pools_data_size = []
    pools_data_ixs = []

    pool_gen = all_vs_all(sizes, max_pool_depth, skip_pairs=skip_pairs)

    if max_pools is not None:
        pool_gen = itertools.islice(pool_gen, max_pools)

    for pool_ixs, pool_size in pool_gen:
        pool_ids_subset = sorted([protein_ids[ix] for ix in pool_ixs])
        pools_data_id.append("_".join(pool_ids_subset))
        pools_data_size.append(pool_size)
        pools_data_ixs.append(list(pool_ixs))

    if not pools_data_id:
        eprint("No pools generated")
        return

    eprint(len(pools_data_id), "pools generated")

    # Redundancy and completeness sanity checks (optimised)
    all_interactions_count = len(sizes) * (len(sizes) - 1) // 2
    all_sum = (np.sum(sizes) ** 2 - np.sum(sizes**2)) / 2

    # Check if all possible interactions are covered across all pools
    unique_pairs_count = 0
    covered_check = np.zeros((len(sizes), len(sizes)), dtype=np.uint8)
    actual_gen_sum = 0
    for pool_ixs in pools_data_ixs:
        for i, p1 in enumerate(pool_ixs):
            for p2 in pool_ixs[i + 1 :]:
                if p1 < p2:
                    low, high = p1, p2
                else:
                    low, high = p2, p1
                if not covered_check[low, high]:
                    covered_check[low, high] = 1
                    unique_pairs_count += 1
                actual_gen_sum += sizes[p1] * sizes[p2]

    eprint(all_interactions_count, "interactions expected")
    eprint(unique_pairs_count, "interactions across all pools generated")
    eprint(
        unique_pairs_count == all_interactions_count,
        "pools include all possible interactions",
    )
    eprint(
        actual_gen_sum / all_sum, "length-weighted redundancy factor across all pools"
    )

    # Export pools to the specified output directory in the chosen format
    # Build a .tsv file with three columns: pool_id, protein_ids (underscore-separated), pool_size (sum of lengths of all proteins in the pool)
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(f"{output_dir}/pools.tsv", "w") as f:
        f.write("pool_id\tprotein_ids\tpool_size\n")
        for i, (pool_id, pool_size) in enumerate(zip(pools_data_id, pools_data_size)):
            f.write(f"{i}\t{pool_id}\t{pool_size}\n")
            export_pool(
                {
                    protein_ids[ix]: proteins[protein_ids[ix]]
                    for ix in pools_data_ixs[i]
                },
                i,
                output_dir,
                mode=export_mode,
            )


def main():

    parser = argparse.ArgumentParser(
        description="FASTA file of pooled proteins to a selected pooled operation"
    )
    parser.add_argument(
        "--mode",
        "-m",
        required=True,
        choices=["bait_vs_all", "all_vs_all"],
        help="Mode of operation: 'bait_vs_all' or 'all_vs_all'",
    )
    parser.add_argument(
        "--pool-fasta",
        "-p",
        required=True,
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
        type=int,
        help="Set the maximum depth of each pool (default: 5000)",
    )
    parser.add_argument(
        "--seed",
        "-s",
        default=None,
        type=int,
        help="Random seed, for reproducible pools (default: non-deterministic)",
    )
    parser.add_argument(
        "--export-mode",
        "-e",
        default="boltz",
        choices=["boltz", "alphafold3", "colabfold"],
        help="Export mode for pool files (boltz, alphafold3, or colabfold)",
    )
    # Bait vs all argument
    parser.add_argument(
        "--bait-fasta",
        "-b",
        required=False,
        help="FASTA file containing the bait protein",
    )
    # All vs all arguments
    parser.add_argument(
        "--init-pools",
        "-i",
        required=False,
        help="Initial pools to omit from the all_vs_all generation",
    )
    parser.add_argument(
        "--max-pools",
        "-n",
        required=False,
        type=int,
        help="Maximum number of pools to sample in all_vs_all mode",
    )

    args = parser.parse_args()

    output_dir = Path(args.output)

    # Optionally sets seed for reproducibility of the pools
    if args.seed is not None:
        random.seed(args.seed)

    # Check if input files exist
    if not os.path.isfile(args.pool_fasta):
        print(f"Error: Pool FASTA file '{args.pool_fasta}' does not exist")
        sys.exit(1)

    # Split based on mode
    if args.mode == "bait_vs_all":
        if not args.bait_fasta:
            print("Error: --bait-fasta is required for bait_vs_all mode")
            sys.exit(1)
        if not os.path.isfile(args.bait_fasta):
            print(f"Error: Bait FASTA file '{args.bait_fasta}' does not exist")
            sys.exit(1)
        run_bait_vs_all(
            args.bait_fasta,
            args.pool_fasta,
            output_dir,
            args.max_pool_depth,
            args.export_mode,
        )

    elif args.mode == "all_vs_all":
        init_pools = None
        if args.init_pools:
            if not os.path.isfile(args.init_pools):
                print(f"Error: Initial pools file '{args.init_pools}' does not exist")
                sys.exit(1)
            init_pools = args.init_pools

        max_pools = None
        if args.max_pools is not None:
            max_pools = args.max_pools

        run_all_vs_all(
            args.pool_fasta,
            init_pools,
            max_pools,
            output_dir,
            args.max_pool_depth,
            args.export_mode,
        )


if __name__ == "__main__":
    main()
