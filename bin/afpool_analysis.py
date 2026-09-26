#!/usr/bin/env python3
"""
Script to explode pooled Alphafold3 predictions into per-pair summary confidence scores and output a single ranked TSV file.

Reads the *_summary_confidences.json files Alphafold3 writes for each pool (both the top-level
best-ranked file and the per-sample files under seed-*_sample-*/), and uses the pools.tsv written by
afpool.py to recover the protein id of every chain in the pool.

Alphafold3 reports chain_pair_iptm and chain_pair_pae_min as n_chains x n_chains matrices whose rows
and columns follow the order the chains were declared in the input JSON. afpool.py labels those
chains A, B, C, ... in the same order it lists them in the protein_ids column of pools.tsv, so chain
i belongs to the i-th protein id of that pool. Every entry of the upper triangle is therefore one
protein pair, which is the mapping used by explode_iptms in the afpool paper.

The two pair scores are summarised over the pool's samples in opposite ways, on purpose. pair_iptm is
the mean, which smooths sampling noise and recovers resolution that Alphafold3 loses by rounding its
summary values to two decimals; pair_iptm_best is the same score in the single best-ranked model,
which is the model written to <pool>_model.cif. pae_min instead takes the lowest value seen in any
sample and either direction, so that a close approach predicted by even one model still shows up.

Output is one row per protein pair, sorted by pair_iptm_z in descending order.
"""

import csv
import glob
import itertools
import json
import os
import statistics
import string
import sys
from pathlib import Path

SUMMARY_SUFFIX = "_summary_confidences.json"

HEADERS = [
    "pool_name",
    "protein_a",
    "protein_b",
    "chain_a",
    "chain_b",
    "pair_iptm",
    "pair_iptm_best",
    "pair_iptm_z",
    "pool_iptm",
    "pae_min",
    "min_chain_ptm",
    "has_clash",
]


def enumerate_chains() -> list[str]:
    """Generate the chain labels A to ZZ that afpool.py assigns to pooled proteins, in order."""
    # Copied from export_pool in afpool.py so the labels stay in step with the exported pools
    return list(string.ascii_uppercase) + [
        a + b for a in string.ascii_uppercase for b in string.ascii_uppercase
    ]


def read_pools(pools_tsv: str) -> dict[str, list[str]]:
    """Read afpool.py's pools.tsv and return a lookup of pool name to list of protein ids.

    Args:
        pools_tsv: Path to the pools.tsv written by afpool.py.

    Pools are keyed both by their bare pool_id ("0") and by the job name afpool.py gives them
    ("pool_0"), so either naming of the Alphafold3 output resolves.
    """
    pools = {}
    with open(pools_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pool_id = row["pool_id"].strip()
            protein_ids = row["protein_ids"].split("_")
            pools[pool_id] = protein_ids
            pools["pool_" + pool_id] = protein_ids

    if not pools:
        raise ValueError(f"Pool file {pools_tsv} contains no pools.")

    return pools


def collect_summaries(input_dir: str) -> dict[str, dict]:
    """Group every summary confidence file under input_dir by the pool it belongs to.

    Args:
        input_dir: Directory to search recursively for Alphafold3 summary confidence files.

    Alphafold3 writes one summary file per pool holding the best-ranked sample, plus one per sample
    in seed-<seed>_sample-<sample>/. The sample files carry the seed and sample in their filename,
    which is how the two are told apart here.
    """
    pattern = os.path.join(input_dir, "**", "*" + SUMMARY_SUFFIX)
    summary_files = sorted(glob.glob(pattern, recursive=True))

    pools = {}
    for summary_file in summary_files:
        stem = Path(summary_file).name.removesuffix(SUMMARY_SUFFIX)
        pool_name, _, seed_sample = stem.partition("_seed-")

        pool = pools.setdefault(pool_name, {"best": None, "samples": []})
        if seed_sample:
            pool["samples"].append(summary_file)
        else:
            pool["best"] = summary_file

    return pools


def load_json(path: str) -> dict:
    with open(path, "r") as f:
        return json.load(f)


def mean_matrix(matrices: list[list[list[float]]]) -> list[list[float]]:
    """Average a list of equally shaped n x n matrices elementwise."""
    return [
        [sum(values) / len(matrices) for values in zip(*rows)]
        for rows in zip(*matrices)
    ]


def robust_z(values: list[float]) -> list:
    """Score each value against the median and median absolute deviation of the whole list.

    Args:
        values: Scores to compare against one another.

    The median and MAD are used instead of the mean and standard deviation so that a genuine hit
    does not inflate the spread it is being measured against. Returns empty strings when the MAD is
    zero, which happens when more than half the values are identical and no spread can be estimated.
    """
    median = statistics.median(values)
    mad = statistics.median(abs(value - median) for value in values)

    if mad == 0:
        return [""] * len(values)

    # 0.6745 rescales the MAD so the score reads on the same scale as a standard deviation
    return [round(0.6745 * (value - median) / mad, 2) for value in values]


def pool_rows(pool_name: str, summaries: dict, protein_ids: list[str]) -> list[dict]:
    """Explode one pool's summary confidences into one row per protein pair.

    Args:
        pool_name: Alphafold3 job name of the pool, e.g. "pool_0".
        summaries: Mapping with the pool's best-ranked file and its per-sample files.
        protein_ids: Protein ids of the pool, in the chain order used by Alphafold3.
    """
    samples = [load_json(path) for path in summaries["samples"]]

    if summaries["best"] is not None:
        best = load_json(summaries["best"])
    elif samples:
        # Fall back to the highest ranking_score sample, which is what the top-level file holds
        best = max(samples, key=lambda d: d.get("ranking_score", -float("inf")))
    else:
        raise ValueError(f"Pool {pool_name} has no summary confidence files.")

    chain_pair_iptm_best = best["chain_pair_iptm"]
    n_chains = len(chain_pair_iptm_best)

    if len(protein_ids) != n_chains:
        raise ValueError(
            f"Pool {pool_name} has {n_chains} chains in the Alphafold3 output but "
            f"{len(protein_ids)} protein ids in the pool file."
        )

    # Averaged over every sample, as in agg_chain_pair_iptms in the afpool paper
    iptm_matrices = [d["chain_pair_iptm"] for d in samples] or [chain_pair_iptm_best]
    chain_pair_iptm_mean = mean_matrix(iptm_matrices)

    pae_matrices = [d["chain_pair_pae_min"] for d in samples] or [
        best["chain_pair_pae_min"]
    ]

    chains = enumerate_chains()[:n_chains]
    chain_ptm = best.get("chain_ptm", [None] * n_chains)

    rows = []
    # combinations() walks the upper triangle of the matrices row by row, so the pair of chain
    # indices, the pair of protein ids and the pair of chain labels all stay in step
    pairs = zip(
        itertools.combinations(range(n_chains), 2),
        itertools.combinations(protein_ids, 2),
        itertools.combinations(chains, 2),
    )
    for (i, j), (protein_a, protein_b), (chain_a, chain_b) in pairs:
        rows.append(
            {
                "pool_name": pool_name,
                "protein_a": protein_a,
                "protein_b": protein_b,
                "chain_a": chain_a,
                "chain_b": chain_b,
                # Rounded to drop the float noise the averages pick up
                "pair_iptm": round(chain_pair_iptm_mean[i][j], 4),
                "pair_iptm_best": chain_pair_iptm_best[i][j],
                # Filled in below, once every pair of the pool is known
                "pair_iptm_z": "",
                "pool_iptm": best.get("iptm", ""),
                # The closest approach seen between the two chains in any sample. chain_pair_pae_min
                # is not symmetric, so both directions are considered
                "pae_min": min(
                    min(matrix[i][j], matrix[j][i]) for matrix in pae_matrices
                ),
                # The weaker of the two chains: a pair score cannot be trusted past the confidence
                # Alphafold3 has in the individual folds it is scoring
                "min_chain_ptm": min(chain_ptm[i], chain_ptm[j])
                if None not in (chain_ptm[i], chain_ptm[j])
                else "",
                "has_clash": best.get("has_clash", ""),
            }
        )

    # Rank each pair against the other pairs of its own pool, so that pools of different size and
    # composition stay comparable to one another
    for row, z in zip(rows, robust_z([row["pair_iptm"] for row in rows])):
        row["pair_iptm_z"] = z

    return rows


def process_pools(input_dir: str, pools_tsv: str, output_file: str) -> None:
    """Explode every pool under input_dir into per-pair scores and write them to a TSV file.

    Args:
        input_dir: Directory holding the Alphafold3 output for one or more pools.
        pools_tsv: pools.tsv written by afpool.py, used to name the chains of each pool.
        output_file: Output TSV filename.
    """
    pools = read_pools(pools_tsv)
    summaries = collect_summaries(input_dir)

    if not summaries:
        print(f"No {SUMMARY_SUFFIX} files found in {input_dir}")
        return

    print(f"Found {len(summaries)} pools in {input_dir}")

    data_rows = []
    for pool_name, pool_summaries in sorted(summaries.items()):
        protein_ids = pools.get(pool_name)
        if protein_ids is None:
            print(f"Skipping {pool_name}: no matching pool in {pools_tsv}")
            continue

        try:
            data_rows.extend(pool_rows(pool_name, pool_summaries, protein_ids))
        except (KeyError, ValueError, IOError, json.JSONDecodeError) as e:
            print(f"Error processing {pool_name}: {e}")
            continue

    if not data_rows:
        print("No valid pools processed")
        return

    # Sort on the within-pool score so the best candidate of every pool rises to the top, falling
    # back on the raw score for pools too uniform to give a spread
    data_rows.sort(
        key=lambda row: (
            row["pair_iptm_z"] if row["pair_iptm_z"] != "" else -float("inf"),
            row["pair_iptm"],
        ),
        reverse=True,
    )

    try:
        with open(output_file, "w") as f:
            # Write header
            f.write("\t".join(HEADERS) + "\n")

            # Write data rows
            for row in data_rows:
                f.write("\t".join(str(row[header]) for header in HEADERS) + "\n")

        print(f"Successfully wrote {len(data_rows)} pairs to {output_file}")
        print("Results sorted by pair_iptm_z (highest to lowest)")

    except IOError as e:
        print(f"Error writing output file {output_file}: {e}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Explode pooled Alphafold3 predictions into per-pair summary scores"
    )
    parser.add_argument(
        "--input-dir",
        "-i",
        default=".",
        help="Directory holding the Alphafold3 pool output (default: current directory)",
    )
    parser.add_argument(
        "--pools",
        "-p",
        default="pools.tsv",
        help="Location of pools.tsv file written by afpool.py (default: pools.tsv)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="pair_results.tsv",
        help="Output TSV filename (default: pair_results.tsv)",
    )

    args = parser.parse_args()

    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    if not os.path.isfile(args.pools):
        print(f"Error: Pool file '{args.pools}' does not exist")
        sys.exit(1)

    process_pools(args.input_dir, args.pools, args.output)


if __name__ == "__main__":
    main()
