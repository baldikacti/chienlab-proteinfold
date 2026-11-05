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

Output is one row per protein pair, sorted by pair_iptm in descending order. With --bait, only pairs
that include the bait protein (the first protein id of the first pool in pools.tsv) are kept. A
self-contained HTML report (<output>_report.html) summarising the TSV files is written alongside.
"""

import csv
import glob
import html
import itertools
import json
import math
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

# Thresholds a pair must exceed to count as a potential interaction
IPTM_THRESHOLD = 0.5
Z_THRESHOLD = 1.5


def enumerate_chains() -> list[str]:
    """Generate the chain labels A to ZZ that afpool.py assigns to pooled proteins, in order."""
    # Copied from export_pool in afpool.py so the labels stay in step with the exported pools
    return list(string.ascii_uppercase) + [
        a + b for a in string.ascii_uppercase for b in string.ascii_uppercase
    ]


def read_pools(pools_tsv: str) -> tuple[dict[str, list[str]], str]:
    """Read afpool.py's pools.tsv and return a lookup of pool name to list of protein ids.

    Args:
        pools_tsv: Path to the pools.tsv written by afpool.py.

    Pools are keyed both by their bare pool_id ("0") and by the job name afpool.py gives them
    ("pool0"), so either naming of the Alphafold3 output resolves. Also returns the first protein
    id of the first pool, which is the bait in bait vs all runs.
    """
    pools = {}
    first_protein = None
    with open(pools_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pool_id = row["pool_id"].strip()
            protein_ids = row["protein_ids"].split("_")
            pools[pool_id] = protein_ids
            pools["pool" + pool_id] = protein_ids
            if first_protein is None:
                first_protein = protein_ids[0]

    if not pools:
        raise ValueError(f"Pool file {pools_tsv} contains no pools.")

    return pools, first_protein


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
        pool_name: Alphafold3 job name of the pool, e.g. "pool0".
        summaries: Mapping with the pool's best-ranked file and its per-sample files.
        protein_ids: Protein ids of the pool, in the chain order used by Alphafold3.
    """
    samples = [load_json(path) for path in summaries["samples"]]
    best = load_json(summaries["best"])

    chain_pair_iptm_best = best["chain_pair_iptm"]
    n_chains = len(chain_pair_iptm_best)

    if len(protein_ids) != n_chains:
        raise ValueError(
            f"Pool {pool_name} has {n_chains} chains in the Alphafold3 output but "
            f"{len(protein_ids)} protein ids in the pool file."
        )

    # Averaged over every sample, as in agg_chain_pair_iptms in the afpool paper. Without sample files
    # the mean falls back to the best-ranked model's value
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


def write_thresholded_pairs(
    data_rows: list[dict],
    output_file: str,
    iptm_threshold: float = IPTM_THRESHOLD,
    z_threshold: float = Z_THRESHOLD,
) -> None:
    """Write potential_pairs.tsv and inconclusive_pairs.tsv next to output_file.

    Args:
        data_rows: Per-pair rows, already sorted.
        output_file: Path of the full results TSV.
        iptm_threshold: pair_iptm a potential pair must exceed.
        z_threshold: pair_iptm_z a potential pair must exceed.

    Rows are thresholded first. Pairs from models with a clash are inconclusive. Pairs without a
    pair_iptm_z are dropped. Within each output, a pair found in several pools is then reduced to
    its row with the highest pair_iptm. A pair with a potential row is left out of inconclusive.
    """

    def pair_key(row: dict) -> tuple:
        return tuple(sorted((row["protein_a"], row["protein_b"])))

    def best_per_pair(rows: list[dict]) -> list[dict]:
        best = {}
        for row in rows:
            key = pair_key(row)
            if key not in best or row["pair_iptm"] > best[key]["pair_iptm"]:
                best[key] = row
        return sorted(best.values(), key=lambda row: row["pair_iptm"], reverse=True)

    potential = best_per_pair([
        row
        for row in data_rows
        if not row["has_clash"]
        and row["pair_iptm"] > iptm_threshold
        and row["pair_iptm_z"] != ""
        and row["pair_iptm_z"] > z_threshold
    ])
    potential_keys = {pair_key(row) for row in potential}
    inconclusive = best_per_pair([
        row
        for row in data_rows
        if row["has_clash"] and pair_key(row) not in potential_keys
    ])

    output_dir = os.path.dirname(output_file)
    for name, rows in (
        ("potential_pairs.tsv", potential),
        ("inconclusive_pairs.tsv", inconclusive),
    ):
        path = os.path.join(output_dir, name)
        with open(path, "w") as f:
            f.write("\t".join(HEADERS) + "\n")
            for row in rows:
                f.write("\t".join(str(row[header]) for header in HEADERS) + "\n")
        print(f"Successfully wrote {len(rows)} pairs to {path}")


# --- HTML report -----------------------------------------------------------------------------------

STATUS_LABELS = {
    "pass": "Passes both thresholds",
    "clash": "Clash (inconclusive)",
    "below": "Below threshold",
}

TABLE_COLUMNS = [
    "pool_name",
    "pair_iptm",
    "pair_iptm_best",
    "pair_iptm_z",
    "pae_min",
    "min_chain_ptm",
    "pool_iptm",
]

CHART_W = 720
CHART_H = 340
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 52, 20, 16, 56

REPORT_CSS = """
:root { color-scheme: light; --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink2: #52514e;
  --muted: #898781; --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,.10);
  --pass: #2a78d6; --clash: #eb6834; --below: #b5b3ab; }
@media (prefers-color-scheme: dark) {
  :root { color-scheme: dark; --surface: #1a1a19; --page: #0d0d0d; --ink: #f0efec; --ink2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,.10);
    --pass: #3987e5; --clash: #d95926; --below: #5f5e5a; }
}
body { margin: 0; background: var(--page); color: var(--ink);
  font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 800px; margin: 0 auto; padding: 32px 20px; }
h1 { font-size: 22px; margin: 0 0 2px; }
h2 { font-size: 16px; margin: 0 0 2px; }
.sub { color: var(--ink2); margin: 0 0 10px; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
  padding: 20px; margin: 16px 0; }
.tiles { display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px; }
.tile .label { color: var(--ink2); font-size: 13px; }
.tile .value { font-size: 24px; font-weight: 600; }
.legend { display: flex; flex-wrap: wrap; gap: 16px; color: var(--ink2); font-size: 13px; margin: 0 0 6px; }
.key { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; }
.key.pass { background: var(--pass); } .key.clash { background: var(--clash); }
.key.below { background: var(--below); }
svg { width: 100%; height: auto; display: block; }
svg text { fill: var(--muted); font-size: 12px; }
.grid { stroke: var(--grid); } .axis { stroke: var(--axis); }
.thr { stroke: var(--ink2); stroke-width: 1.5; stroke-dasharray: 4 4; }
svg text.thr-label, svg text.direct { fill: var(--ink2); }
.zone { fill: var(--pass); opacity: .07; }
.dot { stroke: var(--surface); stroke-width: 2; }
.pt.pass .dot { fill: var(--pass); } .pt.clash .dot { fill: var(--clash); }
.pt.below .dot { fill: var(--below); }
.hit { fill: transparent; }
.pt:hover .dot, .pt:focus .dot { stroke: var(--ink); }
.pt:focus { outline: none; }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
th, td { text-align: right; padding: 6px 8px; border-bottom: 1px solid var(--grid); }
th:first-child, td:first-child { text-align: left; }
th { color: var(--ink2); font-weight: 600; }
details summary { cursor: pointer; color: var(--ink2); }
details summary .note { display: block; font-size: 13px; margin: 2px 0 0 18px; }
.help { margin: 0; padding-left: 18px; } .help li { margin: 3px 0; }
#tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink);
  border: 1px solid var(--border); border-radius: 6px; padding: 6px 8px; font-size: 12px;
  white-space: pre; box-shadow: 0 2px 8px rgba(0,0,0,.15); }
"""

REPORT_JS = """
const tip = document.getElementById("tip");
document.querySelectorAll(".pt").forEach((el) => {
  const show = (x, y) => {
    tip.textContent = el.dataset.tip;
    tip.hidden = false;
    const w = tip.offsetWidth;
    tip.style.left = (x + 14 + w > innerWidth ? x - 14 - w : x + 14) + "px";
    tip.style.top = y + 14 + "px";
  };
  el.addEventListener("mousemove", (e) => show(e.clientX, e.clientY));
  el.addEventListener("mouseleave", () => { tip.hidden = true; });
  el.addEventListener("focus", () => { const r = el.getBoundingClientRect(); show(r.right, r.top); });
  el.addEventListener("blur", () => { tip.hidden = true; });
});
"""


def _number(value):
    """Parse a TSV cell to a float, or None when the cell is empty."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _is_clash(value) -> bool:
    return str(value).strip().lower() in ("1", "1.0", "true")


def _status(row: dict, iptm_threshold: float, z_threshold: float) -> str:
    """Classify a row the same way write_thresholded_pairs does."""
    if _is_clash(row["has_clash"]):
        return "clash"
    iptm, z = _number(row["pair_iptm"]), _number(row["pair_iptm_z"])
    if iptm is not None and z is not None and iptm > iptm_threshold and z > z_threshold:
        return "pass"
    return "below"


def _partner(row: dict, bait_id) -> str:
    """Name the non-bait protein of a pair, or both proteins when there is no bait."""
    a, b = row["protein_a"], row["protein_b"]
    if bait_id in (a, b):
        return b if a == bait_id else a
    return f"{a} – {b}"


def _ticks(low: float, high: float, count: int = 5) -> list[float]:
    """Round tick values (steps of 1, 2, 2.5 or 5 x 10^n) covering low to high."""
    if high <= low:
        high = low + 1
    raw = (high - low) / count
    magnitude = 10 ** math.floor(math.log10(raw))
    step = next(m * magnitude for m in (1, 2, 2.5, 5, 10) if m * magnitude >= raw)
    start = math.floor(low / step) * step
    end = math.ceil(high / step) * step
    return [round(start + i * step, 10) for i in range(round((end - start) / step) + 1)]


def _point(x: float, y: float, status: str, tip: str) -> str:
    """A hoverable dot: a visible marker plus a larger transparent hit target."""
    tip = html.escape(tip, quote=True).replace("\n", "&#10;")
    return (
        f'<g class="pt {status}" tabindex="0" data-tip="{tip}" aria-label="{tip}">'
        f'<circle class="hit" cx="{x:.1f}" cy="{y:.1f}" r="11"/>'
        f'<circle class="dot" cx="{x:.1f}" cy="{y:.1f}" r="5"/></g>'
    )


def _tip(row: dict, bait_id, status: str) -> str:
    return (
        f"{_partner(row, bait_id)}\npool: {row['pool_name']}\npair_iptm: {row['pair_iptm']}\n"
        f"pair_iptm_z: {row['pair_iptm_z'] or 'n/a'}\npae_min: {row['pae_min']}\n"
        f"{STATUS_LABELS[status]}"
    )


def _iptm_axis(height: float, label: str = "pair_iptm"):
    """Gridlines and labels for an ipTM y-axis fixed at 0 to 1. Returns (svg parts, y scale)."""
    top, bottom = MARGIN_T, height - MARGIN_B

    def y_of(value: float) -> float:
        return bottom - min(max(value, 0), 1) * (bottom - top)

    parts = []
    for t in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
        y = y_of(t)
        parts.append(
            f'<line class="{"axis" if t == 0 else "grid"}" x1="{MARGIN_L}" '
            f'x2="{CHART_W - MARGIN_R}" y1="{y:.1f}" y2="{y:.1f}"/>'
            f'<text x="{MARGIN_L - 8}" y="{y + 4:.1f}" text-anchor="end">{t:g}</text>'
        )
    parts.append(
        f'<text transform="translate(14 {(top + bottom) / 2:.1f}) rotate(-90)" '
        f'text-anchor="middle">{label}</text>'
    )
    return parts, y_of


def _iptm_threshold_line(y: float, iptm_threshold: float) -> str:
    return (
        f'<line class="thr" x1="{MARGIN_L}" x2="{CHART_W - MARGIN_R}" y1="{y:.1f}" y2="{y:.1f}"/>'
        f'<text class="thr-label" x="{CHART_W - MARGIN_R}" y="{y - 6:.1f}" '
        f'text-anchor="end">ipTM {iptm_threshold:g}</text>'
    )


def _svg(height: float, parts: list[str], label: str) -> str:
    return (
        f'<svg viewBox="0 0 {CHART_W} {height:.0f}" role="img" aria-label="{html.escape(label)}">'
        + "".join(parts)
        + "</svg>"
    )


def _top_chart(rows: list[dict], bait_id, iptm_threshold: float, z_threshold: float) -> str:
    """Scatter of the 10 proteins with the highest pair_iptm, one dot per protein.

    Rows from models with a clash are uninformative and left out.
    """
    best = {}
    for row in sorted(rows, key=lambda row: _number(row["pair_iptm"]) or 0, reverse=True):
        if not _is_clash(row["has_clash"]):
            best.setdefault(_partner(row, bait_id), row)
    top = list(best.items())[:10]

    plot_w = CHART_W - MARGIN_L - MARGIN_R
    band = plot_w / max(len(top), 1)
    # Rotate the protein names when they would not fit side by side
    rotate = any(len(name) * 7 > band - 6 for name, _ in top)
    height = CHART_H + (50 if rotate else 0)
    baseline = height - MARGIN_B

    parts, y_of = _iptm_axis(height)
    parts.append(_iptm_threshold_line(y_of(iptm_threshold), iptm_threshold))
    for i, (name, row) in enumerate(top):
        x = MARGIN_L + band * (i + 0.5)
        name = html.escape(name)
        if rotate:
            parts.append(
                f'<text text-anchor="end" transform="translate({x + 4:.1f} {baseline + 14}) '
                f'rotate(-35)">{name}</text>'
            )
        else:
            parts.append(f'<text x="{x:.1f}" y="{baseline + 18}" text-anchor="middle">{name}</text>')
        status = _status(row, iptm_threshold, z_threshold)
        parts.append(_point(x, y_of(_number(row["pair_iptm"]) or 0), status, _tip(row, bait_id, status)))

    return _svg(height, parts, "Top 10 proteins by pair_iptm")


def _scatter_chart(
    rows: list[dict],
    x_key: str,
    x_label: str,
    bait_id,
    iptm_threshold: float,
    z_threshold: float,
    x_threshold: float = None,
    show_iptm_threshold: bool = True,
) -> str:
    """Scatter of pair_iptm against the x_key column, one dot per row.

    Rows with an empty x_key are skipped. With x_threshold, a dotted line is drawn there and the
    area passing both thresholds is shaded. The 10 highest passing dots are labelled by name.
    """
    points = [
        (_number(row[x_key]), _number(row["pair_iptm"]) or 0, row, _status(row, iptm_threshold, z_threshold))
        for row in rows
        if _number(row[x_key]) is not None
    ]
    xs = [x for x, _, _, _ in points] + [0] + ([x_threshold] if x_threshold is not None else [])
    x_ticks = _ticks(min(xs), max(xs))
    plot_w = CHART_W - MARGIN_L - MARGIN_R

    def x_of(value: float) -> float:
        return MARGIN_L + (value - x_ticks[0]) / (x_ticks[-1] - x_ticks[0]) * plot_w

    height = CHART_H
    baseline = height - MARGIN_B
    parts, y_of = _iptm_axis(height)
    iy = y_of(iptm_threshold)

    for t in x_ticks:
        x = x_of(t)
        if t != x_ticks[0]:
            parts.append(f'<line class="grid" x1="{x:.1f}" x2="{x:.1f}" y1="{MARGIN_T}" y2="{baseline}"/>')
        parts.append(f'<text x="{x:.1f}" y="{baseline + 18}" text-anchor="middle">{t:g}</text>')
    parts.append(
        f'<text x="{MARGIN_L + plot_w / 2:.1f}" y="{height - 10}" text-anchor="middle">{x_label}</text>'
    )
    if show_iptm_threshold:
        parts.append(_iptm_threshold_line(iy, iptm_threshold))

    if x_threshold is not None:
        # Shade the region a potential pair must fall in
        tx = x_of(x_threshold)
        parts.insert(
            0,
            f'<rect class="zone" x="{tx:.1f}" y="{MARGIN_T}" width="{CHART_W - MARGIN_R - tx:.1f}" '
            f'height="{iy - MARGIN_T:.1f}"/>',
        )
        parts.append(
            f'<line class="thr" x1="{tx:.1f}" x2="{tx:.1f}" y1="{MARGIN_T}" y2="{baseline}"/>'
            f'<text class="thr-label" x="{tx + 6:.1f}" y="{MARGIN_T + 12}">z {x_threshold:g}</text>'
        )

    # Draw passing dots last, highest ipTM on top, and label the 10 highest by name
    points.sort(key=lambda p: (p[3] == "pass", p[1]))
    labelled = {id(p[2]) for p in [p for p in points if p[3] == "pass"][-10:]}
    for xv, iptm, row, status in points:
        x, y = x_of(xv), y_of(iptm)
        parts.append(_point(x, y, status, _tip(row, bait_id, status)))
        if id(row) in labelled:
            anchor, dx = ("end", -9) if x > CHART_W - 140 else ("start", 9)
            parts.append(
                f'<text class="direct" x="{x + dx:.1f}" y="{y + 4:.1f}" text-anchor="{anchor}">'
                f"{html.escape(_partner(row, bait_id))}</text>"
            )

    return _svg(height, parts, f"pair_iptm against {x_key}")


def _table(rows: list[dict], bait_id, statuses: list = None) -> str:
    """HTML table of rows with the bait column removed."""
    first = "Protein" if bait_id else "Pair"
    header = [first] + TABLE_COLUMNS + (["status"] if statuses else [])
    lines = ["<table><thead><tr>" + "".join(f"<th>{h}</th>" for h in header) + "</tr></thead><tbody>"]
    for i, row in enumerate(rows):
        cells = [_partner(row, bait_id)] + [row[c] for c in TABLE_COLUMNS]
        if statuses:
            cells.append(STATUS_LABELS[statuses[i]])
        lines.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in cells) + "</tr>")
    lines.append("</tbody></table>")
    return "".join(lines)


def _legend(statuses=tuple(STATUS_LABELS)) -> str:
    return '<div class="legend">' + "".join(
        f'<span><span class="key {s}"></span>{STATUS_LABELS[s]}</span>' for s in statuses
    ) + "</div>"


def write_report(
    output_file: str,
    bait_id=None,
    iptm_threshold: float = IPTM_THRESHOLD,
    z_threshold: float = Z_THRESHOLD,
) -> None:
    """Write <output>_report.html from the full, potential and inconclusive TSV files.

    Args:
        output_file: Path of the full results TSV; the thresholded TSVs are read from beside it.
        bait_id: Bait protein, left out of pair names. None for all vs all runs.
        iptm_threshold: pair_iptm threshold used for potential pairs.
        z_threshold: pair_iptm_z threshold used for potential pairs.
    """

    def read(path: str) -> list[dict]:
        with open(path, newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    output_dir = os.path.dirname(output_file)
    rows = read(output_file)
    potential = read(os.path.join(output_dir, "potential_pairs.tsv"))
    inconclusive = read(os.path.join(output_dir, "inconclusive_pairs.tsv"))

    n_pools = len({row["pool_name"] for row in rows})
    n_proteins = len({_partner(row, bait_id) for row in rows})
    # Models with a clash are uninformative, so the scatter plots leave them out
    clean = [row for row in rows if not _is_clash(row["has_clash"])]
    n_no_z = sum(1 for row in clean if _number(row["pair_iptm_z"]) is None)
    statuses = [_status(row, iptm_threshold, z_threshold) for row in rows]
    one = "protein" if bait_id else "pair"
    noun = "proteins" if bait_id else "pairs"
    rule = f"no clash, pair_iptm &gt; {iptm_threshold:g} and pair_iptm_z &gt; {z_threshold:g}"

    tiles = "".join(
        f'<div class="tile"><div class="label">{label}</div><div class="value">{value}</div></div>'
        for label, value in (
            ("Pools", n_pools),
            (f"{noun.capitalize()} screened", n_proteins),
            ("Potential pairs", len(potential)),
            ("Clashing pairs", len(inconclusive)),
        )
    )
    no_z_note = (
        f" {n_no_z} rows without a z-score (no spread within their pool) are not shown."
        if n_no_z
        else ""
    )
    potential_table = (
        _table(potential, bait_id) if potential else '<p class="sub">No pairs passed both thresholds.</p>'
    )
    clash_table = _table(inconclusive, bait_id) if inconclusive else '<p class="sub">No clashing pairs.</p>'
    potential_pae_chart = (
        _scatter_chart(
            potential, "pae_min", "pae_min (Å)", bait_id, iptm_threshold, z_threshold,
            show_iptm_threshold=False,
        )
        if potential
        else '<p class="sub">No pairs passed both thresholds.</p>'
    )
    subject = f"Bait {html.escape(bait_id)}" if bait_id else "All vs all"

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pool analysis report</title><style>{REPORT_CSS}</style></head>
<body><main>
<h1>Pool analysis report</h1>
<p class="sub">{subject} · Alphafold3 pooled predictions</p>
<div class="tiles">{tiles}</div>

<section class="card">
<h2>Metric definitions</h2>
<ul class="help">
<li><b>pair_iptm</b>: Confidence that the two proteins interact based on Alphafold3 predictions (0 to 1), averaged over samples.
Higher means a more likely interaction. <b>This is the most important metric for identifying potential interactions.</b></li>
<li><b>pair_iptm_best</b> is the highest achieved ipTM score from any Alphafold3 model.</li>
<li><b>pair_iptm_z</b>: z-score, or standard deviation, of ipTM compared to other protein pairs in the same pool.
This is used to identify and threshold out proteins in pools with artificially high ipTM.
Apart from thresholding, you can ignore this metric - the actual values are only relevant in relation to the threshold.</li>
<li><b>pae_min</b>: Lowest predicted aligned error between the two chains (Å). Lower means higher confidence in the predicted interaction.
<b>This is a second useful metric for interaction.</b></li>
<li><b>Potential pair</b>: {rule}. Proteins that pass a very loose ipTM threshold, showing a POTENTIAL interaction.</li>
<li><b>Clash</b>: Flag for when the generated structure is low-confidence and unreliable. If a protein is flagged, it means we cannot make conclusions about it.</li>
<li><b>Files</b>: <code>pool_analysis.tsv</code> has every pair in every pool.
<code>potential_pairs.tsv</code> has potential pairs that passed initial thresholding.
<code>inconclusive_pairs.tsv</code> has clash pairs that never passed in another pool.</li>
</ul>
</section>

<section class="card">
<h2>Top 10 {noun} by ipTM</h2>
<p class="sub">Highest pair_iptm across pools for each {"protein paired with the bait" if bait_id else "pair"}.
Models with a clash are left out. Dotted line: ipTM threshold.</p>
{_legend(("pass", "below"))}
{_top_chart(rows, bait_id, iptm_threshold, z_threshold)}
</section>

<section class="card">
<h2>ipTM against z-score thresholds</h2>
<p class="sub">Each dot is one {one} in one pool. Models with a clash are left out. Passing dots fall
in the shaded area ({rule}).{no_z_note}</p>
{_legend(("pass", "below"))}
{_scatter_chart(clean, "pair_iptm_z", "pair_iptm_z", bait_id, iptm_threshold, z_threshold, z_threshold)}
</section>

<section class="card">
<h2>ipTM against minimum PAE</h2>
<p class="sub">Each dot is one {one} in one pool. Models with a clash are left out. Higher ipTM and
lower pae_min (upper left) mean a more confident interaction. Dotted line: ipTM threshold.</p>
{_legend(("pass", "below"))}
{_scatter_chart(clean, "pae_min", "pae_min (Å)", bait_id, iptm_threshold, z_threshold)}
</section>

<section class="card">
<h2>ipTM against minimum PAE, potential pairs</h2>
<p class="sub">Each dot is one {one} from <code>potential_pairs.tsv</code>, at its best pool. Higher ipTM
and lower pae_min (upper left) mean a more confident interaction.</p>
{potential_pae_chart}
</section>

<section class="card">
<h2>Potential interacting pairs</h2>
<p class="sub">From <code>potential_pairs.tsv</code>, sorted by pair_iptm.</p>
{potential_table}
</section>

<section class="card">
<details><summary>Clashing pairs ({len(inconclusive)} pair{"" if len(inconclusive) == 1 else "s"}, from <code>inconclusive_pairs.tsv</code>)
<span class="note">Proteins for which we cannot make conclusions about interaction with the bait.
</span></summary>
{clash_table}
</details>
</section>

<section class="card">
<details><summary>All pairs ({len(rows)} rows, from <code>{html.escape(os.path.basename(output_file))}</code>)</summary>
{_table(rows, bait_id, statuses)}
</details>
</section>
</main>
<div id="tip" hidden></div>
<script>{REPORT_JS}</script>
</body></html>
"""

    path = os.path.splitext(output_file)[0] + "_report.html"
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"Successfully wrote report to {path}")


def process_pools(
    input_dir: str, pools_tsv: str, output_file: str, bait: bool = False
) -> None:
    """Explode every pool under input_dir into per-pair scores and write them to a TSV file.

    Args:
        input_dir: Directory holding the Alphafold3 output for one or more pools.
        pools_tsv: pools.tsv written by afpool.py, used to name the chains of each pool.
        output_file: Output TSV filename, without the .tsv extension.
        bait: Only keep pairs including the bait, taken as the first protein id of the first pool.
    """
    output_file = output_file + ".tsv"
    pools, bait_id = read_pools(pools_tsv)
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

    if bait:
        print(f"Keeping only pairs including bait {bait_id}")
        data_rows = [
            row for row in data_rows if bait_id in (row["protein_a"], row["protein_b"])
        ]

    if not data_rows:
        print("No valid pools processed")
        return

    data_rows.sort(key=lambda row: row["pair_iptm"], reverse=True)

    try:
        with open(output_file, "w") as f:
            # Write header
            f.write("\t".join(HEADERS) + "\n")

            # Write data rows
            for row in data_rows:
                f.write("\t".join(str(row[header]) for header in HEADERS) + "\n")

        print(f"Successfully wrote {len(data_rows)} pairs to {output_file}")
        print("Results sorted by pair_iptm (highest to lowest)")

        write_thresholded_pairs(data_rows, output_file)
        write_report(output_file, bait_id if bait else None)

    except IOError as e:
        print(f"Error writing output file {output_file}: {e}")


def main():
    # Note: currently assumes 1 pool passed in. Does not work for multiple pools; we can change this later if needed.
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
        default="pair_results",
        help="Output TSV filename without the .tsv extension (default: pair_results)",
    )
    # Optional argument for JUST bait vs. all mode
    parser.add_argument(
        "--bait",
        action="store_true",
        help="Only keep pairs including bait protein (the first protein id of the first pool)",
    )

    args = parser.parse_args()

    # Check if input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    # Check if pools.tsv exists
    if not os.path.isfile(args.pools):
        print(f"Error: Pool file '{args.pools}' does not exist")
        sys.exit(1)
    # Check that every pool has its main (best-ranked) summary confidence file
    for pool_name, pool_summaries in collect_summaries(args.input_dir).items():
        if pool_summaries["best"] is None:
            print(f"Error: Pool {pool_name} has no {pool_name}{SUMMARY_SUFFIX} file")
            sys.exit(1)

    process_pools(args.input_dir, args.pools, args.output, args.bait)


if __name__ == "__main__":
    main()
