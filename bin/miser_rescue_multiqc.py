#!/usr/bin/env python3
"""
miser_rescue_multiqc.py — build a MultiQC multi-dataset bargraph from the
merged MisER rescue-metrics TSV.

Instead of one bargraph with 13 columns on 13 different scales (read counts in
the millions next to a 0-100 pass rate next to exon sizes in bp), this writes a
single self-describing MultiQC custom-content file with a *dataset switcher*:
each view is one metric on its own scale, samples on the x-axis — so a sample
that behaved differently in MisER stands out instead of being crushed flat.

Input : the merged TSV written by MISER_QC_MERGE (one header row, one row per
        sample, tab-separated). Column names come from bin/miser_qc_summary.py.
Output: <name>_mqc.json — auto-detected by MultiQC as custom content.

Usage:
    python3 miser_rescue_multiqc.py \\
        --input  all_samples_rescue_metrics_mqc.tsv \\
        --output miser_rescue_mqc.json
"""

import argparse
import csv
import json
import sys

ID_COL = "sample"

# Switcher views: (display name, {output category label: source TSV column}).
# Each view is on its own scale. Source columns must exist in the merged TSV
# (see compute_metrics() in bin/miser_qc_summary.py).
VIEWS = [
    ("Pass rate (%)",      {"Pass rate (%)": "pass_rate_pct"}),
    ("Read events",        {"Passed": "passed_read_events", "Failed": "failed_read_events"}),
    ("Unique micro-exons", {"Unique micro-exons": "unique_microexons"}),
    ("Exon support",       {"Single-read": "exons_single_read", "Multi-read": "exons_multi_read"}),
    ("Median exon size",   {"Median exon size (bp)": "median_exon_size"}),
]

# Explicit, meaningful colour per category (hex). Avoids the flat/monochrome
# default: greens = good, reds = failed, blue shades = confidence tiers.
COLORS = {
    "Pass rate (%)":         "#2ca02c",  # green
    "Passed":                "#2ca02c",  # green
    "Failed":                "#d62728",  # red
    "Unique micro-exons":    "#9467bd",  # purple
    "Single-read":           "#aec7e8",  # light blue (lower confidence)
    "Multi-read":            "#1f77b4",  # dark blue  (higher confidence)
    "Median exon size (bp)": "#ff7f0e",  # orange
}


def as_number(cell):
    """Parse a numeric TSV cell; return None if empty or non-numeric."""
    if cell is None or cell == "":
        return None
    try:
        value = float(cell)
    except ValueError:
        return None
    return int(value) if value.is_integer() else value


def build_dataset(rows, mapping):
    """Return {sample: {category_label: value}} for one switcher view."""
    dataset = {}
    for row in rows:
        sample = row[ID_COL]
        values = {}
        for label, column in mapping.items():
            number = as_number(row.get(column))
            if number is not None:
                values[label] = number
        if values:
            dataset[sample] = values
    return dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="merged rescue-metrics TSV")
    parser.add_argument("--output", required=True, help="output *_mqc.json")
    args = parser.parse_args()

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        sys.exit(f"ERROR: no data rows in {args.input}")
    if ID_COL not in rows[0]:
        sys.exit(f"ERROR: '{ID_COL}' column not found in {args.input}. "
                 f"Columns present: {list(rows[0].keys())}")

    # Only keep views whose source columns are actually present, so a future
    # metric rename degrades gracefully instead of emitting empty datasets.
    available = set(rows[0].keys())
    views = [(name, m) for name, m in VIEWS if all(c in available for c in m.values())]
    if not views:
        sys.exit(f"ERROR: none of the expected metric columns are present in {args.input}")

    data = [build_dataset(rows, mapping) for _, mapping in views]
    data_labels = [{"name": name, "ylab": name} for name, _ in views]

    # Colour every category that appears in any view. MultiQC keys bar colours
    # on category name, so one map covers all datasets.
    categories = {cat for _, mapping in views for cat in mapping}
    colours = {cat: COLORS[cat] for cat in categories if cat in COLORS}

    report = {
        "id": "miser_rescue",
        "section_name": "MisER Micro-exon Rescue",
        "description": (
            "Per-sample MisER splice-correction QC. Use the dropdown above the plot "
            "to switch metric; each view is on its own scale, so an outlier sample "
            "(e.g. low pass rate or few micro-exons) is easy to spot."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "miser_rescue_bargraph",
            "title": "MisER: per-sample rescue metrics",
            "cpswitch": False,          # metrics are on different scales; no counts/% toggle
            "data_labels": data_labels,  # the switcher
            "colors": colours,           # per-category colours (see note in module docstring)
        },
        "data": data,                    # list of datasets, one per data_label
    }

    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)

    print(f"Wrote {args.output}: {len(rows)} samples, {len(views)} metric views "
          f"({', '.join(name for name, _ in views)})")


if __name__ == "__main__":
    main()
