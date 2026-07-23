#!/usr/bin/env python3
"""
miser_rescue_multiqc.py — split the merged MisER rescue-metrics TSV into
scale-coherent, colour-able MultiQC custom-content tables.

Why not one switcher: MultiQC 1.25.x only applies custom colours to
config-defined custom_content sections (the mechanism the isoquant_assignment
block in multiqc_config.yml uses). It ignores colours on multi-dataset switcher
JSON. So instead of one dropdown we emit one small TSV per metric group, and
multiqc_config.yml defines a coloured bargraph section for each — each on its
own scale, samples on the x-axis.

Input : merged rescue-metrics TSV from MISER_QC_MERGE. Column names come from
        bin/miser_qc_summary.py (compute_metrics).
Output (in --outdir): miser_passrate_mqc.tsv, miser_events_mqc.tsv,
        miser_support_mqc.tsv, miser_exonsize_mqc.tsv

Usage:
    python3 miser_rescue_multiqc.py \\
        --input  all_samples_rescue_metrics_mqc.tsv \\
        --outdir .
"""

import argparse
import csv
import os
import sys

ID_COL = "sample"

# (output filename, [(TSV column header, source column in merged TSV)])
# Each group is one coloured bargraph section; columns within a group share a
# scale so they can be grouped/stacked sensibly. Colours are set per column in
# multiqc_config.yml (headers: <column>: colour:).
GROUPS = [
    ("miser_passrate_mqc.tsv", [("pass_rate",   "pass_rate_pct")]),
    ("miser_events_mqc.tsv",   [("passed",      "passed_read_events"),
                                ("failed",      "failed_read_events")]),
    ("miser_support_mqc.tsv",  [("single_read", "exons_single_read"),
                                ("multi_read",  "exons_multi_read")]),
    ("miser_exonsize_mqc.tsv", [("median_size", "median_exon_size")]),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="merged rescue-metrics TSV")
    parser.add_argument("--outdir", default=".",   help="where to write the *_mqc.tsv files")
    args = parser.parse_args()

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        sys.exit(f"ERROR: no data rows in {args.input}")
    if ID_COL not in rows[0]:
        sys.exit(f"ERROR: '{ID_COL}' column not found in {args.input}. "
                 f"Columns present: {list(rows[0].keys())}")

    os.makedirs(args.outdir, exist_ok=True)
    available = set(rows[0].keys())
    written = []

    for fname, columns in GROUPS:
        # keep only columns whose source exists (graceful if a metric is renamed)
        columns = [(header, source) for header, source in columns if source in available]
        if not columns:
            continue
        path = os.path.join(args.outdir, fname)
        with open(path, "w") as out:
            out.write("Sample\t" + "\t".join(h for h, _ in columns) + "\n")
            for row in rows:
                values = [row.get(source, "") for _, source in columns]
                out.write(row[ID_COL] + "\t" + "\t".join(values) + "\n")
        written.append(fname)

    if not written:
        sys.exit(f"ERROR: none of the expected metric columns are present in {args.input}")

    print(f"Wrote {len(written)} MisER MultiQC tables: {', '.join(written)}")


if __name__ == "__main__":
    main()
