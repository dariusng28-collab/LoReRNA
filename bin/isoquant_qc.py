#!/usr/bin/env python3
"""
isoquant_qc.py — per-sample IsoQuant read assignment QC

Parses IsoQuant's read_assignments.tsv.gz to produce:
  1. {sample}_isoquant_qc.tsv     — human-readable per-sample stats
  2. {sample}_isoquant_mqc.tsv    — MultiQC custom_content bargraph file

IsoQuant read_assignments.tsv.gz column layout (v3.x):
  #read_id | chr | strand | isoform_id | gene_id | assignment_type | assignment_events | exons

Assignment types (v3.x):
  unique                    — reads assigned to exactly one isoform
  unique_minor_difference   — unique with minor alignment differences
  ambiguous                 — consistent with multiple isoforms
  inconsistent              — mapped but inconsistent with annotation
  not_assigned              — not assigned (includes noninformative)

Usage:
  python3 isoquant_qc.py \\
      --input  sample.read_assignments.tsv.gz \\
      --sample sample_id \\
      --outdir .
"""

import argparse
import gzip
import os
import sys
from collections import Counter

ASSIGNMENT_ORDER = [
    "unique",
    "unique_minor_difference",
    "ambiguous",
    "inconsistent",
    "not_assigned",
]

# IsoQuant uses slightly different names in some versions — normalise them
NORMALISE = {
    "noninformative": "not_assigned",
    "not_assigned":   "not_assigned",
    "unique":         "unique",
    "unique_minor_difference": "unique_minor_difference",
    "ambiguous":      "ambiguous",
    "inconsistent":   "inconsistent",
    "supplementary":  "not_assigned",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True,
                   help="IsoQuant read_assignments.tsv.gz file")
    p.add_argument("--sample", default=None,
                   help="Sample ID label (default: derived from filename)")
    p.add_argument("--outdir", default=".",
                   help="Output directory (default: current dir)")
    return p.parse_args()


def parse_assignments(tsv_gz: str) -> tuple[Counter, int]:
    """Return (counts_by_type, total_reads)."""
    counts: Counter = Counter()
    total = 0

    opener = gzip.open if tsv_gz.endswith(".gz") else open

    with opener(tsv_gz, "rt") as fh:
        header_seen = False
        col_idx = None

        for line in fh:
            line = line.rstrip("\n")

            # Skip comment / header lines
            if line.startswith("#"):
                # Try to detect column index from header comment
                if "assignment_type" in line:
                    fields = line.lstrip("#").split("\t")
                    try:
                        col_idx = fields.index("assignment_type")
                    except ValueError:
                        pass
                continue

            # First non-comment line might be a header without '#'
            if not header_seen:
                fields = line.split("\t")
                if "assignment_type" in fields:
                    col_idx = fields.index("assignment_type")
                    header_seen = True
                    continue
                # No header line — use default column position (index 5)
                header_seen = True
                if col_idx is None:
                    col_idx = 5

            fields = line.split("\t")
            if len(fields) <= col_idx:
                continue

            raw_type = fields[col_idx].strip().lower()
            norm_type = NORMALISE.get(raw_type, "not_assigned")
            counts[norm_type] += 1
            total += 1

    return counts, total


def write_qc_tsv(sample: str, counts: Counter, total: int, outdir: str):
    """Human-readable per-sample QC TSV."""
    out_path = os.path.join(outdir, f"{sample}_isoquant_qc.tsv")
    with open(out_path, "w") as fh:
        fh.write("sample\tassignment_type\tcount\tpercent\n")
        for atype in ASSIGNMENT_ORDER:
            n = counts.get(atype, 0)
            pct = (n / total * 100) if total > 0 else 0.0
            fh.write(f"{sample}\t{atype}\t{n}\t{pct:.2f}\n")
        fh.write(f"{sample}\ttotal\t{total}\t100.00\n")
    print(f"  Written: {out_path}")
    return out_path


def write_multiqc_tsv(sample: str, counts: Counter, outdir: str):
    """
    MultiQC custom_content bargraph file.

    Format expected by multiqc_config.yml section 'isoquant_assignment':
      Header row: Sample <TAB> unique <TAB> unique_minor_difference <TAB> ...
      Data row:   <sample_id> <TAB> <count> <TAB> ...
    """
    out_path = os.path.join(outdir, f"{sample}_isoquant_mqc.tsv")
    with open(out_path, "w") as fh:
        # MultiQC custom content header
        fh.write("# id: 'isoquant_assignment'\n")
        fh.write("# section_name: 'IsoQuant Read Assignments'\n")
        fh.write("# description: 'Read assignment type breakdown from IsoQuant'\n")
        fh.write("# plot_type: 'bargraph'\n")
        fh.write("# pconfig:\n")
        fh.write("#   id: 'isoquant_assignment_plot'\n")
        fh.write("#   title: 'IsoQuant: Read Assignment Types'\n")
        fh.write("#   ylab: 'Number of reads'\n")
        # Column header
        fh.write("Sample\t" + "\t".join(ASSIGNMENT_ORDER) + "\n")
        # Data row
        row = [sample] + [str(counts.get(t, 0)) for t in ASSIGNMENT_ORDER]
        fh.write("\t".join(row) + "\n")
    print(f"  Written: {out_path}")
    return out_path


def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # Derive sample ID from filename if not provided
    sample = args.sample
    if sample is None:
        basename = os.path.basename(args.input)
        # Strip known suffixes
        for suffix in (".read_assignments.tsv.gz", ".tsv.gz", ".tsv"):
            if basename.endswith(suffix):
                sample = basename[: -len(suffix)]
                break
        if sample is None:
            sample = basename.split(".")[0]

    print(f"IsoQuant QC: {sample}")
    print(f"  Input: {args.input}")

    counts, total = parse_assignments(args.input)

    if total == 0:
        print(f"WARNING: No reads found in {args.input} — output will be all zeros",
              file=sys.stderr)

    # Print summary to stdout
    print(f"  Total reads: {total:,}")
    for atype in ASSIGNMENT_ORDER:
        n = counts.get(atype, 0)
        pct = n / total * 100 if total > 0 else 0.0
        print(f"    {atype:<32}: {n:>8,}  ({pct:5.1f}%)")

    write_qc_tsv(sample, counts, total, args.outdir)
    write_multiqc_tsv(sample, counts, args.outdir)

    print(f"Done: {sample}")


if __name__ == "__main__":
    main()
