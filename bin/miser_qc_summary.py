#!/usr/bin/env python3
"""
miser_qc_summary.py — MisER micro-exon rescue QC summarisation (V1)

Parses the MisER missed_small.bed event table, collapses read-level rows into
unique rescued exon events, and produces three outputs:

    <sample>.rescued_microexons.tsv  — unique exon catalog with support counts
    <sample>.rescued_microexons.bed  — IGV/UCSC/JBrowse-compatible BED track
    <sample>.rescue_metrics.tsv      — per-sample QC metrics

Usage:
    python miser_qc_summary.py \\
        --input  Control1.missed_small.bed \\
        --outdir results/miser_qc

Options:
    --input         Path to MisER missed_small.bed file (required)
    --outdir        Output directory (default: current directory)
    --min-support   Minimum supporting reads to include an exon (default: 1)
    --all-events    Include events that did not pass realn_flag filter.
                    Default: only passed events are used for the exon catalog
                    and BED track. Metrics always report both.
    --sep           Column separator in input file (default: auto-detect)
    --inspect       Print column names and first 3 rows then exit.
                    Use this first if you get a column not found error.

Notes on output naming:
    Output filenames are derived automatically from the input filename.
    Control1.missed_small.bed → Control1.rescued_microexons.tsv etc.
    Use --outdir to organise outputs per sample:
        --outdir results/miser_qc/Control1

Notes:
    - MisER stores coordinates and sizes as stringified Python lists e.g.
      [19234179]. These are parsed safely with ast.literal_eval().
    - Multi-exon rescue events (list length > 1) are exploded to one row
      per rescued exon before collapsing.
    - BED score = supporting read count, capped at 1000 (BED spec limit).
    - realn_flag == 1 means the rescue passed MisER's internal filter.

Author: LoReRNA pipeline — VYP Lab, UCL
Version: 1.0.0
"""

import argparse
import ast
import sys
from pathlib import Path

import pandas as pd


# ── Column name constants ──────────────────────────────────────────────────────
# These are the expected MisER output column names.
# If your MisER version uses different names, pass --inspect to see what
# columns are actually present, then update these constants accordingly.

COL_CHROM   = 'chrom'
COL_STARTS  = 'small_exon_starts'
COL_ENDS    = 'small_exon_ends'
COL_SIZES   = 'exon_sizes'
COL_STRAND  = 'read_strand'
COL_SCORE   = 'realn_increase_score'
COL_FLAG    = 'realn_flag'      # 'True'/'False' string in MisER output

# Additional columns present in MisER output (not used in collapsing
# but available for future extension):
#   read_i, tx_i, sum_exon_size, realign_start, realign_end,
#   delta_ratio, delta_ratio_flag, margin_len

REQUIRED_COLS = [COL_CHROM, COL_STARTS, COL_ENDS, COL_SIZES,
                 COL_STRAND, COL_SCORE, COL_FLAG]

BED_SCORE_MAX = 1000


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description='MisER micro-exon rescue QC summarisation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage:')[0].strip()
    )
    parser.add_argument('--input',       required=True,  help='MisER missed_small.bed file')
    parser.add_argument('--outdir',      default='.',    help='Output directory (default: current dir)')
    parser.add_argument('--min-support', type=int, default=1,
                        help='Minimum supporting reads to report an exon (default: 1)')
    parser.add_argument('--all-events',  action='store_true',
                        help='Include events that did not pass realn_flag. '
                             'Default: only passed events used for catalog + BED')
    parser.add_argument('--sep',         default=None,
                        help='Column separator (default: auto-detect tab or whitespace)')
    parser.add_argument('--inspect',     action='store_true',
                        help='Print columns and first 3 rows then exit')
    return parser.parse_args()


# ── File loading ───────────────────────────────────────────────────────────────

def load_miser_file(path: Path, sep=None) -> pd.DataFrame:
    """Load MisER output file, auto-detecting separator if not specified."""
    if not path.exists():
        sys.exit(f"ERROR: Input file not found: {path}")
    if path.stat().st_size == 0:
        sys.exit(f"ERROR: Input file is empty: {path}")

    # Try tab first, then whitespace
    separators = [sep] if sep else ['\t', r'\s+']
    df = None
    for s in separators:
        try:
            df = pd.read_csv(path, sep=s, engine='python', comment='#')
            if len(df.columns) > 3:
                break
        except Exception:
            continue

    if df is None or df.empty:
        sys.exit(f"ERROR: Could not parse input file: {path}")

    print(f"  Loaded {len(df):,} rows, {len(df.columns)} columns")
    return df


def validate_columns(df: pd.DataFrame, path: Path):
    """Check required columns are present. Print helpful error if not."""
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        print(f"\nERROR: Missing required columns: {missing}")
        print(f"\nColumns found in {path.name}:")
        for i, col in enumerate(df.columns, 1):
            print(f"  {i:2d}. {col}")
        print("\nTip: Run with --inspect to see the first few rows.")
        print("     If column names differ, update the COL_* constants at the top of this script.")
        sys.exit(1)


# ── Parsing stringified lists ──────────────────────────────────────────────────

def safe_parse_list(value):
    """
    Parse a stringified Python list from MisER output.
    e.g. '[19234179]' → [19234179]
    Uses ast.literal_eval — safe, no code execution.
    Returns None if parsing fails.
    """
    if pd.isna(value):
        return None
    try:
        result = ast.literal_eval(str(value).strip())
        return result if isinstance(result, list) else [result]
    except (ValueError, SyntaxError):
        return None


def parse_list_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Parse all stringified list columns and drop unparseable rows."""
    list_cols = [COL_STARTS, COL_ENDS, COL_SIZES]

    n_before = len(df)
    for col in list_cols:
        df[col] = df[col].apply(safe_parse_list)

    # Drop rows where any list column failed to parse
    mask_valid = df[list_cols].notna().all(axis=1)
    n_dropped = (~mask_valid).sum()
    if n_dropped > 0:
        print(f"  WARNING: Dropped {n_dropped:,} rows with unparseable list columns")
    df = df[mask_valid].copy()

    if df.empty:
        sys.exit("ERROR: No parseable rows remaining after list column parsing. "
                 "Check input file format with --inspect.")

    print(f"  Parsed list columns: {n_before - n_dropped:,} of {n_before:,} rows valid")
    return df


# ── Explode multi-exon events ──────────────────────────────────────────────────

def explode_events(df: pd.DataFrame) -> pd.DataFrame:
    """
    Explode list columns so each row represents a single rescued exon.

    MisER can rescue multiple micro-exons in a single read — in these cases
    small_exon_starts, small_exon_ends and exon_sizes contain lists with
    more than one element. Exploding gives one row per exon per read.

    All three list columns must have the same length per row — rows where
    they differ are dropped with a warning (indicates a malformed event).
    """
    # Check list lengths are consistent within each row
    list_cols = [COL_STARTS, COL_ENDS, COL_SIZES]
    lengths = df[list_cols].apply(lambda col: col.apply(len))
    consistent = (lengths[COL_STARTS] == lengths[COL_ENDS]) & \
                 (lengths[COL_STARTS] == lengths[COL_SIZES])

    n_inconsistent = (~consistent).sum()
    if n_inconsistent > 0:
        print(f"  WARNING: Dropped {n_inconsistent:,} rows with mismatched list lengths")
        df = df[consistent].copy()

    n_before = len(df)
    df = df.explode(list_cols)
    n_after = len(df)

    multi_exon = n_after - n_before
    if multi_exon > 0:
        print(f"  Exploded {multi_exon:,} additional exon entries from multi-exon read events")

    # Convert exploded coordinates to integers
    for col in list_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Drop rows where coordinates became NaN after explosion
    df = df.dropna(subset=list_cols)
    df[list_cols] = df[list_cols].astype(int)

    return df.reset_index(drop=True)


# ── Filter on realn_flag ───────────────────────────────────────────────────────

def apply_flag_filter(df: pd.DataFrame, use_all: bool):
    """
    Filter events by realn_flag.

    realn_flag == 1: rescue passed MisER's internal confidence filter
    realn_flag == 0: rescue did not pass

    Returns:
        df_passed  — events that passed (always computed for metrics)
        df_catalog — events used for exon catalog and BED (passed only,
                     or all events if --all-events is set)
    """
    # MisER writes realn_flag as Python string 'True'/'False', not integer.
    # Normalise to boolean: accept 'True', 'true', '1', True
    def parse_flag(val):
        if isinstance(val, bool):
            return val
        return str(val).strip().lower() in ('true', '1')

    df['_flag_bool'] = df[COL_FLAG].apply(parse_flag)

    df_passed = df[df['_flag_bool']].copy()
    df_failed = df[~df['_flag_bool']].copy()

    n_total  = len(df)
    n_passed = len(df_passed)
    n_failed = len(df_failed)
    pass_rate = (n_passed / n_total * 100) if n_total > 0 else 0.0

    print(f"  realn_flag: {n_passed:,} passed / {n_total:,} total ({pass_rate:.1f}% pass rate)")
    if n_failed > 0 and not use_all:
        print(f"  {n_failed:,} events excluded (realn_flag == 0). Use --all-events to include them.")

    df_catalog = df.drop(columns=['_flag_bool']) if use_all else df_passed.drop(columns=['_flag_bool'])
    return df_passed, df_catalog


# ── Collapse to unique rescued exons ──────────────────────────────────────────

def collapse_exons(df: pd.DataFrame, min_support: int) -> pd.DataFrame:
    """
    Collapse read-level events to unique rescued exon events.

    Grouping key: (chrom, start, end, strand, exon_size)
    This is the biological unit of interest — one row per unique rescued exon.

    Aggregations per unique exon:
        support_reads  — number of reads supporting the rescue
        mean_score     — mean realignment improvement score
        max_score      — maximum realignment improvement score
    """
    group_key = [COL_CHROM, COL_STARTS, COL_ENDS, COL_STRAND, COL_SIZES]

    df[COL_SCORE] = pd.to_numeric(df[COL_SCORE], errors='coerce')

    collapsed = df.groupby(group_key, as_index=False).agg(
        support_reads  = (COL_SCORE, 'count'),
        mean_score     = (COL_SCORE, 'mean'),
        max_score      = (COL_SCORE, 'max')
    ).rename(columns={
        COL_CHROM:  'chrom',
        COL_STARTS: 'start',
        COL_ENDS:   'end',
        COL_STRAND: 'strand',
        COL_SIZES:  'exon_size'
    })

    collapsed['mean_score'] = collapsed['mean_score'].round(4)
    collapsed['max_score']  = collapsed['max_score'].round(4)

    # Apply minimum support filter
    n_before = len(collapsed)
    collapsed = collapsed[collapsed['support_reads'] >= min_support].copy()
    n_filtered = n_before - len(collapsed)
    if n_filtered > 0:
        print(f"  Removed {n_filtered:,} exons with support_reads < {min_support}")

    # Sort by chromosome then start position
    collapsed = collapsed.sort_values(['chrom', 'start']).reset_index(drop=True)

    return collapsed[['chrom', 'start', 'end', 'strand', 'exon_size',
                       'support_reads', 'mean_score', 'max_score']]


# ── Build BED track ────────────────────────────────────────────────────────────

def build_bed(exons: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """
    Build a BED6 track from the unique exon catalog.

    Columns: chrom, start, end, name, score, strand
    - name  = zero-padded rescued exon ID e.g. rescued_exon_00001
    - score = support_reads capped at 1000 (BED spec: 0–1000)
    """
    n = len(exons)
    bed = pd.DataFrame({
        'chrom':  exons['chrom'],
        'start':  exons['start'],
        'end':    exons['end'],
        'name':   [f"rescued_exon_{i+1:05d}" for i in range(n)],
        'score':  exons['support_reads'].clip(upper=BED_SCORE_MAX),
        'strand': exons['strand']
    })
    return bed


# ── Compute QC metrics ─────────────────────────────────────────────────────────

def compute_metrics(df_all: pd.DataFrame,
                    df_passed: pd.DataFrame,
                    exons: pd.DataFrame,
                    prefix: str) -> pd.DataFrame:
    """
    Compute per-sample QC metrics.

    Reports raw event counts (pre-collapse) and unique exon counts (post-collapse).
    Pass rate is computed from read-level events before collapsing — this is the
    most meaningful indicator of MisER run quality.
    """
    n_total   = len(df_all)
    n_passed  = len(df_passed)
    pass_rate = round(n_passed / n_total * 100, 2) if n_total > 0 else 0.0

    metrics = {
        'sample':                  prefix,
        'total_read_events':       n_total,
        'passed_read_events':      n_passed,
        'failed_read_events':      n_total - n_passed,
        'pass_rate_pct':           pass_rate,
        'unique_microexons':       len(exons),
        'median_exon_size':        int(exons['exon_size'].median()) if len(exons) > 0 else 0,
        'mean_exon_size':          round(exons['exon_size'].mean(), 2) if len(exons) > 0 else 0.0,
        'min_exon_size':           int(exons['exon_size'].min()) if len(exons) > 0 else 0,
        'max_exon_size':           int(exons['exon_size'].max()) if len(exons) > 0 else 0,
        'mean_support_reads':      round(exons['support_reads'].mean(), 2) if len(exons) > 0 else 0.0,
        'max_support_reads':       int(exons['support_reads'].max()) if len(exons) > 0 else 0,
        'exons_single_read':       int((exons['support_reads'] == 1).sum()),
        'exons_multi_read':        int((exons['support_reads'] > 1).sum()),
    }

    return pd.DataFrame([metrics])


# ── Write outputs ──────────────────────────────────────────────────────────────

def write_outputs(outdir: Path, prefix: str,
                  exons: pd.DataFrame,
                  bed: pd.DataFrame,
                  metrics: pd.DataFrame):
    """Write all three output files."""
    outdir.mkdir(parents=True, exist_ok=True)

    path_tsv     = outdir / f"{prefix}.rescued_microexons.tsv"
    path_bed     = outdir / f"{prefix}.rescued_microexons.bed"
    path_metrics = outdir / f"{prefix}.rescue_metrics.tsv"

    exons.to_csv(path_tsv, sep='\t', index=False)
    print(f"  Written: {path_tsv}")

    # BED track — no header, tab-separated (IGV/UCSC requirement)
    bed.to_csv(path_bed, sep='\t', index=False, header=False)
    print(f"  Written: {path_bed}")

    metrics.to_csv(path_metrics, sep='\t', index=False)
    print(f"  Written: {path_metrics}")

    return path_tsv, path_bed, path_metrics


# ── Inspect mode ───────────────────────────────────────────────────────────────

def inspect_file(path: Path, sep=None):
    """Print column names and first 3 rows then exit. Use when debugging."""
    df = load_miser_file(path, sep)
    print(f"\nFile: {path}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"\nColumns ({len(df.columns)}):")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i:2d}. {col!r}")
    print(f"\nFirst 3 rows:")
    print(df.head(3).to_string())
    sys.exit(0)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    input_path = Path(args.input)
    outdir     = Path(args.outdir)

    # Inspect mode — print columns and exit
    if args.inspect:
        inspect_file(input_path, args.sep)

    # Derive prefix from input filename
    # Control1.missed_small.bed → Control1
    prefix = input_path.name.replace('.missed_small.bed', '')
    if prefix == input_path.name:
        # Fallback: strip everything from first dot
        prefix = input_path.stem.split('.')[0]

    print(f"\n{'='*60}")
    print(f" MisER QC Summary — {prefix}")
    print(f"{'='*60}")
    print(f"\n[1/6] Loading input file: {input_path}")

    df = load_miser_file(input_path, args.sep)
    validate_columns(df, input_path)

    print(f"\n[2/6] Parsing stringified list columns")
    df = parse_list_columns(df)

    print(f"\n[3/6] Exploding multi-exon events")
    df = explode_events(df)
    print(f"  Total exon-level events: {len(df):,}")

    print(f"\n[4/6] Filtering on realn_flag")
    df_passed, df_catalog = apply_flag_filter(df, args.all_events)

    print(f"\n[5/6] Collapsing to unique rescued exons")
    exons = collapse_exons(df_catalog, args.min_support)
    print(f"  Unique rescued exons: {len(exons):,}")

    # BED track
    bed = build_bed(exons, prefix)

    # Metrics — always uses df (all events) and df_passed for pass rate
    metrics = compute_metrics(df.drop(columns=['_flag_bool'], errors='ignore'), df_passed, exons, prefix)

    print(f"\n[6/6] Writing outputs to: {outdir}")
    write_outputs(outdir, prefix, exons, bed, metrics)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f" Summary: {prefix}")
    print(f"{'='*60}")
    print(f"  Total read events    : {metrics['total_read_events'].iloc[0]:,}")
    print(f"  Passed events        : {metrics['passed_read_events'].iloc[0]:,} "
          f"({metrics['pass_rate_pct'].iloc[0]}%)")
    print(f"  Unique micro-exons   : {metrics['unique_microexons'].iloc[0]:,}")
    print(f"  Median exon size     : {metrics['median_exon_size'].iloc[0]} bp")
    print(f"  Mean support reads   : {metrics['mean_support_reads'].iloc[0]}")
    print(f"  Multi-read supported : {metrics['exons_multi_read'].iloc[0]:,} exons")
    print(f"  Single-read only     : {metrics['exons_single_read'].iloc[0]:,} exons")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
