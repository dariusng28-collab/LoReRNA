#!/usr/bin/env python

import argparse
import csv
import os
from collections import defaultdict


ID_COLUMNS = ("transcript_id", "target_id", "target_name", "name", "Name")
COUNT_COLUMNS = ("num_reads", "est_counts", "count", "counts")


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames or []


def find_column(candidates, fieldnames):
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    return fieldnames[0] if fieldnames else None


def parse_samplesheet(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    required = {"sample", "condition", "bam"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise ValueError(f"Samplesheet is missing required columns: {', '.join(sorted(missing))}")
    return rows


def parse_tx2gene(path):
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        mapping = {}
        for row in reader:
            mapping[row["transcript_id"]] = {
                "gene_id": row.get("gene_id", ""),
                "gene_name": row.get("gene_name", row.get("gene_id", "")),
            }
    return mapping


def collect_oarfish_outputs(paths):
    grouped = defaultdict(dict)
    for path in paths:
        basename = os.path.basename(path)
        if basename.endswith(".quant"):
            sample = basename[:-6]
            grouped[sample]["quant"] = path
        elif basename.endswith(".infreps.pq"):
            sample = basename[:-11]
            grouped[sample]["infreps"] = path
        elif basename.endswith(".meta_info.json"):
            sample = basename[:-15]
            grouped[sample]["meta_info"] = path
    return grouped


def main():
    parser = argparse.ArgumentParser(description="Merge per-sample oarfish outputs into count and usage matrices.")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("paths", nargs="+")
    args = parser.parse_args()

    samples = parse_samplesheet(args.samplesheet)
    tx2gene = parse_tx2gene(args.tx2gene)
    grouped_outputs = collect_oarfish_outputs(args.paths)

    sample_names = [row["sample"] for row in samples]
    counts_by_transcript = defaultdict(dict)

    for sample in sample_names:
        if sample not in grouped_outputs or "quant" not in grouped_outputs[sample]:
            raise FileNotFoundError(f"Missing oarfish quant output for sample '{sample}'")
        quant_rows, fieldnames = read_tsv(grouped_outputs[sample]["quant"])
        id_col = find_column(ID_COLUMNS, fieldnames)
        count_col = find_column(COUNT_COLUMNS, fieldnames)
        for row in quant_rows:
            transcript_id = row[id_col]
            raw_count = row.get(count_col, "0") if count_col else "0"
            try:
                value = float(raw_count)
            except ValueError:
                value = 0.0
            counts_by_transcript[transcript_id][sample] = value

    ordered_transcripts = list(tx2gene.keys())
    ordered_transcripts.extend(sorted(set(counts_by_transcript.keys()) - set(ordered_transcripts)))

    counts_path = os.path.join(args.outdir, "oarfish.transcript_counts.tsv")
    usage_path = os.path.join(args.outdir, "oarfish.isoform_usage.tsv")
    manifest_path = os.path.join(args.outdir, "oarfish.infreps_manifest.tsv")

    gene_totals = {sample: defaultdict(float) for sample in sample_names}
    for transcript_id in ordered_transcripts:
        gene_id = tx2gene.get(transcript_id, {}).get("gene_id", "")
        for sample in sample_names:
            gene_totals[sample][gene_id] += counts_by_transcript[transcript_id].get(sample, 0.0)

    with open(counts_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_name", *sample_names])
        for transcript_id in ordered_transcripts:
            gene_id = tx2gene.get(transcript_id, {}).get("gene_id", "")
            gene_name = tx2gene.get(transcript_id, {}).get("gene_name", gene_id)
            row = [transcript_id, gene_id, gene_name]
            row.extend(counts_by_transcript[transcript_id].get(sample, 0.0) for sample in sample_names)
            writer.writerow(row)

    with open(usage_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_name", *sample_names])
        for transcript_id in ordered_transcripts:
            gene_id = tx2gene.get(transcript_id, {}).get("gene_id", "")
            gene_name = tx2gene.get(transcript_id, {}).get("gene_name", gene_id)
            row = [transcript_id, gene_id, gene_name]
            for sample in sample_names:
                denom = gene_totals[sample].get(gene_id, 0.0)
                value = counts_by_transcript[transcript_id].get(sample, 0.0)
                row.append(0.0 if denom == 0 else value / denom)
            writer.writerow(row)

    with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "quant", "infreps", "meta_info"])
        for sample in sample_names:
            writer.writerow([
                sample,
                grouped_outputs[sample].get("quant", ""),
                grouped_outputs[sample].get("infreps", ""),
                grouped_outputs[sample].get("meta_info", ""),
            ])


if __name__ == "__main__":
    main()
