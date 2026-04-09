#!/usr/bin/env python

import argparse
from collections import defaultdict


def parse_attributes(field):
    attrs = {}
    for raw in field.strip().split(";"):
        raw = raw.strip()
        if not raw or " " not in raw:
            continue
        key, value = raw.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def block_starts_and_sizes(exons, chrom_start):
    starts = []
    sizes = []
    for start, end in exons:
        starts.append(start - chrom_start)
        sizes.append(end - start)
    return starts, sizes


def main():
    parser = argparse.ArgumentParser(description="Convert transcript/exon features in a GTF file to BED12.")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    transcripts = {}
    exons = defaultdict(list)

    with open(args.gtf, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            start = int(fields[3]) - 1
            end = int(fields[4])
            chrom = fields[0]
            strand = fields[6]

            if feature == "transcript":
                transcripts[transcript_id] = {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
            elif feature == "exon":
                exons[transcript_id].append((start, end))

    with open(args.output, "w", encoding="utf-8") as out_handle:
        for transcript_id, tx in transcripts.items():
            transcript_exons = sorted(exons.get(transcript_id, []))
            if not transcript_exons:
                transcript_exons = [(tx["start"], tx["end"])]

            chrom_start = transcript_exons[0][0]
            chrom_end = transcript_exons[-1][1]
            block_starts, block_sizes = block_starts_and_sizes(transcript_exons, chrom_start)
            block_count = len(transcript_exons)

            out_handle.write(
                "\t".join(
                    [
                        tx["chrom"],
                        str(chrom_start),
                        str(chrom_end),
                        transcript_id,
                        "0",
                        tx["strand"],
                        str(chrom_start),
                        str(chrom_end),
                        "0",
                        str(block_count),
                        ",".join(str(size) for size in block_sizes),
                        ",".join(str(start) for start in block_starts),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
