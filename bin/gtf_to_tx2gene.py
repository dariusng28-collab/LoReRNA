#!/usr/bin/env python

import argparse


def parse_attributes(field):
    attrs = {}
    for raw in field.strip().split(";"):
        raw = raw.strip()
        if not raw:
            continue
        if " " not in raw:
            continue
        key, value = raw.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def main():
    parser = argparse.ArgumentParser(description="Extract transcript-to-gene mapping from a GTF file.")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    seen = set()
    with open(args.output, "w", encoding="utf-8") as out_handle:
        out_handle.write("transcript_id\tgene_id\tgene_name\n")
        with open(args.gtf, "r", encoding="utf-8") as in_handle:
            for line in in_handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue
                attrs = parse_attributes(fields[8])
                transcript_id = attrs.get("transcript_id")
                gene_id = attrs.get("gene_id")
                gene_name = attrs.get("gene_name", gene_id)
                if not transcript_id or not gene_id:
                    continue
                key = (transcript_id, gene_id, gene_name)
                if key in seen:
                    continue
                seen.add(key)
                out_handle.write(f"{transcript_id}\t{gene_id}\t{gene_name}\n")


if __name__ == "__main__":
    main()
