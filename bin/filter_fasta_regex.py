#!/usr/bin/env python

import argparse
import re


import screed


def write_records_to_fasta(records, fasta):
    with open(fasta, "w") as f:
        for record in records:
            f.write(f'>{record["name"]}\n{record["sequence"]}\n')


def filter_records(fasta, pattern):
    filtered_records = []
    with screed.open(fasta) as records:
        for record in records:
            name = record["name"]
            if re.findall(pattern, name, flags=re.I):
                filtered_records.append(record)
    return filtered_records


def filter_fasta_with_regex(fasta_to_filter, out_fasta, regex):
    record_subset = filter_records(fasta_to_filter, regex)
    write_records_to_fasta(record_subset, out_fasta)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Extract sequences whose names match a pattern"""
    )
    parser.add_argument("--input-fasta", type=str, help="Sequence file to filter")
    parser.add_argument("--output-fasta", type=str, help="File to write")
    parser.add_argument(
        "--regex-pattern",
        type=str,
        help="Regular expression pattern to match for the names of seuqences in the file",
    )
    args = parser.parse_args()

    filter_fasta_with_regex(args.input_fasta, args.output_fasta, args.regex_pattern)
