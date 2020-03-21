#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
import logging
from itertools import groupby
import argparse

import pandas as pd
from tqdm import tqdm
import screed

from count_umis_per_cell import get_cell_barcode, get_molecular_barcode, \
    CELL_BARCODE_PATTERN

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def get_good_cell_barcode_records(reads, good_barcodes, cell_barcode_pattern):
    good_cell_barcode_records = defaultdict(list)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record, cell_barcode_pattern)
            if cell_barcode in good_barcodes:
                good_cell_barcode_records[cell_barcode].append(record)
    return good_cell_barcode_records


def record_to_fastq_string(record):
    """Convert screed record into valid fastq output"""
    return f"@{record['name']}\n{record['sequence']}\n+\n{record['quality']}\n"


def write_records(records, filename):
    """Write screed records to a file (compression auto-detected)"""
    if filename.endswith('gz'):
        import gzip
        opener = gzip.open
        mode = 'wb'
    else:
        opener = open
        mode = 'w'

    with opener(filename, mode) as f:
        f.writelines([record_to_fastq_string(r) for r in records])


def read_barcodes(filename):
    """Read a barcodes.tsv filename, output set of unique barcodes

    They should already be unique.. the "frozenset" datastructure is just for
    quick checking of set membership
    """
    with open(filename) as f:
        barcodes = frozenset(x.strip() for x in f.readlines())
    return barcodes


def main(reads, good_barcodes_filename, outdir, channel_id=None,
         cell_barcode_pattern=CELL_BARCODE_PATTERN):
    if channel_id is not None:
        prefix = f"{channel_id}_"
    else:
        prefix = ''

    good_barcodes = read_barcodes(good_barcodes_filename)

    good_cell_barcode_records = get_good_cell_barcode_records(
        reads, good_barcodes, cell_barcode_pattern)

    for cell_barcode, records in good_cell_barcode_records.items():
        filename = f"{outdir}/{prefix}_{cell_barcode}.fastq.gz"
        write_records(records, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Extract reads for each cell barcode with sufficient 
        UMIs""")
    parser.add_argument("--reads", type=str,
                        help="Reads file (fastq/fasta, can be gzipped or not)")
    parser.add_argument('-b', "--good-barcodes", type=str,
                        default='barcodes.csv',
                        help="QC-passing barcodes from 10x")
    parser.add_argument("--cell-barcode-pattern", type=str,
                        help="Regular expressions for cell barcodes. Default "
                             "is 10x Genomics 'CB:Z' tag",
                        default=CELL_BARCODE_PATTERN)
    parser.add_argument("--outdir", type=str,
                        default='.',
                        help="Output directory for fastqs. Default is current"
                             " directory")
    parser.add_argument("--channel-id", type=str,
                        help="Output prefix for fastqs")

    args = parser.parse_args()
    main(args.reads, args.good_barcodes, args.outdir, args.channel_id,
         args.cell_barcode_pattern)
