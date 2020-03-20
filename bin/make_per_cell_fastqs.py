#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
import logging
from itertools import groupby
import argparse

import pandas as pd
from tqdm import tqdm
import screed

from count_umis_per_cell import get_cell_barcode, get_molecular_barcode

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def get_good_cell_barcode_records(reads, cells_with_minimum_n_umi):
    good_cell_barcode_records = defaultdict(list)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record)
            if cell_barcode in cells_with_minimum_n_umi:
                good_cell_barcode_records[cell_barcode].append(record)
    return good_cell_barcode_records


def record_to_fastq_string(record):
    """Convert screed record into valid fastq output"""
    return f"@{record['name']}\n{record['sequence']}\n+\n{record['quality']}\n"


def write_records(records, filename):
    with open(filename, 'w') as f:
        f.writelines([record_to_fastq_string(r) for r in records])


def main(reads, csv, outdir, prefix):
    cells_with_minimum_n_umi = pd.read_csv(csv, header=False, index_col=0,
                                           squeeze=True)
    good_cell_barcode_records = get_good_cell_barcode_records(
        reads, cells_with_minimum_n_umi)

    for cell_barcode, records in good_cell_barcode_records.items():
        filename = f"{outdir}/{prefix}_{cell_barcode}.fastq"
        write_records(records, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Count number of unique molecular identifiers per cell""")
    parser.add_argument("--reads", type=str,
                        help="Reads file (fastq/fasta, can be gzipped or not)")
    parser.add_argument("--min-umi-per-cell", type=int,
                        default=MIN_UMI_PER_CELL,
                        help="Minimum number of unique molecular identifiers "
                             "(barcodes) per cell")
    parser.add_argument("-o", "--csv", type=str,
                        default='n_umis_per_cell_barcode.csv',
                        help="Number of UMIs counted per cell")

    args = parser.parse_args()
    main(args.reads, args.csv, args.cell_barcode_pattern,
         args.molecular_barcode_pattern,
         args.min_umi_per_cell)
