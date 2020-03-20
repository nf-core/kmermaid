#!/usr/bin/env python
from __future__ import print_function
from collections import defaultdict
import logging
from itertools import groupby
import argparse

import pandas as pd
from tqdm import tqdm
import screed

# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)

CELL_BARCODE_PATTERN = '(CB|XC):Z:([ACGT]+)\-1'
MOLECULAR_BARCODE_PATTERN = '(UB|XB):Z:([ACGT]+)'
MIN_UMI_PER_CELL = 100


def get_cell_barcode(record, cell_barcode_pattern=CELL_BARCODE_PATTERN):
    """Extract cell barcode sequence"""
    found_cell_barcode = re.findall(cell_barcode_pattern, record['name'])
    if found_cell_barcode:
        return found_cell_barcode[0][1]


def get_molecular_barcode(record,
                          molecular_barcode_pattern=MOLECULAR_BARCODE_PATTERN):
    """Extract molecular barcode (UMI) sequence"""
    found_molecular_barcode = re.findall(molecular_barcode_pattern,
                                         record['name'])
    if found_molecular_barcode:
        return found_molecular_barcode[0][1]


def get_cell_barcode_umi_counts(reads,
                                cell_barcode_pattern=CELL_BARCODE_PATTERN,
                                molecular_barcode_pattern=MOLECULAR_BARCODE_PATTERN):
    """Count number of unique molecules per cell"""
    barcode_counter = defaultdict(set)

    with screed.open(reads) as f:
        for record in tqdm(f):
            cell_barcode = get_cell_barcode(record, cell_barcode_pattern)
            if cell_barcode is not None:
                molecular_barcode = get_molecular_barcode(record,
                                                          molecular_barcode_pattern)
                barcode_counter[cell_barcode].add(molecular_barcode)
    return barcode_counter


def main(reads, csv, cell_barcode_pattern=CELL_BARCODE_PATTERN,
         molecular_barcode_pattern=MOLECULAR_BARCODE_PATTERN,
         min_umi_per_cell=MIN_UMI_PER_CELL):
    barcode_counter = get_cell_barcode_umi_counts(reads, cell_barcode_pattern,
                                                  molecular_barcode_pattern)
    umi_per_barcode = {k: len(v) for k, v in barcode_counter.items()}
    series = pd.Series(umi_per_barcode)
    cells_with_minimum_n_umi = series[series >= min_umi_per_cell]
    cells_with_minimum_n_umi.to_csv(csv, header=False, index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Count number of unique molecular identifiers per cell""")
    parser.add_argument("--reads", type=str,
                        help="Reads file (fastq/fasta, can be gzipped or not)")
    parser.add_argument("--cell-barcode-pattern", type=str,
                        help="Regular expressions for cell barcodes. Default is" \
                             " 10x Genomics 'CB:Z' tag",
                        default=CELL_BARCODE_PATTERN)
    parser.add_argument("--molecular-barcode-pattern", type=str,
                        help="Regular expressions for molecular barcodes. " \
                             "Default is 10x Genomics 'CB:Z' tag",
                        default=MOLECULAR_BARCODE_PATTERN)
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
