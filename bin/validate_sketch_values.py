#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict, Counter
import logging
import argparse
import glob
import os
import sys


# Create a logger
logging.basicConfig(format='%(name)s - %(asctime)s %(levelname)s: %(message)s')
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def get_sketch_values(values, values_log2):
    if values:
        sketch_values = [int(x) for x in values.split(',')]
    else:
        sketch_values = [2 ** int(x) for x in values_log2.split(',')]
    return sketch_values


def value_or_bool(value):
    if value == 'false':
        return False
    if value == 'true':
        logger.exception("Must set a value for the --sketch_num_hashes, "
                        "--sketch_num_hashes_log2, --sketch_scaled, "
                        "--sketch_scaled_log2 options! Cannot simply set the "
                        "flag. E.g. '--sketch_num_hashes 5' is valid but "
                        "'--sketch_num_hashes' on its own is not")
        sys.exit(1)
    else:
        return value


def main(sketch_num_hashes, sketch_num_hashes_log2, sketch_scaled, sketch_scaled_log2, out,
         sketch_style):
    sketch_num_hashes = value_or_bool(sketch_num_hashes)
    sketch_num_hashes_log2 = value_or_bool(sketch_num_hashes_log2)
    sketch_scaled = value_or_bool(sketch_scaled)
    sketch_scaled_log2 = value_or_bool(sketch_scaled_log2)

    using_size = sketch_num_hashes or sketch_num_hashes_log2
    using_scaled = sketch_scaled or sketch_scaled_log2

    if using_size and using_scaled:
        logger.exception("Cannot specify both sketch scales and sizes! Can only"
                         " use one of --sketch_num_hashes, --sketch_num_hashes_log2, --sketch_scaled, "
                         "--sketch_scaled_log2. Exiting.")
        sys.exit(1)


    if using_size:
        if sketch_num_hashes and sketch_num_hashes_log2:
            logger.exception("Cannot specify both --sketch_num_hashes and --sketch_num_hashes_log2! Exiting.")
            sys.exit(1)
        sketch_values = get_sketch_values(sketch_num_hashes, sketch_num_hashes_log2)
        with open(sketch_style, 'w') as f:
            f.write('size')
    elif using_scaled:
        if sketch_scaled and sketch_scaled_log2:
            logger.exception("Cannot specify both --sketch_scaled and --sketch_scaled_log2! Exiting.")
            sys.exit(1)
        sketch_values = get_sketch_values(sketch_scaled, sketch_scaled_log2)
        with open(sketch_style, 'w') as f:
            f.write('scaled')

    else:
        logger.info("Did not specify a sketch size or scale with any of "
                    "--sketch_num_hashes, --sketch_num_hashes_log2, --sketch_scaled, --sketch_scaled_log2! "
                    "Falling back on sourmash's default of --sketch_scaled 500")
        sketch_values = [500]

    with open(out, 'w') as f:
        for number in sketch_values:
            f.write(f"{number}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Ensure that sketch sizes/scaleds provided are valid""")
    parser.add_argument("--sketch_num_hashes", type=str, help="Flat size of the sketches")
    parser.add_argument("--sketch_num_hashes_log2", type=str, help="Flat size of the sketches, log2")
    parser.add_argument("--sketch_scaled", type=str, help="Fraction of total observed hashes to observe")
    parser.add_argument("--sketch_scaled_log2", type=str, help="Fraction of total observed hashes to observe, log2")
    parser.add_argument("-o", "--output", dest='output', default='sketch_values.txt',
        type=str, help="file with output")
    parser.add_argument("--sketch_style", default='sketch_style.txt', type=str, help="file indicating 'size' or 'scaled'")

    args = parser.parse_args()
    main(args.sketch_num_hashes, args.sketch_num_hashes_log2, args.sketch_scaled,
         args.sketch_scaled_log2, args.output, args.sketch_style)
