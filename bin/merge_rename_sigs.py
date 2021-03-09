#!/usr/bin/env python
import argparse
import itertools
import shutil

import sourmash


def merge(filenames, ksize, moltype, name=None, outsig=None):
    """
    merge one or more signatures.

    Adapted from 'sourmash sig merge' command
    """

    first_sig = None
    mh = None
    total_loaded = 0

    for sigfile in filenames:
        this_n = 0

        for sigobj in sourmash.load_file_as_signatures(
            sigfile, ksize=ksize, select_moltype=moltype
        ):

            if sigobj is None:
                error(
                    "No signature in file {}",
                    sigfile,
                )

            # first signature? initialize a bunch of stuff
            if first_sig is None:
                first_sig = sigobj
                mh = first_sig.minhash.copy_and_clear()

            try:
                sigobj_mh = sigobj.minhash

                mh.merge(sigobj_mh)
            except:
                error(
                    "ERROR when merging signature '{}' ({}) from file {}",
                    sigobj.name(),
                    sigobj.md5sum()[:8],
                    sigfile,
                )
                raise

            this_n += 1
            total_loaded += 1

    merged_sigobj = sourmash.SourmashSignature(mh)
    if name is not None:
        merged_sigobj._name = name

    if outsig is not None:
        with open(outsig, "wt") as f:
            sourmash.save_signatures([merged_sigobj], fp=f)

    return merged_sigobj


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Merge signatures with same ksize and moltype"""
    )
    # Signature files
    parser.add_argument("sigfiles", nargs="+")
    parser.add_argument(
        "--moltypes",
        type=str,
        help="Molecule types, comma-separated. e,g. 'protein,dayhoff'",
        required=True,
    )
    parser.add_argument(
        "--ksizes", type=str, help="K-mer sizes to combine", required=True
    )

    parser.add_argument(
        "-o", "--outsig", type=str, help="Signature file to output to", required=True
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="Name of the signature to use",
    )

    args = parser.parse_args()

    # Only iterate and read over the sigfiles if there is really something to merge
    # "something to merge" = there is more than one sigfile. otherwise there's no point
    # in reading in the files only to make the same file again
    if len(args.sigfiles) > 1:
        ksizes = map(int, args.ksizes.split(","))
        moltypes = args.moltypes.split(",")

        merged_sigobjs = []
        for moltype, ksize in itertools.product(moltypes, ksizes):
            merged_sigobj = merge(
                args.sigfiles, moltype=moltype, ksize=ksize, name=args.name
            )
            merged_sigobjs.append(merged_sigobj)

        with open(args.outsig, "wt") as f:
            sourmash.save_signatures(merged_sigobjs, fp=f)
    else:
        # Otherwise, nothing to merge. Simply copy the file to the
        # output signature location
        shutil.copyfile(args.sigfiles[0], args.outsig)
