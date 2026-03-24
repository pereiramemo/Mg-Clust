#!/usr/bin/env python3
"""
mg-clust module 3: Concatenation and formatting of ORFs fasta files and coverage tables.

Assumes execution inside the conda environment "mg-clust-module-3" (or equivalent)
where dependencies are available on PATH.

- Concatenates ORF protein sequences from all samples, prefixing each header with the sample name
- Filters ORFs by minimum length using bbduk
- Creates an MMseqs2 sequence database from the filtered ORFs
- Concatenates per-sample mean coverage tables, prefixing each ORF ID with the sample name
- Concatenates per-sample reads coverage tables, prefixing each ORF ID with the sample name
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import sys, os
import subprocess
import shutil
sys.path.insert(0, os.path.dirname(__file__))
from utils import run, check_tools

bbduk = "bbduk.sh"
mmseqs = "mmseqs"

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description="mg-clust module 3", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_dir", dest="input_dir", required=True,
        help="input dir where to find the ORF fasta files and coverage tables")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--min_orf_length", dest="min_orf_length", type=int, default=60,
        help="minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default: 60)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    check_tools([bbduk, mmseqs])
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    # Check that input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Input directory {args.input_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.2. Check output directory
    ###########################################################################

    if os.path.isdir(args.output_dir):
        if not args.overwrite:
            print(f"{args.output_dir} already exists; use --overwrite to overwrite")
            sys.exit(0)
        try:
            shutil.rmtree(args.output_dir)
        except Exception:
            print(f"rm -r output directory {args.output_dir} failed", file=sys.stderr)
            sys.exit(1)
    
    ###########################################################################
    # 3.3. Create output directory
    ###########################################################################
    
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.4 Concat ORF files
    ###########################################################################

    concat_orfs = os.path.join(args.output_dir, "orfs.faa")

    # Find and concat all *_orfs.faa files in subdirectories
    orf_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs.faa"))

    if not orf_files:
        print(f"No ORF files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    open(concat_orfs, "w").close()  # create/clear the output file

    for orf_file in orf_files:
        # sample_name is the filename without the "_orfs.faa" suffix
        # set when running FragGeneScanRs in mg-clust-module-2 with the -o argument
        sample_name = os.path.basename(orf_file).replace("_orfs.faa", "")

        try:
            with open(concat_orfs, "a") as out_fh:
                # here sed is used instead of python to reduce computation time
                proc = subprocess.run(
                    ["sed", f"s/^>/>{ sample_name}|/", orf_file],
                    stdout=out_fh, check=False)
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, proc.args)
            
        except subprocess.CalledProcessError:
            print(f"prefixing headers in {orf_file} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.5 Filter ORFs by length
    ###########################################################################

    concat_orfs_filt = concat_orfs.replace(".faa", "_filt.faa")

    try:
        run(
            [
                bbduk,
                f"in={concat_orfs}",
                f"out={concat_orfs_filt}",
                "overwrite=t",
                f"threads={args.nslots}",
                f"minlength={args.min_orf_length}",
                "amino=t"
            ]
        )
    except subprocess.CalledProcessError:
        print("bbduk failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6 Create mmseqs database
    ###########################################################################

    orfs_db = concat_orfs_filt.replace(".faa", "_db")

    if not os.path.exists(orfs_db):
        try:
            run(
                [
                    mmseqs,
                    "createdb",
                    concat_orfs_filt,
                    orfs_db,
                    "--dbtype", "1"
                ]
            )
        except subprocess.CalledProcessError:
            print("mmseqs createdb failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.7 Concat and format mean coverage table
    ###########################################################################

    meancov_table = os.path.join(args.output_dir, "orfs_meancov.tsv")

    # Find and concat all *_orfs_meancov.tsv files in subdirectories
    meancov_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs_meancov.tsv"))

    if not meancov_files:
        print(f"No coverage files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(meancov_table, "w", encoding="utf-8") as out_fh:
            for meancov_file in meancov_files:

                # sample_name is the filename without the "_orfs_meancov.tsv" suffix
                # set when running bedtools in mg-clust-module-2.
                sample_name = os.path.basename(meancov_file).replace("_orfs_meancov.tsv", "")

                with open(meancov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            orf_id = parts[4]
                            coverage = parts[5]

                            # Format: sample_name<tab>orf_id<tab>mean_coverage
                            out_fh.write(f"{sample_name}\t{orf_id}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8 Concat and format reads coverage table
    ###########################################################################

    readscov_table = os.path.join(args.output_dir, "orfs_readscov.tsv")

    # Find and concat all *_orfs_readscov.tsv files in subdirectories
    readscov_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs_readscov.tsv"))

    if not readscov_files:
        print(f"No reads coverage files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(readscov_table, "w", encoding="utf-8") as out_fh:
            for readscov_file in readscov_files:
                sample_name = os.path.basename(readscov_file).replace("_orfs_readscov.tsv", "")

                with open(readscov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            orf_id = parts[4]
                            coverage = parts[5]

                            # Format: sample_name<tab>orf_id<tab>reads_coverage
                            out_fh.write(f"{sample_name}\t{orf_id}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format reads coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-3.py exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()
