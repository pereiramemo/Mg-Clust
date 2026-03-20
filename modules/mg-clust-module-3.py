#!/usr/bin/env python3
"""
Python rewrite of mg-clust_module-3.bash

Replicates the original behavior assuming execution inside the conda environment
"mg-clust-module-3" (or equivalent) where dependencies are available on PATH.

- Accepts the same CLI options
- Executes bbduk, mmseqs with equivalent parameters
- Performs the same filesystem operations and concatenation
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import glob
import os
import shutil
import subprocess
import sys
from typing import List, Optional

bbduk = "bbduk.sh"
mmseqs = "mmseqs"

# Test values for development/debugging using test/data
TEST_VALUES = {
    "input_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output",
    "output_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-3",
    "nslots": 4,
    "min_orf_length": 60,
}

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Run a command and check return code; optionally redirect stdout to file
###############################################################################

def run(cmd: List[str], stdout_path: Optional[str] = None) -> None:
    """Run a command, stream output or redirect to file, and fail on non-zero return code."""
    try:
        if stdout_path:
            with open(stdout_path, "w") as out_f:
                proc = subprocess.run(cmd, stdout=out_f, check=False)
        else:
            proc = subprocess.run(cmd, check=False)
    except FileNotFoundError as exc:
        print(f"Command not found: {cmd[0]} ({exc})", file=sys.stderr)
        sys.exit(1)

    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)

###############################################################################
# 2.2 Check that required tools are available on PATH
###############################################################################

def check_tools(tools: List[str]) -> None:
    missing = [t for t in tools if subprocess.run(["which", t], capture_output=True).returncode != 0]
    if missing:
        print(f"Missing tools: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.3 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description="mg-clust module 3", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_dir", dest="input_dir", default=TEST_VALUES["input_dir"],
        help="input dir where to find the ORF fasta files and coverage tables")
    
    parser.add_argument("--nslots", dest="nslots", type=int, default=TEST_VALUES["nslots"],
        help="number of threads used (default: 4)")
    
    parser.add_argument("--min_orf_length", dest="min_orf_length", type=int, default=TEST_VALUES["min_orf_length"],
        help="minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default: 60)")
    
    parser.add_argument("--output_dir", dest="output_dir", default=TEST_VALUES["output_dir"],
        help="directory to output generated data")
    
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=TEST_VALUES["overwrite"],
        help="overwrite previous folder if present (default: True)")

    return parser.parse_args()

###############################################################################
# 2.4 Ensure a file exists; exit with error if not
###############################################################################

def ensure_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"{label} is not a real file", file=sys.stderr)
        sys.exit(1)

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

    # Find all *_orfs.faa files in subdirectories
    orf_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs.faa"))

    if not orf_files:
        print(f"No ORF files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(concat_orfs, "w", encoding="utf-8") as out_fh:
            for orf_file in orf_files:
                sample = os.path.basename(orf_file).replace("_orfs.faa", "")

                # Read file and prefix headers with sample name
                with open(orf_file, "r", encoding="utf-8", errors="ignore") as in_fh:
                    for line in in_fh:
                        if line.startswith(">"):
                            # Add sample prefix to header
                            out_fh.write(f">{sample}-{line[1:]}")
                        else:
                            out_fh.write(line)
    except Exception as exc:
        print(f"concat orf files failed: {exc}", file=sys.stderr)
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

    # Find all *_orfs_meancov.tsv files in subdirectories
    meancov_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs_meancov.tsv"))

    if not meancov_files:
        print(f"No coverage files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(meancov_table, "w", encoding="utf-8") as out_fh:
            for meancov_file in meancov_files:
                sample = os.path.basename(meancov_file).replace("_orfs_meancov.tsv", "")

                # Read coverage file and reformat
                with open(meancov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            orf_id = parts[4]
                            coverage = parts[5]

                            # Format: SAMPLE-contig_id_start_end<tab>mean_coverage
                            out_fh.write(f"{sample}-{orf_id}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8 Concat and format reads coverage table
    ###########################################################################

    readscov_table = os.path.join(args.output_dir, "orfs_readscov.tsv")

    # Find all *_orfs_readscov.tsv files in subdirectories
    readscov_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs_readscov.tsv"))

    if not readscov_files:
        print(f"No reads coverage files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(readscov_table, "w", encoding="utf-8") as out_fh:
            for readscov_file in readscov_files:
                sample = os.path.basename(readscov_file).replace("_orfs_readscov.tsv", "")

                # Read coverage file and reformat
                with open(readscov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            orf_id = parts[4]
                            coverage = parts[5]

                            # Format: SAMPLE-contig_id_start_end<tab>reads_coverage
                            out_fh.write(f"{sample}-{orf_id}\t{coverage}\n")

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
