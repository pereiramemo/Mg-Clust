#!/usr/bin/env python3
"""
Python rewrite of mg-clust_module-3.bash

Replicates the original behavior assuming execution inside the conda environment
"mg-clust-module-3" (or equivalent) where dependencies are available on PATH.

- Accepts the same CLI options
- Executes bbduk, mmseqs with equivalent parameters
- Performs the same filesystem operations and concatenation
"""

import argparse
import glob
import os
import subprocess
import sys
from typing import List, Optional

# Test values for development/debugging using test/data
TEST_VALUES = {
    "input_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output",
    "output_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-3",
    "nslots": 4,
    "min_orf_length": 60,
}


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="mg-clust module 3 (Python)", add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")
    parser.add_argument("--input_dir", dest="input_dir", default=TEST_VALUES["input_dir"],
        help="input dir where to find the ORF fasta files and coverage tables")
    parser.add_argument("--nslots", dest="nslots", type=int, default=TEST_VALUES["nslots"],
        help="number of threads used (default: 4)")
    parser.add_argument("--min_orf_length", dest="min_orf_length", type=int, default=TEST_VALUES["min_orf_length"],
        help="minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default: 60)")
    parser.add_argument("--output_dir", dest="output_dir", default=TEST_VALUES["output_dir"],
        help="directory to output generated data")

    return parser.parse_args()


def main() -> None:
    # Fail early on command failures by checking return codes
    # Tools are expected to be available on PATH via the active conda env
    bbduk = "bbduk.sh"
    mmseqs = "mmseqs"

    args = parse_args()

    # Create output directory if it doesn't exist
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir output directory {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    # 4) Concat ORF files
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

    # 5) Filter ORFs by length
    concat_orfs_filtered = concat_orfs.replace(".faa", "_filtered.faa")

    try:
        run(
            [
                bbduk,
                f"in={concat_orfs}",
                f"out={concat_orfs_filtered}",
                "overwrite=t",
                f"threads={args.nslots}",
                f"minlength={args.min_orf_length}",
                "amino=t"
            ]
        )
    except subprocess.CalledProcessError:
        print("bbduk failed", file=sys.stderr)
        sys.exit(1)

    # 6) Create mmseqs database
    orfs_db = concat_orfs_filtered.replace(".faa", "_db")

    if not os.path.exists(orfs_db):
        try:
            run(
                [
                    mmseqs,
                    "createdb",
                    concat_orfs_filtered,
                    orfs_db,
                    "--dbtype", "1"
                ]
            )
        except subprocess.CalledProcessError:
            print("mmseqs createdb failed", file=sys.stderr)
            sys.exit(1)

    # 7) Concat and format coverage table
    cov_table = os.path.join(args.output_dir, "orfs_cov.tsv")

    # Find all *_orfs.cov files in subdirectories
    cov_files = glob.glob(os.path.join(args.input_dir, "*", "*_orfs.cov"))

    if not cov_files:
        print(f"No coverage files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        with open(cov_table, "w", encoding="utf-8") as out_fh:
            for cov_file in cov_files:
                sample = os.path.basename(cov_file).replace("_orfs.cov", "")

                # Read coverage file and reformat
                with open(cov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            contig_id = parts[0]
                            start = parts[1]
                            end = parts[2]
                            coverage = parts[3]

                            # Format: SAMPLE-contig_id_start_end<tab>coverage
                            out_fh.write(f"{sample}-{contig_id}_{start}_{end}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-3.py exited successfully")
    sys.exit(0)


if __name__ == "__main__":
    main()
