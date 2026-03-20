#!/usr/bin/env python3
"""
Python rewrite of mg-clust_module-2.bash

Replicates the original behavior assuming execution inside the conda environment
"mg-clust-module-2" (or equivalent) where dependencies are available on PATH.

- Accepts the same CLI options
- Executes FragGeneScan, bedtools with equivalent parameters
- Performs the same filesystem operations and cleanup
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import re
import shutil
import subprocess
import sys
from typing import List, Optional

# Test values for development/debugging using test/data
TEST_VALUES = {
    "assembly_file": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-1/assembly/test_sample.contigs.fa",
    "bam_file": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-1/test_sample_sorted_markdup.bam",
    "sample_name": "test_sample",
    "output_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-2",
    "nslots": 4,
    "train_file_name": "illumina_1",
    "overwrite": True
}

###############################################################################
# 2. Define utility functions
###############################################################################

# 2.1 Run a command and check return code; optionally redirect stdout to file
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

# 2.2 Parse command-line arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="mg-clust module 2 (Python)", add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--assembly_file", dest="assembly_file", default=TEST_VALUES["assembly_file"],
        help="input assembled metagenome (fasta format)")
    
    parser.add_argument("--bam_file", dest="bam_file", default=TEST_VALUES["bam_file"],
        help="bam input file (reads mapped to contigs)")
    
    parser.add_argument("--nslots", dest="nslots", type=int, default=TEST_VALUES["nslots"],
        help="number of threads used (default: 4)")
    
    parser.add_argument("--output_dir", dest="output_dir", default=TEST_VALUES["output_dir"],
        help="directory to output generated data (default: test_output)")
    
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=TEST_VALUES["overwrite"],
        help="overwrite previous folder if present (default: True)")
    
    parser.add_argument("--sample_name", dest="sample_name", default=TEST_VALUES["sample_name"],
        help="sample name used to name the files")
    
    parser.add_argument("--train_file_name", dest="train_file_name", default=TEST_VALUES["train_file_name"],
        help="train file name used to run FragGeneScan (default: illumina_1)")

    return parser.parse_args()

# 2.3 Ensure a file exists; exit with error if not
def ensure_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"{label} is not a real file", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:
    # Fail early on command failures by checking return codes
    # Tools are expected to be available on PATH via the active conda env
    fraggenescan = "FragGeneScan"
    bedtools = "bedtools"

    args = parse_args()

    # 4) Check mandatory files
    ensure_file(args.assembly_file, "input assembly file")
    ensure_file(args.bam_file, "input bam file")

    # 5) Create output folder
    if os.path.isdir(args.output_dir):
        if not args.overwrite:
            print(f"{args.output_dir} already exists; use --overwrite to overwrite")
            sys.exit(0)
        try:
            shutil.rmtree(args.output_dir)
        except Exception:
            print(f"rm -r output directory {args.output_dir} failed", file=sys.stderr)
            sys.exit(1)

    # Create directory (either new or after removal)
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir output directory {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    # 6) Predict ORFs
    # Note: FragGeneScan requires FRAGGENESCAN_TRAIN_FILE_DIR env variable
    # This should be set in the conda environment
    fraggenescan_train_dir = os.environ.get("FRAGGENESCAN_TRAIN_FILE_DIR", "")

    try:
        run(
            [
                fraggenescan,
                "-s", args.assembly_file,
                "-o", os.path.join(args.output_dir, f"{args.sample_name}_orfs"),
                "-w", "0",
                "--unordered",
                "-p", str(args.nslots),
                "-t", args.train_file_name,
                "-r", fraggenescan_train_dir,
            ]
        )
    except subprocess.CalledProcessError:
        print("FragGeneScan failed", file=sys.stderr)
        sys.exit(1)

    # 7) Create BED file
    ffn_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.ffn")
    bed_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.bed")

    # Parse FFN headers and create BED file
    # Format: >contig_prefix_contig_id_start_end_strand
    # BED format: contig_id, start-1 (zero-based), end (one-based)
    try:
        with open(ffn_file, "r", encoding="utf-8", errors="ignore") as fh_in, \
             open(bed_file, "w", encoding="utf-8") as fh_out:
            for line in fh_in:
                if line.startswith(">"):
                    # Remove '>' and split by '_'
                    parts = line[1:].strip().split("_")
                    if len(parts) >= 4:
                        contig_prefix = parts[0]
                        contig_id = parts[1]
                        start = int(parts[2])
                        end = int(parts[3])

                        # BED format is zero-based for start, one-based for end
                        fh_out.write(f"{contig_prefix}_{contig_id}\t{start - 1}\t{end}\n")
    except Exception as exc:
        print(f"Creating bed file failed: {exc}", file=sys.stderr)
        sys.exit(1)

    # 8) Get ORFs coverage
    tmp_cov_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs_tmp.cov")

    try:
        run(
            [
                bedtools,
                "coverage",
                "-a", bed_file,
                "-b", args.bam_file,
                "-mean"
            ],
            stdout_path=tmp_cov_file
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute coverage failed", file=sys.stderr)
        sys.exit(1)

    # 9) Fix ORFs start coordinates (back to one-based)
    cov_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.cov")

    try:
        with open(tmp_cov_file, "r", encoding="utf-8") as fh_in, \
             open(cov_file, "w", encoding="utf-8") as fh_out:
            for line in fh_in:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    contig_id = parts[0]
                    start_coord = int(parts[1]) + 1
                    end_coord = parts[2]
                    abund = parts[3]

                    fh_out.write(f"{contig_id}\t{start_coord}\t{end_coord}\t{abund}\n")
    except Exception as exc:
        print(f"Formatting start coords failed: {exc}", file=sys.stderr)
        sys.exit(1)

    # Remove temporary file
    try:
        os.remove(tmp_cov_file)
    except Exception:
        print(f"removing {args.sample_name}_orfs_tmp.cov failed", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-2.py exited successfully")
    sys.exit(0)


if __name__ == "__main__":
    main()
