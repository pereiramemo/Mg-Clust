#!/usr/bin/env python3
"""
mg-clust module 2: ORF prediction and coverage estimation.

Assumes execution inside the conda environment "mg-clust-module-2" (or equivalent)
where dependencies are available on PATH.

- Predicts ORFs from an assembled metagenome using FragGeneScanRs
- Converts ORF coordinates from the FragGeneScanRs .out file to BED format,
  converting 1-based start coordinates to 0-based as required by bedtools
- Computes the number of reads per ORF using bedtools coverage
- Computes mean depth per ORF using bedtools coverage -mean
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

fraggenescan = "FragGeneScanRs"
bedtools = "bedtools"

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
        description="mg-clust module 2", add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--assembly_file", dest="assembly_file", required=True,
        help="input assembled metagenome (fasta format)")

    parser.add_argument("--bam_file", dest="bam_file", required=True,
        help="bam input file (reads mapped to contigs)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name used to name the files")

    parser.add_argument("--train_file_name", dest="train_file_name", default="illumina_1",
        help="train file name used to run FragGeneScan (default: illumina_1)")

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

    check_tools([fraggenescan, bedtools])
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    ensure_file(args.assembly_file, "input assembly file")
    ensure_file(args.bam_file, "input bam file")

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
        print(f"mkdir output directory {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.4. Predict ORFs
    ###########################################################################

    try:
        run(
            [
                fraggenescan,
                "-s", args.assembly_file,
                "-o", os.path.join(args.output_dir, f"{args.sample_name}_orfs"),
                "-w", "0",
                "--unordered",  
                "-p", str(args.nslots),
                "-t", args.train_file_name
            ]
        )
    except subprocess.CalledProcessError:
        print("FragGeneScan failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Create BED file
    ###########################################################################

    out_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.out")
    bed_file = os.path.join(args.output_dir, f"{args.sample_name}_orfs.bed")

    # check that .out file was created
    ensure_file(out_file, "FragGeneScan output .out file")

    # Parse .out file and create BED file
    # .out format: contig_id    start    end    strand    ...
    # BED format: contig_id, start-1 (zero-based), end (one-based)

    try:
        with open(out_file, "r", encoding="utf-8") as fh_in, \
             open(bed_file, "w", encoding="utf-8") as fh_out:
            first_data_line = True
            for line in fh_in:
                if line.startswith("#") or not line.strip():
                    continue
                if first_data_line:
                    if not line.startswith(">"):
                        print(f"Unexpected format in {out_file}: \n"
                              f"first line does not start with '>'", file=sys.stderr)
                        sys.exit(1)
                    first_data_line = False
                parts = line.strip().split("\t")
                if re.match(r"^>(.+)$", parts[0]):
                    contig_id = parts[0][1:]  # remove leading '>'
                if not re.match(r"^>(.+)$", parts[0]) and len(parts) >= 3:
                    start = int(parts[0])
                    end = int(parts[1])
                    strand = parts[2]
                    orf_id = f"{contig_id}_{start}_{end}_{strand}"
                    fh_out.write(f"{contig_id}\t{start - 1}\t{end}\t{strand}\t{orf_id}\n")
    except Exception as exc:
        print(f"Creating bed file failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Get number of reads per ORF
    ###########################################################################

    # Run bedtools coverage -counts to get read counts per ORF
    bedtool_reads = os.path.join(args.output_dir, f"{args.sample_name}_orfs_readscov.tsv")
    try:
        run(
            [bedtools, "coverage", 
             "-a", bed_file, 
             "-b", args.bam_file,
             "-counts"],
            stdout_path=bedtool_reads
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute read counts failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Get mean coverage per ORF
    ###########################################################################

    # Run bedtools coverage -mean to get mean depth per ORF
    bedtools_mean = os.path.join(args.output_dir, f"{args.sample_name}_orfs_meancov.tsv")
    try:
        run(
            [bedtools, "coverage", "-a", bed_file, "-b", args.bam_file, "-mean"],
            stdout_path=bedtools_mean
        )
    except subprocess.CalledProcessError:
        print("bedtools to compute mean coverage failed", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-2.py exited successfully")
    sys.exit(0)

###########################################################################
# 4. Execute main function
###########################################################################

if __name__ == "__main__":
    main()
