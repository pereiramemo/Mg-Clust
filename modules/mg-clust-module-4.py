#!/usr/bin/env python3
"""
Python rewrite of mg-clust_module-4.bash

Replicates the original behavior assuming execution inside the conda environment
"mg-clust-module-4" (or equivalent) where dependencies are available on PATH.

- Accepts the same CLI options
- Executes mmseqs for clustering with equivalent parameters
- Performs the same data processing and formatting
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from typing import List, Optional

# Test values for development/debugging using test/data
TEST_VALUES = {
    "orfs_db": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-3/orfs_filtered_db",
    "cov_table": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-3/orfs_cov.tsv",
    "output_dir": "/home/epereira/workspace/repositories/tools/Mg-Clust/test/test_output/output-4",
    "sample_name": "test_sample",
    "nslots": 4,
    "clust_thres": 0.7,
    "min_opu_occup": 2,
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
        description="mg-clust module 4 (Python)", add_help=False
    )

    parser.add_argument("--help", action="help", help="print this help")
    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=TEST_VALUES["clust_thres"],
        help="clustering threshold (default: 0.7)")
    parser.add_argument("--cov_table", dest="cov_table", default=TEST_VALUES["cov_table"],
        help="ORFs' coverage table")
    parser.add_argument("--min_opu_occup", dest="min_opu_occup", type=int, default=TEST_VALUES["min_opu_occup"],
        help="minimum OPU occupancy (smaller than this will be discarded; default: 2)")
    parser.add_argument("--nslots", dest="nslots", type=int, default=TEST_VALUES["nslots"],
        help="number of threads used (default: 4)")
    parser.add_argument("--orfs_db", dest="orfs_db", default=TEST_VALUES["orfs_db"],
        help="mmseqs orfs db")
    parser.add_argument("--output_dir", dest="output_dir", default=TEST_VALUES["output_dir"],
        help="directory to output generated data (default: mg-clust_output-2)")
    parser.add_argument("--sample_name", dest="sample_name", default=TEST_VALUES["sample_name"],
        help="sample name used to name the files")

    return parser.parse_args()


def ensure_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"{label} is not a real file", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    # Fail early on command failures by checking return codes
    # Tools are expected to be available on PATH via the active conda env
    mmseqs = "mmseqs"

    args = parse_args()

    # Check mandatory files
    ensure_file(args.cov_table, "coverage table")
    if not os.path.exists(args.orfs_db):
        print(f"ORFs database {args.orfs_db} does not exist", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir output directory {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    # 4) Run orfs clustering
    # Create subdirectory with threshold in name (remove decimal point)
    clust_thres_str = str(args.clust_thres).replace("0.", "")
    clust_dir = os.path.join(args.output_dir, f"clust_orfs_id{clust_thres_str}")

    try:
        os.makedirs(clust_dir, exist_ok=False)
    except Exception:
        print(f"mkdir {clust_dir} failed", file=sys.stderr)
        sys.exit(1)

    tmp_dir = os.path.join(clust_dir, "tmp")
    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            print(f"rm -r {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

    clust_db = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}")

    try:
        run(
            [
                mmseqs,
                "cluster",
                args.orfs_db,
                clust_db,
                tmp_dir,
                "--min-seq-id", str(args.clust_thres),
                "--threads", str(args.nslots),
                "--cov-mode", "0",
                "-c", "0.85"
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs2 cluster failed", file=sys.stderr)
        sys.exit(1)

    # Remove tmp directory
    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            print(f"rm {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

    # 5) Convert to tab tables
    clust_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}.tsv")

    try:
        run(
            [
                mmseqs,
                "createtsv",
                args.orfs_db,
                args.orfs_db,
                clust_db,
                clust_table
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createtsv failed", file=sys.stderr)
        sys.exit(1)

    # 6) Cross tables
    clust2abund_table = os.path.join(clust_dir, f"orfs_clust2abund_id{clust_thres_str}.tsv")

    # Build cluster mapping (seq_id -> clust_id)
    array_clust = {}
    duplicated_list_path = os.path.join(clust_dir, f"duplicated_{clust_thres_str}.list")

    try:
        with open(clust_table, "r", encoding="utf-8") as fh, \
             open(duplicated_list_path, "w", encoding="utf-8") as dup_fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    clust_id = re.sub(r"_[+,-]$", "", parts[0])
                    seq_id = re.sub(r"_[+,-]$", "", parts[1])

                    if seq_id not in array_clust:
                        array_clust[seq_id] = clust_id
                    else:
                        # Duplicated seq_id (ORF in exact same location, opposite strand)
                        dup_fh.write(f"{parts[0]}\n")
    except Exception as exc:
        print(f"Processing cluster table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    # Map coverage to clusters
    not_found_list_path = os.path.join(clust_dir, f"not_found_{clust_thres_str}.list")

    try:
        with open(args.cov_table, "r", encoding="utf-8") as cov_fh, \
             open(clust2abund_table, "w", encoding="utf-8") as out_fh, \
             open(not_found_list_path, "w", encoding="utf-8") as nf_fh:
            for line in cov_fh:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    seq_id = parts[0]
                    abund = parts[1]

                    if seq_id in array_clust:
                        clust_id = array_clust[seq_id]
                        # Extract sample name from seq_id (format: SAMPLE-k###_###...)
                        sample_match = re.match(r"(.*)-k[0-9]+_[0-9]+.*", seq_id)
                        if sample_match:
                            sample = sample_match.group(1)
                            out_fh.write(f"{sample}\t{clust_id}\t{seq_id}\t{abund}\n")
                    else:
                        # seq_ids filtered by length in module 3 will not be found
                        nf_fh.write(f"{seq_id}\n")
    except Exception as exc:
        print(f"Map cluster to abund failed: {exc}", file=sys.stderr)
        sys.exit(1)

    # 7) Create workables: collapsed and filtered
    parent_dir = os.path.dirname(args.output_dir)
    clust2abund_table_workable = os.path.join(parent_dir, f"orfs_clust2abund_id{clust_thres_str}_workable.tsv")

    # Process with Python (equivalent to gawk script)
    clust2abund_array = defaultdict(lambda: defaultdict(float))
    clust2occup_array = defaultdict(lambda: defaultdict(int))

    try:
        with open(clust2abund_table, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    sample = parts[0]
                    opu_id = parts[1]
                    abund = float(parts[3])

                    clust2abund_array[sample][opu_id] += abund
                    clust2occup_array[opu_id][sample] = 1

        # Write filtered output
        with open(clust2abund_table_workable, "w", encoding="utf-8") as out_fh:
            for sample in clust2abund_array:
                for opu_id in clust2abund_array[sample]:
                    occup_total = len(clust2occup_array[opu_id])
                    abund_total = clust2abund_array[sample][opu_id]

                    if occup_total >= args.min_opu_occup:
                        out_fh.write(f"{sample}\t{opu_id}\t{abund_total}\n")

    except Exception as exc:
        print(f"collapse and filter opus failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-4.py exited successfully")
    sys.exit(0)


if __name__ == "__main__":
    main()
