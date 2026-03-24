#!/usr/bin/env python3
"""
Python rewrite of mg-clust_module-4.bash

Replicates the original behavior assuming execution inside the conda environment
"mg-clust-module-4" (or equivalent) where dependencies are available on PATH.

- Accepts the same CLI options
- Executes mmseqs for clustering with equivalent parameters
- Performs the same data processing and formatting
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import shutil
import subprocess
import sys
import pandas as pd
from typing import List, Optional 

mmseqs = "mmseqs"

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
        description="mg-clust module 4", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=TEST_VALUES["clust_thres"],
        help="clustering threshold, passed as -t to mmseqs cluster (default: 0.7)")

    parser.add_argument("--clust_cov_len", dest="clust_cov_len", type=float, default=0.85,
        help="minimum fraction of aligned residues for clustering, passed as -c to mmseqs cluster (default: 0.85)")
    
    parser.add_argument("--meancov_table", dest="meancov_table", default=TEST_VALUES["cov_table"],
        help="ORFs' mean coverage table")
    
    parser.add_argument("--readscov_table", dest="readscov_table", default=TEST_VALUES["cov_table"],
        help="ORFs' reads coverage table")
    
    parser.add_argument("--nslots", dest="nslots", type=int, default=TEST_VALUES["nslots"],
        help="number of threads used (default: 4)")
    
    parser.add_argument("--orfs_db", dest="orfs_db", default=TEST_VALUES["orfs_db"],
        help="mmseqs orfs db")
    
    parser.add_argument("--output_dir", dest="output_dir", default=TEST_VALUES["output_dir"],
        help="directory to output generated data (default: mg-clust_output-2)")
    
    parser.add_argument("--sample_name", dest="sample_name", default=TEST_VALUES["sample_name"],
        help="sample name used to name the files")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

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

    check_tools([mmseqs])
    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    ensure_file(args.meancov_table, "mean coverage table")
    ensure_file(args.readscov_table, "reads coverage table")
    ensure_file(args.orfs_db, "ORFs database")

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
    # 3.4. Run orfs clustering
    ###########################################################################

    # Create subdirectory with threshold in name (remove decimal point)
    clust_thres_str = str(args.clust_thres * 100).rstrip("0").rstrip(".") + "perc"
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
                "-c", str(args.clust_cov_len)
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

    ###########################################################################
    # 3.5. Convert to tab tables
    ###########################################################################

    mmseqs_clust_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}.tsv")

    try:
        run(
            [
                mmseqs,
                "createtsv",
                args.orfs_db,
                args.orfs_db,
                clust_db,
                mmseqs_clust_table
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createtsv failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Map mean coverage to clusters
    ###########################################################################

    clust2meancov_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}2meancov.tsv")
    not_found_list_path = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}_not_found.list")

    try:
        meancov_table_df = pd.read_csv(args.meancov_table, sep="\t", header=None,
                             names=["sample_name", "seq_id", "abund"])

        mmseqs_clust_table_df = pd.read_csv(mmseqs_clust_table, sep="\t", header=None,
                                            names=["clust_id", "orf_id"])

        # Construct orf_id to match the format used in mg-clust-module-3.py
        # when ORFs fasta files are concatenated and formatted with sample name prefix
        meancov_table_df["orf_id"] = meancov_table_df["sample_name"] + "|" + meancov_table_df["seq_id"]

        # Left join to map cluster IDs; ORFs filtered by length in module 3 will have NaN clust_id
        merged_df = meancov_table_df.merge(
            mmseqs_clust_table_df,
            on="orf_id",
            how="left"
        )

        # Split found / not found
        clust2meancov_table_df = merged_df[merged_df["clust_id"].notna()]
        notfound_df = merged_df[merged_df["clust_id"].isna()]

        clust2meancov_table_df[["sample_name", "clust_id", "orf_id", "abund"]].to_csv(
            clust2meancov_table, sep="\t", header=False, index=False)

        notfound_df["orf_id"].to_csv(not_found_list_path, header=False, index=False)

    except Exception as exc:
        print(f"Map cluster to mean coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Map reads coverage to clusters
    ###########################################################################

    clust2readscov_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}2readscov.tsv")

    try:
        readscov_table_df = pd.read_csv(args.readscov_table, sep="\t", header=None,
                             names=["sample_name", "seq_id", "abund"])

        # Construct orf_id to match the format used in mg-clust-module-3.py
        # when ORFs fasta files are concatenated and formatted with sample name prefix
        readscov_table_df["orf_id"] = readscov_table_df["sample_name"] + "|" + readscov_table_df["seq_id"]

        # Left join to map cluster IDs; ORFs filtered by length in module 3 will have NaN clust_id
        merged_df = readscov_table_df.merge(
            mmseqs_clust_table_df,
            on="orf_id",
            how="left"
        )

        # Split found / not found
        clust2readscov_table_df = merged_df[merged_df["clust_id"].notna()]

        clust2readscov_table_df[["sample_name", "clust_id", "orf_id", "abund"]].to_csv(
            clust2readscov_table, sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Map cluster to reads coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)
        
    ###########################################################################
    # 3.7. Create workable: Collapsed mean coverage table
    ###########################################################################

    clust2meancov_table_workable = os.path.join(args.output_dir, f"orfs_clust_id{clust_thres_str}_meancov_workable.tsv")

    try:
        # Sum mean coverage per sample + cluster
        clust2meancov_table_workable_df = clust2meancov_table_df.groupby(["sample_name", "clust_id"], 
                                                                       as_index=False)["abund"].sum()
        clust2meancov_table_workable_df.to_csv(clust2meancov_table_workable, 
                                             sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Collapse and filter opus failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8. Create workable: Collapsed reads coverage table
    ###########################################################################

    clust2readscov_table_workable = os.path.join(args.output_dir, f"orfs_clust_id{clust_thres_str}_readscov_workable.tsv")

    try:
        # Sum reads coverage per sample + cluster
        clust2readscov_table_workable_df = clust2readscov_table_df.groupby(["sample_name", "clust_id"], 
                                                                       as_index=False)["abund"].sum()
        clust2readscov_table_workable_df.to_csv(clust2readscov_table_workable, 
                                             sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Collapse and filter opus failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-4.py exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()

