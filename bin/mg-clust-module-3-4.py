#!/usr/bin/env python3
"""
mg-clust module 3-4: Concatenation, filtering, MMseqs2 DB creation, ORF clustering,
and abundance table generation.

Assumes execution inside the conda environment "mg-clust-module-3-4" (or equivalent)
where dependencies are available on PATH.

- Concatenates ORF protein sequences from all samples, prefixing each header with the sample name
- Filters ORFs by minimum length using bbduk
- Creates an MMseqs2 sequence database from the filtered ORFs
- Concatenates per-sample mean coverage tables, prefixing each ORF ID with the sample name
- Concatenates per-sample reads coverage tables, prefixing each ORF ID with the sample name
- Clusters filtered ORFs using MMseqs2 at a given sequence identity threshold
- Maps mean coverage and reads coverage tables to cluster IDs
- Collapses per-ORF abundance to per-cluster abundance summed across samples
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import sys, os
import subprocess
import shutil
import pandas as pd
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
        description="mg-clust module 3-4", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--orf_files", dest="orf_files", required=True, nargs="+",
        help="list of per-sample ORF fasta files (*_orfs.faa)")

    parser.add_argument("--meancov_files", dest="meancov_files", required=True, nargs="+",
        help="list of per-sample mean coverage tables (*_orfs_meancov.tsv)")

    parser.add_argument("--readscov_files", dest="readscov_files", required=True, nargs="+",
        help="list of per-sample reads coverage tables (*_orfs_readscov.tsv)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads used (default: 4)")

    parser.add_argument("--min_orf_length", dest="min_orf_length", type=int, default=60,
        help="minimum length of ORFs (amino acids); ORFs shorter than this will be discarded (default: 60)")

    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=0.7,
        help="clustering threshold, passed as --min-seq-id to mmseqs cluster (default: 0.7)")

    parser.add_argument("--clust_cov_len", dest="clust_cov_len", type=float, default=0.85,
        help="minimum fraction of aligned residues for clustering, passed as -c to mmseqs cluster (default: 0.85)")

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
    # 3.1. Check output directory
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
    # 3.2. Create output directory
    ###########################################################################

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except Exception:
        print(f"mkdir {args.output_dir} failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.3. Concat ORF files
    ###########################################################################

    concat_orfs = os.path.join(args.output_dir, "orfs.faa")

    open(concat_orfs, "w").close()

    for orf_file in args.orf_files:
        # sample_name is the filename without the "_orfs.faa" suffix
        # set when running FragGeneScanRs in mg-clust-module-2 with the -o argument
        sample_name = os.path.basename(orf_file).replace("_orfs.faa", "")

        try:
            with open(concat_orfs, "a") as out_fh:
                # here sed is used instead of python to reduce computation time
                proc = subprocess.run(
                    ["sed", f"s/^/>{ sample_name}|/", orf_file],
                    stdout=out_fh, check=False)
            if proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, proc.args)

        except subprocess.CalledProcessError:
            print(f"prefixing headers in {orf_file} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.4. Filter ORFs by length
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
    # 3.5. Create mmseqs database
    ###########################################################################

    orfs_db = concat_orfs_filt.replace(".faa", "_db")

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
    # 3.6. Concat and format mean coverage table
    ###########################################################################

    meancov_table = os.path.join(args.output_dir, "orfs_meancov.tsv")

    try:
        with open(meancov_table, "w", encoding="utf-8") as out_fh:
            for meancov_file in args.meancov_files:

                # sample_name is the filename without the "_orfs_meancov.tsv" suffix
                # set when running bedtools in mg-clust-module-2
                sample_name = os.path.basename(meancov_file).replace("_orfs_meancov.tsv", "")

                with open(meancov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            orf_id = parts[4]
                            coverage = parts[5]
                            out_fh.write(f"{sample_name}\t{orf_id}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format mean coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.7. Concat and format reads coverage table
    ###########################################################################

    readscov_table = os.path.join(args.output_dir, "orfs_readscov.tsv")

    try:
        with open(readscov_table, "w", encoding="utf-8") as out_fh:
            for readscov_file in args.readscov_files:
                sample_name = os.path.basename(readscov_file).replace("_orfs_readscov.tsv", "")

                with open(readscov_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 6:
                            orf_id = parts[4]
                            coverage = parts[5]
                            out_fh.write(f"{sample_name}\t{orf_id}\t{coverage}\n")
    except Exception as exc:
        print(f"concat and format reads coverage table failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.8. Run ORF clustering
    ###########################################################################

    clust_thres_str = str(args.clust_thres * 100).rstrip("0").rstrip(".") + "perc"
    clust_dir = os.path.join(args.output_dir, f"clust_orfs_id{clust_thres_str}")

    try:
        os.makedirs(clust_dir, exist_ok=False)
    except Exception:
        print(f"mkdir {clust_dir} failed", file=sys.stderr)
        sys.exit(1)

    tmp_dir = os.path.join(clust_dir, "tmp")
    clust_db = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}")

    try:
        run(
            [
                mmseqs,
                "cluster",
                orfs_db,
                clust_db,
                tmp_dir,
                "--min-seq-id", str(args.clust_thres),
                "--threads", str(args.nslots),
                "--cov-mode", "0",
                "-c", str(args.clust_cov_len)
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs cluster failed", file=sys.stderr)
        sys.exit(1)

    if os.path.isdir(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            print(f"rm {tmp_dir} failed", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.9. Convert clustering results to TSV
    ###########################################################################

    mmseqs_clust_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}.tsv")

    try:
        run(
            [
                mmseqs,
                "createtsv",
                orfs_db,
                orfs_db,
                clust_db,
                mmseqs_clust_table
            ]
        )
    except subprocess.CalledProcessError:
        print("mmseqs createtsv failed", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.10. Map mean coverage to clusters
    ###########################################################################

    clust2meancov_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}2meancov.tsv")
    not_found_list_path = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}_not_found.list")

    try:
        meancov_table_df = pd.read_csv(meancov_table, sep="\t", header=None,
                             names=["sample_name", "seq_id", "abund"])

        mmseqs_clust_table_df = pd.read_csv(mmseqs_clust_table, sep="\t", header=None,
                                            names=["clust_id", "orf_id"])

        # Construct orf_id to match the prefixed format created in section 3.3
        meancov_table_df["orf_id"] = meancov_table_df["sample_name"] + "|" + meancov_table_df["seq_id"]

        # Left join; ORFs filtered by length will have NaN clust_id
        merged_df = meancov_table_df.merge(
            mmseqs_clust_table_df,
            on="orf_id",
            how="left"
        )

        clust2meancov_table_df = merged_df[merged_df["clust_id"].notna()]
        notfound_df = merged_df[merged_df["clust_id"].isna()]

        clust2meancov_table_df[["sample_name", "clust_id", "orf_id", "abund"]].to_csv(
            clust2meancov_table, sep="\t", header=False, index=False)

        notfound_df["orf_id"].to_csv(not_found_list_path, header=False, index=False)

    except Exception as exc:
        print(f"Map cluster to mean coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.11. Map reads coverage to clusters
    ###########################################################################

    clust2readscov_table = os.path.join(clust_dir, f"orfs_clust_id{clust_thres_str}2readscov.tsv")

    try:
        readscov_table_df = pd.read_csv(readscov_table, sep="\t", header=None,
                             names=["sample_name", "seq_id", "abund"])

        # Construct orf_id to match the prefixed format created in section 3.3
        readscov_table_df["orf_id"] = readscov_table_df["sample_name"] + "|" + readscov_table_df["seq_id"]

        # Left join; ORFs filtered by length will have NaN clust_id
        merged_df = readscov_table_df.merge(
            mmseqs_clust_table_df,
            on="orf_id",
            how="left"
        )

        clust2readscov_table_df = merged_df[merged_df["clust_id"].notna()]

        clust2readscov_table_df[["sample_name", "clust_id", "orf_id", "abund"]].to_csv(
            clust2readscov_table, sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Map cluster to reads coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.12. Create workable: Collapsed mean coverage table
    ###########################################################################

    clust2meancov_table_workable = os.path.join(
        args.output_dir, f"orfs_clust_id{clust_thres_str}_meancov_workable.tsv")

    try:
        clust2meancov_table_workable_df = clust2meancov_table_df.groupby(
            ["sample_name", "clust_id"], as_index=False)["abund"].sum()
        clust2meancov_table_workable_df.to_csv(
            clust2meancov_table_workable, sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Collapse mean coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.13. Create workable: Collapsed reads coverage table
    ###########################################################################

    clust2readscov_table_workable = os.path.join(
        args.output_dir, f"orfs_clust_id{clust_thres_str}_readscov_workable.tsv")

    try:
        clust2readscov_table_workable_df = clust2readscov_table_df.groupby(
            ["sample_name", "clust_id"], as_index=False)["abund"].sum()
        clust2readscov_table_workable_df.to_csv(
            clust2readscov_table_workable, sep="\t", header=False, index=False)

    except Exception as exc:
        print(f"Collapse reads coverage failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print("mg-clust_module-3-4.py exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()
