#!/usr/bin/env python3
"""
mg-clust module 6: Merge taxonomy, function, and abundance tables.

Assumes execution inside the conda environment "mg-clust-module-6" (or equivalent)
where duckdb is importable.

- Concatenates per-sample ORF taxonomy annotation files (*_orf_tax_annot_workable.tsv)
  into a single orf_tax_annot_workable.tsv
- Concatenates per-sample ORF function annotation files (*_orf_fun_annot_workable.tsv)
  into a single orf_fun_annot_workable.tsv
- Merges the clustering abundance tables (meancov, readscov) from module 4 with the
  concatenated taxonomy and function tables using DuckDB, joining on cluster representative
  ORF ID (clust_id = orf_id), to produce orfs_clust_id<N>perc_meancov2taxa2fun_workable.tsv
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import sys, os
import shutil
import duckdb
sys.path.insert(0, os.path.dirname(__file__))
from utils import check_file

###############################################################################
# 2. Define utility functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Merge taxonomy, function, and OPU abundance tables", add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--tax_files", dest="tax_files", required=True, nargs="+",
        help="per-sample ORF taxonomy annotation files (*_orf_tax_annot_workable.tsv)")

    parser.add_argument("--fun_files", dest="fun_files", required=True, nargs="+",
        help="per-sample ORF function annotation files (*_orf_fun_annot_workable.tsv)")

    parser.add_argument("--meancov_table", dest="meancov_table", required=True,
        help="collapsed mean coverage workable table from module 4 "
             "(orfs_clust_id*_meancov_workable.tsv)")

    parser.add_argument("--readscov_table", dest="readscov_table", required=True,
        help="collapsed reads coverage workable table from module 4 "
             "(orfs_clust_id*_readscov_workable.tsv)")

    parser.add_argument("--clust_thres", dest="clust_thres", type=float, default=0.7,
        help="clustering threshold used in module 4; used to name the output file (default: 0.7)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 3. Define the main function
###############################################################################

def main() -> None:

    args = parse_args()

    ###########################################################################
    # 3.1. Check mandatory files
    ###########################################################################

    for f in args.tax_files:
        check_file(f, "taxonomy annotation file")
    for f in args.fun_files:
        check_file(f, "function annotation file")
    check_file(args.meancov_table, "mean coverage workable table")
    check_file(args.readscov_table, "reads coverage workable table")

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
    # 3.4. Concatenate taxonomy annotation files
    ###########################################################################

    orf_tax_annot = os.path.join(args.output_dir, "orf_tax_annot_workable.tsv")

    try:
        with open(orf_tax_annot, "w", encoding="utf-8") as out_fh:
            for tax_file in args.tax_files:
                with open(tax_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        out_fh.write(line)
    except Exception as exc:
        print(f"concatenating taxonomy files failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.5. Concatenate function annotation files
    ###########################################################################

    orf_fun_annot = os.path.join(args.output_dir, "orf_fun_annot_workable.tsv")

    try:
        with open(orf_fun_annot, "w", encoding="utf-8") as out_fh:
            for fun_file in args.fun_files:
                with open(fun_file, "r", encoding="utf-8") as in_fh:
                    for line in in_fh:
                        out_fh.write(line)
    except Exception as exc:
        print(f"concatenating function files failed: {exc}", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.6. Merge tables with DuckDB
    ###########################################################################

    clust_thres_str = str(args.clust_thres * 100).rstrip("0").rstrip(".") + "perc"
    merged_table = os.path.join(
        args.output_dir,
        f"orfs_clust_id{clust_thres_str}_cov2taxa2fun_workable.tsv"
    )

    try:
        con = duckdb.connect()

        con.execute("""
            CREATE TABLE meancov AS
            SELECT column0 AS sample_name, column1 AS clust_id, column2::DOUBLE AS abund
            FROM read_csv(?, delim='\t', header=false)
        """, [args.meancov_table])

        con.execute("""
            CREATE TABLE readscov AS
            SELECT column0 AS sample_name, column1 AS clust_id, column2::DOUBLE AS abund
            FROM read_csv(?, delim='\t', header=false)
        """, [args.readscov_table])

        # tax columns: sample_name, orf_id, taxid, rank, name, lineage
        con.execute("""
            CREATE TABLE tax AS
            SELECT column0 AS sample_name, column1 AS orf_id,
                   column2 AS taxid, column3 AS rank, column4 AS name, column5 AS lineage
            FROM read_csv(?, delim='\t', header=false)
        """, [orf_tax_annot])

        # fun columns: sample_name, orf_id, ko_id, score, evalue
        con.execute("""
            CREATE TABLE fun AS
            SELECT column0 AS sample_name, column1 AS orf_id,
                   column2 AS ko_id, column3::DOUBLE AS score, column4::DOUBLE AS evalue
            FROM read_csv(?, delim='\t', header=false)
        """, [orf_fun_annot])

        # Left join on clust_id = orf_id: taxonomy and function are looked up
        # from the cluster representative ORF (clust_id encodes the representative's
        # sample and ORF ID in sample_name|orf_id format, matching the orf_id column
        # in the tax and fun tables).
        con.execute(f"""
            COPY (
                SELECT
                    meancov.sample_name,
                    meancov.clust_id,
                    meancov.abund           AS meancov,
                    readscov.abund          AS readscov,
                    tax.taxid,
                    tax.rank,
                    tax.name,
                    tax.lineage,
                    fun.ko_id,
                    fun.score               AS fun_score,
                    fun.evalue              AS fun_evalue
                FROM meancov
                LEFT JOIN readscov
                    ON  meancov.sample_name = readscov.sample_name
                    AND meancov.clust_id    = readscov.clust_id
                LEFT JOIN tax ON meancov.clust_id = tax.orf_id
                LEFT JOIN fun ON meancov.clust_id = fun.orf_id
                ORDER BY meancov.sample_name, meancov.clust_id
            ) TO '{merged_table}' (DELIMITER '\t', HEADER true)
        """)

        con.close()

    except Exception as exc:
        print(f"DuckDB merge failed: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"{os.path.basename(__file__)} exited successfully")
    sys.exit(0)

###########################################################################
# 4. Run the main function
###########################################################################

if __name__ == "__main__":
    main()