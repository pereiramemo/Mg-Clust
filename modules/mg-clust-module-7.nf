// ─────────────────────────────────────────────────────────────────────────────
// MODULE 7: Merge taxonomy, function, and abundance tables
// Input:  per-sample tax/fun workable TSVs (collected) + module 4 workable tables
// Output: merged orfs_clust_id<N>perc_meancov2taxa2fun_workable.tsv
// ─────────────────────────────────────────────────────────────────────────────

process MODULE7 {
    container "ghcr.io/epereira/mg-clust/module-7:latest"
    publishDir "${params.output_dir}/",
                mode: "copy"

    input:
    path(tax_files)
    path(fun_files)
    path(meancov_table)
    path(readscov_table)

    output:
    path("module7")

    script:
    """
    mg-clust-module-7.py \
        --tax_files      ${tax_files} \
        --fun_files      ${fun_files} \
        --meancov_table  ${meancov_table} \
        --readscov_table ${readscov_table} \
        --clust_thres    ${params.clust_thres} \
        --output_dir     module7 \
        --overwrite
    """
}