// ─────────────────────────────────────────────────────────────────────────────
// MODULE 6: Merge taxonomy, function, and abundance tables
// Input:  per-sample tax/fun workable TSVs (collected) + module 4 workable tables
// Output: merged orfs_clust_id<N>perc_meancov2taxa2fun_workable.tsv
// ─────────────────────────────────────────────────────────────────────────────

process MODULE6 {

    container "ghcr.io/epereira/mg-clust/${task.process.toLowerCase().replaceFirst('module', 'module-')}:latest"
    publishDir "${params.output_dir}/${task.process.toLowerCase().replaceFirst('module', 'module-')}/",
           mode: "copy",
           enabled: params.full_output || params.stop_at_module == 6            

    input:
    path(tax_files)
    path(fun_files)
    path(meancov_table)
    path(readscov_table)

    output:
    path("${task.process.toLowerCase().replaceFirst('module', 'module-')}")

    script:
    """
    mg-clust-module-6.py \
        --tax_files      ${tax_files} \
        --fun_files      ${fun_files} \
        --meancov_table  ${meancov_table} \
        --readscov_table ${readscov_table} \
        --clust_thres    ${params.clust_thres} \
        --output_dir     ${task.process.toLowerCase().replaceFirst('module', 'module-')} \
        --overwrite
    """
}