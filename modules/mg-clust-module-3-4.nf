// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3-4: concatenate all samples, filter ORFs, create MMseqs2 DB,
//             cluster ORFs, and generate per-cluster abundance tables
// Input:  all per-sample ORF FASTAs + coverage tables (collected)
// Output: per-cluster abundance tables (workable)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE3_4 {
    container "ghcr.io/epereira/mg-clust/module-3-4:latest"
    publishDir "${params.output_dir}/",
                mode: "copy"

    input:
    path(orf_files)
    path(meancov_files)
    path(readscov_files)

    output:
    path("module3-4/orfs_clust_id*_meancov_workable.tsv"),   emit: meancov_workable
    path("module3-4/orfs_clust_id*_readscov_workable.tsv"),  emit: readscov_workable

    script:
    """
    mg-clust-module-3-4.py \
        --orf_files      ${orf_files} \
        --meancov_files  ${meancov_files} \
        --readscov_files ${readscov_files} \
        --output_dir     module3-4 \
        --nslots         ${params.nslots} \
        --min_orf_length ${params.min_orf_length} \
        --clust_thres    ${params.clust_thres} \
        --clust_cov_len  ${params.clust_cov_len} \
        --overwrite
    """
}
