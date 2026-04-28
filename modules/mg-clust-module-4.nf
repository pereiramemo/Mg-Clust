// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: cluster ORFs + generate abundance tables
// Input:  MMseqs2 DB + merged coverage tables
// Output: per-cluster abundance tables
// Note:  orf_filt_db is a literal prefix — Nextflow stages all DB companion files under this name
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4 {
    container "ghcr.io/epereira/mg-clust/module-4:latest"
    publishDir "${params.output_dir}/", 
                mode: "copy"

    input:
    path(orfs_filt_db)          
    path(concat_meancov_table)
    path(concat_readscov_table)

    output:
    path("module4"),                                            emit: module4_dir
    path("module4/orfs_clust_id*_meancov_workable.tsv"),       emit: meancov_workable
    path("module4/orfs_clust_id*_readscov_workable.tsv"),      emit: readscov_workable

    script:
    """
    mg-clust-module-4.py \
        --orfs_db        orfs_filt_db \
        --meancov_table  ${concat_meancov_table} \
        --readscov_table ${concat_readscov_table} \
        --output_dir     module4 \
        --nslots         ${params.nslots} \
        --clust_thres    ${params.clust_thres} \
        --clust_cov_len  ${params.clust_cov_len} \
        --overwrite
    """
}
