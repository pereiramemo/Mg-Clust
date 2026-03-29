// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: concatenate all samples, filter ORFs, create MMseqs2 DB
// Input:  all per-sample ORF FASTAs + coverage tables (collected)
// Output: filtered ORF MMseqs2 DB + merged coverage tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE3 {
    container "ghcr.io/epereira/mg-clust/module-3:latest"
    publishDir "${params.output_dir}/module3/", mode: "copy"

    input:
    path(orf_files)
    path(meancov_files)
    path(readscov_files)

    output:
    path("module3/orfs_filt_db*"),       emit: orfs_filt_db
    path("module3/orfs_meancov.tsv"),    emit: concat_meancov
    path("module3/orfs_readscov.tsv"),   emit: concat_readscov

    script:
    """
    mg-clust-module-3.py \
        --orf_files      ${orf_files} \
        --meancov_files  ${meancov_files} \
        --readscov_files ${readscov_files} \
        --output_dir     module3 \
        --nslots         ${params.nslots} \
        --min_orf_length ${params.min_orf_length} \
        --overwrite
    """
}

