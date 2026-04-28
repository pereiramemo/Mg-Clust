// ─────────────────────────────────────────────────────────────────────────────
// MODULE 6: Functional annotation of ORFs against KO HMMs using pyHMMER
// Input:  ORF protein sequences from module 2 (per sample)
// Output: per-ORF functional annotations
// ─────────────────────────────────────────────────────────────────────────────

process MODULE6 {
    container "ghcr.io/epereira/mg-clust/module-6:latest"
    publishDir "${params.output_dir}/module6/",
                mode: "copy",
                enabled: params.full_output || params.stop_at_module == 6

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(orfs_faa)

    output:
    path("${sample_name}"),                                                           emit: sample_dir
    path("${sample_name}/${sample_name}_orf_fun_annot_workable.tsv"),                 emit: fun_workable

    script:
    """
    mg-clust-module-6.py \
        --orfs_faa      ${orfs_faa} \
        --hmm_db        ${params.hmm_db} \
        --sample_name   ${sample_name} \
        --evalue_thres  ${params.evalue_thres} \
        ${params.cut_ga ? '--cut_ga' : '--no-cut_ga'} \
        --nslots        ${params.nslots} \
        --output_dir    ${sample_name} \
        --overwrite
    """
}
