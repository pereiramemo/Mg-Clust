// ─────────────────────────────────────────────────────────────────────────────
// MODULE 4: Taxonomic annotation of contigs against GTDB using MMseqs2
// Input:  contig sequences from module 1 (per sample)
// Output: per-contig taxonomy annotations
// ─────────────────────────────────────────────────────────────────────────────

process MODULE4 {
    
    container "ghcr.io/epereira/mg-clust/${task.process.toLowerCase().replaceFirst('module', 'module-')}:latest"
    publishDir "${params.output_dir}/${task.process.toLowerCase().replaceFirst('module', 'module-')}/",
           mode: "copy",
           enabled: params.full_output || params.stop_at_module == 4           


    tag "${sample_name}"

    input:
    tuple val(sample_name), path(assembly), path(orf_bed)

    output:
    path("${sample_name}"),                                                           emit: sample_dir
    path("${sample_name}/${sample_name}_orf_tax_annot_workable.tsv"),                 emit: tax_workable

    script:
    """
    mg-clust-module-4.py \
        --contigs       ${assembly} \
        --bed_file      ${orf_bed} \
        --gtdb          ${params.gtdb} \
        --sample_name   ${sample_name} \
        --lca_mode      ${params.lca_mode} \
        --sensitivity   ${params.sensitivity} \
        --tax_lineage   ${params.tax_lineage} \
        --nslots        ${params.nslots} \
        --output_dir    ${sample_name} \
        --overwrite
    """
}
