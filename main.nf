#!/usr/bin/env nextflow

// Include modules
include { MODULE1 } from './modules/mg-clust-module-1.nf'
include { MODULE2 } from './modules/mg-clust-module-2.nf'
include { MODULE3 } from './modules/mg-clust-module-3.nf'
include { MODULE4 } from './modules/mg-clust-module-4.nf'
include { MODULE5 } from './modules/mg-clust-module-5.nf'
include { MODULE6 } from './modules/mg-clust-module-6.nf'

workflow {

    main:
    // Read paired-end samples from reads_dir using the pattern defined in nextflow.config
    reads_ch = channel.fromFilePairs(
        "${params.input_dir}/${params.reads_pattern}",
        checkIfExists: true
    )

    log.info "Pipeline will run up to module ${params.stop_at_module}"

    // MODULE1: Assembly + read mapping (always runs)
    module1_out = MODULE1(reads_ch)

    // MODULE2: ORF prediction + coverage estimation (per sample)
    if (params.stop_at_module >= 2) {
        module2_out = MODULE2(module1_out)
    }

    // MODULE3: Concatenate samples, filter ORFs, create MMseqs2 DB, cluster ORFs, generate abundance tables
    if (params.stop_at_module >= 3) {
        orf_files_ch      = module2_out.faa.map     { _sn, faa -> faa }.collect()
        meancov_files_ch  = module2_out.meancov.map  { _sn, cov -> cov }.collect()
        readscov_files_ch = module2_out.readscov.map { _sn, cov -> cov }.collect()
        module3_out       = MODULE3(orf_files_ch, meancov_files_ch, readscov_files_ch)
    }

    // MODULE4: Taxonomic annotation of contigs against GTDB using MMseqs2
    if (params.stop_at_module >= 4) {
        contigs_ch = module1_out
                    .map { sample_name, assembly, _bam -> tuple(sample_name, assembly) }
                    .join(module2_out.bed, by : 0) // joins on sample_name, produces tuple(sample_name, assembly, orf_bed)
            
        module4_out = MODULE4(contigs_ch)
    }

    // MODULE5: Functional annotation of ORFs against KO HMMs using pyHMMER
    if (params.stop_at_module >= 5) {
        faa = module2_out.faa.map { sample_name, faa_path -> tuple(sample_name, faa_path) }
        module5_out = MODULE5(faa)
    }

    // MODULE6: Merge taxonomy, function, and abundance tables
    if (params.stop_at_module >= 6) {
        tax_files_ch = module4_out.tax_workable.collect()
        fun_files_ch = module5_out.fun_workable.collect()
        MODULE6(
            tax_files_ch,
            fun_files_ch,
            module3_out.meancov_workable,
            module3_out.readscov_workable
        )
    }

}
