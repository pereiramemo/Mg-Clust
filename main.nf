#!/usr/bin/env nextflow

// Include modules
include { MODULE1 } from './modules/mg-clust-module-1.nf'
include { MODULE2 } from './modules/mg-clust-module-2.nf'
include { MODULE3 } from './modules/mg-clust-module-3.nf'
include { MODULE4 } from './modules/mg-clust-module-4.nf'


workflow {

    main:
    // Read paired-end samples from reads_dir using the pattern defined in nextflow.config
    reads_ch = channel.fromFilePairs(
        "${params.input_dir}/${params.reads_pattern}",
        checkIfExists: true
    )

    if (params.stop_at_module < 4) {
        log.warn "Workflow will stop at module ${params.stop_at_module}. "
    }

    // MODULE1: Assembly + read mapping
    module1_out = MODULE1(reads_ch)

    if (params.stop_at_module >= 2) {
        // MODULE2: ORF prediction + coverage estimation (per sample)
        module2_out = MODULE2(module1_out)

        if (params.stop_at_module >= 3) {
            // MODULE3: Concatenate samples, filter ORFs, create MMseqs2 DB
            orf_files_ch      = module2_out.faa.collect()
            meancov_files_ch  = module2_out.meancov.collect()
            readscov_files_ch = module2_out.readscov.collect()
            module3_out = MODULE3(orf_files_ch, meancov_files_ch, readscov_files_ch)

            if (params.stop_at_module >= 4) {
                // MODULE4: ORF clustering + abundance tables
                MODULE4(module3_out.orfs_filt_db, module3_out.concat_meancov, module3_out.concat_readscov)
            }
        }
    }

}
