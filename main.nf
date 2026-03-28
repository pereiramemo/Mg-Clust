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
    
    module1_out = MODULE1(reads_ch)
    module2_out = MODULE2(module1_out)
    // module2_out.collect().view()

    // Collect per-sample outputs into flat lists before running module 3
    orf_files_ch      = module2_out.faa.collect()
    meancov_files_ch  = module2_out.meancov.collect()
    readscov_files_ch = module2_out.readscov.collect()

    module3_out = MODULE3(orf_files_ch, meancov_files_ch, readscov_files_ch)
    orfs_filt_db_ch = module3_out.orfs_filt_db
    concat_meancov_ch = module3_out.concat_meancov
    concat_readscov_ch = module3_out.concat_readscov
    
    MODULE4(orfs_filt_db_ch, concat_meancov_ch, concat_readscov_ch)
}