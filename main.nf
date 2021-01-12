#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
assert params.input: "Input samplesheet not specified!" 

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main.nf" addParams( options: modules['SCQC'])
// include { SCNORM } from "./modules/local/scnorm/main.nf" addParams( options: modules['SCNORM'])
include { SOLO_SPLIT_BATCHES } from "./modules/local/solo/main.nf" addParams( options: modules['SOLO_SPLIT_BATCHES'])
include { SOLO } from "./modules/local/solo/main.nf" addParams( options: modules['SOLO'])
include { SOLO_MERGE_BATCHES } from "./modules/local/solo/main.nf" addParams( options: modules['SOLO_MERGE_BATCHES'])


workflow { 
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples)

    SOLO_MERGE_BATCHES(
        SOLO (
            SOLO_SPLIT_BATCHES(SCQC.out.adata).adata.transpose()
        ).is_doublet.groupTuple()
    )

    // SCNORM(SCQC.out.adata)
}