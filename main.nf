#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
assert params.input: "Input samplesheet not specified!" 

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main.nf" addParams( options: modules['SCQC'])
include { SCNORM } from "./modules/local/scnorm/main.nf" addParams( options: modules['SCNORM'])

workflow { 
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples) 

    SCNORM(SCQC.out.adata)
}