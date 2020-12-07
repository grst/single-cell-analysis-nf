#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { CHECK_SAMPLESHEET; check_samplesheet_paths } from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main.nf" addParams( options: modules['SCQC'])

workflow { 
    ch_samples = CHECK_SAMPLESHEET(ch_input)
        .splitCsv(header:true, sep:',')
        .map { check_samplesheet_paths(it) }

    SCQC(ch_samples)

}