nextflow.enable.dsl = 2

def modules = params.modules.clone()
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

import { SCQC } from "./modules/local/scqc/main.nf"

workflow { 
    ch_samples = CHECK_SAMPLESHEET(ch_input)
        .splitCsv(header:true, sep:',')
        .map { check_samplesheet_paths(it) }


}