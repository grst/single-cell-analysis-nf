#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
assert params.input: "Input samplesheet not specified!" 

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main.nf" addParams( 
    options: modules['SCQC']
)
include { SOLO_SPLIT_BATCHES } from "./modules/local/solo/main.nf" addParams( 
    options: modules['SOLO_SPLIT_BATCHES']
)
include { SOLO } from "./modules/local/solo/main.nf" addParams( 
    options: modules['SOLO']
)
include { SOLO_MERGE_BATCHES } from "./modules/local/solo/main.nf" addParams( 
    options: modules['SOLO_MERGE_BATCHES']
)

process MERGE_STATS {
    publishDir {params.outdir}, mode:params.publish_dir_mode
    input:
        path(stats_tsv)

    output:
        path("qc_stats_all.tsv")

    script:
    """
    mkdir out
    head -n1 ${stats_tsv[0]} > out/qc_stats_all.tsv
    for f in *.tsv; do
        tail -n+2 \$f >> out/qc_stats_all.tsv
    done
    mv out/qc_stats_all.tsv . 
    """
}


workflow { 
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples)

    if(!params.skip_solo) {
        SOLO_MERGE_BATCHES(
            SOLO (
                SOLO_SPLIT_BATCHES(SCQC.out.adata.filter{ it["run_solo"] == "true" }).adata.transpose()
            ).is_doublet.groupTuple()
        )
    }

    MERGE_STATS(SCQC.out.qc_stats.collect())

    // SCNORM(SCQC.out.adata)
}