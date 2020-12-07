/*
 * Reformat input samplesheet and check validity
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/bactanti/bin/
    //TODO meaningful check
    """
    cp $samplesheet samplesheet.valid.csv
    # check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ sample, single_end?, [ fastq_1, fastq_2 ] ]
def check_samplesheet_paths(LinkedHashMap row) {
    def meta = row.findAll { key, val -> key != "input_adata" }

    def array = [ meta, file(row.input_adata, checkIfExists: true)]
    return array
}
