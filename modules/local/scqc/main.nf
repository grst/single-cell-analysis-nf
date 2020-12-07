include { nxfVars } from "../nxfvars.nf"
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process SCQC {
    tag = { meta.id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta), path(output_adata), emit: adata
    path("*.html"), emit: notebook

    script:
    output_adata = "output_adata.h5ad"
    min_genes = meta.min_genes
    min_counts = meta.min_counts
    max_pct_mito = meta.max_pct_mito
    max_counts = meta.max_counts
    """
    ${nxfVars(task)}
    export PYTHONPATH="${moduleDir}"
    jupytext --execute -o scqc-notebook.ipynb ${moduleDir}/scqc-notebook.py
    jupyter nbconvert scqc-notebook.ipynb --to html
    """    
}