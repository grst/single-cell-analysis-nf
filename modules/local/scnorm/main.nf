include { nxfVars } from "../nxfvars.nf"
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process SCNORM {
    tag = { meta.id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda "/home/sturm/.conda/envs/single-cell-analysis-nf-norm"
    cpus 22

    input:
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta), path(output_adata), emit: adata
    path("*.html"), emit: notebook

    script:
    output_adata = "output_adata.h5ad"
    cell_cycle_genes = "${moduleDir}/Macosko_cell_cycle_genes.txt"
    """
    ${nxfVars(task)}
    jupytext --execute -o \$(pwd)/scnorm-notebook.ipynb ${moduleDir}/scnorm-notebook.py
    jupyter nbconvert scnorm-notebook.ipynb --to html
    """    
}