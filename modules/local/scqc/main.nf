import { nxfVars } from "../nxfvars.nf"


process SCQC {
    input:
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta), path(output_adata), emit: adata
    path("*.html"), emit: notebook

    script:
    """
    ${nxfVars(task)}

    jupyter nbconvert ${moduleDir}/scqc-notebook.py --execute --to html
    """    
}