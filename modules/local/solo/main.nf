include { nxfVars } from "../nxfvars.nf"
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SOLO_SPLIT_BATCHES {
    tag = { meta.id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda "/home/sturm/.conda/envs/single-cell-analysis-nf-solo"

    input:
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta), path("*.h5ad"), emit: adata 

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc

    adata = sc.read_h5ad("${input_adata}") 
    batch_key = "${meta.batch_key}"

    # select highly variable genes to speedup solo
    # salmon counts are not integers. However this breaks highly variable genes. 
    adata.X.data = adata.X.data.astype("int")
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000, subset=True)

    # write out adata separated by batches
    if batch_key == "":
        adata.write_h5ad("${meta.id}.solo.h5ad", compression="lzf")
    else:
        for batch in adata.obs[batch_key].unique():
            tmp_adata = adata[adata.obs[batch_key] == batch, :].copy()
            tmp_adata.write_h5ad(f"${meta.id}_{batch}.solo.h5ad", compression="lzf")
    """
}

process SOLO {
    tag = { id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    errorStrategy 'ignore'
    cpus 4
    conda "/home/sturm/.conda/envs/single-cell-analysis-nf-solo"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

    input:
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta), path("${id}_is_doublet.csv"), emit: is_doublet

    script:
    id = input_adata.baseName
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}                         
    solo --set-reproducible-seed 42 -o out ${moduleDir}/model.json ${input_adata}
    python ${moduleDir}/make_solo_csv.py out/is_doublet.npy "${input_adata}" > ${id}_is_doublet.csv
    """
}

process SOLO_MERGE_BATCHES {
    tag = { meta.id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda "/home/sturm/.conda/envs/single-cell-analysis-nf-solo"

    input:
    tuple val(meta), path(is_doublet)

    output:
    tuple val(meta), path("${meta.id}.is_doublet.csv"), emit: is_doublet

    script:
    """
    cat *.csv > "${meta.id}.is_doublet.csv"
    """
}