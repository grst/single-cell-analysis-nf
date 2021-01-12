# +
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

import rpy2.rinterface_lib.callbacks
import logging

from rpy2.robjects import pandas2ri
import anndata2ri
import scanpy as sc
from nxfvars import nxf
from threadpoolctl import threadpool_limits
import matplotlib
from pathlib import Path
# -


cpus = int(nxf.task("cpus", 16))
cell_cycle_markers = nxf.input("cell_cycle_genes", "Macosko_cell_cycle_genes.txt")
input_adata = nxf.input("input_adata", "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_qc_norm_scrnaseq/01_qc_and_filtering/Adams_Kaminski_2020_COPD/output_adata.h5ad")
output_adata = nxf.input("output_adata", "/tmp/adata_norm.h5ad")

threadpool_limits(cpus)

# +
# Ignore R warning messages
# Note: this can be commented out to get more verbose R output
# rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

pandas2ri.activate()
anndata2ri.activate()

# Automatically convert rpy2 outputs to pandas dataframes
# %load_ext rpy2.ipython

# + magic_args="-i cpus" language="R"
#
# # Load all the R libraries we will be using in the notebook
# library(scran)
# library(BiocParallel)
# BPPARAM = MulticoreParam(workers=cpus)
# -

adata = sc.read_h5ad(input_adata)

# ### Scran normalization
#
# Scran requires a clustering before normalization

adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups", resolution=0.5)

adata_scran = adata.copy()

# Preprocess variables for scran normalization
input_groups = adata_pp.obs["groups"]
data_mat = adata_scran.X.T.tocsr()

# + magic_args="-i data_mat -i input_groups -o size_factors" language="R"
#
# size_factors = computeSumFactors(
#     data_mat, clusters=input_groups, min.mean=0.1, BPPARAM=BPPARAM
# )
# -
sns.distplot(size_factors)


adata_scran.obs["size_factors"] = size_factors

sc.pl.scatter(adata_scran, "size_factors", "n_genes_by_counts")

adata_scran.layers["counts"] = adata_scran.X.copy()

sparse_size_factors = sp.sparse.diags(1 / adata_scran.obs["size_factors"])

adata_scran.X = sparse_size_factors @ adata_scran.X

sc.pp.log1p(adata_scran)

adata_scran.raw = adata_scran

adata = adata_scran

# ### Cell-cycle scoring

# Prepare genes for human and mouse
cc_genes = pd.read_table(cell_cycle_markers, delimiter="\t")
s_genes = cc_genes["S"].dropna()
g2m_genes = cc_genes["G2.M"].dropna()

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# ### Visualize

batch_key = None if "sample" not in adata.obs.columns else "sample"

try:
    sc.pp.highly_variable_genes(
        adata_scran, flavor="cell_ranger", n_top_genes=5000, batch_key=batch_key
    )
except (ValueError, IndexError):
    # This tends to fail with a batch key... if it does re-do it without a batch key. 
    sc.pp.highly_variable_genes(
        adata_scran, flavor="cell_ranger", n_top_genes=5000, batch_key=None
    )

sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.tl.leiden(adata)

markers = {
    "CD8A",
    "CD4",
    "FOXP3",
    "phase",
    "leiden",
    "sample",
    "origin",
    "tissue",
    "condition",
} & (set(adata.var_names) | set(adata.obs.columns))

sc.set_figure_params(figsize=(5,5))

sc.pl.umap(adata, color=markers, ncols=3)

adata.write_h5ad(output_adata, compression="lzf")
