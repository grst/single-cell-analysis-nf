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
import sys
sys.path.insert(0, "/home/sturm/projects/2020/scanpy/")
import scanpy as sc
import scanpy.external as sce
from sctransform import sctransform


# +
# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
# rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
# %load_ext rpy2.ipython

# + language="R"
#
# # Load all the R libraries we will be using in the notebook
# library(scran)
# -

adata = sc.read_h5ad("/tmp/adata.h5ad")

# +
# adata = sc.datasets.pbmc3k()
# -

adata_raw = adata.copy()
sc.pp.log1p(adata_raw)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# ### CPM normalization

adata_cpm = adata.copy()
adata_cpm.raw = adata_raw
sc.pp.normalize_per_cell(adata_cpm, counts_per_cell_after=1e6)
sc.pp.log1p(adata_cpm)
sc.pp.highly_variable_genes(adata_cpm, flavor="cell_ranger", n_top_genes=4000)
sc.pp.pca(adata_cpm)
sc.pp.neighbors(adata_cpm)
sc.tl.leiden(adata_cpm)
sc.tl.umap(adata_cpm)

# ### Scran normalization

adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.5)

adata_scran = adata.copy()

#Preprocess variables for scran normalization
input_groups = adata_pp.obs['groups']
data_mat = adata_scran.X.T.toarray()

pandas2ri.activate()
anndata2ri.activate()

# + magic_args="-i data_mat -i input_groups -o size_factors" language="R"
#
# size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)
# -
pandas2ri.deactivate()
anndata2ri.deactivate()

sns.distplot(size_factors)


adata_scran.obs['size_factors'] = size_factors

sc.pl.scatter(adata_scran, 'size_factors', 'n_genes')

adata_scran.layers["counts"] = adata_scran.X.copy()

sparse_size_factors = sp.sparse.diags(1/adata_scran.obs['size_factors'])

adata_scran.X = sparse_size_factors @ adata_scran.X

sc.pp.log1p(adata_scran)

adata_scran.raw = adata_scran

sc.pp.highly_variable_genes(adata_scran, flavor="cell_ranger", n_top_genes=4000)

sc.pp.pca(adata_scran, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_scran)
sc.tl.umap(adata_scran)

sc.tl.leiden(adata_scran)

# ### sctransform normalization

adata_sct = adata.copy()
adata_sct.raw = adata_raw

# +
# adata_sct.X = adata_sct.X.toarray()
# -

sctransform(adata_sct, n_top_genes=4000)

# +
# sc.pp.highly_variable_genes(adata_sct, flavor="cell_ranger", n_top_genes=4000)
# -

sc.pp.log1p(adata_sct)

sc.pp.pca(adata_sct, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_sct)
sc.tl.umap(adata_sct)

# ### scran + sctransform

adata_scran_sct = adata.copy()

adata_scran_sct.X = sparse_size_factors @ adata_scran_sct.X

sctransform(adata_scran_sct, n_top_genes=4000)

sc.pp.log1p(adata_scran_sct)

sc.pp.pca(adata_scran_sct, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_scran_sct)
sc.tl.umap(adata_scran_sct)

# ### Visualize

markers = ["CD8A", "CD4", "FOXP3", "sample"]

sc.pl.umap(adata_cpm, color=markers)

sc.pl.umap(adata_scran, color=markers)

sc.pl.umap(adata_sct, color=markers)

sc.pl.umap(adata_scran_sct, color=markers)


