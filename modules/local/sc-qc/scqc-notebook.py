# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
from nxfvars import nxf
import matplotlib.pyplot as plt
import seaborn as sns
from qc_plots import plot_qc_metrics

# %%
input_adata = nxf.input("input_adata", "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/10_public_datasets/Laughney_Massague_2020_NSCLC/h5ad/laughney_massague_2020_nsclc.h5ad")
meta = nxf.input("meta", {"min_genes": 500, "min_counts": 1800, "max_pct_mito": 20, "max_counts": 50000 })

# %%
adata = sc.read_h5ad(input_adata)

# %%
if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

# %%
sc.pp.calculate_qc_metrics(adata, qc_vars=("mito", ), log1p = False, inplace=True, percent_top=None)

# %%
adata.obs.columns

# %%
sc.pl.violin(adata, 'total_counts', groupby="sample", rotation=90, log=True, cut=0)

# %%
sc.pl.violin(adata, 'pct_counts_mito', groupby="sample", rotation=90)

# %%
plot_qc_metrics(adata, **meta)

# %%
plot_qc_metrics(adata, cumulative=True, **meta)

# %%
