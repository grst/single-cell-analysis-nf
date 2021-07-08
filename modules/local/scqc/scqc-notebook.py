# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
from nxfvars import nxf
import matplotlib.pyplot as plt
import seaborn as sns
from qc_plots import plot_qc_metrics
import pandas as pd
import numpy as np

# %%
dataset_id = nxf.input("dataset_id", "lambrechts_2018_luad_6653")
input_adata = nxf.input(
    "input_adata",
    "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/10_public_datasets/Lambrechts_2018_LUAD/E-MTAB-6653/h5ad_raw/lambrechts_2018_luad_6653.h5ad",
)
output_adata = nxf.input("output_adata", "/tmp/adata.h5ad")
output_stats = nxf.input("output_stats", "/tmp/qc_stats.tsv")
thresholds = {
    key: int(nxf.input(key, default_value))
    for key, default_value in {
        "min_genes": 500,
        "min_counts": 1800,
        "max_pct_mito": 20,
        "max_counts": 50000,
    }.items()
}

# %%
adata = sc.read_h5ad(input_adata)

# %%
# Add fake sample if its not in obs
if "sample" not in adata.obs.columns:
    no_sample = True
    adata.obs["sample"] = "1"
else:
    no_sample = False

# %%
if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

# %%
sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# %%
def get_stats_df(adata):
    return pd.DataFrame().assign(
        dataset_id=[dataset_id],
        min_genes=[np.min(adata.obs["n_genes_by_counts"])],
        max_genes=[np.max(adata.obs["n_genes_by_counts"])],
        min_counts=[np.min(adata.obs["total_counts"])],
        max_counts=[np.max(adata.obs["total_counts"])],
        min_pct_mito=[np.min(adata.obs["pct_counts_mito"])],
        max_pct_mito=[np.max(adata.obs["pct_counts_mito"])],
        n_obs=len(adata.obs_names),
        n_var=len(adata.var_names),
    )


stats_before = get_stats_df(adata).assign(status=["before_qc"])


# %%
adata.obs.columns

# %%
figwidth = min(max(adata.obs["sample"].unique().size * 0.5, 2), 20)

# %%
fig, ax = plt.subplots(1, 1, figsize=(figwidth, 5))
sc.pl.violin(
    adata, "total_counts", groupby="sample", rotation=90, log=True, cut=0, ax=ax
)

# %%
fig, ax = plt.subplots(1, 1, figsize=(figwidth, 5))
sc.pl.violin(adata, "pct_counts_mito", groupby="sample", rotation=90, ax=ax)

# %%
plot_qc_metrics(adata, **thresholds)

# %%
plot_qc_metrics(adata, cumulative=True, **thresholds)

# %%
# very basic gene filtering - genes with 0 cells cause some downstream processes to fail.
print("Filtering genes")
print(f"    Before: {adata.shape[1]}")
sc.pp.filter_genes(adata, min_counts=3)
print(f"    After: {adata.shape[1]}")

# %%
# Apply thresholds
print("Filter by min_counts")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, min_counts=thresholds["min_counts"])
print(f"    After: {adata.shape[0]}")


print("Filter by max_counts")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, max_counts=thresholds["max_counts"])
print(f"    After: {adata.shape[0]}")


print("Filter by min_genes")
print(f"    Before: {adata.shape[0]}")
sc.pp.filter_cells(adata, min_genes=thresholds["min_genes"])
print(f"    After: {adata.shape[0]}")


print("Filter by max_pct_mito")
print(f"    Before: {adata.shape[0]}")
adata = adata[adata.obs["pct_counts_mito"] < thresholds["max_pct_mito"]].copy()
print(f"    After: {adata.shape[0]}")

# %% [markdown]
# ## After filtering

# %%
plot_qc_metrics(adata, **thresholds)

# %%
plot_qc_metrics(adata, cumulative=True, **thresholds)

# %% [markdown]
# ### Save AnnData object

# %%
if no_sample:
    del adata.obs["sample"]

# %%
stats_after = get_stats_df(adata).assign(status=["after_qc"])
pd.concat([stats_before, stats_after]).to_csv(output_stats, sep="\t")
adata.write_h5ad(output_adata)

# %%
