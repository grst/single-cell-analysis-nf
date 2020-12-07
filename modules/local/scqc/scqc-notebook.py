# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
from nxfvars import nxf
import matplotlib.pyplot as plt
import seaborn as sns
from qc_plots import plot_qc_metrics

# %%
input_adata = nxf.input(
    "input_adata",
    "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/10_public_datasets/Laughney_Massague_2020_NSCLC/h5ad/laughney_massague_2020_nsclc.h5ad",
)
output_adata = nxf.input("output_adata", "/tmp/adata.h5ad")
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
if "mito" not in adata.var.columns:
    adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")

# %%
sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# %%
adata.obs.columns

# %%
sc.pl.violin(adata, "total_counts", groupby="sample", rotation=90, log=True, cut=0)

# %%
sc.pl.violin(adata, "pct_counts_mito", groupby="sample", rotation=90)

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

# %%
adata.write_h5ad(output_adata)
