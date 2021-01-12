#!/usr/bin/env python

import sys

import scanpy as sc
import pandas as pd
import numpy as np

adata_path = sys.argv[2]
is_doublet_path = sys.argv[1]

is_doublet = np.load(is_doublet_path)
adata = sc.read_h5ad(adata_path)

df = (
    pd.DataFrame()
    .assign(var_names=adata.obs_names, is_doublet=is_doublet)
    .set_index("var_names")
)
df.to_csv(sys.stdout, header=False)