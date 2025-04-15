import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad
import os

sc.settings.verbosity = 3

adata = sc.read_h5ad('/home/enno/code/scAnalysis/data/merged/combined_filtered.h5ad')

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 25, :]

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, inplace=True, flavor="seurat_v3", layer="counts")

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts'])
