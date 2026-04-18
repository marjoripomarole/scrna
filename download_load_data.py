import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.figdir = './figures/'

# Download and load PBMC 3k dataset via scanpy's built-in downloader
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

print(adata)
# AnnData object with n_obs x n_vars = 2700 x 32738
