import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.figdir = './figures/'

# --- Load data ---
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()
print(adata)

# --- Quality control ---
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                            log1p=False, inplace=True)
sc.pl.violin(adata,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
print(f"Cells after QC: {adata.n_obs}")
print(f"Genes after QC: {adata.n_vars}")

# --- Normalization & feature selection ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"Highly variable genes: {adata.var.highly_variable.sum()}")
sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

# --- Dimensionality reduction ---
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pl.pca(adata, color='CST3')

# --- Clustering ---
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata, resolution=0.5)
print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
           title='PBMC clusters', frameon=False)

# --- Marker genes ---
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head(20))
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

# --- Cell type annotation ---
marker_genes = ['CD3D', 'CD4', 'CD8A', 'MS4A1', 'GNLY',
                'CD14', 'LYZ', 'FCGR3A', 'FCER1A', 'PPBP']
sc.pl.umap(adata, color=marker_genes, ncols=5)

cell_type_map = {
    '0': 'CD4+ T cells',
    '1': 'CD14+ Monocytes',
    '2': 'CD4+ T cells',
    '3': 'B cells',
    '4': 'CD8+ T cells',
    '5': 'NK cells',
    '6': 'CD14+ Monocytes',
    '7': 'Dendritic cells',
    '8': 'CD16+ Monocytes',
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)

sc.pl.umap(adata, color='cell_type', legend_loc='on data',
           title='PBMC Cell Types', frameon=False,
           palette='tab10')

print("\nCell type counts:")
print(adata.obs['cell_type'].value_counts())

adata.write('pbmc3k_processed.h5ad')
