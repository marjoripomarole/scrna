# Visualize known markers on UMAP
marker_genes = ['CD3D', 'CD4', 'CD8A', 'MS4A1', 'GNLY',
                'CD14', 'LYZ', 'FCGR3A', 'FCER1A', 'PPBP']
sc.pl.umap(adata, color=marker_genes, ncols=5)

# Assign cell type labels to clusters
# (adjust based on YOUR cluster numbers)
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

# Final annotated UMAP
sc.pl.umap(adata, color='cell_type', legend_loc='on data',
           title='PBMC Cell Types', frameon=False,
           palette='tab10')
