# Build neighborhood graph using top 10 PCs
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Cluster using Leiden algorithm (better than Louvain)
sc.tl.leiden(adata, resolution=0.5)
print(f"Number of clusters: {adata.obs['leiden'].nunique()}")

# Run UMAP for visualization
sc.tl.umap(adata)

# Plot UMAP colored by cluster
sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
           title='PBMC clusters', frameon=False)
