# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Get top markers as a dataframe
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head(20))

# Dot plot of top markers per cluster
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
