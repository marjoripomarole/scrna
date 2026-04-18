# Run PCA
sc.tl.pca(adata, svd_solver='arpack')

# Plot variance explained by each PC
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
# Look for the "elbow" — where adding more PCs stops helping

# Visualize cells in PCA space
sc.pl.pca(adata, color='CST3')  # CST3 = monocyte marker
