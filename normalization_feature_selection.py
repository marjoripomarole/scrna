# Normalize: scale each cell to 10,000 total counts, then log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save raw normalized counts for later (needed for DE analysis)
adata.raw = adata

# Find highly variable genes (top 2000)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                             min_disp=0.5)
print(f"Highly variable genes: {adata.var.highly_variable.sum()}")

# Plot variable genes
sc.pl.highly_variable_genes(adata)

# Keep only highly variable genes for downstream analysis
adata = adata[:, adata.var.highly_variable]

# Scale data (zero mean, unit variance) — needed for PCA
sc.pp.scale(adata, max_value=10)
