# Calculate QC metrics
# Mitochondrial genes start with "MT-"
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                            log1p=False, inplace=True)

# Visualize QC metrics
sc.pl.violin(adata,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Filter cells:
# - At least 200 genes detected
# - No more than 2500 genes (doublet filter)
# - Less than 5% mitochondrial reads (dead cell filter)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print(f"Cells after QC: {adata.n_obs}")
print(f"Genes after QC: {adata.n_vars}")
