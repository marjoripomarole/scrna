import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.figdir = './figures/'

adata = sc.read_h5ad('pbmc3k_processed.h5ad')
print(adata)
print("\nCell types:", adata.obs['cell_type'].unique().tolist())

# --- DE: marker genes per cell type (Wilcoxon) ---
sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon',
                         use_raw=True)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, use_raw=True)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, use_raw=True,
                                 show_gene_labels=True)

# --- Top markers per cell type as a dataframe ---
markers_df = sc.get.rank_genes_groups_df(adata, group=None)
print("\nTop DE genes per cell type:")
print(markers_df.groupby('group').head(5).to_string())

# --- Pairwise DE: CD4+ T cells vs CD8+ T cells ---
print("\n--- CD4+ vs CD8+ T cells ---")
sc.tl.rank_genes_groups(adata, groupby='cell_type', groups=['CD4+ T cells'],
                         reference='CD8+ T cells', method='wilcoxon',
                         use_raw=True, key_added='de_cd4_vs_cd8')
sc.pl.rank_genes_groups(adata, key='de_cd4_vs_cd8', n_genes=20)
de_cd4_cd8 = sc.get.rank_genes_groups_df(adata, group='CD4+ T cells',
                                          key='de_cd4_vs_cd8')
print(de_cd4_cd8.head(20).to_string())

# --- Pairwise DE: CD14+ Monocytes vs NK cells ---
print("\n--- CD14+ Monocytes vs NK cells ---")
sc.tl.rank_genes_groups(adata, groupby='cell_type',
                         groups=['CD14+ Monocytes'],
                         reference='NK cells', method='wilcoxon',
                         use_raw=True, key_added='de_mono_vs_nk')
sc.pl.rank_genes_groups(adata, key='de_mono_vs_nk', n_genes=20)

# --- Volcano plot: CD4+ vs CD8+ ---
de = de_cd4_cd8.copy()
de = de[de['pvals_adj'] > 0]
de['neg_log10_padj'] = -np.log10(de['pvals_adj'])
top_up = de[de['logfoldchanges'] > 0].head(10)
top_dn = de[de['logfoldchanges'] < 0].head(10)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(de['logfoldchanges'], de['neg_log10_padj'],
           alpha=0.4, s=10, color='grey')
ax.scatter(top_up['logfoldchanges'], top_up['neg_log10_padj'],
           color='red', s=20, label='Up in CD4+')
ax.scatter(top_dn['logfoldchanges'], top_dn['neg_log10_padj'],
           color='blue', s=20, label='Up in CD8+')
for _, row in pd.concat([top_up, top_dn]).iterrows():
    ax.annotate(row['names'], (row['logfoldchanges'], row['neg_log10_padj']),
                fontsize=7)
ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
ax.axvline(0, color='black', linestyle='--', linewidth=0.8)
ax.set_xlabel('log2 Fold Change (CD4+ vs CD8+)')
ax.set_ylabel('-log10(adjusted p-value)')
ax.set_title('Volcano: CD4+ T cells vs CD8+ T cells')
ax.legend()
plt.tight_layout()
plt.savefig('./figures/volcano_cd4_vs_cd8.png', dpi=150)
plt.show()

print("\nDone. Figures saved to ./figures/")
