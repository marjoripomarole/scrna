# scRNA-seq Analysis — PBMC 3k

Single-cell RNA sequencing analysis of the classic **PBMC 3k** dataset (10x Genomics, 2,700 peripheral blood mononuclear cells).

## Dataset

- **Source:** 10x Genomics public PBMC 3k dataset
- **Raw data:** 2,700 cells × 32,738 genes
- Downloaded via `scanpy.datasets.pbmc3k()`

## Pipeline

Run the full pipeline end-to-end with:

```bash
python run_pipeline.py
```

### 1. Download & Load (`download_load_data.py`)

Downloads the PBMC 3k dataset using Scanpy's built-in downloader.
Also includes `load_data.R` to load the same data via Seurat from the 10X matrix files in `filtered_gene_bc_matrices/hg19/`.

### 2. Quality Control (`quality_control.py`)

QC metrics computed per cell:
- `n_genes_by_counts` — number of genes detected
- `total_counts` — total UMI counts
- `pct_counts_mt` — percentage of mitochondrial reads (genes starting with `MT-`)

**Filters applied:**
| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| Min genes per cell | ≥ 200 | Remove empty droplets |
| Max genes per cell | < 2,500 | Remove likely doublets |
| Mitochondrial reads | < 5% | Remove dead/damaged cells |
| Min cells per gene | ≥ 3 | Remove very lowly detected genes |

**Results after QC:** 2,700 → **2,638 cells**, 32,738 → **13,714 genes**

### 3. Normalization & Feature Selection (`normalization_feature_selection.py`)

- Normalize each cell to 10,000 total counts, then log1p-transform
- Select top **1,838 highly variable genes** (mean: 0.0125–3, dispersion ≥ 0.5)
- Scale to zero mean, unit variance (capped at 10)

### 4. Dimensionality Reduction (`dimensionality_reduction.py`)

- PCA (50 components, arpack solver)
- UMAP for 2D visualization

### 5. Clustering (`clustering.py`)

- k-NN graph: 10 neighbors, top 40 PCs
- Leiden algorithm (resolution = 0.5) → **6 clusters**

### 6. Cell Type Annotation (`cell_type_annotation.py`)

Clusters annotated using canonical PBMC marker genes:

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0, 2 | CD4+ T cells | CD3D, CD4, IL32 |
| 1, 6 | CD14+ Monocytes | CD14, LYZ, CST3 |
| 3 | B cells | MS4A1, CD79A |
| 4 | CD8+ T cells | CD8A, GNLY |
| 5 | NK cells | GNLY, FCGR3A |

**Cell type counts:**
| Cell Type | Cells |
|-----------|-------|
| CD4+ T cells | 1,610 |
| CD14+ Monocytes | 639 |
| B cells | 340 |
| CD8+ T cells | 36 |
| NK cells | 13 |

### 7. Marker Genes (`find_marker_genes.py`)

Wilcoxon rank-sum test per cluster. Top 5 marker genes visualized as dot plot and heatmap.

### 8. Differential Expression (`differential_expression.py`)

Runs after the main pipeline (requires `pbmc3k_processed.h5ad`):

```bash
python differential_expression.py
```

- **One-vs-rest DE** per cell type (Wilcoxon, using raw counts)
- **Pairwise DE:** CD4+ T cells vs CD8+ T cells
- **Pairwise DE:** CD14+ Monocytes vs NK cells
- **Volcano plot** for CD4+ vs CD8+ (`figures/volcano_cd4_vs_cd8.png`)

Top hits (CD4+ vs CD8+): MALAT1, CD3D, IL32 enriched in CD4+; HLA-DPA1, HLA-DPB1 enriched in CD8+.

## Requirements

```
scanpy
anndata
pandas
numpy
matplotlib
leidenalg
```

For R:
```
Seurat
ggplot2
dplyr
```

## Project Structure

```
scrna/
├── run_pipeline.py                 # Full pipeline (steps 1–7)
├── differential_expression.py      # DE analysis (step 8)
├── download_load_data.py
├── quality_control.py
├── normalization_feature_selection.py
├── dimensionality_reduction.py
├── clustering.py
├── cell_type_annotation.py
├── find_marker_genes.py
├── load_data.R
├── filtered_gene_bc_matrices/      # 10X sparse matrix files
│   └── hg19/
└── figures/
    ├── violin_qc.png
    ├── filter_genes_dispersion_hvg.png
    ├── pca_variance_ratio_variance_ratio.png
    ├── pca_CST3.png
    ├── umap_leiden.png
    ├── rank_genes_groups_leiden_markers.png
    ├── dotplot__markers_dotplot.png
    └── volcano_cd4_vs_cd8.png
```
