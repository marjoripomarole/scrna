# scRNA-seq Analysis — PBMC 3k

Single-cell RNA sequencing analysis of the classic **PBMC 3k** dataset (10x Genomics, 2,700 peripheral blood mononuclear cells).

## Dataset

- **Source:** 10x Genomics public PBMC 3k dataset
- **Raw data:** 2,700 cells × 32,738 genes
- Downloaded via `scanpy.datasets.pbmc3k()`

## Pipeline

### 1. Download & Load (`download_load_data.py`)

Downloads the PBMC 3k dataset using Scanpy's built-in downloader and saves the raw AnnData object to `data/pbmc3k_raw.h5ad`.

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

**Results after QC:**
- Cells: 2,700 → **2,638**
- Genes: 32,738 → **13,714**

QC violin plots are saved to `figures/violin_qc.png`.

## Requirements

```
scanpy
anndata
pandas
matplotlib
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
├── download_load_data.py       # Download PBMC3k and save raw data
├── quality_control.py          # QC metrics, filtering, violin plots
├── load_data.R                 # Load 10X data into Seurat
├── data/
│   └── pbmc3k_raw.h5ad         # Raw AnnData object
├── filtered_gene_bc_matrices/  # 10X sparse matrix files
│   └── hg19/
└── figures/
    └── violin_qc.png           # QC violin plots
```
