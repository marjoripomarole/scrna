library(Seurat)
library(ggplot2)
library(dplyr)

# Load data (already downloaded via Python)
pbmc.data <- Read10X(data.dir = 'filtered_gene_bc_matrices/hg19/')

# Create Seurat object (keep genes in ≥3 cells, cells with ≥200 genes)
pbmc <- CreateSeuratObject(counts = pbmc.data,
                           project = 'PBMC3k',
                           min.cells = 3,
                           min.features = 200)
pbmc
# 2700 samples x 13714 features
