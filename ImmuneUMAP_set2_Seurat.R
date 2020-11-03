library(dplyr)
library(Seurat)
library(patchwork)

print("loaded packages")

matrix = readRDS(file=path.expand("~/GSM3722075_PBMC_Rep3_fragments.rds"))
print("loaded matrix")

colors <- scan(path.expand("~/set2_top10k_colors.txt"), what=character(), sep='\n')
print("loaded colors")

set2umap <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 2000, min.features = 3)
print("created seurat object")

set2umap <- NormalizeData(set2umap, normalization.method = "LogNormalize", scale.factor = 10000)
print("normalized data")

all_cells <- rownames(set2umap)
set2umap <- ScaleData(set2umap, features = all_cells)

set2umap <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

set2umap <- RunUMAP(set2umap, dims = 1:10)

DimPlot(set2umap, reduction = "umap", do.return=TRUE)
ggsave(path.expand("~/umap_colored_set2_seurat_default.svg"), device=)



