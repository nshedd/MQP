library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

print("loaded packages")

matrix = readRDS(file=path.expand("~/GSM3722075_PBMC_Rep3_fragments.rds"))
print("loaded matrix")
print(matrix[1:10, 1:10, drop=FALSE])

colors <- scan(path.expand("~/set2_top10k_colors.txt"), what=character(), sep='\n')
print("loaded colors")

set2umap <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

set2umap <- NormalizeData(set2umap, normalization.method = "LogNormalize", scale.factor = 10000)

all_cells <- rownames(set2umap)
set2umap <- ScaleData(set2umap, features = all_cells)

set2umap <- FindVariableFeatures(object = set2umap)
set2umap <- RunPCA(set2umap, features = VariableFeatures(object = set2umap))

set2umap <- FindNeighbors(set2umap, dims = 1:10)
set2umap <- FindClusters(set2umap, resolution = 0.5)

set2umap <- RunUMAP(set2umap, dims = 1:10)

plot = DimPlot(set2umap, reduction = "umap")
ggsave(path.expand("~/umap_colored_set2_seurat_default.png"), device=)

set2umap.markers <- FindAllMarkers(set2umap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
set2umap.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
