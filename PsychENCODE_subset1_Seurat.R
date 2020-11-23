library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

print("loaded packages")

matrix = readRDS(file=path.expand("~/GRCh38-rDHSs_GW17_Cortex_50k_fragments_subset1.rds"))
print("loaded matrix")
matrix = t(matrix)
print("transposed matrix")
print(matrix[1:10, 1:10, drop=FALSE])
print("Number of rows:")
print(nrow(matrix))
print("Number of columns:")
print(ncol(matrix))
print(matrix)

set2umap <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

set2umap <- NormalizeData(set2umap, normalization.method = "LogNormalize", scale.factor = 10000)

all_cells <- rownames(set2umap)
set2umap <- ScaleData(set2umap, features = all_cells)

set2umap <- FindVariableFeatures(object = set2umap)
set2umap <- RunPCA(set2umap, features = VariableFeatures(object = set2umap))

set2umap <- FindNeighbors(set2umap, dims = 1:10)
set2umap <- FindClusters(set2umap, resolution = 0.5)

set2umap <- RunUMAP(set2umap, dims = 1:10)

saveRDS(set2umap, file = path.expand("~/GRCh38-rDHSs_GW17_Cortex_50k_fragments_subset1_Seurat.rds"))

plot = DimPlot(set2umap, reduction = "umap")
ggtitle("GRCh38-rDHSs_GW17_Cortex_50k_fragments_subset1")
ggsave(path.expand("~/PsychENCODE_subset1_seurat_default.png"), device=)

set2umap.markers <- FindAllMarkers(set2umap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = set2umap.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

write.table(diff_expressed, file = path.expand("~/set2_differentiallyexpressed.txt"), sep="\t")
