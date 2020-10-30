install.packages("dplyr")
install.packages("patchwork")

library(dplyr)
library(Seurat)
library(patchwork)


matrix = read.table("/data/zusers/pratth/sc/atac/GSM3722075_PBMC_Rep3_fragments.tsv.gz.rDHS.matrix.tsv", header=TRUE, sep="\t", row.names=1)
print("loaded matrix")

colors <- scan(path.expand("~/set2_top10k_colors.txt"), what=character, sep='\n')
print("loaded colors")

set2umap <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)
print("created seurat object")

set2umap <- NormalizeData(set2umap, normalization.method = "LogNormalize", scale.factor = 10000)
print("normalized data")

all.cells <- rownames(set2umap)
set2umap <- ScaleData(set2umap, features = all.cells)

set2umap <- RunUMAP(set2umap, dims = 1:10)

DimPlot(set2umap, reduction = "umap", do.return=TRUE)
ggsave(path.expand("~/umap_colored_set2_seurat_default.svg"), device=)



