remotes::install_github("satijalab/seurat", ref = "release/4.0.0")

library(dplyr)
library(Seurat)
library(patchwork)


matrix = read.table("/data/zusers/pratth/sc/atac/GSM3722075_PBMC_Rep3_fragments.tsv.gz.rDHS.matrix.tsv", header=TRUE, sep="/t", row.names=1)

#sums= rowSums(matrix)

#n_matrix <-apply(matrix,2,function(x) matrix[,x]/rowsum(matrix))

colors <- scan(path.expand("~/set2_top10k_colors.txt"), character(), quote = "")

set2umap <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

matrix <- NormalizeData(matrix, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(matrix)
matrix <- ScaleData(matrix, features = all.genes)

matrix <- RunUMAP(matrix, dims = 1:10)

DimPlot(matrix, reduction = "umap", do.return=TRUE)
ggsave(path.expand("~/umap_colored_set2_seurat_default.svg"), device=)



