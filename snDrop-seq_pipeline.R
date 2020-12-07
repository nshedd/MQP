library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

cerebellarhem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

cerebellarhem <- NormalizeData(cerebellarhem, normalization.method = "LogNormalize", scale.factor = 10000)

cerebellarhem <- FindVariableFeatures(object = cerebellarhem)

all_cells <- rownames(cerebellarhem)
cerebellarhem <- ScaleData(cerebellarhem, features = all_cells)

cerebellarhem <- RunPCA(cerebellarhem, features = VariableFeatures(object = cerebellarhem))

cerebellarhem <- FindNeighbors(cerebellarhem, dims = 1:10)
cerebellarhem <- FindClusters(cerebellarhem, resolution = 0.5)

cerebellarhem <- RunUMAP(cerebellarhem, dims = 1:10)

saveRDS(set2umap, file = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(set2umap, reduction = "umap")
ggsave(path.expand("~/umap_GSE97930_CerebellarHem_Seurat_default.png"), device=)

set2umap.markers <- FindAllMarkers(set2umap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = set2umap.markers %>% group_by(cluster)

path2 = path.expand("~/Zlab single-cell marker genes - Brain.tsv")

brain_genes = read.table(path2, header=TRUE)
blood_genes = 

celltypes <- character()
for (gene in diff_expressed$gene) {
  if (gene %in% brain_genes$Human_Gene) {
    celltypes <- c(celltypes, T_cells$Human_Gene)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/GSE97930_CerebellarHem_differentiallyexpressed.txt"), sep="\t")