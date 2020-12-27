library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:100)
#CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:100, metric="euclidean")

#new.cluster.ids <- c("Gran1", "Gran2", "Purk", "Ast", "Gran3", "OPC1", "Oli", "Mic", "End/Per", "OPC2")
#names(new.cluster.ids) <- levels(CerebellarHem)
#CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_100pc_oglabels.png"), device=)
