library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:200)
#CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:200, metric="euclidean")

#new.cluster.ids <- c("Gran", "Purk", "Ast", "OPC1", "Oli", "Mic", "End/Per", "OPC2")
#names(new.cluster.ids) <- levels(CerebellarHem)
#CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_200pc_oglabels.png"), device=)
