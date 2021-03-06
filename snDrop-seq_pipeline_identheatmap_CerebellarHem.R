devtools::install_github("nshedd/scclusteval")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scclusteval)

path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

## Paper labels

CerebellarHem2 <- readRDS(file=path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem2 <- RunUMAP(CerebellarHem2, dims = 1:20, metric="euclidean")

plot = DimPlot(CerebellarHem2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_findct_oglabels.png"), device=)

## My labels
CerebellarHem <- readRDS(file=path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:20)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Gran1","Gran2","Gran3","Purk1","Ast1","Purk2","OPC1","Oli","Mic","End/Per","OPC2","Ast2")
names(new.cluster.ids) <- levels(CerebellarHem)
CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_findct.png"), device=)

## Heatmap comparing labels

oglabels <- Idents(CerebellarHem2)
oglabels <- factor(oglabels)

newlabels <- Idents(CerebellarHem)
newlabels <- factor(newlabels)


png(file=path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_JaccardHeatmap.png"))
PairWiseJaccardSetsHeatmap(ident1 = oglabels, ident2 = newlabels)
dev.off()
