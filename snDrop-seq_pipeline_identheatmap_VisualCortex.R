devtools::install_github("nshedd/scclusteval")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scclusteval)

path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

## Paper labels

VisualCortex2 <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))	

VisualCortex2 <- RunUMAP(VisualCortex2, dims = 1:20, metric="euclidean")

plot = DimPlot(VisualCortex2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_findct_oglabels.png"), device=)

## My labels
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))	

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:20)
VisualCortex <- FindClusters(VisualCortex, resolution = 1)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Oli","Ex1","Ex2","In1","Ex3","Ex4","Ast1","Ex5","Ex6","Ex7",
                     "In2","In3","OPC","In4","Ex8","In5","Mic","Ex9","End/Per","Ex10","Ast")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_findct.png"), device=)


## Heatmap comparing labels

oglabels <- Idents(VisualCortex2)
oglabels <- factor(oglabels)

newlabels <- Idents(VisualCortex)
newlabels <- factor(newlabels)


png(file=path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_JaccardHeatmap.png"))
PairWiseJaccardSetsHeatmap(ident1 = oglabels, ident2 = newlabels)
dev.off()
