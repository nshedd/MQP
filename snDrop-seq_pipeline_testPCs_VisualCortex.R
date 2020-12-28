library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


## 2 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:2)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("In1","Ex1","Ex2","End","Oli1","Ex3","Oli1","Ex4","Ex5","Per/Mic","OPC","Ex6","Ex7","In2","In3","Ast","Oli2")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_002pc.png"), device=)

## 5 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:5)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:5, metric="euclidean")

#new.cluster.ids <- c("Mic/Per/End","Oli","Ex1","Ex2","In1","Ex3","In2","Ast","OPC")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_005pc.png"), device=)

## 10 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:10)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:10, metric="euclidean")

#new.cluster.ids <- c("Ex1/End","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","Ex5","OPC","Per/Mic")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_010pc.png"), device=)

## 20 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:20)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ast","In2","In3","Ex3","Ex4","OPC","Ex5","Ex6","In4","Mic","Ex7","Ex8","Ex9","Ex10","Per")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_020pc.png"), device=)

## 50 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:50)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:50, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ex3","Ast","In2","Ex4","In3","OPC","Ex5","In4","Mic","Ex6","Ex7","In5","Ex8","Per")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_050pc.png"), device=)

## 75 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:75)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:75, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","Ex6","In4","Ex7","Per")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_075pc.png"), device=)

## 100 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:100)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:100, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","In4","Ex6","Per")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_100pc.png"), device=)

## 150 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:150)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:150, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","OPC","Ex5","Mic/Per","Ex6")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_150pc.png"), device=)

## 200 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:200)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:200, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","OPC","Ex4","Mic/Per","Ex5")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_200pc.png"), device=)
