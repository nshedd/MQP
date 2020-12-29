library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


## 2 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:2)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("In1","Ex1","Ex2","End","Oli1","Ex3","Oli1","Ex4","Ex5","Per/Mic","OPC","Ex6","Ex7","In2","In3","Ast","Oli2")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_002pc_oglabels.png"), device=)

## 5 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:5)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:5, metric="euclidean")

#new.cluster.ids <- c("Mic/Per/End","Oli","Ex1","Ex2","In1","Ex3","In2","Ast","OPC")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_005pc_oglabels.png"), device=)

## 10 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:10)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:10, metric="euclidean")

#new.cluster.ids <- c("Ex1/End","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","Ex5","OPC","Per/Mic")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_010pc_oglabels.png"), device=)

## 20 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:20)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:20, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ast","In2","In3","Ex3","Ex4","OPC","Ex5","Ex6","In4","Mic","Ex7","Ex8","Ex9","Ex10","Per")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_020pc_oglabels.png"), device=)

## 50 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:50)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:50, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ex3","Ast","In2","Ex4","In3","OPC","Ex5","In4","Mic","Ex6","Ex7","In5","Ex8","Per")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_050pc_oglabels.png"), device=)

## 75 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:75)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:75, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","Ex6","In4","Ex7","Per")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_075pc_oglabels.png"), device=)

## 100 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:100)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:100, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","In4","Ex6","Per")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_100pc_oglabels.png"), device=)

## 150 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:150)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:150, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","OPC","Ex5","Mic/Per","Ex6")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_150pc_oglabels.png"), device=)

## 200 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#All <- FindNeighbors(All, dims = 1:200)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:200, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","OPC","Ex4","Mic/Per","Ex5")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_200pc_oglabels.png"), device=)
