library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


## 2 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:2)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:2, metric="euclidean")

new.cluster.ids <- c("In1","Ex1","Ex2","End","Oli1","Ex3","Oli1","Ex4","Ex5","Per/Mic","OPC","Ex6","Ex7","In2","In3","Ast","Oli2")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_002pc.png"), device=)

## 5 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:5)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:5, metric="euclidean")

new.cluster.ids <- c("Mic/Per/End","Oli","Ex1","Ex2","In1","Ex3","In2","Ast","OPC")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_005pc.png"), device=)

## 10 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:10)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:10, metric="euclidean")

new.cluster.ids <- c("Ex1/End","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","Ex5","OPC","Per/Mic")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_010pc.png"), device=)

## 20 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ast","In2","In3","Ex3","Ex4","OPC","Ex5","Ex6","In4","Mic","Ex7","Ex8","Ex9","Ex10","Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_020pc.png"), device=)

## 50 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:50)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:50, metric="euclidean")

new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ex3","Ast","In2","Ex4","In3","OPC","Ex5","In4","Mic","Ex6","Ex7","In5","Ex8","Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_050pc.png"), device=)

## 75 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:75)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:75, metric="euclidean")

new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","Ex6","In4","Ex7","Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_075pc.png"), device=)

## 100 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:100)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:100, metric="euclidean")

new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","In4","Ex6","Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_100pc.png"), device=)

## 150 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:150)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:150, metric="euclidean")

new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","OPC","Ex5","Mic/Per","Ex6")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_150pc.png"), device=)

## 200 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:200)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:200, metric="euclidean")

new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","OPC","Ex4","Mic/Per","Ex5")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_200pc.png"), device=)

## 4 PCs -- optimal
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:4)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:4, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Ex3","In1","End/Mic/Per","Oli1","In2","Ast","Oli2","OPC")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_004pc.png"), device=)