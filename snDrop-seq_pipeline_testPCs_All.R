library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


## 2 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:2)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:2, metric="euclidean")

new.cluster.ids <- c("In1","Ex1","In2","Ex2","In3","Ex3","Ex4","Ex5","Ast/Mic/Gran","Gran/OPC","Ex6","Ex7","Gran/Purk","Ex8",
                    "Oli1","Oli2","Ex9","Mic/End/Per","Purk","Gran/Ast","Oli3","Ex10","Oli4","Ex11","In4","Ast","Ex12")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_002pc.png"), device=)

## 5 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:5)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:5, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli","Ex3","In1","Gran","In2","Ast","OPC/End","Purk","Ast/Per/Mic","Ex4","Ast/Oli")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_005pc.png"), device=)

## 10 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:10)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:10, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli","Ex3","In1","In2","Gran1","Ast1","Ex4","OPC","Purk","Gran2","Mic","Ast2","End/Per","Ast/Mic")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_010pc.png"), device=)

## 20 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:20)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Ex1","Oli","Ex2","Ex3","Gran1","In1","Ast1","In2","In3","OPC","Ex4",
                     "Purk","Ex5","Ex6","Gran2","Mic","In4","Ast2","End/Per","Ex7")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_020pc.png"), device=)

## 50 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:50)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:50, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli1","Ex3","Gran1","Ast1","In1","In2","Ex4","OPC1","In3","Purk","Ex5",
                     "In4","Mic","Gran2","Ex6","Ast2","Ex7","Per/End","In5","Ex8","OPC2","Oli2")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_050pc.png"), device=)

## 75 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:75)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:75, metric="euclidean")

new.cluster.ids <- c("Ex1","Oli","Ex2","Gran","Ex3","Ast1","In1","In2","Ex4","In3","OPC1","Purk",
                    "Ex5","In4","Mic","Ex6","Ast2","Ex7","Per/End","In5","Ex8","Ex9","OPC2")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_075pc_oglabels.png"), device=)

## 100 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:100)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:100, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli1","Gran","Ex3","Ast1","In1","In2","In3","OPC1","Ex4","Purk",
                    "Ex5","In4","Mic","Ex6","Ast2","Ex7","End/Per","Ex8","Ex9","OPC2","Oli2")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_100pc.png"), device=)

## 150 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:150)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:150, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli","In1","Gran","Ex3","Ast1","In2","OPC1","Purk",
                    "Ex4","Ex5","In3","Mic","Ex6","Ast2","Ex7","End/Per","Ex8","OPC2")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_150pc.png"), device=)

## 200 PCs
All <- readRDS(file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:200)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:200, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Gran/Purk","Oli","Ex3","In1","In2","Ast1","OPC",
                     "Ex4","Mic","Ex5","Ast2","End/Per","Ex6","Ex7","In3")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_200pc.png"), device=)
