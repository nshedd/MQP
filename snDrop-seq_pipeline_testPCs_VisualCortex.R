library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


## 2 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:2)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:2, metric="euclidean")

new.cluster.ids <- c("Ast","Ex1","Oli1","In1","Ex2","In2","Ex3","Ex4","Ex5","Oli2","In2","OPC",
                    "In3","Ex6","End","Ex7","Ex8","Mic","Ex9","Ex10","Ex11","Oli3","Ex12")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_002pc.png"), device=)

## 5 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:5)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:5, metric="euclidean")

new.cluster.ids <- c("Ex1/In1","Ex2","Oli","In2","In3","Ex2","Ex3","Ex4","Ast","OPC/Per","End/Mic")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_005pc.png"), device=)

## 10 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:10)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:10, metric="euclidean")

new.cluster.ids <- c("Ex2","Ex1/In1","Ex3","Oli","In2","Ex4","In3","Ast1","OPC","In4","Mic1","Per/End","Ast2/Mic2")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_010pc.png"), device=)

## 20 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:20)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Ex1/In1","Ex2","Oli","Ex3","In2","In3","Ast","Ex4","Ex5","In4","OPC","Ex6","In5","Mic","Ex7","End/Per","Ex8")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_020pc.png"), device=)

## 50 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:50)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:50, metric="euclidean")

new.cluster.ids <- c("Ex2","Ex1/In1","Oli","In2","Ex3","Ex4","Ast","In3","Ex5","In4","OPC","In5","Ex6","Mic","Ex7","End/Per","In6","Ex8")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_050pc.png"), device=)

## 75 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:75)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:75, metric="euclidean")

new.cluster.ids <- c("Ex2","Ex1/In1","Oli","Ex3","In2","In3","Ast","Ex4","In4","OPC","In5","Ex5","Mic","Ex6","End/Per","In6","Ex7")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_075pc.png"), device=)

## 100 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:100)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:100, metric="euclidean")

new.cluster.ids <- c("Ex2","Ex1/In1","Ex3","Oli","In2","Ast","In3","Ex4","In4","OPC","In5","Ex5","Mic","Ex6","End/Per","In6")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_100pc.png"), device=)

## 150 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:150)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:150, metric="euclidean")

new.cluster.ids <- c("Ex1/In1","Ex2","Oli","Ex3","In2","In3","Ast","Ex4","OPC","In4","Ex5","Mic","Ex6","Per/End","Ex7")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_150pc.png"), device=)

## 200 PCs
VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:200)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:200, metric="euclidean")

new.cluster.ids <- c("Ex1/In1","Ex2","Oli","Ex3","In2","In3","Ast","Ex4","OPC","Ex5","Mic","End/Per","Ex6","Ex7")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_200pc.png"), device=)
