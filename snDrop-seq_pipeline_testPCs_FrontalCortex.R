library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

## 2 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:2)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("In1","Ex1","Ex2","End","Oli1","Ex3","Oli1","Ex4","Ex5","Per/Mic","OPC","Ex6","Ex7","In2","In3","Ast","Oli2")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_002pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_002pc.png"), device=)


##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_002.txt"), sep="\t")

## 5 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:5)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:5, metric="euclidean")

#new.cluster.ids <- c("Mic/Per/End","Oli","Ex1","Ex2","In1","Ex3","In2","Ast","OPC")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_005pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_005pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_005.txt"), sep="\t")

## 10 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:10)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:10, metric="euclidean")

#new.cluster.ids <- c("Ex1/End","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","Ex5","OPC","Per/Mic")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_010pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_010pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_010.txt"), sep="\t")

## 20 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","End/Ex2","In1","Ast","In2","In3","Ex3","Ex4","OPC","Ex5","Ex6","In4","Mic","Ex7","Ex8","Ex9","Ex10","Per")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_020pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_020pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_020.txt"), sep="\t")

## 50 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:50)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:50, metric="euclidean")

#new.cluster.ids <- c("Ex1","Oli","Ex2","In1","Ex3","Ast","In2","Ex4","In3","OPC","Ex5","In4","Mic","Ex6","Ex7","In5","Ex8","Per")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_050pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_050pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_050.txt"), sep="\t")

## 75 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:75)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:75, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","Ex6","In4","Ex7","Per")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_075pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_050pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_075.txt"), sep="\t")

## 100 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:100)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:100, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","Ex4","OPC","Ex5","In3","Mic","In4","Ex6","Per")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_100pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_100pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_100.txt"), sep="\t")

## 150 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:150)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:150, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","In2","Ex3","Ast","Ex4","OPC","Ex5","Mic/Per","Ex6")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_150pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_150pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_150.txt"), sep="\t")

## 200 PCs
FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:200)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:200, metric="euclidean")

#new.cluster.ids <- c("End/Ex1","Oli","Ex2","In1","Ex3","In2","Ast","OPC","Ex4","Mic/Per","Ex5")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/umap_200pc.png"), device=)

## Dotplot
FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("/data/rusers/sheddn/datavis4/dotplot_200pc.png"), device=)

##Save embedding
embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings_200.txt"), sep="\t")
