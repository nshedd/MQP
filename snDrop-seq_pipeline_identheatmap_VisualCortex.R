devtools::install_github("nshedd/scclusteval")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scclusteval)

path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

## Paper labels

VisualCortex2 <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

VisualCortex2[["percent.mt"]] <- PercentageFeatureSet(VisualCortex2, pattern = "^MT-")		

VisualCortex2 <- subset(VisualCortex2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

VisualCortex2 <- NormalizeData(VisualCortex2, normalization.method = "LogNormalize", scale.factor = 10000)		

VisualCortex2 <- FindVariableFeatures(VisualCortex2, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(VisualCortex2)		
VisualCortex2 <- ScaleData(VisualCortex2, features = all_cells)		

VisualCortex2 <- RunPCA(VisualCortex2, features = VariableFeatures(object = VisualCortex2))	

VisualCortex2 <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

plot = DimPlot(VisualCortex2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_findct_oglabels.png"), device=)

## My labels
VisualCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT-")		

VisualCortex <- subset(VisualCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

VisualCortex <- NormalizeData(VisualCortex, normalization.method = "LogNormalize", scale.factor = 10000)		

VisualCortex <- FindVariableFeatures(VisualCortex, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(VisualCortex)		
VisualCortex <- ScaleData(VisualCortex, features = all_cells)		

VisualCortex <- RunPCA(VisualCortex, features = VariableFeatures(object = VisualCortex2))	

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:20)
VisualCortex <- FindClusters(VisualCortex, resolution = 1)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Oli","Neuron","Ex1","Ex2","Ex3","In1","Ex4","Ast","Ex5","In2","In3",
                    "Ex6","OPC","In4","Ex7","Ex8","Mic","In5","Ex9","End/Per","Ex10","Ast")
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
