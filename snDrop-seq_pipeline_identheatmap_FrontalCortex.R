library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scclusteval)

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")		

FrontalCortex <- subset(FrontalCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)		

FrontalCortex <- FindVariableFeatures(FrontalCortex, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(FrontalCortex)		
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)		

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))		

saveRDS(FrontalCortex, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))


## My labels
FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Ex1","Ex2","Oli1","Ast","Ex3","In1","Oli2","In2","In3","Ex4","OPC","Ex5",
                    "Ex6","Ex7","In4","In5","Mic","Ex8","Ex9","Ex10","Ex11","Ex12","End/Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct.png"), device=)


## Paper labels

FrontalCortex2 <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex2 <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct_oglabels.png"), device=)


## Heatmap comparing labels

oglabels <- Idents(FrontalCortex2)
oglabels <- factor(oglabels)

newlabels <- Idents(FrontalCortex)
newlabels <- factor(newlabels)


png(file=path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_Seurat_JaccardHeatmap.png"))
heatmap <- PairWiseJaccardSetsHeatmap(oglabels, newlabels)
dev.off()




