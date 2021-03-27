library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)

## Load SingleR reference dataset
path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)


BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

BA4.6_SCE <- as.SingleCellExperiment(BA4.6)
BA4.6_clust <- Idents(BA4.6)

print('Running SingleR...')
BA4.6_SingleR <- SingleR(test=BA4.6_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=BA4.6_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
BA4.6$SingleR.pruned.calls <- BA4.6_SingleR$pruned.labels
BA4.6$SingleR.calls <- BA4.6_SingleR$labels

DimPlot(BA4.6, group.by="SingleR.cluster", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

BA9_SCE <- as.SingleCellExperiment(BA9)
BA9_clust <- Idents(BA9)

print('Running SingleR...')
BA9_SingleR <- SingleR(test=BA9_SCE,
                         ref=Lake_SCE,
                         labels=Lake_labels,
                         clusters=BA9_clust,
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts")

print('Plotting...')
BA9$SingleR.pruned.calls <- BA9_SingleR$pruned.labels
BA9$SingleR.calls <- BA9_SingleR$labels

DimPlot(BA9, group.by="SingleR.cluster", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)
