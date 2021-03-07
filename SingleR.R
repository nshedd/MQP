library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(SingleR)
library(SingleCellExperiment)

## Load reference dataset
path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)

## SingleR on BA4/6 data
print('Loading UCLA-ASD BA4.6 data...')
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

CTL_SCE <- as.SingleCellExperiment(CTL)
CTL_clust <- Idents(CTL)

print('Running SingleR...')
CTL_SingleR <- SingleR(test=CTL_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=CTL_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
CTL$SingleR.pruned.calls <- CTL_SingleR$pruned.labels
CTL$SingleR.calls <- CTL_SingleR$labels

DimPlot(CTL, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_SingleRlabel.png', width = 8, height = 7)

## SingleR on BA9 data
print('Loading UCLA-ASD BA9 data...')
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')

CTL_SCE <- as.SingleCellExperiment(CTL)
CTL_clust <- Idents(CTL)

print('Running SingleR...')
CTL_SingleR <- SingleR(test=CTL_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=CTL_SingleR,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
CTL$SingleR.pruned.calls <- CTL_SingleR$pruned.labels
CTL$SingleR.calls <- CTL_SingleR$labels

DimPlot(CTL, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_SingleRlabel.png', width = 8, height = 7)
