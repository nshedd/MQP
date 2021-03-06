if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")

library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(SingleR)
library(SingleCellExperiment)

matrix = readRDS(path.expand("~/GSE97930_All.RDS"))
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)



CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

CTL_SCE <- as.SingleCellExperiment(CTL)
CTL_clust <- Idents(CTL)

CTL_SingleR <- SingleR(test=CTL_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=CTL_SingleR,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

CTL$SingleR.pruned.calls <- CTL_SingleR$pruned.labels
CTL$SingleR.calls <- CTL_SingleR$labels

DimPlot(CTL, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_SingleRlabel.png', width = 8, height = 7)
