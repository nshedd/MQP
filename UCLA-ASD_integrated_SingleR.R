library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)

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
