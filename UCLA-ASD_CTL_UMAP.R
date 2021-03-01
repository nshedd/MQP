library(dplyr)
library(Seurat)
library(ggplot2)

CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')

CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

VlnPlot(CTL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, )
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-VlnPlot.png', width = 14, height = 7)

plot1 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/FtrSctr.png", width = 14, height = 7)

CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CTL)
CTL <- ScaleData(CTL, features = all.genes)

CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))

VizDimLoadings(CTL, dims = 1:20, reduction = "pca")
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/VizDimLoadings.png", width = 14, height = 7)

DimPlot(CTL, reduction = "pca")
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/PCADimPlot.png"), width = 14, height = 7)

DimHeatmap(CTL, dims = 1:20, balanced = TRUE)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/DimHeatmap.png", width = 14, height = 7)

CTL <- JackStraw(CTL, num.replicate = 100)
CTL <- ScoreJackStraw(CTL, dims = 1:20)
JackStrawPlot(CTL, dims = 1:20)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/JackstrawPlot.png", device=, width = 14, height = 7)

ElbowPlot(CTL)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/ElbowPlot.png", device=, width = 14, height = 7)

CTL <- FindNeighbors(CTL, dims = 1:20)
CTL <- FindClusters(CTL, resolution = 0.5)

CTL <- RunUMAP(CTL, dims = 1:20)

DimPlot(CTL, reduction = "umap")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP.png')
