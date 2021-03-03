library(dplyr)
library(Seurat)
library(ggplot2)

ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_SampleLabels')

ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")

VlnPlot(ASD, features = c("nFeature_RNA"))
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-VlnPlot-nFeatureRNA.png', width = 14, height = 7)
VlnPlot(ASD, features = c("nCount_RNA"))
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-VlnPlot-nCountRNA.png', width = 14, height = 7)
VlnPlot(ASD, features = c("percent.mt"))
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-VlnPlot-percentmt.png', width = 14, height = 7)

plot1 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/FtrSctr.png", width = 14, height = 7)

ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)

ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ASD)
ASD <- ScaleData(ASD, features = all.genes)

ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))

saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BySample.RDS')

VizDimLoadings(ASD, dims = 1:20, reduction = "pca")
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/VizDimLoadings.png", width = 14, height = 7)

DimPlot(ASD, reduction = "pca")
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/PCADimPlot.png", width = 14, height = 7)

DimHeatmap(ASD, dims = 1:20, balanced = TRUE)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/DimHeatmap.png", width = 14, height = 7)

ASD <- JackStraw(ASD, num.replicate = 100)
ASD <- ScoreJackStraw(ASD, dims = 1:20)
JackStrawPlot(ASD, dims = 1:20)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/JackstrawPlot.png", device=, width = 14, height = 7)

ElbowPlot(ASD)
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/ElbowPlot.png", device=, width = 14, height = 7)

# ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BySample.RDS')

# ASD <- FindNeighbors(ASD, dims = 1:20)
# ASD <- FindClusters(ASD, resolution = 0.5)

ASD <- RunUMAP(ASD, dims = 1:20)

DimPlot(ASD, reduction = "umap")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_samples.png', width = 8, height = 7)

# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample.RDS')
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_OrigIdent.RDS')

ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BySample.RDS')

ASD <- FindNeighbors(ASD, dims = 1:20)
ASD <- FindClusters(ASD, resolution = 0.5)

ASD <- RunUMAP(ASD, dims = 1:20)

DimPlot(ASD, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP.png', width = 8, height = 7)

# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample.RDS')
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample.RDS')
