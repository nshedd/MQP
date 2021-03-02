library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')

# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

# VlnPlot(CTL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, )
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-VlnPlot.png', width = 14, height = 7)

# plot1 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot <- plot1 + plot2
# ggsave("/data/rusers/sheddn/UCLA-ASD/plots/FtrSctr.png", width = 14, height = 7)

# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_ScaledData.RDS')

CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_PCAprocessed_BySample.RDS')

CTL <- RunHarmony(CTL, "orig.ident")

saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony.png')

saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony.RDS')
