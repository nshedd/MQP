library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

## BA4/6
ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')

ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")

ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)

ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ASD)
ASD <- ScaleData(ASD, features = all.genes)

print("Saving Scaled data BA4.6...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_ScaledData_BA4.6.RDS')

ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))

print("Saving PCA data BA4.6...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BA4.6.RDS')

ASD <- RunHarmony(ASD, "orig.ident")

print("Saving Harmony data BA4.6...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BA4.6.RDS')

ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:30)

ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:30) %>% FindClusters()

DimPlot(ASD, group.by="ident")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6.png', width = 8, height = 7)

print("Saving UMAP data BA 46...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6.RDS')


## BA9
ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA9')

ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")

ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)

ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(ASD)
ASD <- ScaleData(ASD, features = all.genes)

print("Saving Scaled data BA9...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_ScaledData_BA9.RDS')

ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))

print("Saving PCA data BA9...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BA9.RDS')

ASD <- RunHarmony(ASD, "orig.ident")

print("Saving Harmony data BA9...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BA9.RDS')

ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:30)

ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:30) %>% FindClusters()

DimPlot(ASD, group.by="ident")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9.png', width = 8, height = 7)

print("Saving UMAP data BA 46...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9.RDS')

# # ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_SampleLabels')
# 
# # ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")
# 
# # ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)
# 
# # all.genes <- rownames(ASD)
# # ASD <- ScaleData(ASD, features = all.genes)
# 
# # saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_ScaledData.RDS')
# 
# ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_PCAprocessed_BySample.RDS')
# 
# ASD <- RunHarmony(ASD, "orig.ident")
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BySample.RDS')
# 
# # ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BySample.RDS')
# 
# ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:30)
# 
# # ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:30) %>% FindClusters()
# 
# DimPlot(ASD, group.by="ident")
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_samples.png', width = 8, height = 7)
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_OrigIdent.RDS')
# 
# ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BySample.RDS')
# 
# ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:30)
# 
# ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:30) %>% FindClusters()
# 
# DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony.png', width = 8, height = 7)
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony.RDS')
