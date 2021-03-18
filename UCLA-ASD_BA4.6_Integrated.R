library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)

ASD_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')

new.cluster.ids <- c('ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD')
names(new.cluster.ids) <- levels(ASD_BA4.6)
ASD_BA4.6 <- RenameIdents(ASD_BA4.6, new.cluster.ids)


CTL_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')

new.cluster.ids <- c('CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL')
names(new.cluster.ids) <- levels(CTL_BA4.6)
CTL_BA4.6 <- RenameIdents(CTL_BA4.6, new.cluster.ids)


BA4.6 <- merge(CTL_BA4.6, y=ASD_BA4.6, add.cell.ids=c('CTL','ASD'), project='UCLA-ASD')


saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6')

BA4.6[["percent.mt"]] <- PercentageFeatureSet(BA4.6, pattern = "^MT-")

BA4.6 <- NormalizeData(BA4.6, normalization.method = "LogNormalize", scale.factor = 10000)

BA4.6 <- FindVariableFeatures(BA4.6, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(BA4.6)
BA4.6 <- ScaleData(BA4.6, features = all.genes)

BA4.6 <- RunPCA(BA4.6, features = VariableFeatures(object = BA4.6))

BA4.6 <- RunHarmony(BA4.6, "orig.ident")

BA4.6 <- RunUMAP(CTL_BA4.6, reduction = "harmony", dims = 1:20)

BA4.6 <- FindNeighbors(CTL_BA4.6, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(BA4.6, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- BA4.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(BA4.6@meta.data))  ## Assuming 15% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

BA4.6 <- doubletFinder_v3(BA4.6, PCs = 1:20, pN = 0.15, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

p1 <- DimPlot(BA4.6, reduction = "umap", group.by = "ident")
p2 <- DimPlot(BA4.6, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_combined.png', width = 8, height = 7)

DimPlot(BA4.6, reduction = "umap", split.by = "ident")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_separate.png', width = 8, height = 7)





