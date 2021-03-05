library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)

## BA4/6
# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')
# 
# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")
# 
# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)
# 
# CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))
# 
# CTL <- RunHarmony(CTL, "orig.ident")
# 
# print("Saving Harmony data BA4.6...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA4.6.RDS')
# 
# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)
# 
# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()
# 
# DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6.png', width = 8, height = 7)
# 
# print("Saving UMAP data BA4.6...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

print("Loading UMAP data BA4.6...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','Ex4','Oli','Ex5','In1','In2','In3','OPC','Ex6',
                    'Ex7','Ast1','End?','Ex8','Mic','Ex9','In4','Ex10','In5','Ex11','In6','Per?',
                    'In6','Dop?','In7','Ex12','Ast2','Ex13','Ast3','Ex14','Ex15','Ex16','End/Per')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_labeled.png', width = 8, height = 7)

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(CTL, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- CTL@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sobj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CTL <- doubletFinder_v3(CTL, PCs = 1:30, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(CTL,group.by = colnames(CTL@meta.data)[grep("DF", colnames(CTL@meta.data))])
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_doublets.png', width = 7, height = 7)

CTL <- SubsetData(CTL, cells=rownames(CTL@meta.data)[which(CTL@meta.data$DF.classification == "Singlet")])
DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_doubletsremoved.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA4.6...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

# CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# CTL.markers %>% group_by(cluster)
# 
# marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# all_known_marker_genes = marker_gene_table$Human.Gene
# 
# intersection = intersect(CTL.markers$gene, all_known_marker_genes)
# 
# dotplot <- DotPlot(CTL, features = intersection) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
# ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA4.6_dotplot.png", width = 14, height = 7)


## BA9
# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')
# 
# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")
# 
# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)
# 
# CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))
# 
# CTL <- RunHarmony(CTL, "orig.ident", max.iter.harmony=20)
# 
# print("Saving Harmony data BA9...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA9.RDS')
# 
# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)
# 
# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution=0.5)
# 
# DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9.png', width = 8, height = 7)
# 
# print("Saving UMAP data BA9...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9.RDS')

print("Loading UMAP data BA9...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','Ast1','OPC','Ex4','In1','In2','In3','?','Ex5','Ex6','Ex7',
                    'Ast2','Mic','Ex8','End','Oli1','In4','In5','Ex9','Ex10/Ast3','Ex11','Dop?','End/Per?',
                    'Ex12','Ex13','In6','Ex14','In7','Ex15','Ex16','Ex17/Ast4','Ex18','In8','Ex19','In9','Oli2')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_labeled.png', width = 8, height = 7)

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(CTL, PCs = 1:30, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- CTL@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(sobj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CTL <- doubletFinder_v3(CTL, PCs = 1:30, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(CTL,group.by = colnames(CTL@meta.data)[grep("DF", colnames(CTL@meta.data))])
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_doublets.png', width = 7, height = 7)

CTL <- SubsetData(CTL, cells=rownames(CTL@meta.data)[which(CTL@meta.data$DF.classification == "Singlet")])
DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_doubletsremoved.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA9...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')
# 
# CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# CTL.markers %>% group_by(cluster)
# 
# marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# all_known_marker_genes = marker_gene_table$Human.Gene
# 
# intersection = intersect(CTL.markers$gene, all_known_marker_genes)
# 
# dotplot <- DotPlot(CTL, features = intersection) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
# ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA9_dotplot.png", width = 14, height = 7)



# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')

# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_ScaledData.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_PCAprocessed_BySample.RDS')

# CTL <- RunHarmony(CTL, "orig.ident")

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()

# DimPlot(CTL, group.by="ident")
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_samples.png', width = 8, height = 7)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_OrigIdent.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()

# DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony.png', width = 8, height = 7)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony.RDS')
