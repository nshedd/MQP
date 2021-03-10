library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)

# ## Load reference dataset
# path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")
# 
# matrix = read.table(path1, header=TRUE, row.names=1)
# Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)
# 
# Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")
# 
# Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
# 
# Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# Lake_SCE <- as.SingleCellExperiment(Lake)
# Lake_labels <- Idents(Lake)

## BA4/6
# ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')
# 
# ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")
# 
# ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(ASD)
# ASD <- ScaleData(ASD, features = all.genes)
# 
# ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))
# 
# ASD <- RunHarmony(ASD, "orig.ident")
# 
# print("Saving Harmony data BA4.6...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BA4.6.RDS')
# 
# ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:20)
# 
# ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)
# 
# print("Saving UMAP data BA 4.6...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

print("Loading UMAP data BA4.6...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6.png', width = 8, height = 7)

# ## Doublet removal
# ## pK indetification
# sweep.res.list_pbmc <- paramSweep_v3(ASD, PCs = 1:20, sct = FALSE)
# sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
# bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
# bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))
# 
# ## Doublet proportion estimate
# annotations <- ASD@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
# nExp_poi <- round(0.15*nrow(ASD@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# ASD <- doubletFinder_v3(ASD, PCs = 1:20, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# 
# DimPlot(ASD,group.by = colnames(ASD@meta.data)[grep("DF", colnames(ASD@meta.data))])
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6_doublets.png', width = 7, height = 7)
# 
# ASD <- SubsetData(ASD, cells=rownames(ASD@meta.data)[which(ASD@meta.data$DF.classification == "Singlet")])
# DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6_doubletsremoved.png', width = 8, height = 7)
# 
# print("Saving UMAP data w/o Doublets BA4.6...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

# print("Loading UMAP data w/o Doublets BA4.6...")
# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')
# 
# new.cluster.ids <- c('Ex1','Ex2','Ex3','Ex4','Oli','Ex5','In1','In2','In3','OPC','Ex6',
#                     'Ex7','Ast1','End?','Ex8','Mic','Ex9','In4','Ex10','In5','Ex11','In6','Per?',
#                     'In6','Dop?','In7','Ex12','Ast2','Ex13','Ast3','Ex14','Ex15','Ex16','End/Per')
# names(new.cluster.ids) <- levels(ASD)
# ASD <- RenameIdents(ASD, new.cluster.ids)
# 
# DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6_labeled.png', width = 8, height = 7)

# ASD.markers <- FindAllMarkers(ASD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ASD.markers %>% group_by(cluster)
# 
# marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# all_known_marker_genes = marker_gene_table$Human.Gene
# 
# intersection = intersect(ASD.markers$gene, all_known_marker_genes)
# 
# dotplot <- DotPlot(ASD, features = intersection) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   scale_y_discrete(limits = rev(levels(ASD$seurat_clusters)))
# ggsave("/data/rusers/sheddn/UCLA-ASD/plots/ASD_BA4.6_dotplot.png", width = 14, height = 7)
# 
# ASD_SCE <- as.SingleCellExperiment(ASD)
# ASD_clust <- Idents(ASD)
# 
# print('Running SingleR...')
# ASD_SingleR <- SingleR(test=ASD_SCE,
#                        ref=Lake_SCE,
#                        labels=Lake_labels,
#                        clusters=ASD_clust,
#                        assay.type.test = "logcounts",
#                        assay.type.ref = "logcounts")
# 
# print('Plotting...')
# ASD$SingleR.pruned.calls <- ASD_SingleR$pruned.labels
# ASD$SingleR.calls <- ASD_SingleR$labels
# 
# DimPlot(ASD, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA4.6_SingleRlabel.png', width = 8, height = 7)


## BA4/6
# ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA9')
# 
# ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")
# 
# ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(ASD)
# ASD <- ScaleData(ASD, features = all.genes)
# 
# ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))
# 
# ASD <- RunHarmony(ASD, "orig.ident")
# 
# print("Saving Harmony data BA9...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_Harmonyprocessed_BA9.RDS')
# 
# ASD <- RunUMAP(ASD, reduction = "harmony", dims = 1:20)
# 
# ASD <- FindNeighbors(ASD, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)
# 
# print("Saving UMAP data BA 9...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9.RDS')

print("Loading UMAP data BA9...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9.RDS')

DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9.png', width = 8, height = 7)

# ## Doublet removal
# ## pK indetification
# sweep.res.list_pbmc <- paramSweep_v3(ASD, PCs = 1:20, sct = FALSE)
# sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
# bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
# bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))
# 
# ## Doublet proportion estimate
# annotations <- ASD@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
# nExp_poi <- round(0.15*nrow(ASD@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# ASD <- doubletFinder_v3(ASD, PCs = 1:20, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# 
# DimPlot(ASD,group.by = colnames(ASD@meta.data)[grep("DF", colnames(ASD@meta.data))])
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9_doublets.png', width = 7, height = 7)
# 
# ASD <- SubsetData(ASD, cells=rownames(ASD@meta.data)[which(ASD@meta.data$DF.classification == "Singlet")])
# DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9_doubletsremoved.png', width = 8, height = 7)
# 
# print("Saving UMAP data w/o Doublets BA9...")
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')

# print("Loading UMAP data w/o Doublets BA9...")
# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')
# 
# new.cluster.ids <- c('Ex1','Ex2','Ex3','Ex4','Oli','Ex5','In1','In2','In3','OPC','Ex6',
#                     'Ex7','Ast1','End?','Ex8','Mic','Ex9','In4','Ex10','In5','Ex11','In6','Per?',
#                     'In6','Dop?','In7','Ex12','Ast2','Ex13','Ast3','Ex14','Ex15','Ex16','End/Per')
# names(new.cluster.ids) <- levels(ASD)
# ASD <- RenameIdents(ASD, new.cluster.ids)
# 
# DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9_labeled.png', width = 8, height = 7)

# ASD.markers <- FindAllMarkers(ASD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ASD.markers %>% group_by(cluster)
# 
# marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# all_known_marker_genes = marker_gene_table$Human.Gene
# 
# intersection = intersect(ASD.markers$gene, all_known_marker_genes)
# 
# dotplot <- DotPlot(ASD, features = intersection) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   scale_y_discrete(limits = rev(levels(ASD$seurat_clusters)))
# ggsave("/data/rusers/sheddn/UCLA-ASD/plots/ASD_BA9_dotplot.png", width = 14, height = 7)
# 
# ASD_SCE <- as.SingleCellExperiment(ASD)
# ASD_clust <- Idents(ASD)
# 
# print('Running SingleR...')
# ASD_SingleR <- SingleR(test=ASD_SCE,
#                        ref=Lake_SCE,
#                        labels=Lake_labels,
#                        clusters=ASD_clust,
#                        assay.type.test = "logcounts",
#                        assay.type.ref = "logcounts")
# 
# print('Plotting...')
# ASD$SingleR.pruned.calls <- ASD_SingleR$pruned.labels
# ASD$SingleR.calls <- ASD_SingleR$labels
# 
# DimPlot(ASD, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9_SingleRlabel.png', width = 8, height = 7)
