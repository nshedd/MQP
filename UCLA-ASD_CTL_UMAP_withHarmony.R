library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
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

## BA4/6
CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')

CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

CTL <- subset(CTL, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CTL)
CTL <- ScaleData(CTL, features = all.genes)

CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))

CTL <- RunHarmony(CTL, "orig.ident")

# print("Saving Harmony data BA4.6...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA4.6.RDS')

CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:20)

CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6.png', width = 8, height = 7)

# print("Saving UMAP data BA4.6...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(CTL, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- CTL@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(CTL@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CTL <- doubletFinder_v3(CTL, PCs = 1:20, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(CTL,group.by = colnames(CTL@meta.data)[grep("DF", colnames(CTL@meta.data))])
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_doublets.png', width = 7, height = 7)

CTL <- SubsetData(CTL, cells=rownames(CTL@meta.data)[which(CTL@meta.data$DF.classification == "Singlet")])
DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_doubletsremoved.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA4.6...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

# new.cluster.ids <- c('Ex1','In1','Ex2','Ex3','Ex4','Ex5','Oli','In2','End/Per',
#                     'Ex6','In3','Ex7','OPC','Ast1','In4','Ex8','Mic','Ex9','Ex10',
#                     'Ex11','In5','In6','Ex12','Ex13','In7','Ex14','Ast2','Ex15')
# names(new.cluster.ids) <- levels(CTL)
# CTL <- RenameIdents(CTL, new.cluster.ids)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_labeled.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA4.6...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTL.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA4.6_dotplot.png", width = 14, height = 7)

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


## BA9
CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')

CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

CTL <- subset(CTL, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CTL)
CTL <- ScaleData(CTL, features = all.genes)

CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))

CTL <- RunHarmony(CTL, "orig.ident", max.iter.harmony=20)

# print("Saving Harmony data BA9...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA9.RDS')

CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:20)

CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9.png', width = 8, height = 7)

# print("Saving UMAP data BA9...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9.RDS')

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(CTL, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- CTL@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(CTL@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

CTL <- doubletFinder_v3(CTL, PCs = 1:20, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(CTL,group.by = colnames(CTL@meta.data)[grep("DF", colnames(CTL@meta.data))])
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_doublets.png', width = 7, height = 7)

CTL <- SubsetData(CTL, cells=rownames(CTL@meta.data)[which(CTL@meta.data$DF.classification == "Singlet")])
DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_doubletsremoved.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA9...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','Ast1','OPC','In1','Ex4','Ex5','In2', 
                     'Ex6','In3','Ex7','Ast2','Mic','In4','End/Per1','Ex8','Oli1',
                     'Ex9','Ex10','In5','In6','End/Per2','Ex11','In7','Ex12','Oli2')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_labeled.png', width = 8, height = 7)

print("Saving UMAP data w/o Doublets BA9...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved_Relabeled.RDS')

CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTL.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA9_dotplot.png", width = 14, height = 7)

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
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_SingleRlabel.png', width = 8, height = 7)
