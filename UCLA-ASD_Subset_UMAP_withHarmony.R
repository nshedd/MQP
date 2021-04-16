library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)


## Load SingleR reference dataset
path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)



ASD_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/ASD_BA4.6')
ASD_BA4.6$Group <- "ASD"


CTL_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/CTL_BA4.6')
CTL_BA4.6$Group <- "CTL"

BA4.6 <- merge(CTL_BA4.6, y=ASD_BA4.6,add.cell.ids=c('CTL','ASD'),project = "UCLA-ASD")

BA4.6 <- NormalizeData(BA4.6, normalization.method = "LogNormalize", scale.factor = 10000)
BA4.6 <- FindVariableFeatures(BA4.6, selection.method = "vst", nfeatures = 2000)


saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6')

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6')

all.genes <- rownames(BA4.6)
BA4.6 <- ScaleData(BA4.6, features = all.genes)

BA4.6 <- RunPCA(BA4.6, features = VariableFeatures(object = BA4.6))

BA4.6_noharmony <- FindNeighbors(BA4.6, dims = 1:20) %>% FindClusters(resolution=1)
BA4.6_noharmony <- RunUMAP(BA4.6, dims = 1:20, metric="euclidean")


BA4.6 <- RunHarmony(BA4.6, "orig.ident")

BA4.6 <- RunUMAP(BA4.6, reduction = "harmony", dims = 1:20, metric="euclidean")

BA4.6 <- FindNeighbors(BA4.6, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=1)

before <- DimPlot(BA4.6_noharmony, reduction = "umap", group.by="Sample")
after <- DimPlot(BA4.6, reduction = "umap", group.by="Sample")
before + after
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_before_and_after_Harmony.png', width = 16, height = 7)

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

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_DoubletsRemoved')

# BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_DoubletsRemoved')



p1 <- DimPlot(BA4.6, reduction = "umap", group.by = "Group")
p2 <- DimPlot(BA4.6, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA4.6_integrated_combined.png', width = 16, height = 7)

DimPlot(BA4.6, reduction = "umap", split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA4.6_integrated_separate.png', width = 16, height = 7)

BA4.6.markers <- FindAllMarkers(BA4.6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BA4.6.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(BA4.6.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(BA4.6, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(BA4.6$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/subset/plots/combined_BA4.6_dotplot.png", width = 14, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_WithDEGs.RDS')


print('Running SingleR...')
BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
                       ref=Lake_SCE,
                       method='cluster',
                       labels=Lake_labels,
                       clusters=Idents(BA4.6),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print(BA4.6_SingleR$labels)
write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/subset/data/BA4.6_cluster_ids.txt")

new.cluster.ids <- BA4.6_SingleR$labels
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

print('Plotting...')

DimPlot(BA4.6, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_SingleR.RDS')








ASD_BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/ASD_BA9')
ASD_BA9$Group <- "ASD"


CTL_BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/CTL_BA9')
CTL_BA9$Group <- "CTL"

BA9 <- merge(CTL_BA9, y=ASD_BA9,add.cell.ids=c('CTL','ASD'),project = "UCLA-ASD")

BA9 <- NormalizeData(BA9, normalization.method = "LogNormalize", scale.factor = 10000)
BA9 <- FindVariableFeatures(BA9, selection.method = "vst", nfeatures = 2000)


saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9')

BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9')

all.genes <- rownames(BA9)
BA9 <- ScaleData(BA9, features = all.genes)

BA9 <- RunPCA(BA9, features = VariableFeatures(object = BA9))

BA9_noharmony <- FindNeighbors(BA9, dims = 1:20) %>% FindClusters(resolution=1)
BA9_noharmony <- RunUMAP(BA9, dims = 1:20, metric="euclidean")


BA9 <- RunHarmony(BA9, "orig.ident")

BA9 <- RunUMAP(BA9, reduction = "harmony", dims = 1:20, metric="euclidean")

BA9 <- FindNeighbors(BA9, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=1)

before <- DimPlot(BA9_noharmony, reduction = "umap", group.by="Sample")
after <- DimPlot(BA9, reduction = "umap", group.by="Sample")
before + after
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_before_and_after_Harmony.png', width = 16, height = 7)

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(BA9, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- BA9@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(BA9@meta.data))  ## Assuming 15% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

BA9 <- doubletFinder_v3(BA9, PCs = 1:20, pN = 0.15, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_DoubletsRemoved')

# BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_DoubletsRemoved')



p1 <- DimPlot(BA9, reduction = "umap", group.by = "Group")
p2 <- DimPlot(BA9, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA9_integrated_combined.png', width = 16, height = 7)

DimPlot(BA9, reduction = "umap", split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA9_integrated_separate.png', width = 16, height = 7)

BA9.markers <- FindAllMarkers(BA9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BA9.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(BA9.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(BA9, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(BA9$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/subset/plots/combined_BA9_dotplot.png", width = 14, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_WithDEGs.RDS')


print('Running SingleR...')
BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
                         ref=Lake_SCE,
                         method='cluster',
                         labels=Lake_labels,
                         clusters=Idents(BA9),
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts")

print(BA9_SingleR$labels)
write.table(BA9_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/subset/data/BA9_cluster_ids.txt")

new.cluster.ids <- BA9_SingleR$labels
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

print('Plotting...')

DimPlot(BA9, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_SingleR.RDS')
