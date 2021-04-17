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

## BA4.6
# CTL

print('23BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/23BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW23 = CreateSeuratObject(counts = expression_matrix, project='23BW')
BW23[["percent.mt"]] <- PercentageFeatureSet(BW23, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW23$nCount_RNA)>=minCov){
  countLOW=min(BW23$nCount_RNA)
}else{
  countLOW=quantile(BW23$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW23$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW23$nFeature_RNA, prob=0.01)
BW23 <- subset(BW23, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW23$Sample <- "BW23"
BW23 <- NormalizeData(BW23, verbose = FALSE)
BW23 <- FindVariableFeatures(BW23, selection.method = "vst", nfeatures = 2000)
BW23$Group <- "CTL"

print('25BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/25BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW25 = CreateSeuratObject(counts = expression_matrix, project='25BW')
BW25[["percent.mt"]] <- PercentageFeatureSet(BW25, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW25$nCount_RNA)>=minCov){
  countLOW=min(BW25$nCount_RNA)
}else{
  countLOW=quantile(BW25$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW25$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW25$nFeature_RNA, prob=0.01)
BW25 <- subset(BW25, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW25$Sample <- "BW25"
BW25 <- NormalizeData(BW25, verbose = FALSE)
BW25 <- FindVariableFeatures(BW25, selection.method = "vst", nfeatures = 2000)
BW25$Group <- "CTL"

print('30BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/30BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW30 = CreateSeuratObject(counts = expression_matrix, project='30BW')
BW30[["percent.mt"]] <- PercentageFeatureSet(BW30, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW30$nCount_RNA)>=minCov){
  countLOW=min(BW30$nCount_RNA)
}else{
  countLOW=quantile(BW30$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW30$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW30$nFeature_RNA, prob=0.01)
BW30 <- subset(BW30, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW30$Sample <- "BW30"
BW30 <- NormalizeData(BW30, verbose = FALSE)
BW30 <- FindVariableFeatures(BW30, selection.method = "vst", nfeatures = 2000)
BW30$Group <- "CTL"

# ASD

print('60BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/60BW_S30/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW60 = CreateSeuratObject(counts = expression_matrix, project='60BW')
BW60[["percent.mt"]] <- PercentageFeatureSet(BW60, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW60$nCount_RNA)>=minCov){
  countLOW=min(BW60$nCount_RNA)
}else{
  countLOW=quantile(BW60$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW60$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW60$nFeature_RNA, prob=0.01)
BW60 <- subset(BW60, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW60$Sample <- "BW60"
BW60 <- NormalizeData(BW60, verbose = FALSE)
BW60 <- FindVariableFeatures(BW60, selection.method = "vst", nfeatures = 2000)
BW60$Group <- "ASD"

print('31BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/31BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW31 = CreateSeuratObject(counts = expression_matrix, project='31BW')
BW31[["percent.mt"]] <- PercentageFeatureSet(BW31, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW31$nCount_RNA)>=minCov){
  countLOW=min(BW31$nCount_RNA)
}else{
  countLOW=quantile(BW31$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW31$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW31$nFeature_RNA, prob=0.01)
BW31 <- subset(BW31, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW31$Sample <- "BW31"
BW31 <- NormalizeData(BW31, verbose = FALSE)
BW31 <- FindVariableFeatures(BW31, selection.method = "vst", nfeatures = 2000)
BW31$Group <- "ASD"

print('32BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/32BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW32 = CreateSeuratObject(counts = expression_matrix, project='32BW')
BW32[["percent.mt"]] <- PercentageFeatureSet(BW32, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW32$nCount_RNA)>=minCov){
  countLOW=min(BW32$nCount_RNA)
}else{
  countLOW=quantile(BW32$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW32$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW32$nFeature_RNA, prob=0.01)
BW32 <- subset(BW32, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW32$Sample <- "BW32"
BW32 <- NormalizeData(BW32, verbose = FALSE)
BW32 <- FindVariableFeatures(BW32, selection.method = "vst", nfeatures = 2000)
BW32$Group <- "ASD"


anchors <- FindIntegrationAnchors(object.list = list(BW31,BW32,BW60, BW23, BW25, BW30), dims = 1:20)
BA4.6 <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(BA4.6) <- "integrated"

# Run the standard workflow for visualization and clustering
BA4.6 <- ScaleData(BA4.6, verbose = FALSE)
BA4.6 <- RunPCA(BA4.6, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
BA4.6 <- RunUMAP(BA4.6, reduction = "pca", dims = 1:20, metric="euclidean")
BA4.6 <- FindNeighbors(BA4.6, reduction = "pca", dims = 1:20)
BA4.6 <- FindClusters(BA4.6, resolution = 0.5)

p1 <- DimPlot(BA4.6, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(BA4.6, reduction = "umap", label = TRUE)
p1+p2

ggsave('/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/BA4.6_UMAP_bysample.png', width = 15, height = 7)


# ## pK indetification
# sweep.res.list_pbmc <- paramSweep_v3(BA4.6, PCs = 1:20, sct = FALSE)
# sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
# bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
# bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))
# 
# ## Doublet proportion estimate
# annotations <- BA4.6@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
# nExp_poi <- round(0.15*nrow(BA4.6@meta.data))  ## Assuming 15% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# BA4.6 <- doubletFinder_v3(BA4.6, PCs = 1:20, pN = 0.15, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


DimPlot(BA4.6, reduction = "umap", split.by = "Group")

# saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_DoubletsRemoved.RDS')

BA4.6.markers <- FindAllMarkers(BA4.6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BA4.6.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(BA4.6.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(BA4.6, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(BA4.6$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/combined_BA4.6_dotplot.png", width = 14, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_WithDEGs.RDS')


print('Running SingleR...')
BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
                         ref=Lake_SCE,
                         method='cluster',
                         labels=Lake_labels,
                         clusters=Idents(BA4.6),
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts")

print(BA4.6_SingleR$labels)
write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_cluster_ids.txt")

new.cluster.ids <- BA4.6_SingleR$labels
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

print('Plotting...')

DimPlot(BA4.6, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/UMAP_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_SingleR.RDS')

DefaultAssay(BA4.6) <- "RNA"
nk.markers <- FindConservedMarkers(BA4.6, ident.1 = 0, grouping.var = "Group", verbose = FALSE)
head(nk.markers)

write.table(nk.markers,  '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_conserved_markers_cluster0.txt')


## BA9
#CTL

print('55BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/55BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW55 = CreateSeuratObject(counts = expression_matrix, project='55BW')
BW55[["percent.mt"]] <- PercentageFeatureSet(BW55, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW55$nCount_RNA)>=minCov){
  countLOW=min(BW55$nCount_RNA)
}else{
  countLOW=quantile(BW55$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW55$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW55$nFeature_RNA, prob=0.01)
BW55 <- subset(BW55, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW55$Sample <- "BW55"
BW55 <- NormalizeData(BW55, verbose = FALSE)
BW55 <- FindVariableFeatures(BW55, selection.method = "vst", nfeatures = 2000)
BW55$Group <- "CTL"

print('64BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/64BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW64 = CreateSeuratObject(counts = expression_matrix, project='64BW')
BW64[["percent.mt"]] <- PercentageFeatureSet(BW64, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW64$nCount_RNA)>=minCov){
  countLOW=min(BW64$nCount_RNA)
}else{
  countLOW=quantile(BW64$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW64$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW64$nFeature_RNA, prob=0.01)
BW64 <- subset(BW64, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW64$Sample <- "BW64"
BW64 <- NormalizeData(BW64, verbose = FALSE)
BW64 <- FindVariableFeatures(BW64, selection.method = "vst", nfeatures = 2000)
BW64$Group <- "CTL"

print('69BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/69BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW69 = CreateSeuratObject(counts = expression_matrix, project='69BW')
BW69[["percent.mt"]] <- PercentageFeatureSet(BW69, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW69$nCount_RNA)>=minCov){
  countLOW=min(BW69$nCount_RNA)
}else{
  countLOW=quantile(BW69$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW69$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW69$nFeature_RNA, prob=0.01)
BW69 <- subset(BW69, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW69$Sample <- "BW69"
BW69 <- NormalizeData(BW69, verbose = FALSE)
BW69 <- FindVariableFeatures(BW69, selection.method = "vst", nfeatures = 2000)
BW69$Group <- "CTL"

# ASD

print('66BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/66BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW66 = CreateSeuratObject(counts = expression_matrix, project='66BW')
BW66[["percent.mt"]] <- PercentageFeatureSet(BW66, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW66$nCount_RNA)>=minCov){
  countLOW=min(BW66$nCount_RNA)
}else{
  countLOW=quantile(BW66$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW66$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW66$nFeature_RNA, prob=0.01)
BW66 <- subset(BW66, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW66$Sample <- "BW66"
BW66 <- NormalizeData(BW66, verbose = FALSE)
BW66 <- FindVariableFeatures(BW66, selection.method = "vst", nfeatures = 2000)
BW66$Group <- "ASD"

print('54BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/54BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW54 = CreateSeuratObject(counts = expression_matrix, project='54BW')
BW54[["percent.mt"]] <- PercentageFeatureSet(BW54, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW54$nCount_RNA)>=minCov){
  countLOW=min(BW54$nCount_RNA)
}else{
  countLOW=quantile(BW54$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW54$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW54$nFeature_RNA, prob=0.01)
BW54 <- subset(BW54, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW54$Sample <- "BW54"
BW54 <- NormalizeData(BW54, verbose = FALSE)
BW54 <- FindVariableFeatures(BW54, selection.method = "vst", nfeatures = 2000)
BW54$Group <- "ASD"

print('34BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/34BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW34 = CreateSeuratObject(counts = expression_matrix, project='34BW')
BW34[["percent.mt"]] <- PercentageFeatureSet(BW34, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW34$nCount_RNA)>=minCov){
  countLOW=min(BW34$nCount_RNA)
}else{
  countLOW=quantile(BW34$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW34$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW34$nFeature_RNA, prob=0.01)
BW34 <- subset(BW34, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW34$Sample <- "BW34"
BW34 <- NormalizeData(BW34, verbose = FALSE)
BW34 <- FindVariableFeatures(BW34, selection.method = "vst", nfeatures = 2000)
BW34$Group <- "ASD"


anchors <- FindIntegrationAnchors(object.list = list(BW34,BW54,BW66, BW69, BW64, BW55), dims = 1:20)
BA9 <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(BA9) <- "integrated"

# Run the standard workflow for visualization and clustering
BA9 <- ScaleData(BA9, verbose = FALSE)
BA9 <- RunPCA(BA9, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
BA9 <- RunUMAP(BA9, reduction = "pca", dims = 1:20, metric="euclidean")
BA9 <- FindNeighbors(BA9, reduction = "pca", dims = 1:20)
BA9 <- FindClusters(BA9, resolution = 0.5)

p1 <- DimPlot(BA9, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(BA9, reduction = "umap", label = TRUE)
p1+p2

ggsave('/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/BA9_UMAP_bysample.png', width = 15, height = 7)


# ## pK indetification
# sweep.res.list_pbmc <- paramSweep_v3(BA9, PCs = 1:20, sct = FALSE)
# sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
# bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
# bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))
# 
# ## Doublet proportion estimate
# annotations <- BA9@meta.data$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
# nExp_poi <- round(0.15*nrow(BA9@meta.data))  ## Assuming 15% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
# BA9 <- doubletFinder_v3(BA9, PCs = 1:20, pN = 0.15, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# 

DimPlot(BA9, reduction = "umap", split.by = "Group")

# saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA9_DoubletsRemoved.RDS')

BA9.markers <- FindAllMarkers(BA9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BA9.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(BA9.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(BA9, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(BA9$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/combined_BA9_dotplot.png", width = 14, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA9_WithDEGs.RDS')


print('Running SingleR...')
BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
                       ref=Lake_SCE,
                       method='cluster',
                       labels=Lake_labels,
                       clusters=Idents(BA9),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print(BA9_SingleR$labels)
write.table(BA9_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA9_cluster_ids.txt")

new.cluster.ids <- BA9_SingleR$labels
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

print('Plotting...')

DimPlot(BA9, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/UMAP_BA9_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA9_SingleR.RDS')

DefaultAssay(BA9) <- "RNA"
nk.markers <- FindConservedMarkers(BA9, ident.1 = 0, grouping.var = "Group", verbose = FALSE)
head(nk.markers)

write.table(nk.markers,  '/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA9_conserved_markers_cluster0.txt')


