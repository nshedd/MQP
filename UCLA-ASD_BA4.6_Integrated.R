library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)


ASD_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')
ASD_BA4.6$Group <- "ASD"


CTL_BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')
CTL_BA4.6$Group <- "CTL"

BA4.6 <- merge(CTL_BA4.6, y=ASD_BA4.6,add.cell.ids=c('CTL','ASD'),project = "UCLA-ASD")

BA4.6 <- NormalizeData(BA4.6, normalization.method = "LogNormalize", scale.factor = 10000)
BA4.6 <- FindVariableFeatures(BA4.6, selection.method = "vst", nfeatures = 2000)


saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6')

# BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6')

all.genes <- rownames(BA4.6)
BA4.6 <- ScaleData(BA4.6, features = all.genes)

BA4.6 <- RunPCA(BA4.6, features = VariableFeatures(object = BA4.6))

BA4.6 <- RunHarmony(BA4.6, "orig.ident")

BA4.6 <- RunUMAP(BA4.6, reduction = "harmony", dims = 1:20)

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

p1 <- DimPlot(BA4.6, reduction = "umap", group.by = "Group")
p2 <- DimPlot(BA4.6, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_combined.png', width = 8, height = 7)

DimPlot(BA4.6, reduction = "umap", split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_separate.png', width = 8, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_DoubletsRemoved')

CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTL.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA9_dotplot.png", width = 14, height = 7)
