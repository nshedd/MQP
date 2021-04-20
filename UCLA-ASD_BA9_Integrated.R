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



ASD_BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA9')
ASD_BA9$Group <- "ASD"


CTL_BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')
CTL_BA9$Group <- "CTL"

BA9 <- merge(CTL_BA9, y=ASD_BA9,add.cell.ids=c('CTL','ASD'),project = "UCLA-ASD")

BA9 <- NormalizeData(BA9, normalization.method = "LogNormalize", scale.factor = 10000)
BA9 <- FindVariableFeatures(BA9, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(BA9)
BA9 <- ScaleData(BA9, features = all.genes)

BA9 <- RunPCA(BA9, features = VariableFeatures(object = BA9))

BA9 <- RunHarmony(BA9, "orig.ident")

BA9 <- RunUMAP(BA9, reduction = "harmony", dims = 1:20)

BA9 <- FindNeighbors(BA9, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution=0.5)

## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(BA9, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- BA9@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.07*nrow(BA9@meta.data))  ## Assuming 7% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

BA9 <- doubletFinder_v3(BA9, PCs = 1:20, pN = 0.15, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


p1 <- DimPlot(BA9, reduction = "umap", group.by = "Group")
p2 <- DimPlot(BA9, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_combined.png', width = 16, height = 7)

DimPlot(BA9, reduction = "umap", split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_separate.png', width = 16, height = 7)

BA9.markers <- FindAllMarkers(BA9, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BA9.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(BA9.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(BA9, features = intersection) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits = rev(levels(BA9$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/combined_BA9_dotplot.png", width = 14, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
                       ref=Lake_SCE,
                       method='cluster',
                       labels=Lake_labels,
                       clusters==Idents(BA9),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print(BA9_SingleR$labels)
write.table(BA9_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA9_cluster_ids.txt")

new.cluster.ids <- BA9_SingleR$labels
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

print('Plotting...')

DimPlot(BA9, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)
