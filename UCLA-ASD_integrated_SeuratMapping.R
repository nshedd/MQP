library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)

## Load SingleR reference dataset
path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)



BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

anchors <- FindTransferAnchors(
  reference = Lake,
  query = BA4.6,
  dims = 1:20
)

BA4.6 <- MapQuery(
  anchorset = anchors,
  query = BA4.6,
  reference = Lake,
  refdata = list(celltype = "orig.ident"),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(BA4.6, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, 
    pt.size = 0.5, repel = TRUE)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_tranferredlabel.png', width = 8, height = 7)


saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_TransferLabel.RDS')


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

BA9_SCE <- as.SingleCellExperiment(BA9)
BA9_clust <- Idents(BA9)

print('Running SingleR...')
BA9_SingleR <- SingleR(test=BA9_SCE,
                         ref=Lake_SCE,
                         labels=Lake_labels,
                         clusters=BA9_clust,
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts")

print('Plotting...')
BA9$SingleR.pruned.calls <- BA9_SingleR$pruned.labels
BA9$SingleR.calls <- BA9_SingleR$labels
BA9$SingleR.cluster <- BA9_SingleR$cluster

DimPlot(BA9, group.by="SingleR.cluster", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_SingleR.RDS')
