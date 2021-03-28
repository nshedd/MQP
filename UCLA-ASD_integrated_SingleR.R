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
Lake_labels <- Lake$orig.ident


BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

print('Running SingleR...')
BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
                       ref=Lake_SCE,
                       method='cluster',
                       labels=Lake_labels,
                       clusters=Idents(BA4.6),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")


# BA4.6$SingleR.pruned.calls <- BA4.6_SingleR$pruned.labels
# BA4.6$SingleR.calls <- BA4.6_SingleR$labels

print(BA4.6_SingleR$labels)
write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA4.6_cluster_ids.txt")

new.cluster.ids <- BA4.6_SingleR$labels
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

print('Plotting...')

DimPlot(BA4.6, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_SingleR.RDS')


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

print('Running SingleR...')
BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
                       ref=Lake_SCE,
                       method='cluster',
                       labels=Lake_labels,
                       clusters=Idents(BA9),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

# BA9$SingleR.pruned.calls <- BA9_SingleR$pruned.labels
# BA9$SingleR.calls <- BA9_SingleR$labels

print(BA9_SingleR$labels)
write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA9_cluster_ids.txt")

new.cluster.ids <- BA9_SingleR$labels
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

print('Plotting...')

DimPlot(BA9, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_SingleR.RDS')
