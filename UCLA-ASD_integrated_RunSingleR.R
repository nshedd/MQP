library(Seurat)
library(ggplot2)
library(SingleR)
library(SingleCellExperiment)


## Load SingleR reference dataset
path1 = "/data/rusers/sheddn/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)


BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')


print('Running SingleR...')
BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
                         ref=Lake_SCE,
                         method='cluster',
                         labels=Lake_labels,
                         clusters=Idents(BA4.6),
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts")

print(BA4.6_SingleR$labels)
write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA4.6_cluster_ids.txt")

new.cluster.ids <- BA4.6_SingleR$labels
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

print('Plotting...')

DimPlot(BA4.6, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)
