library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")
FrontalCortex_VlnPlot <- VlnPlot(FrontalCortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("~/Lake/FrontalCortex/qcvlnplot_GSE97930_FrontalCortex_Seurat.png"), device=)

plot1 <- FeatureScatter(FrontalCortex, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(FrontalCortex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("~/Lake/FrontalCortex/featurescatter_GSE97930_FrontalCortex_Seurat.png"), device=)

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

pcs_plot_FrontalCortex <- VizDimLoadings(FrontalCortex, dims = 1:2, reduction = "pca")
ggsave(path.expand("~/Lake/FrontalCortex/qcdimloadings_GSE97930_FrontalCortex_Seurat.png"), device=)

pca_plot_FrontalCortex <- DimPlot(FrontalCortex, reduction = "pca")
ggsave(path.expand("~/Lake/FrontalCortex/pcascatter_GSE97930_FrontalCortex_Seurat.png"), device=)

pca_heatmap_FrontalCortex <- DimHeatmap(FrontalCortex, dims = 1, balanced = TRUE)
ggsave(path.expand("~/Lake/FrontalCortex/pcaheatmap_GSE97930_FrontalCortex_Seurat.png"), plot=pca_heatmap_FrontalCortex, device=)

FrontalCortex <- JackStraw(FrontalCortex, num.replicate = 100)
FrontalCortex <- ScoreJackStraw(FrontalCortex, dims = 1:20)
Jackstraw_FrontalCortex <- JackStrawPlot(FrontalCortex, dims = 1:15)
ggsave(path.expand("~/Lake/FrontalCortex/jackstrawplot_GSE97930_FrontalCortex_Seurat.png"), device=)

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:10)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex@active.ident <- plyr::mapvalues(x = FrontalCortex@active.ident)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = FrontalCortex.markers %>% group_by(cluster)

saveRDS(FrontalCortex, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_default.png"), device=)

featureplot_FrontalCortex < - FeaturePlot(FrontalCortex, features = c("SYT1", "SYT1", "RBFOX3", "GAD2", "SLC1A3", "GRIN2B", "PCDH15", "MBP", "APBB1IP", "SLC6A1", "PCDH15", "SLC1A3"))
ggsave(path.expand("~/Lake/FrontalCortex/featureplot_GSE97930_FrontalCortex_Seurat_default.png"), device=)

write.table(diff_expressed, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_differentiallyexpressed.txt"), sep="\t")
