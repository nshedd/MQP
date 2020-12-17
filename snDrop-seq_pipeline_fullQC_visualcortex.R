library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

VisualCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT-")
VisualCortex_VlnPlot <- VlnPlot(VisualCortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("~/Lake/VisualCortex/qcvlnplot_GSE97930_VisualCortex_Seurat.png"), device=)

plot1 <- FeatureScatter(VisualCortex, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(VisualCortex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("~/Lake/VisualCortex/featurescatter_GSE97930_VisualCortex_Seurat.png"), device=)

VisualCortex <- NormalizeData(VisualCortex, normalization.method = "LogNormalize", scale.factor = 10000)

VisualCortex <- FindVariableFeatures(object = VisualCortex)

all_cells <- rownames(VisualCortex)
VisualCortex <- ScaleData(VisualCortex, features = all_cells)

VisualCortex <- RunPCA(VisualCortex, features = VariableFeatures(object = VisualCortex))

pcs_plot_VisualCortex <- VizDimLoadings(VisualCortex, dims = 1:2, reduction = "pca")
ggsave(path.expand("~/Lake/VisualCortex/qcdimloadings_GSE97930_VisualCortex_Seurat.png"), device=)

pca_plot_VisualCortex <- DimPlot(VisualCortex, reduction = "pca")
ggsave(path.expand("~/Lake/VisualCortex/pcascatter_GSE97930_VisualCortex_Seurat.png"), device=)

pca_heatmap_VisualCortex <- DimHeatmap(VisualCortex, dims = 1, balanced = TRUE)
ggsave(path.expand("~/Lake/VisualCortex/pcaheatmap_GSE97930_VisualCortex_Seurat.png"), plot=pca_heatmap_VisualCortex, device=)

VisualCortex <- JackStraw(VisualCortex, num.replicate = 100)
VisualCortex <- ScoreJackStraw(VisualCortex, dims = 1:20)
Jackstraw_VisualCortex <- JackStrawPlot(VisualCortex, dims = 1:15)
ggsave(path.expand("~/Lake/VisualCortex/jackstrawplot_GSE97930_VisualCortex_Seurat.png"), device=)

#VisualCortex <- FindNeighbors(VisualCortex, dims = 1:10)
#VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

VisualCortex.markers <- FindAllMarkers(VisualCortex, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = VisualCortex.markers %>% group_by(cluster)

saveRDS(VisualCortex, file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_orgident.png"), device=)

featureplot_VisualCortex < - FeaturePlot(VisualCortex, features = c("SYT1"))
ggsave(path.expand("~/Lake/VisualCortex/featureplot_GSE97930_VisualCortex_Seurat_default.png"), device=)

write.table(diff_expressed, file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_differentiallyexpressed_ident.txt"), sep="\t")
