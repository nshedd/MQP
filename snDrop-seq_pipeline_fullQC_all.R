library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


matrix = readRDS(path.expand("~/GSE97930_All.RDS"))

All <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^MT-")
All_VlnPlot <- VlnPlot(All, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("~/Lake/All/qcvlnplot_GSE97930_All_Seurat.png"), device=, width = 21, height = 7)

plot1 <- FeatureScatter(All, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(All, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("~/Lake/All/featurescatter_GSE97930_All_Seurat.png"), device, width = 14, height = 7)

All <- NormalizeData(All, normalization.method = "LogNormalize", scale.factor = 10000)

All <- FindVariableFeatures(object = All)

all_cells <- rownames(All)
All <- ScaleData(All, features = all_cells)

All <- RunPCA(All, features = VariableFeatures(object = All))

pcs_plot_All <- VizDimLoadings(All, dims = 1:2, reduction = "pca")
ggsave(path.expand("~/Lake/All/qcdimloadings_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)

pca_plot_All <- DimPlot(All, reduction = "pca")
ggsave(path.expand("~/Lake/All/pcascatter_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)

pca_heatmap_All <- DimHeatmap(All, dims = 1, balanced = TRUE)
ggsave(path.expand("~/Lake/All/pcaheatmap_GSE97930_All_Seurat.png"), plot=pca_heatmap_All, device=, width = 14, height = 7)

All <- JackStraw(All, num.replicate = 100)
All <- ScoreJackStraw(All, dims = 1:20)
Jackstraw_All <- JackStrawPlot(All, dims = 1:15)
ggsave(path.expand("~/Lake/All/jackstrawplot_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)

Elbow_All<- ElbowPlot(All)
ggsave(path.expand("~/Lake/All/elbowplot_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)

#All <- FindNeighbors(All, dims = 1:10)
#All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(All)
#All <- RenameIdents(All, new.cluster.ids)

All.markers <- FindAllMarkers(All, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = All.markers %>% group_by(cluster)

saveRDS(All, file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_origident.png"), device=)

#featureplot_All < - FeaturePlot(All, features = c("SYT1", "SYT1", "RBFOX3", "GAD2", "SLC1A3", "GRIN2B", "PCDH15", "MBP", "APBB1IP", "SLC6A1", "PCDH15", "SLC1A3"))
#ggsave(path.expand("~/Lake/All/featureplot_GSE97930_All_Seurat_default.png"), device=)

write.table(diff_expressed, file = path.expand("~/Lake/All/GSE97930_All_differentiallyexpressed_origident.txt"), sep="\t")
