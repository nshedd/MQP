library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT-")
CerebellarHem_VlnPlot <- VlnPlot(CerebellarHem, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("~/Lake/CerebellarHem/qcvlnplot_GSE97930_CerebellarHem_Seurat.png"), device=)

plot1 <- FeatureScatter(CerebellarHem, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CerebellarHem, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("~/Lake/CerebellarHem/featurescatter_GSE97930_CerebellarHem_Seurat.png"), device=)

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)

CerebellarHem <- FindVariableFeatures(object = CerebellarHem)

all_cells <- rownames(CerebellarHem)
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem))

pcs_plot_CerebellarHem <- VizDimLoadings(CerebellarHem, dims = 1:2, reduction = "pca")
ggsave(path.expand("~/Lake/CerebellarHem/qcdimloadings_GSE97930_CerebellarHem_Seurat.png"), device=)

pca_plot_CerebellarHem <- DimPlot(CerebellarHem, reduction = "pca")
ggsave(path.expand("~/Lake/CerebellarHem/pcascatter_GSE97930_CerebellarHem_Seurat.png"), device=)

pca_heatmap_CerebellarHem <- DimHeatmap(CerebellarHem, dims = 1, cells = 500, balanced = TRUE)
ggsave(path.expand("~/Lake/CerebellarHem/pcaheatmap_GSE97930_CerebellarHem_Seurat.png"), device=)

CerebellarHem <- JackStraw(CerebellarHem, num.replicate = 100)
CerebellarHem <- ScoreJackStraw(CerebellarHem, dims = 1:20)
Jackstraw_CerebellarHem <- JackStrawPlot(CerebellarHem, dims = 1:15)
ggsave(path.expand("~/Lake/CerebellarHem/jackstrawplot_GSE97930_CerebellarHem_Seurat.png"), device=)

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:10)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(CerebellarHem)
#CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = CerebellarHem.markers %>% group_by(cluster)

saveRDS(CerebellarHem, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_default.png"), device=)

featureplot_CerebellarHem < - FeaturePlot(CerebellarHem, features = c("SYT1", "SYT1", "RBFOX3", "GAD2", "SLC1A3", "GRIN2B", "PCDH15", "MBP", "APBB1IP", "SLC6A1", "PCDH15", "SLC1A3"))
ggsave(path.expand("~/Lake/CerebellarHem/featureplot_GSE97930_CerebellarHem_Seurat_default.png"), device=)
  
write.table(diff_expressed, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_differentiallyexpressed.txt"), sep="\t")
