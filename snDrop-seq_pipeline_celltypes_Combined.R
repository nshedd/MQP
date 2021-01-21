library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Download table
matrix = readRDS(path.expand("~/GSE97930_All.RDS"))

# Use Seurat Pipeline with appropriate parameters
All <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^MT-")

All <- subset(All, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

All <- NormalizeData(All, normalization.method = "LogNormalize", scale.factor = 10000)

All <- FindVariableFeatures(object = All)

all_cells <- rownames(All)
All <- ScaleData(All, features = all_cells)

All <- RunPCA(All, features = VariableFeatures(object = All))

All <- FindNeighbors(All, dims = 1:20)
All <- FindClusters(All, resolution = 1.5)

All <- RunUMAP(All, dims = 1:20, metric="euclidean")

#Save Seurat object to create heatmap later
saveRDS(All, file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#Plot UMAP
plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_findct.png"), device=)

#Find overexpressed genes
All.markers <- FindAllMarkers(All, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All.markers %>% group_by(cluster)

# Compare to lab marker gene list and create/save csv of matches
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = All.markers
celltypelist = brain_genes$Cell.type

celltypes <- character()
for (gene in diff_expressed$gene) {
  if (gene %in% brain_genes$Human.Gene) {
    celltype = brain_genes$Cell.type[brain_genes$Human.Gene == gene]
    celltype = paste(celltype, collapse = ', ')
    celltypes <- c(celltypes, celltype)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

diff_expressed_condensed = diff_expressed[diff_expressed$cell_type != "unknown",]

#List of all overexpressed genes, with potential cell type of "unknown" if there were no matches
write.table(diff_expressed, file = path.expand("~/Lake/All/umap_GSE97930_All_Seurat_diffexpressed.txt"), sep = '\t')
#List of only overexpressed genes that matched marker genes
write.table(diff_expressed_condensed, file = path.expand("~/Lake/All/umap_GSE97930_All_Seurat_diffexpressed_condensed.txt"), sep = '\t')

#Create feature plots of marker genes for each cell type
FeaturePlot(All, features = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_ExcitatoryFeatures.png"), device=)

FeaturePlot(All, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_InhibitoryFeatures.png"), device=)

FeaturePlot(All, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_OligodendrocyteFeatures.png"), device=)

FeaturePlot(All, features = c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_AstrocyteFeatures.png"), device=)

FeaturePlot(All, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_EndothelialFeatures.png"), device=)

FeaturePlot(All, features = c("COBLL1", "PDGFRB", "COBLL1", "PDGFRB"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_PericyteFeatures.png"), device=)

FeaturePlot(All, features = c("PCDH15", "OLIG1", "PCDH15", "OLIG1"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_OPCFeatures.png"), device=)

FeaturePlot(All, features = c("APBB1IP", "P2RY12", "APBB1IP", "P2RY12"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_MicrogliaFeatures.png"), device=)

FeaturePlot(All, features = c("RYR1"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_PurkinjeFeatures.png"), device=)

FeaturePlot(All, features = c("RELN", "GRM4", "RBFOX3"))
ggsave(path.expand("~/Lake/All/GSE97930_All_Seurat_GranuleFeatures.png"), device=)
