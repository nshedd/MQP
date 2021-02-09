library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/matrix.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

UMB5168 <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

UMB5168[["percent.mt"]] <- PercentageFeatureSet(UMB5168, pattern = "^MT-")

UMB5168 <- subset(UMB5168, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

UMB5168 <- NormalizeData(UMB5168, normalization.method = "LogNormalize", scale.factor = 10000)

UMB5168 <- FindVariableFeatures(object = UMB5168)

all_cells <- rownames(UMB5168)
UMB5168 <- ScaleData(UMB5168, features = all_cells)

UMB5168 <- RunPCA(UMB5168, features = VariableFeatures(object = UMB5168))

saveRDS(UMB5168, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/Matrix_Seurat.rds"))

UMB5168 <- FindNeighbors(UMB5168, dims = 1:20)
UMB5168 <- FindClusters(UMB5168, resolution = 1)

UMB5168 <- RunUMAP(UMB5168, dims = 1:20, metric="euclidean")

plot = DimPlot(UMB5168, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/umap_numberedclusters.png"), device=)

UMB5168.markers <- FindAllMarkers(UMB5168, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UMB5168.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(UMB5168.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(UMB5168, features = intersection) + 
   theme(axis.text.x = element_text(angle = 90)) + 
   scale_y_discrete(limits = rev(levels(UMB5168$seurat_clusters)))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/dotplot.png"), device=)



brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = UMB5168.markers

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


write.table(diff_expressed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(UMB5168, features = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/ExcitatoryFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/InhibitoryFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/OligodendrocyteFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/AstrocyteFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/EndothelialFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("PCDH15", "OLIG1", "PCDH15", "OLIG1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/OPCFeatures.png"), device=)

FeaturePlot(UMB5168, features = c("APBB1IP", "P2RY12", "APBB1IP", "P2RY12"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/analysis/MicrogliaFeatures.png"), device=)
