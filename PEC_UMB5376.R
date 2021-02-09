library(Seurat)
library(ggplot2)
library(dplyr)

print("loaded libraries")

path1 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/matrix.tsv")

matrix = read.table(path1, header=TRUE, row.names=1, sep='\t')

UMB5376 <- CreateSeuratObject(counts = matrix, project = "PEC", min.cells = 3, min.features = 200)

UMB5376[["percent.mt"]] <- PercentageFeatureSet(UMB5376, pattern = "^MT-")

UMB5376 <- NormalizeData(UMB5376, normalization.method = "LogNormalize", scale.factor = 10000)

UMB5376 <- FindVariableFeatures(object = UMB5376)

all_cells <- rownames(UMB5376)
UMB5376 <- ScaleData(UMB5376, features = all_cells)

UMB5376 <- RunPCA(UMB5376, features = VariableFeatures(object = UMB5376))

ElbowPlot(UMB5376)
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/ElbowPlot.png"), device=)

#saveRDS(UMB5376, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/Matrix_Seurat.rds"))

UMB5376 <- FindNeighbors(UMB5376, dims = 1:20)
UMB5376 <- FindClusters(UMB5376, resolution = 0.5)

UMB5376 <- RunUMAP(UMB5376, dims = 1:20, metric="euclidean", n.neighbors=10)

plot = DimPlot(UMB5376, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/umap_numberedclusters.png"), device=)

UMB5376.markers <- FindAllMarkers(UMB5376, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UMB5376.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(UMB5376.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(UMB5376, features = intersection) + 
   theme(axis.text.x = element_text(angle = 90)) + 
   scale_y_discrete(limits = rev(levels(UMB5376$seurat_clusters)))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/dotplot.png"), device=)



brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = UMB5376.markers

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


write.table(diff_expressed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(UMB5376, features = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/ExcitatoryFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/InhibitoryFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OligodendrocyteFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/AstrocyteFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/EndothelialFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("PCDH15", "OLIG1", "PCDH15", "OLIG1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OPCFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("APBB1IP", "P2RY12", "APBB1IP", "P2RY12"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/MicrogliaFeatures.png"), device=)
